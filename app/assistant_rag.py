"""
MethScope AI Assistant: RAG over PubMed + Semantic Scholar, generation via Cloudflare Workers AI (REST API).
"""
from __future__ import annotations

import asyncio
import logging
import os
import re
import xml.etree.ElementTree as ET
from typing import Any

import httpx

logger = logging.getLogger(__name__)

# NCBI recommends tool + email for e-utils
EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
S2_SEARCH = "https://api.semanticscholar.org/graph/v1/paper/search"
def _cf_ai_url(account_id: str, model: str) -> str:
    # モデル名はパスにそのまま（例: @cf/meta/llama-3.1-8b-instruct）
    return f"https://api.cloudflare.com/client/v4/accounts/{account_id}/ai/run/{model}"

DEFAULT_MODEL = "@cf/meta/llama-3.1-8b-instruct"
MAX_CONTEXT_CHARS = 12000
# 検索は多めに取り、リランキング後に上位だけをコンテキストに載せる
MAX_PUBMED_FETCH = 10
MAX_S2_FETCH = 10
MAX_PAPERS_AFTER_RERANK = 10

S2_MAX_RETRIES = 4
S2_BACKOFF_BASE_SEC = 1.5


def _cf_credentials() -> tuple[str | None, str | None]:
    account = os.environ.get("CLOUDFLARE_ACCOUNT_ID") or os.environ.get("CF_ACCOUNT_ID")
    token = os.environ.get("CLOUDFLARE_API_TOKEN") or os.environ.get("CF_API_TOKEN")
    return account, token


def cf_workers_ai_configured() -> bool:
    a, t = _cf_credentials()
    return bool(a and t)


def _build_search_terms(gene: str, species: str, sample: str, extra: str = "") -> str:
    g = gene.strip()
    if not g:
        return "DNA methylation"
    # 遺伝子名 + メチル化系キーワード（PubMed は AND を暗黙にしないため明示）
    term = f'({g} OR "{g.upper()}") AND (methylation OR "DNA methylation" OR CpG OR epigenetic OR promoter)'
    sp = species.strip()
    if sp:
        sl = sp.lower()
        if "musculus" in sl or sl in ("mouse", "mouse (mus musculus)"):
            term += " AND (mouse OR \"mus musculus\" OR \"Mus musculus\")"
        elif "sapiens" in sl or sl in ("human", "homo sapiens"):
            term += " AND (human OR \"Homo sapiens\" OR hsa)"
        else:
            term += f" AND ({sp})"
    if sample.strip():
        term += f" AND ({sample.strip()})"
    if extra.strip():
        term += f" AND ({extra.strip()})"
    return term


async def pubmed_search_and_fetch(
    client: httpx.AsyncClient, term: str, retmax: int = MAX_PUBMED_FETCH
) -> list[dict[str, Any]]:
    """PubMed: esearch -> efetch XML abstracts."""
    params = {
        "db": "pubmed",
        "term": term,
        "retmax": str(retmax),
        "retmode": "json",
        "sort": "relevance",
        "tool": "methscope",
        "email": os.environ.get("NCBI_EMAIL", "methscope@local"),
    }
    r = await client.get(f"{EUTILS_BASE}/esearch.fcgi", params=params, timeout=60.0)
    r.raise_for_status()
    data = r.json()
    idlist = data.get("esearchresult", {}).get("idlist") or []
    if not idlist:
        return []

    fetch_params = {
        "db": "pubmed",
        "id": ",".join(idlist),
        "retmode": "xml",
        "tool": "methscope",
        "email": os.environ.get("NCBI_EMAIL", "methscope@local"),
    }
    r2 = await client.get(f"{EUTILS_BASE}/efetch.fcgi", params=fetch_params, timeout=90.0)
    r2.raise_for_status()
    return _parse_pubmed_xml(r2.text)


def _local_tag(tag: str) -> str:
    return tag.split("}")[-1] if "}" in tag else tag


def _parse_pubmed_xml(xml_text: str) -> list[dict[str, Any]]:
    root = ET.fromstring(xml_text)
    out: list[dict[str, Any]] = []
    for article in root.iter():
        if _local_tag(article.tag) != "PubmedArticle":
            continue
        pmid = ""
        title = ""
        abstract_parts: list[str] = []
        medline = None
        for el in article:
            if _local_tag(el.tag) == "MedlineCitation":
                medline = el
                break
        scope = medline if medline is not None else article
        if medline is not None:
            for el in medline:
                if _local_tag(el.tag) == "PMID" and el.text:
                    pmid = el.text.strip()
                    break
        for el in scope.iter():
            ln = _local_tag(el.tag)
            if ln == "ArticleTitle":
                title = "".join(el.itertext()).strip()
            elif ln == "AbstractText":
                label = el.get("Label")
                text = "".join(el.itertext()).strip()
                if label:
                    abstract_parts.append(f"{label}: {text}")
                else:
                    abstract_parts.append(text)
        abstract = " ".join(abstract_parts) if abstract_parts else ""
        if not title and not abstract:
            continue
        out.append(
            {
                "source": "PubMed",
                "id": pmid,
                "title": title,
                "abstract": abstract[:4000],
                "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/" if pmid else "",
            }
        )
    return out


def _s2_parse_response(data: dict[str, Any]) -> list[dict[str, Any]]:
    papers = data.get("data") or []
    out: list[dict[str, Any]] = []
    for p in papers:
        pid = p.get("paperId") or ""
        title = (p.get("title") or "").strip()
        abstract = (p.get("abstract") or "").strip()[:4000]
        if not title and not abstract:
            continue
        authors = p.get("authors") or []
        auth_s = ", ".join(a.get("name", "") for a in authors[:5])
        if len(authors) > 5:
            auth_s += " et al."
        out.append(
            {
                "source": "Semantic Scholar",
                "id": pid,
                "title": title,
                "abstract": abstract,
                "year": p.get("year"),
                "authors": auth_s,
                "url": p.get("url") or (f"https://www.semanticscholar.org/paper/{pid}" if pid else ""),
            }
        )
    return out


async def semantic_scholar_search(
    client: httpx.AsyncClient, query: str, limit: int = MAX_S2_FETCH
) -> list[dict[str, Any]]:
    """Semantic Scholar: 429 / 5xx / タイムアウト時は指数バックオフで再試行。失敗時は空リスト（PubMed のみで続行）。"""
    params = {
        "query": query,
        "limit": str(limit),
        "fields": "title,abstract,url,year,authors,paperId",
    }
    headers: dict[str, str] = {}
    key = os.environ.get("SEMANTIC_SCHOLAR_API_KEY")
    if key:
        headers["x-api-key"] = key

    delay = S2_BACKOFF_BASE_SEC
    last_err: str | None = None
    for attempt in range(S2_MAX_RETRIES):
        try:
            r = await client.get(S2_SEARCH, params=params, headers=headers, timeout=60.0)
            if r.status_code == 429:
                logger.warning(
                    "[assistant] Semantic Scholar 429 (attempt %d/%d), retry in %.1fs",
                    attempt + 1,
                    S2_MAX_RETRIES,
                    delay,
                )
                await asyncio.sleep(delay)
                delay *= 2
                continue
            if r.status_code >= 500:
                logger.warning(
                    "[assistant] Semantic Scholar %s (attempt %d/%d)",
                    r.status_code,
                    attempt + 1,
                    S2_MAX_RETRIES,
                )
                await asyncio.sleep(delay)
                delay *= 2
                continue
            r.raise_for_status()
            return _s2_parse_response(r.json())
        except httpx.TimeoutException as e:
            last_err = str(e)
            logger.warning("[assistant] Semantic Scholar timeout (attempt %d/%d)", attempt + 1, S2_MAX_RETRIES)
            await asyncio.sleep(delay)
            delay *= 2
        except httpx.HTTPError as e:
            last_err = str(e)
            logger.warning("[assistant] Semantic Scholar HTTP error: %s", e)
            await asyncio.sleep(delay)
            delay *= 2

    if last_err:
        logger.warning("[assistant] Semantic Scholar gave up: %s", last_err)
    return []


def _dedupe_by_title(items: list[dict[str, Any]]) -> list[dict[str, Any]]:
    seen: set[str] = set()
    out: list[dict[str, Any]] = []
    for it in items:
        key = re.sub(r"\s+", " ", (it.get("title") or "").lower().strip())
        if len(key) < 8:
            key = f"{key}:{it.get('id')}"
        if key in seen:
            continue
        seen.add(key)
        out.append(it)
    return out


def _norm_gene_token(gene: str) -> str:
    g = gene.strip()
    if g.lower().startswith("gene-"):
        g = g[5:]
    return g.lower()


def _query_terms_for_rerank(
    gene: str, species: str, sample: str, extra_query: str = ""
) -> set[str]:
    """遺伝子・種・サンプル・メチル化関連語をトークン化（リランキング用）。"""
    q = f"{gene} {species} {sample} {extra_query} methylation DNA epigenetic CpG promoter"
    terms: set[str] = set()
    for w in re.split(r"[^\w]+", q.lower()):
        if len(w) >= 2:
            terms.add(w)
    gn = _norm_gene_token(gene)
    if gn:
        terms.add(gn)
        # Bdnf など先頭のみ大文字の別表記
        if gn.isalpha():
            terms.add(gn.capitalize())
    return terms


def _score_paper_relevance(p: dict[str, Any], terms: set[str], gene_norm: str) -> float:
    t = (p.get("title") or "").lower()
    a = (p.get("abstract") or "").lower()
    score = 0.0
    for term in terms:
        if len(term) < 2:
            continue
        if term in t:
            score += 3.0
        elif term in a:
            score += 1.0
    if gene_norm:
        if gene_norm in t:
            score += 10.0
        elif gene_norm in a:
            score += 5.0
    return score


def rerank_papers(
    papers: list[dict[str, Any]],
    gene: str,
    species: str,
    sample: str,
    extra_query: str = "",
    top_n: int = MAX_PAPERS_AFTER_RERANK,
) -> list[dict[str, Any]]:
    """キーワード重なりで関連度スコアを付け、上位のみを RAG に渡す。"""
    if not papers:
        return []
    terms = _query_terms_for_rerank(gene, species, sample, extra_query)
    gn = _norm_gene_token(gene)
    scored = [( _score_paper_relevance(p, terms, gn), p) for p in papers]
    scored.sort(key=lambda x: -x[0])
    return [p for _, p in scored[:top_n]]


def _collect_allowed_ids(papers: list[dict[str, Any]]) -> tuple[set[str], set[str]]:
    pmids: set[str] = set()
    s2_ids: set[str] = set()
    for p in papers:
        src = p.get("source") or ""
        pid = str(p.get("id") or "").strip()
        if not pid:
            continue
        if src == "PubMed":
            pmids.add(pid)
        elif src == "Semantic Scholar":
            s2_ids.add(pid)
    return pmids, s2_ids


def _citation_whitelist_prompt_block(pmids: set[str], s2_ids: set[str]) -> str:
    """プロンプトに埋め込み、モデルが提示文献外の ID を引用しないようにする。"""
    pl = ", ".join(sorted(pmids)) if pmids else "（なし）"
    sl = ", ".join(sorted(s2_ids)) if s2_ids else "（なし）"
    return f"""【引用してよい ID の一覧（これ以外の PubMed PMID や Semantic Scholar paperId を本文では引用しない）】
- PubMed PMID: {pl}
- Semantic Scholar paperId: {sl}
引用形式は [PubMed:PMID] または [S2:paperId] のみ。上記に無い数字や ID は書かない。"""


def sanitize_citations(text: str, pmids: set[str], s2_ids: set[str]) -> tuple[str, list[str]]:
    """生成文から、提示文献に無い引用タグを除去し注意メッセージを付与。"""
    notes: list[str] = []

    def bad_pmid(m: re.Match) -> str:
        pid = m.group(1).strip()
        if pid in pmids:
            return m.group(0)
        notes.append(f"PubMed:{pid}")
        return ""

    def bad_s2(m: re.Match) -> str:
        sid = m.group(1).strip()
        if sid in s2_ids:
            return m.group(0)
        notes.append(f"S2:{sid[:24]}…" if len(sid) > 24 else f"S2:{sid}")
        return ""

    out = text
    out = re.sub(r"\[PubMed:\s*(\d+)\s*\]", bad_pmid, out)
    out = re.sub(r"\[S2:\s*([^\]]+?)\s*\]", bad_s2, out)
    # 空行の連続を軽く整理
    out = re.sub(r"\n{3,}", "\n\n", out).strip()
    if notes:
        uniq = sorted(set(notes))
        out += "\n\n※ 提示文献に含まれない引用 " + str(len(uniq)) + " 件を本文から除去しました。"
    return out, notes


def _format_rag_context(papers: list[dict[str, Any]]) -> str:
    chunks: list[str] = []
    total = 0
    for i, p in enumerate(papers, 1):
        block = (
            f"[{i}] {p.get('source')} ID={p.get('id')}\n"
            f"Title: {p.get('title')}\n"
        )
        if p.get("authors"):
            block += f"Authors: {p.get('authors')}\n"
        if p.get("year"):
            block += f"Year: {p.get('year')}\n"
        ab = p.get("abstract") or ""
        if ab:
            block += f"Abstract: {ab}\n"
        if total + len(block) > MAX_CONTEXT_CHARS:
            break
        chunks.append(block)
        total += len(block)
    return "\n---\n".join(chunks)


async def retrieve_literature(
    gene: str,
    species: str,
    sample: str,
    extra_query: str = "",
) -> list[dict[str, Any]]:
    term = _build_search_terms(gene, species, sample, extra_query)
    s2_query = term
    async with httpx.AsyncClient() as client:
        pubmed_papers, s2_papers = await asyncio.gather(
            pubmed_search_and_fetch(client, term, MAX_PUBMED_FETCH),
            semantic_scholar_search(client, s2_query, MAX_S2_FETCH),
        )

    merged = _dedupe_by_title(pubmed_papers + s2_papers)
    if not merged:
        return []
    return rerank_papers(merged, gene, species, sample, extra_query, MAX_PAPERS_AFTER_RERANK)


async def workers_ai_run(
    system_prompt: str,
    user_prompt: str,
) -> str:
    account, token = _cf_credentials()
    if not account or not token:
        raise RuntimeError(
            "Cloudflare Workers AI が未設定です。環境変数 CLOUDFLARE_ACCOUNT_ID と "
            "CLOUDFLARE_API_TOKEN（Workers AI 権限）を設定してください。"
        )

    model = os.environ.get("CLOUDFLARE_WORKERS_AI_MODEL", DEFAULT_MODEL)
    url = _cf_ai_url(account, model)

    # REST API: prompt 形式（公式ドキュメント準拠）
    # max_tokens 未指定だと出力が短く途中で切れやすい（モデル既定の生成上限）
    combined = f"[System]\n{system_prompt}\n\n[User]\n{user_prompt}"
    max_tokens = int(os.environ.get("CLOUDFLARE_WORKERS_AI_MAX_TOKENS", "3072"))
    payload: dict[str, Any] = {"prompt": combined, "max_tokens": max(256, min(max_tokens, 8192))}

    headers = {
        "Authorization": f"Bearer {token}",
        "Content-Type": "application/json",
    }

    async with httpx.AsyncClient() as client:
        r = await client.post(url, json=payload, headers=headers, timeout=120.0)

    if r.status_code >= 400:
        logger.error("[assistant] Workers AI error %s: %s", r.status_code, r.text[:500])
        r.raise_for_status()

    body = r.json()
    if not body.get("success"):
        errs = body.get("errors") or []
        raise RuntimeError(f"Workers AI: {errs}")

    result = body.get("result") or {}
    text = result.get("response")
    if not text:
        raise RuntimeError("Workers AI: empty response")
    usage = result.get("usage") or body.get("usage")
    if usage:
        logger.info("[assistant] Workers AI usage: %s", usage)
    return text.strip()


OVERVIEW_SYSTEM = """あなたは MethScope AI Assistant です。入力で与えられた文献スニペットのみに基づき、
指定された遺伝子・種・サンプル文脈での DNA メチル化パターンの生物学的解釈のたたき台を述べてください。

厳守:
- 文献に無い断定は避け、「一般的には」「本スニペットでは」など区別する。
- 引用はユーザーメッセージに記載された「引用してよい ID の一覧」に含まれるものだけを使う。
  形式は [PubMed:PMID] または [S2:paperId]。一覧に無い PMID や paperId は絶対に書かない。
- 出力は日本語で、簡潔な段落＋最後に「参考文献（提示）」として一覧（タイトル・出典・URL が分かるように）。
- メチル化が疾患・発達・転写制御とどう結びつきうかを文献ベースで述べる。
- 途中で文や箇条書きを切らず、段落の区切りをはっきりさせ、最後まで完結させる。"""


CHAT_SYSTEM = """あなたは MethScope AI Assistant です。与えられた文献コンテキストと会話履歴に基づき回答してください。
文献にない内容は推測と明示する。引用する場合はメッセージ内の「引用してよい ID の一覧」に含まれる ID のみを
[PubMed:PMID] または [S2:paperId] で示す。一覧に無い ID は書かない。日本語で答える。"""


async def generate_overview(
    species: str,
    sample: str,
    gene_name: str,
) -> dict[str, Any]:
    papers = await retrieve_literature(gene_name, species, sample, "")
    if not papers:
        return {
            "gene_name": gene_name,
            "species": species,
            "sample": sample,
            "interpretation": "関連文献が見つかりませんでした。検索語を変えるか、別の遺伝子で試してください。",
            "papers": [],
        }

    ctx = _format_rag_context(papers)
    pmids, s2_ids = _collect_allowed_ids(papers)
    whitelist = _citation_whitelist_prompt_block(pmids, s2_ids)
    user = f"""種: {species}
サンプル/条件: {sample or "（未指定）"}
対象遺伝子: {gene_name}

{whitelist}

以下は検索で得た文献スニペットです（関連度の高い順）:

{ctx}
"""
    interpretation = await workers_ai_run(OVERVIEW_SYSTEM, user)
    interpretation, _ = sanitize_citations(interpretation, pmids, s2_ids)
    return {
        "gene_name": gene_name,
        "species": species,
        "sample": sample,
        "interpretation": interpretation,
        "papers": [
            {
                "source": p["source"],
                "id": p["id"],
                "title": p.get("title"),
                "url": p.get("url"),
            }
            for p in papers
        ],
    }


async def generate_chat_reply(
    species: str,
    sample: str,
    gene_name: str,
    user_message: str,
    history: list[dict[str, str]] | None,
) -> dict[str, Any]:
    extra = user_message[:500]
    papers = await retrieve_literature(gene_name, species, sample, extra_query=extra)
    pmids, s2_ids = _collect_allowed_ids(papers) if papers else (set(), set())
    whitelist = _citation_whitelist_prompt_block(pmids, s2_ids) if papers else "（文献スニペットなしのため引用 ID なし）"
    ctx = _format_rag_context(papers) if papers else "（文献スニペットなし）"

    hist_lines = []
    if history:
        for turn in history[-10:]:
            role = turn.get("role", "user")
            content = (turn.get("content") or "").strip()
            if content:
                hist_lines.append(f"{role}: {content}")
    hist_block = "\n".join(hist_lines) if hist_lines else "（なし）"

    user = f"""種: {species}
サンプル: {sample or "（未指定）"}
遺伝子: {gene_name}

{whitelist}

直近の会話:
{hist_block}

文献コンテキスト:
{ctx}

ユーザーの質問:
{user_message}
"""
    reply = await workers_ai_run(CHAT_SYSTEM, user)
    reply, _ = sanitize_citations(reply, pmids, s2_ids)
    return {
        "reply": reply,
        "papers": [
            {
                "source": p["source"],
                "id": p["id"],
                "title": p.get("title"),
                "url": p.get("url"),
            }
            for p in papers
        ],
    }
