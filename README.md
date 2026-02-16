# Nanopore Gene Methylation Viz

**Short description:** Web app to visualize nanopore methylation (modkit BED) along gene structure—promoter, exon, intron, CDS, downstream—with built-in GFF annotation. Upload BED, search by gene name, view interactive plots with strand and exon/CDS counts.

modkit 出力 **BED のみ**をアップロードし、遺伝子名で検索して**プロモーター（上流 2k）・エクソン・イントロン・CDS・下流**のメチル化をインタラクティブグラフで表示する Web アプリ。**アノテーション（GFF）はアプリに内蔵**し、事前処理キャッシュで起動が速い。

## 必要な環境

- Python 3.8+
- pandas, numpy, pyranges, fastapi, uvicorn

## セットアップ

**推奨**: Jupyter / spacy 等と Pydantic のバージョンが衝突するため、**このプロジェクト専用の仮想環境**を作成してから入れると安全です。

```bash
cd meth_viz_nanopore
python3 -m venv meth_viz_nanopore_env
source meth_viz_nanopore_env/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

仮想環境を使わない場合:

```bash
cd meth_viz_nanopore
pip install -r requirements.txt
```

### GFF の内蔵（初回のみ）

アプリは **data/gene_regions.pkl** から遺伝子領域を読み込みます。このファイルは **genomic.gff** を事前に処理して作成します。

1. **genomic.gff** を `data/` に置く（プロジェクトルートの `genomic.gff` をコピーしても可）。
2. 次のコマンドを実行（数分かかります。約 4 万遺伝子分の領域を抽出します）。

```bash
python scripts/build_gene_regions.py
# または GFF のパスを指定:
python scripts/build_gene_regions.py /path/to/genomic.gff
```

これで `data/gene_regions.pkl` が作成され、アプリ起動時に自動で読み込まれます。**genomic.gff はそのまま data/ に保管**しておけば、パラメータ（プロモーター長など）を変えて再実行するときに使えます。

## 起動方法

**仮想環境を使っている場合は、必ず activate してから、その Python で uvicorn を実行してください**（そうしないと conda base の Python が使われて ImportError になることがあります）。

```bash
source meth_viz_nanopore_env/bin/activate 
python -m uvicorn app.main:app --reload --port 8000
```

ブラウザで **http://localhost:8000** を開く。

## 使い方

1. **ファイルアップロード**  
   - **modkit BED** のみ選択: 例 `r0081_m.bed` または `r0081_m_chr.bed`（染色体名は `chr1` でも `NC_000067.7` など RefSeq でも可）  
   - 「アップロード」をクリック。（アノテーションはアプリ内蔵のため GFF のアップロードは不要）

2. **遺伝子検索**  
   - 検索ボックスに遺伝子名（例: `Xkr4`, `Bdnf`）を入力。  
   - 候補から選択するか、そのまま「グラフ表示」をクリック。

3. **グラフ**  
   - X 軸: ゲノム位置  
   - Y 軸: メチル化率（%）  
   - 色付きの縦帯: プロモーター / エクソン / イントロン / 下流領域（凡例あり）  
   - ホバーで各サイトの位置・メチル化率・カバレッジを表示。

## プロジェクト構成

- `PROJECT_PLAN.md` — 開発計画（Step 1〜4）
- `gene_methylation.py` — GFF から遺伝子領域抽出、1 遺伝子分のメチル化サイト取得（プロット用）
- `scripts/build_gene_regions.py` — genomic.gff を処理して `data/gene_regions.pkl` を生成（初回・更新時のみ実行）
- `data/genomic.gff` — アノテーション GFF（手元に保管。上記スクリプトの入力）
- `data/gene_regions.pkl` — 事前処理した遺伝子領域キャッシュ（アプリ起動時に読み込み）
- `app/main.py` — FastAPI: BED アップロード・遺伝子一覧・遺伝子メチル化 API
- `app/static/index.html` — フロント（BED アップロード UI、遺伝子検索、Plotly グラフ）
- `requirements.txt` — Python 依存関係

## API（参考）

- `POST /api/upload` — BED のみ multipart でアップロード
- `GET /api/genes` — 遺伝子名一覧（内蔵 gene_regions から）
- `GET /api/gene/{gene_id}` — 指定遺伝子のメチル化サイトと領域情報（JSON）

## 注意

- BED の染色体名が `chr1`, `chr2` … で、GFF が RefSeq (`NC_000067.7` 等) の場合は、`gene_methylation.py` 内の `CHROM_GFF_TO_CHR` で対応しています（GRCm39 想定）。
- アップロードしたファイルはサーバーの一時ディレクトリに保存され、再起動や再アップロードで上書きされます。
