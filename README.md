# Nanopore Gene Methylation Viz (methscope)

**English:** Web app to visualize nanopore methylation (modkit BED) along gene structure—promoter, exon, intron, CDS, downstream—with **built-in GFF annotation** (pre-processed cache). Upload BED, search by gene, view interactive Plotly charts.

**日本語:** modkit 出力 **BED のみ**をアップロードし、遺伝子名で検索してメチル化を表示。**GFF はアプリ内蔵**（事前ビルドした `gene_list.json` + `gene_regions.db`）。

---

## 必要な環境

- Python 3.11 推奨（3.8+）
- 依存: `requirements.txt` を参照（FastAPI, pandas, pyranges, boto3 など）

## セットアップ

```bash
cd methscope
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

**Pydantic 2 が必要**です。conda base と混ぜると `IncEx` エラーになることがあるので、**必ずこの venv の Python で** `python -m uvicorn ...` を実行してください。

## 参照データのビルド（初回）

`data/genomic.gff` を置き、次を実行（数分かかります）。

```bash
python scripts/build_gene_regions.py
```

出力: `data/gene_list.json`, `data/gene_regions.db`, `data/gene_regions.pkl`  
起動時は **gene_list.json + gene_regions.db** があれば遅延読み込み（低メモリ・Render 無料枠向け）。無い場合は `gene_regions.pkl` をまとめて読み込み。

## ローカル起動

```bash
source .venv/bin/activate
python -m uvicorn app.main:app --reload --port 8000
```

→ **http://localhost:8000**

## Docker

```bash
docker build -t methscope .
docker run -p 8000:8000 -v "$(pwd)/data:/app/data" methscope
```

`data/` に `gene_list.json` と `gene_regions.db` をマウントする。コンテナ内 `PORT` は環境変数に対応。

---

## デプロイ（Render など）

**向かない:** Vercel 等のサーバーレス（常時プロセス・大きな参照データが前提）。

**向いている:** Render / Railway / Fly.io など常時起動の Web サービス。

| 項目 | 例 |
|------|-----|
| Build | `pip install -r requirements.txt` または `bash scripts/render_build.sh`（GFF をビルドで取る場合） |
| Start | `uvicorn app.main:app --host 0.0.0.0 --port $PORT` |

**参照データの渡し方（どちらか）**

1. **Cloudflare R2 または AWS S3**（推奨・ビルド短縮）  
   環境変数で `R2_*` または `S3_*` を設定。起動時に `methscope-data/reference/` から `gene_list.json` と `gene_regions.db` を取得。`BUILD_GFF_URL` は不要。

2. **ビルド時に GFF から生成**  
   `BUILD_GFF_URL` に genomic.gff の直接ダウンロード URL を設定し、`bash scripts/render_build.sh` でビルド（10〜15 分程度のことあり）。

無料枠（512 MB）では遅延読み込み + 非同期 BED ジョブを利用。スリープ・初回応答遅延は Render 無料枠の仕様です。

---

## オブジェクトストレージ（任意・R2 / S3）

アップロードした BED をクラウドに置き、参照データもクラウドから取得できます。

**キー構成（共通）**

```
methscope-data/
  reference/   gene_list.json, gene_regions.db（必須）, genomic.gff（任意）
  uploads/     {job_id}.bed（アプリが保存）
```

**アップロード（ローカル → バケット）**

```bash
python scripts/build_gene_regions.py
export R2_BUCKET=... R2_ACCOUNT_ID=... R2_ACCESS_KEY_ID=... R2_SECRET_ACCESS_KEY=...
# または S3: S3_BUCKET, AWS_ACCESS_KEY_ID, AWS_SECRET_ACCESS_KEY, AWS_REGION
python scripts/upload_reference_to_s3.py
```

**アプリ側の環境変数**

| R2 | S3 |
|----|-----|
| `R2_BUCKET`, `R2_ACCOUNT_ID`, `R2_ACCESS_KEY_ID`, `R2_SECRET_ACCESS_KEY` | `S3_BUCKET`, `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY`, `AWS_REGION`（任意） |

R2 の詳細手順は Cloudflare ダッシュボードで **R2 → バケット作成 → Manage R2 API Tokens（Read & Write）** を参照。

---

## CI/CD

`.github/workflows/deploy.yml`: push/PR で依存インストール・インポート確認・`docker build`。`main` への push 後、Secret `RENDER_DEPLOY_HOOK_URL` があれば Render の Deploy Hook を叩く。

---

## 使い方（UI）

1. modkit **BED** をアップロード（処理完了までポーリング）。
2. 遺伝子名を検索 → **グラフ表示**。

---

## プロジェクト構成

```
methscope/
├── app/
│   ├── main.py           # FastAPI
│   └── static/index.html
├── gene_methylation.py
├── requirements.txt
├── Dockerfile
├── scripts/
│   ├── build_gene_regions.py
│   ├── render_build.sh
│   └── upload_reference_to_s3.py
├── data/                 # ローカル用（大きいファイルは .gitignore）
└── .github/workflows/deploy.yml
```

---

## API（参考）

| メソッド | パス | 説明 |
|----------|------|------|
| GET | `/health` | ヘルスチェック |
| GET | `/api/status` | gene_regions / BED / current_job_id |
| POST | `/api/upload` | BED → `{ job_id, status: "processing" }` |
| GET | `/api/jobs/{job_id}` | `processing` / `ready` / `failed` |
| GET | `/api/genes` | 遺伝子一覧、`?q=` でフィルタ |
| GET | `/api/gene/{gene_id}` | プロット用 JSON |

Swagger: `/docs`

---

## 注意

- BED の `chr1` と GFF の RefSeq（例 `NC_000067.7`）は `gene_methylation.py` の `CHROM_GFF_TO_CHR`（GRCm39 想定）で対応。
- BED はローカル `data/uploads/` または **R2/S3 設定時はクラウド**に保存。ジョブ完了後はメモリ上の DataFrame で参照。
