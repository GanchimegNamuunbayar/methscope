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

## デプロイ（本番環境）

このアプリは **常時起動の Python サーバー**と**メモリ上の状態**（起動時に読み込む `gene_regions.pkl`、アップロードした BED）に依存しているため、**Vercel のようなサーバーレス環境には向いていません**（インスタンスごとにメモリが分かれる・制限が厳しい・大容量 pkl のコールドスタートが重い）。

**推奨:** 次のような「常時稼働コンテナ / サーバー」向けのサービスでデプロイしてください。

| サービス | 手順の目安 |
|----------|------------|
| **[Railway](https://railway.app)** | リポジトリを連携 → ルートに `Dockerfile` または `railway.toml` + 起動コマンド `uvicorn app.main:app --host 0.0.0.0 --port $PORT` を指定。`data/gene_regions.pkl` はビルド時に生成するか、Volume で配置。 |
| **[Render](https://render.com)** | New → Web Service → リポジトリ選択。Build: `pip install -r requirements.txt`、Start: `uvicorn app.main:app --host 0.0.0.0 --port $PORT`。`data/` は Build 時に `build_gene_regions.py` で生成するか、永続ディスクに置く。 |
| **[Fly.io](https://fly.io)** | `fly launch` 後、`fly.toml` で `cmd = "uvicorn app.main:app --host 0.0.0.0 --port 8080"`。`data/gene_regions.pkl` は Docker イメージに含めるか、Volume でマウント。 |

共通のポイント:

1. **ポート**  
   本番では環境変数 `PORT` を使うことが多いので、起動コマンドを `uvicorn app.main:app --host 0.0.0.0 --port ${PORT:-8000}` のようにする。
2. **data/gene_regions.pkl**  
   GitHub には載せない前提なので、デプロイ先で用意する必要がある。  
   - ビルド時に `genomic.gff` を置き、`python scripts/build_gene_regions.py` を実行して pkl を生成する（ビルド時間・メモリに注意）、  
   - または別ストレージから取得して `data/` に配置する。

### 無料枠で使うなら（Render / Railway / Fly.io）

このアプリは **起動時に約 100 MB の gene_regions.pkl をメモリに載せる**ため、256 MB だけのインスタンスでは不足しがちです。無料で運用する場合の目安です。

| サービス | 無料枠の要点 | このプロジェクト向きか |
|----------|----------------|-------------------------|
| **[Render](https://render.com)** | 無料 Web サービスあり。**月 750 時間**まで。**15 分アクセスがないとスリープ**し、次のリクエストで復帰（数十秒かかることがある）。メモリは無料枠でも 512 MB 程度が一般的。 | **◎ 無料で使うならここ。** スリープはあるが、デモ・個人利用なら許容しやすい。ビルド時に `gene_regions.pkl` を生成してイメージに含めれば、起動時は pkl を読むだけ。 |
| **[Railway](https://railway.app)** | **常設の無料枠はなし**。$5 のトライアル後は Hobby で **$5/月**（クレジット込み）の従量課金。 | △ 無料で長く使いたい場合は不向き。有料なら常時起動・設定しやすい。 |
| **[Fly.io](https://fly.io)** | 無料で **256 MB の VM が 3 台**程度。請求が **$5 未満の月は未請求**になる運用だが、公式保証ではない。256 MB だと pkl + BED で足りない可能性。512 MB にするとわずかに課金の可能性。 | ○ 常時起動にしたい場合の候補。256 MB では厳しいので、512 MB で「ほぼ無料」になるか要確認。 |

**結論:**  
- **完全無料で試したい → Render**（スリープあり・起動遅延ありで我慢）。  
- **常時起動を優先し、月数ドル出せる → Railway（Hobby）** か **Fly.io（512 MB）**。

---

### Render へのデプロイ（手順）

#### 1. 前提：genomic.gff をダウンロード用 URL で用意する

Render のビルドでは `gene_regions.pkl` を生成するために **genomic.gff** が必要です。GitHub には載せないので、次のいずれかで「直接ダウンロードできる URL」を用意します。

- **Dropbox**: ファイルをアップロード → 共有リンク作成 → リンク末尾の `?dl=0` を **`?dl=1`** に変えた URL（例: `https://www.dropbox.com/s/xxxxx/genomic.gff?dl=1`）
- **Google Drive**: 共有リンクの ID を使って `https://drive.google.com/uc?export=download&id=ファイルID` 形式の URL（大きいファイルは要確認）
- **AWS S3 / 他のクラウド**: 公開読み取り可能な URL

この URL をあとで環境変数 **`BUILD_GFF_URL`** に設定します。

#### 2. Render でアカウント作成とリポジトリ連携

1. [Render](https://render.com) にアクセスし、**Sign Up**（GitHub で登録すると連携が簡単です）。
2. ダッシュボードで **New +** → **Web Service** を選択。
3. **Connect a repository** で、このプロジェクトの GitHub リポジトリを選び **Connect**。
4. リポジトリが表示されたら **Connect** をクリック。

#### 3. Web Service の設定

次のように入力します。

| 項目 | 値 |
|------|-----|
| **Name** | 任意（例: `meth-viz-nanopore`） |
| **Region** | お好みのリージョン（例: Singapore または Oregon） |
| **Branch** | `main`（デフォルトのまま） |
| **Runtime** | **Python 3** |
| **Build Command** | `bash scripts/render_build.sh` |
| **Start Command** | `uvicorn app.main:app --host 0.0.0.0 --port $PORT` |
| **Instance Type** | **Free**（無料枠） |

#### 4. 環境変数を追加

**Environment** セクションで **Add Environment Variable** をクリックし、1 つ追加します。

| Key | Value |
|-----|--------|
| `BUILD_GFF_URL` | 手順 1 で用意した genomic.gff の**直接ダウンロード URL** |

これでビルド時に GFF が取得され、`data/gene_regions.pkl` が生成されます。

#### 5. デプロイ実行

1. 画面下部の **Create Web Service** をクリック。
2. ビルドが始まります。GFF のダウンロードと `build_gene_regions.py` の実行で **おおむね 10〜15 分**かかることがあります。
3. ビルドが成功すると自動で起動し、**Your service is live at …** の URL（例: `https://meth-viz-nanopore.onrender.com`）でアクセスできます。

#### 6. 動作確認

- ブラウザで上記 URL を開く。
- BED ファイルをアップロード → 遺伝子名で検索 → グラフ表示まで試す。
- 無料枠では **15 分間アクセスがないとスリープ**します。再度開くと復帰まで数十秒かかることがあります。

#### トラブルシューティング

- **ビルドが "Error" になる**  
  - **Logs** タブでエラー内容を確認。  
  - `BUILD_GFF_URL` が正しいか、その URL をブラウザで開いたときに GFF がダウンロードされるか確認。
- **「Gene regions not loaded」と表示される**  
  - ビルド時に `gene_regions.pkl` が作られていません。`BUILD_GFF_URL` を設定したうえで **Manual Deploy** → **Clear build cache & deploy** で再デプロイ。
- **起動が遅い**  
  - 無料枠のスリープからの復帰です。そのまま待つか、有料プランでスリープを無効にできます。

---

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

```
meth_viz_nanopore/
├── .gitignore
├── README.md
├── PROJECT_PLAN.md
├── requirements.txt
├── gene_methylation.py          # 遺伝子領域抽出・メチル化サイト取得
├── app/
│   ├── main.py                  # FastAPI（BED アップロード・API）
│   └── static/
│       └── index.html            # フロント（アップロード UI・Plotly グラフ）
├── data/
│   └── .gitkeep                 # genomic.gff / gene_regions.pkl はここに配置（git では無視）
└── scripts/
    ├── build_gene_regions.py    # GFF → gene_regions.pkl 事前処理
    └── render_build.sh          # Render 用ビルド（BUILD_GFF_URL で GFF 取得 → pkl 生成）
```

## API（参考）

| メソッド | パス | 説明 |
|----------|------|------|
| GET | `/health` | ヘルスチェック（デプロイ・LB 用）。常に 200 で `{"status":"ok"}` を返す。 |
| GET | `/api/status` | アプリ状態（gene_regions 読み込み有無・遺伝子数・BED 読み込み有無・行数）。 |
| POST | `/api/upload` | BED のみ multipart でアップロード。 |
| GET | `/api/genes` | 遺伝子名一覧。クエリ `?q=xxx` で部分一致検索（省略可）。 |
| GET | `/api/gene/{gene_id}` | 指定遺伝子のメチル化サイトと領域情報（JSON）。 |

Swagger UI: **http://localhost:8000/docs**

## 注意

- BED の染色体名が `chr1`, `chr2` … で、GFF が RefSeq (`NC_000067.7` 等) の場合は、`gene_methylation.py` 内の `CHROM_GFF_TO_CHR` で対応しています（GRCm39 想定）。
- アップロードしたファイルはサーバーの一時ディレクトリに保存され、再起動や再アップロードで上書きされます。
