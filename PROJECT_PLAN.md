# Nanopore メチル化可視化 Web アプリ 開発計画

## 概要
modkit 出力 BED と GFF を入力し、遺伝子名で検索してプロモーター（上流 2kbp・下流 2kbp）を含むメチル化のインタラクティブグラフを表示する Web アプリを開発する。

---

## Step 1: バックエンドデータモジュール

**目的**: 1 遺伝子分の「位置別メチル化」と「領域境界」を返す処理を追加する。

- **1.1** `gene_methylation.py` の `extract_regions` を拡張  
  - プロモーターを「上流 2kbp・下流 2kbp」に変更可能にする（`PROMOTER_DOWN=2000` のオプション）。
- **1.2** 染色体名の対応  
  - GFF: `NC_000067.7`、BED: `chr1` など表記差があるため、マッピング辞書を用意し、BED 読み込み時に GFF の chrom に合わせて検索できるようにする。
- **1.3** 新関数 `get_gene_methylation_for_plot(gff_path, bed_path, gene_name, promoter_up=2000, promoter_down=2000)`  
  - 指定遺伝子の gene_regions を取得。  
  - その遺伝子の染色体・範囲（プロモーター開始 ～ 下流終了）で BED をフィルタ（type は `m` = m5C を想定）。  
  - 返り値:  
    - `sites`: `[{ position, methylation_ratio, coverage }, ...]`（プロットの点）  
    - `regions`: `[{ region_type, start, end }, ...]`（promoter / exon / intron / downstream の範囲。プロットの縦帯に使用）

**成果物**: `gene_methylation.py` の拡張（および必要なら `app/plot_data.py` のような薄いラッパー）。

---

## Step 2: FastAPI バックエンド

**目的**: ファイルアップロードと遺伝子メチル化データ取得の API を提供する。

- **2.1** プロジェクト構成  
  - `app/` に FastAPI アプリを配置。  
  - 依存: `fastapi`, `uvicorn`, `python-multipart`, 既存の `pandas`, `pyranges`, `numpy`。
- **2.2** エンドポイント  
  - `POST /api/upload`  
    - BED ファイルと GFF ファイルを multipart で受信。  
    - 一時ディレクトリに保存し、セッション（またはメモリ）で「現在の BED パス」「現在の GFF パス」を保持。  
  - `GET /api/genes`  
    - 現在の GFF から遺伝子名一覧を返す（検索・オートコンプリート用）。  
  - `GET /api/gene/{gene_id}` または `POST /api/gene` (body: `{ "gene_name": "Xkr4" }`)  
    - 現在の BED + GFF を使い、Step 1 の `get_gene_methylation_for_plot` を呼ぶ。  
    - 返り値: `sites`, `regions`, および遺伝子の chrom/strand などのメタ情報。
- **2.3** エラーハンドリング  
  - 未アップロード時に 400 でメッセージを返す。  
  - 遺伝子が見つからない場合は 404。

**成果物**: `app/main.py`, `requirements.txt` の更新。

---

## Step 3: フロントエンド

**目的**: アップロード UI・遺伝子検索・インタラクティブなメチル化グラフを提供する。

- **3.1** 単一 HTML ページ（+ JS）  
  - FastAPI の静的ファイルまたはテンプレートで配信。  
  - 構成: ファイルアップロード（BED, GFF）→ 「解析」または GFF 解析後に遺伝子一覧表示 → 遺伝子名検索（入力 + 候補表示）→ 「表示」でグラフ取得。
- **3.2** 遺伝子検索  
  - `/api/genes` で取得した一覧をクライアントに保持。  
  - 入力に応じてフィルタして候補を表示し、選択で遺伝子名を確定。
- **3.3** インタラクティブグラフ（Plotly.js）  
  - X 軸: ゲノム位置。  
  - Y 軸: メチル化率（0–100%）。  
  - 点: 各 CpG（または modkit の 1 行）の position と methylation_ratio。  
  - 領域の表示: `regions` に基づき、プロモーター・エクソン・イントロン・下流を異なる色の縦帯（`layout.shapes` の矩形）で表示。  
  - 凡例: どの色がどの領域かを表示（ノートブックと同様の意図）。  
  - ツールチップ: 位置・メチル化率・カバレッジを表示。

**成果物**: `app/static/index.html`（および必要なら `app/static/app.js`）。

---

## Step 4: 統合・動作確認

- **4.1** 起動方法  
  - `uvicorn app.main:app --reload` で API を起動。  
  - ブラウザでフロントを開き、`r0081_m_chr.bed` と `genomic.gff` をアップロードし、遺伝子名（例: Xkr4）で検索してグラフが表示されることを確認。
- **4.2** README 更新  
  - セットアップ手順（`pip install -r requirements.txt`）、起動コマンド、使い方を記載。

---

## ファイル構成（予定）

```
meth_viz_nanopore/
├── PROJECT_PLAN.md          # 本計画書
├── gene_methylation.py       # 拡張（extract_regions オプション + get_gene_methylation_for_plot）
├── requirements.txt         # 依存関係
├── app/
│   ├── main.py               # FastAPI アプリ、エンドポイント
│   └── static/
│       └── index.html        # フロントエンド（アップロード・検索・Plotly グラフ）
├── r0081_m_chr.bed
└── genomic.gff
```

---

## 実装順序

1. ✅ Step 1 のデータモジュール実装（`gene_methylation.py` 拡張）
2. ✅ Step 2 の FastAPI 実装（`app/main.py`）
3. ✅ Step 3 のフロントエンド実装（`app/static/index.html`）
4. ✅ Step 4 の統合・README 更新

以上に従い、各セクションを実装済み。
