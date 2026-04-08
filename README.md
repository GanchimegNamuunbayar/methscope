# methscope

modkit 出力の **BED** をアップロードし、内蔵 GFF 注釈に沿って遺伝子ごとにメチル化を可視化する Web アプリ。

## セットアップ

```bash
cd methscope
python3 -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## 参照データ（初回）

`data/genomic.gff` を置いてから:

```bash
python scripts/build_gene_regions.py
```

## 起動

```bash
source .venv/bin/activate
python -m uvicorn app.main:app --reload --port 8000
```

ブラウザで **http://localhost:8000**

## 動作確認
https://github.com/user-attachments/assets/5a055661-742b-426d-9c7a-dc7ec41f99a6


