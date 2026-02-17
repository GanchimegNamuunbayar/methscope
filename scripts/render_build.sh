#!/usr/bin/env bash
# Render の Build Command で使うスクリプト。
# 環境変数 BUILD_GFF_URL が設定されていれば、GFF を取得して gene_regions.pkl を生成する。
set -e

pip install -r requirements.txt

mkdir -p data
if [ -n "${BUILD_GFF_URL:-}" ]; then
  echo "Downloading GFF from BUILD_GFF_URL..."
  curl -fSL -o data/genomic.gff "$BUILD_GFF_URL"
  echo "Building gene_regions.pkl (this may take several minutes)..."
  python scripts/build_gene_regions.py data/genomic.gff
  echo "Done. gene_regions.pkl is ready."
else
  echo "BUILD_GFF_URL not set. Skip building gene_regions.pkl."
  echo "App will start but /api/genes and /api/upload will return 503 until you add data manually."
fi
