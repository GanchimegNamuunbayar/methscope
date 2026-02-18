#!/usr/bin/env python3
"""
Pre-process genomic.gff into gene_regions cache for fast app startup.

Outputs (for lazy-load / low-memory deploy, e.g. Render free tier):
  - data/gene_list.json   : list of gene IDs (loaded at startup only)
  - data/gene_regions.db  : SQLite, one row per gene (loaded on demand per gene)

Also writes data/gene_regions.pkl for local dev (full in-memory mode).
Run once after placing genomic.gff in data/ (or pass path).

Definitions: promoter = 2k upstream of TSS only (no overlap with gene);
downstream = 2k past TES. Exon, intron, and CDS are extracted from GFF.
"""
import sys
import json
import pickle
import sqlite3
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from gene_methylation import extract_regions


def main():
    data_dir = ROOT / "data"
    data_dir.mkdir(exist_ok=True)

    gff_path = data_dir / "genomic.gff"
    if len(sys.argv) > 1:
        gff_path = Path(sys.argv[1]).resolve()
    elif not gff_path.exists() and (ROOT / "genomic.gff").exists():
        gff_path = ROOT / "genomic.gff"
    if not gff_path.exists():
        print(f"Error: GFF not found: {gff_path}")
        print("  Place genomic.gff in data/ or run: python scripts/build_gene_regions.py /path/to/genomic.gff")
        sys.exit(1)

    print(f"Reading GFF: {gff_path}")
    print("Extracting gene regions (promoter 2k upstream, downstream 2k, exon/intron/CDS)...")
    gene_regions = extract_regions(
        gff_path,
        promoter_up=2000,      # TSS 上流 2 kbp
        downstream_down=2000, # TES 下流 2 kbp
    )
    n = len(gene_regions)
    print(f"  Genes: {n}")

    # 1) Gene list only (small; loaded at app startup)
    gene_list_path = data_dir / "gene_list.json"
    gene_list = sorted(gene_regions.keys())
    with open(gene_list_path, "w", encoding="utf-8") as f:
        json.dump(gene_list, f)
    print(f"Saved: {gene_list_path} ({len(gene_list)} ids)")

    # 2) SQLite: one row per gene, regs as pickle blob (loaded on demand per gene)
    db_path = data_dir / "gene_regions.db"
    with sqlite3.connect(db_path) as conn:
        conn.execute(
            "CREATE TABLE IF NOT EXISTS genes (gene_id TEXT PRIMARY KEY, data BLOB)"
        )
        conn.execute("DELETE FROM genes")
        for gid in gene_list:
            blob = pickle.dumps(gene_regions[gid], protocol=pickle.HIGHEST_PROTOCOL)
            conn.execute("INSERT INTO genes (gene_id, data) VALUES (?, ?)", (gid, blob))
    print(f"Saved: {db_path}")

    # 3) Optional: single pkl for local dev (full in-memory mode)
    out_path = data_dir / "gene_regions.pkl"
    with open(out_path, "wb") as f:
        pickle.dump(gene_regions, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
