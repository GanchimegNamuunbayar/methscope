#!/usr/bin/env python3
"""
Pre-process genomic.gff into gene_regions cache (pickle) for fast app startup.
Run once after placing genomic.gff in data/ (or pass path).

Definitions: promoter = 2k upstream of TSS only (no overlap with gene);
downstream = 2k past TES. Exon, intron, and CDS are extracted from GFF.
"""
import sys
import pickle
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

    out_path = data_dir / "gene_regions.pkl"
    print(f"Reading GFF: {gff_path}")
    print("Extracting gene regions (promoter 2k upstream, downstream 2k, exon/intron/CDS)...")
    gene_regions = extract_regions(
        gff_path,
        promoter_up=2000,
        promoter_down=2000,
        downstream_down=2000,
    )
    print(f"  Genes: {len(gene_regions)}")
    with open(out_path, "wb") as f:
        pickle.dump(gene_regions, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    main()
