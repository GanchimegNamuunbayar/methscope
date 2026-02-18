#!/usr/bin/env python3
import logging
import pandas as pd
import numpy as np
from pathlib import Path
import pyranges as pr
from collections import defaultdict

logger = logging.getLogger(__name__)

############################
PROMOTER_UP = 2000   # TSS 上流 (bp)
DOWNSTREAM_DOWN = 3000  # TES 下流 (bp)
CHUNKSIZE = 500_000

# GRCm39: GFF seqid (RefSeq) <-> BED chrom (chrN)
CHROM_GFF_TO_CHR = {
    "NC_000067.7": "chr1", "NC_000068.8": "chr2", "NC_000069.7": "chr3",
    "NC_000070.7": "chr4", "NC_000071.7": "chr5", "NC_000072.7": "chr6",
    "NC_000073.7": "chr7", "NC_000074.7": "chr8", "NC_000075.7": "chr9",
    "NC_000076.7": "chr10", "NC_000077.7": "chr11", "NC_000078.7": "chr12",
    "NC_000079.7": "chr13", "NC_000080.7": "chr14", "NC_000081.7": "chr15",
    "NC_000082.7": "chr16", "NC_000083.7": "chr17", "NC_000084.7": "chr18",
    "NC_000085.7": "chr19", "NC_000086.8": "chrX", "NC_000087.8": "chrY",
}
CHROM_CHR_TO_GFF = {v: k for k, v in CHROM_GFF_TO_CHR.items()}


def _chrom_for_bed(gff_chrom):
    """Return BED chromosome name (chr1, chr2, ...) for a GFF seqid."""
    return CHROM_GFF_TO_CHR.get(gff_chrom, gff_chrom)

############################
def parse_attr(attr):
    d = {}
    if pd.isna(attr):
        return d
    for x in str(attr).split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            d[k] = v
    return d

############################
def extract_regions(gff_path, promoter_up=None, downstream_down=None):
    """Extract gene regions from GFF. promoter_up = TSS 上流 (bp), downstream_down = TES 下流 (bp)."""
    pu = promoter_up if promoter_up is not None else PROMOTER_UP
    dd = downstream_down if downstream_down is not None else DOWNSTREAM_DOWN

    cols = ["seqid","source","type","start","end","score","strand","phase","attr"]
    gff = pd.read_csv(
        gff_path, sep="\t", comment="#",
        names=cols, low_memory=False
    )
    gff["type"] = gff["type"].astype(str).str.lower()

    genes = gff[gff["type"] == "gene"].copy()

    gene_info = []
    for _, r in genes.iterrows():
        a = parse_attr(r["attr"])
        gid = (
            a.get("ID")
            or a.get("gene_id")
            or a.get("Name")
            or f"{r['seqid']}_{r['start']}_{r['end']}"
        )

        gene_info.append({
            "gene": gid,
            "chrom": r["seqid"],
            "start": int(r["start"]),
            "end": int(r["end"]),
            "strand": r["strand"]
        })

    gene_info = pd.DataFrame(gene_info)
    gene_regions = {}

    for chrom in gene_info.chrom.unique():
        gff_c = gff[gff["seqid"] == chrom]

        for _, g in gene_info[gene_info.chrom == chrom].iterrows():
            regs = {}

            # ----- gene body -----
            regs["gene"] = pd.DataFrame([{
                "seqid": chrom,
                "start": g["start"],
                "end": g["end"]
            }])

            # ----- exon -----
            exon = gff_c[
                (gff_c["type"] == "exon") &
                (gff_c["start"] >= g["start"]) &
                (gff_c["end"] <= g["end"])
            ][["seqid","start","end"]].copy()

            regs["exon"] = exon.reset_index(drop=True)

            # ----- intron = gene - exon -----
            if not exon.empty:
                gene_pr = pr.PyRanges(
                    chromosomes=[chrom],
                    starts=[g["start"]],
                    ends=[g["end"]]
                )
                exon_pr = pr.PyRanges(
                    exon.rename(columns={
                        "seqid":"Chromosome",
                        "start":"Start",
                        "end":"End"
                    })
                )
                intron = gene_pr.subtract(exon_pr).as_df()
                intron = intron.rename(columns={
                    "Chromosome":"seqid",
                    "Start":"start",
                    "End":"end"
                })
            else:
                intron = pd.DataFrame(columns=["seqid","start","end"])

            regs["intron"] = intron.reset_index(drop=True)

            # ----- CDS (coding sequence) -----
            cds = gff_c[
                (gff_c["type"].astype(str).str.lower() == "cds") &
                (gff_c["start"] >= g["start"]) &
                (gff_c["end"] <= g["end"])
            ][["seqid", "start", "end"]].copy()
            regs["cds"] = cds.reset_index(drop=True)

            # ----- promoter: upstream of TSS only (no overlap with gene body) -----
            # ----- downstream: downstream of TES only -----
            if g["strand"] == "+":
                tss = g["start"]
                tes = g["end"]
                regs["promoter"] = pd.DataFrame([{
                    "seqid": chrom,
                    "start": max(0, tss - pu),
                    "end": tss
                }])
                regs["downstream"] = pd.DataFrame([{
                    "seqid": chrom,
                    "start": tes,
                    "end": tes + dd
                }])
            else:
                tss = g["end"]
                tes = g["start"]
                regs["promoter"] = pd.DataFrame([{
                    "seqid": chrom,
                    "start": tss,
                    "end": tss + pu
                }])
                regs["downstream"] = pd.DataFrame([{
                    "seqid": chrom,
                    "start": max(0, tes - dd),
                    "end": tes
                }])

            regs["strand"] = g["strand"]
            regs["exon_count"] = len(regs["exon"])
            regs["cds_count"] = len(regs["cds"])
            gene_regions[g["gene"]] = regs

    return gene_regions


############################
def get_gene_list(gff_path):
    """Return list of gene identifiers for search (from GFF), without building full regions."""
    cols = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attr"]
    gff = pd.read_csv(gff_path, sep="\t", comment="#", names=cols, low_memory=False)
    genes = gff[gff["type"].astype(str).str.lower() == "gene"]
    names = []
    for _, r in genes.iterrows():
        a = parse_attr(r["attr"])
        gid = a.get("ID") or a.get("gene_id") or a.get("Name") or f"{r['seqid']}_{r['start']}_{r['end']}"
        names.append(gid)
    return sorted(set(names))


def load_bed_dataframe(bed_path):
    """Load BED once into a processed DataFrame (type m, numeric). For caching."""
    logger.info("[gene plot] Loading BED file (once per session)...")
    cols = ["chrom", "start", "end", "type", "score", "strand", "ps", "pe", "color", "all", "ratio", "mod", "canonical", "other", "delete", "fail", "diff", "nocall"]
    usecols = ["chrom", "start", "end", "type", "mod", "all"]
    try:
        bed_df = pd.read_csv(bed_path, sep="\t", names=cols, usecols=usecols, comment="#", low_memory=False)
    except Exception:
        bed_df = pd.read_csv(bed_path, sep="\t", names=cols, usecols=usecols, comment="#", low_memory=False, dtype={"chrom": str})
    bed_df = bed_df[bed_df["type"].astype(str) == "m"]
    bed_df["start"] = pd.to_numeric(bed_df["start"], errors="coerce").fillna(0).astype(int)
    bed_df["end"] = pd.to_numeric(bed_df["end"], errors="coerce").fillna(0).astype(int)
    bed_df["mod"] = pd.to_numeric(bed_df["mod"], errors="coerce").fillna(0)
    bed_df["all"] = pd.to_numeric(bed_df["all"], errors="coerce").fillna(0)
    logger.info("[gene plot] BED loaded: %d rows (cached).", len(bed_df))
    return bed_df


def get_gene_methylation_from_cached(gene_regions, bed_df, gene_name):
    """
    Return plot data for one gene using pre-loaded gene_regions and bed_df. No file I/O.
    """
    regs = gene_regions.get(gene_name)
    if regs is None:
        for k in gene_regions:
            if k == gene_name or k.replace("gene-", "") == gene_name or k.endswith("_" + gene_name):
                gene_name = k
                regs = gene_regions[k]
                break
        else:
            gn_lower = gene_name.lower()
            for k in gene_regions:
                if gn_lower in k.lower() or k.replace("gene-", "").lower() == gn_lower:
                    gene_name = k
                    regs = gene_regions[k]
                    break
            else:
                raise KeyError(f"Gene not found: {gene_name}")

    chrom = regs["gene"].iloc[0]["seqid"]
    strand = regs.get("strand", ".")

    span_start = int(min(
        regs["promoter"].iloc[0]["start"],
        regs["gene"].iloc[0]["start"],
        regs["downstream"].iloc[0]["start"],
    ))
    span_end = int(max(
        regs["promoter"].iloc[0]["end"],
        regs["gene"].iloc[0]["end"],
        regs["downstream"].iloc[0]["end"],
    ))

    bed_chrom = _chrom_for_bed(chrom)
    mask_chrom = (bed_df["chrom"] == chrom) | (bed_df["chrom"] == bed_chrom)
    mask_span = (bed_df["end"] > span_start) & (bed_df["start"] < span_end)
    sub = bed_df.loc[mask_chrom & mask_span].copy()
    logger.info("[gene plot] Filtered to %d sites for %s (from cache).", len(sub), gene_name)

    def ratio_or_compute(row):
        a = row["all"]
        if a and a > 0:
            return 100.0 * float(row["mod"]) / float(a)
        return np.nan

    sub["methylation_ratio"] = sub.apply(ratio_or_compute, axis=1)
    sites = []
    for _, row in sub.iterrows():
        pos = int(row["start"])
        sites.append({
            "position": pos,
            "methylation_ratio": float(row["methylation_ratio"]) if not np.isnan(row["methylation_ratio"]) else None,
            "coverage": int(row["all"]),
        })
    sites.sort(key=lambda x: x["position"])

    regions = []
    for rtype in ["promoter", "exon", "intron", "cds", "downstream"]:
        df = regs.get(rtype)
        if df is None or df.empty:
            continue
        for _, row in df.iterrows():
            regions.append({"region_type": rtype, "start": int(row["start"]), "end": int(row["end"])})
    regions.sort(key=lambda r: (r["start"], r["end"]))

    exon_count = regs.get("exon_count", 0)
    cds_count = regs.get("cds_count", 0)

    return {
        "sites": sites,
        "regions": regions,
        "gene": str(gene_name),
        "chrom": str(chrom),
        "strand": str(strand),
        "span_start": span_start,
        "span_end": span_end,
        "exon_count": int(exon_count) if exon_count is not None else 0,
        "cds_count": int(cds_count) if cds_count is not None else 0,
    }


def get_gene_methylation_for_plot(
    gff_path,
    bed_path,
    gene_name,
    promoter_up=2000,
    downstream_down=2000,
):
    """
    Return per-site methylation and region boundaries for one gene (for interactive plot).
    Returns: dict with sites, regions, gene, chrom, strand.
    """
    logger.info("[gene plot] Loading GFF and extracting all gene regions (this can take a while for large GFF)...")
    gene_regions = extract_regions(
        gff_path,
        promoter_up=promoter_up,
        downstream_down=downstream_down,
    )
    logger.info("[gene plot] Extracted %d genes. Looking up gene: %s", len(gene_regions), gene_name)

    # Match gene: exact or strip "gene-" prefix
    regs = gene_regions.get(gene_name)
    if regs is None:
        for k in gene_regions:
            if k == gene_name or k.replace("gene-", "") == gene_name or k.endswith("_" + gene_name):
                gene_name = k
                regs = gene_regions[k]
                break
        else:
            # Try case-insensitive and Name-like match
            gn_lower = gene_name.lower()
            for k in gene_regions:
                if gn_lower in k.lower() or k.replace("gene-", "").lower() == gn_lower:
                    gene_name = k
                    regs = gene_regions[k]
                    break
            else:
                raise KeyError(f"Gene not found: {gene_name}")

    chrom = regs["gene"].iloc[0]["seqid"]
    strand = regs.get("strand", ".")

    span_start = min(
        regs["promoter"].iloc[0]["start"],
        regs["gene"].iloc[0]["start"],
        regs["downstream"].iloc[0]["start"],
    )
    # Use same logic as cached path: load BED once per call, then use shared helper
    bed_df = load_bed_dataframe(bed_path)
    return get_gene_methylation_from_cached(gene_regions, bed_df, gene_name)


############################
def build_master(gene_regions):
    rows = []
    for gene, regs in gene_regions.items():
        for rtype, df in regs.items():
            if df.empty:
                continue
            for i, r in df.iterrows():
                rows.append({
                    "Chromosome": r["seqid"],
                    "Start": int(r["start"]),
                    "End": int(r["end"]),
                    "gene": gene,
                    "region": rtype,
                    "region_id": i + 1
                })
    return pr.PyRanges(pd.DataFrame(rows))

############################
def process_bed(bed, pr_regions, sample):
    cols = [
        "chrom","start","end","type","score","strand",
        "ps","pe","color","all","ratio",
        "mod","canonical","other","delete",
        "fail","diff","nocall"
    ]

    acc = defaultdict(lambda: [0, 0])

    reader = pd.read_csv(
        bed, sep="\t", names=cols,
        usecols=["chrom","start","end","mod","all"],
        chunksize=CHUNKSIZE,
        comment="#"
    )

    for chunk in reader:
        chunk["Start"] = chunk["start"].astype(int)
        chunk["End"] = chunk["end"].astype(int)
        chunk["mod"] = pd.to_numeric(chunk["mod"], errors="coerce").fillna(0)
        chunk["all"] = pd.to_numeric(chunk["all"], errors="coerce").fillna(0)

        chunk_pr = pr.PyRanges(
            chunk.rename(columns={
                "chrom":"Chromosome",
                "Start":"Start",
                "End":"End"
            })[["Chromosome","Start","End","mod","all"]]
        )

        joined = chunk_pr.join(pr_regions).as_df()
        if joined.empty:
            continue

        # 🔥 超重要：category → object
        for c in ["gene", "region", "region_id"]:
            joined[c] = joined[c].astype(str)

        g = joined.groupby(
            ["gene", "region", "region_id"],
            observed=True
        ).agg(
            mod_sum=("mod", "sum"),
            cov_sum=("all", "sum")
        ).reset_index()

        for _, r in g.iterrows():
            key = (r["gene"], r["region"], int(r["region_id"]))
            acc[key][0] += int(r["mod_sum"])
            acc[key][1] += int(r["cov_sum"])

    out = []
    for (gene, region, rid), (m, c) in acc.items():
        out.append({
            "gene": gene,
            "region": region,
            "region_id": rid,
            "condition": sample,
            "methylation": (m / c) if c > 0 else np.nan,
            "meth_reads": m,
            "coverage": c
        })

    return pd.DataFrame(out)

############################
def main():
    gff = "../../data/reference/annotation/genomic.gff"

    beds = {
        "r0081": "../../data/1.0.0/01_preprocessed/modkit/CpG_m/r0081_m.bed",
        "r0082": "../../data/1.0.0/01_preprocessed/modkit/CpG_m/r0082_m.bed",
        "r0095": "../../data/1.0.0/01_preprocessed/modkit/CpG_m/r0095_m.bed",
        "r0096": "../../data/1.0.0/01_preprocessed/modkit/CpG_m/r0096_m.bed",
    }

    out = Path(
        "../../data/1.0.0/04_ora/go_meth/"
        "all_gene_region_methylation_coverageAware.csv"
    )
    out.parent.mkdir(parents=True, exist_ok=True)

    print("[INFO] extracting gene regions")
    gene_regions = extract_regions(gff)
    pr_regions = build_master(gene_regions)
    print(f"[INFO] total annotated regions: {len(pr_regions)}")

    res = []
    for s, bed in beds.items():
        print(f"[INFO] processing {s}")
        df = process_bed(bed, pr_regions, s)
        res.append(df)

    pd.concat(res, ignore_index=True).to_csv(out, index=False)
    print("[DONE]", out)

if __name__ == "__main__":
    main()
