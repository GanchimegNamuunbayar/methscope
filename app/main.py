"""
FastAPI backend for nanopore methylation viz: BED upload only; GFF is bundled (pre-processed).
"""
import logging
import pickle
import tempfile
import shutil
import time
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

from fastapi import FastAPI, File, UploadFile, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse

# Import from project root (gene_methylation.py)
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from gene_methylation import (
    get_gene_methylation_from_cached,
    load_bed_dataframe,
    extract_regions,
)

app = FastAPI(title="Nanopore Gene Methylation Viz")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

ROOT = Path(__file__).resolve().parent.parent
GENE_REGIONS_PATH = ROOT / "data" / "gene_regions.pkl"

# Where to look for genomic.gff if .pkl is missing (in order)
GFF_CANDIDATES = [
    ROOT / "data" / "genomic.gff",
    ROOT / "genomic.gff",
    Path(__file__).resolve().parent / "genomic.gff",
    Path(__file__).resolve().parent / "static" / "genomic.gff",
]


def _load_bundled_gene_regions():
    # 1) Pre-processed cache (fast)
    if GENE_REGIONS_PATH.exists():
        try:
            with open(GENE_REGIONS_PATH, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            logging.warning("[API] Failed to load %s: %s", GENE_REGIONS_PATH, e)

    # 2) Fallback: find genomic.gff and extract regions (slow, first time only)
    for gff_path in GFF_CANDIDATES:
        if gff_path.exists():
            logging.info("[API] No .pkl cache; loading GFF from %s (this may take several minutes)...", gff_path)
            t0 = time.perf_counter()
            try:
                gene_regions = extract_regions(
                    gff_path, promoter_up=2000, promoter_down=2000, downstream_down=2000
                )
                elapsed = round(time.perf_counter() - t0, 1)
                logging.info("[API] GFF done: %d genes in %.1fs", len(gene_regions), elapsed)
                # Save for next startup
                GENE_REGIONS_PATH.parent.mkdir(parents=True, exist_ok=True)
                with open(GENE_REGIONS_PATH, "wb") as f:
                    pickle.dump(gene_regions, f, protocol=pickle.HIGHEST_PROTOCOL)
                logging.info("[API] Saved cache to %s", GENE_REGIONS_PATH)
                return gene_regions
            except Exception as e:
                logging.exception("[API] Failed to load GFF from %s: %s", gff_path, e)
                return None

    return None

# In-memory store: BED path + cached BED data (gene_regions are bundled, loaded at startup)
_state = {
    "bed_path": None,
    "upload_dir": None,
    "gene_regions": None,  # loaded at startup from data/gene_regions.pkl
    "bed_df": None,
}


def _cleanup_upload_dir():
    if _state.get("upload_dir") and Path(_state["upload_dir"]).exists():
        try:
            shutil.rmtree(_state["upload_dir"], ignore_errors=True)
        except Exception:
            pass
        _state["upload_dir"] = None
    _state["bed_path"] = None
    _state["bed_df"] = None


@app.on_event("startup")
def startup_load_gene_regions():
    _state["gene_regions"] = _load_bundled_gene_regions()
    if _state["gene_regions"] is not None:
        logging.info("[API] Loaded gene regions: %d genes", len(_state["gene_regions"]))
    else:
        logging.warning(
            "[API] No gene regions. Place genomic.gff in project root or data/ and restart; "
            "or run: python scripts/build_gene_regions.py"
        )


@app.post("/api/upload")
async def upload_files(
    bed: UploadFile = File(..., description="modkit BED (e.g. *_m_chr.bed)"),
):
    """Upload BED file only. Gene annotation is from bundled pre-processed GFF."""
    if _state.get("gene_regions") is None:
        raise HTTPException(
            status_code=503,
            detail="Gene regions not loaded. Place genomic.gff in project root or data/ and restart the app.",
        )
    _cleanup_upload_dir()
    upload_dir = Path(tempfile.mkdtemp(prefix="meth_viz_"))
    _state["upload_dir"] = str(upload_dir)

    bed_path = upload_dir / (bed.filename or "uploaded.bed")
    try:
        with open(bed_path, "wb") as f:
            content = await bed.read()
            f.write(content)
    except Exception as e:
        _cleanup_upload_dir()
        raise HTTPException(status_code=400, detail=f"Failed to save BED: {e}")

    _state["bed_path"] = str(bed_path)

    load_stats = {"genes_count": len(_state["gene_regions"]), "bed_rows": 0, "time_bed_sec": 0}
    t0 = time.perf_counter()
    try:
        logging.info("[API] Loading BED...")
        _state["bed_df"] = load_bed_dataframe(bed_path)
        load_stats["bed_rows"] = len(_state["bed_df"])
        load_stats["time_bed_sec"] = round(time.perf_counter() - t0, 2)
        logging.info("[API] BED done: %d rows in %.1fs", load_stats["bed_rows"], load_stats["time_bed_sec"])
    except Exception as e:
        logging.exception("[API] BED load failed: %s", e)
        _state["bed_df"] = None
        raise HTTPException(status_code=500, detail=f"Failed to load BED: {e}")

    return {
        "status": "ok",
        "message": "BED uploaded and loaded. You can search genes and view plots.",
        "load_stats": load_stats,
    }


@app.get("/api/genes")
async def list_genes():
    """Return list of gene identifiers from bundled gene regions."""
    if _state.get("gene_regions") is None:
        raise HTTPException(
            status_code=503,
            detail="Gene regions not loaded. Place genomic.gff in project root or data/ and restart the app.",
        )
    return {"genes": sorted(_state["gene_regions"].keys())}


@app.get("/api/gene/{gene_id:path}")
async def get_gene_methylation(gene_id: str):
    """Return methylation sites and region boundaries for the given gene (for interactive plot)."""
    logging.info("[API] Gene plot request: %s", gene_id)
    bed_path = _state.get("bed_path")
    if not bed_path or not Path(bed_path).exists():
        raise HTTPException(status_code=400, detail="No BED uploaded. Please upload BED first.")
    if _state.get("gene_regions") is None:
        raise HTTPException(status_code=503, detail="Gene regions not loaded. Place genomic.gff in project root or data/ and restart.")
    if _state.get("bed_df") is None:
        raise HTTPException(status_code=400, detail="BED not loaded. Please upload BED first.")

    if not gene_id or not gene_id.strip():
        raise HTTPException(status_code=400, detail="gene_id is required.")

    try:

        data = get_gene_methylation_from_cached(
            _state["gene_regions"], _state["bed_df"], gene_id.strip()
        )
        logging.info("[API] Gene plot done: %s (%d sites)", data.get("gene"), len(data.get("sites", [])))
        return data
    except KeyError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


@app.get("/")
async def index():
    """Serve the frontend."""
    static_dir = Path(__file__).parent / "static"
    index_file = static_dir / "index.html"
    if index_file.exists():
        return FileResponse(index_file)
    return {"message": "Static files not found. Put index.html in app/static/."}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
