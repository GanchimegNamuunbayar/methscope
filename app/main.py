"""
FastAPI backend for nanopore methylation viz: BED upload only; GFF is bundled (pre-processed).

Supports two modes:
- Lazy-load (low memory): gene_list.json + gene_regions.db — list at startup, load one gene from DB on demand.
- Legacy: gene_regions.pkl — full dict loaded at startup (local dev).

Upload uses async jobs: save file and return job_id immediately (avoids Render 502 timeout);
BED is loaded in a background task; client polls GET /api/jobs/{job_id} until ready.
"""
from __future__ import annotations

import asyncio
import json
import logging
import pickle
import sqlite3
import shutil
import time
import uuid
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s", datefmt="%H:%M:%S")

from fastapi import FastAPI, File, UploadFile, HTTPException, Query
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
GENE_LIST_PATH = ROOT / "data" / "gene_list.json"
GENE_REGIONS_DB_PATH = ROOT / "data" / "gene_regions.db"
UPLOAD_DIR = ROOT / "data" / "uploads"

# Where to look for genomic.gff if .pkl is missing (in order)
GFF_CANDIDATES = [
    ROOT / "data" / "genomic.gff",
    ROOT / "genomic.gff",
    Path(__file__).resolve().parent / "genomic.gff",
    Path(__file__).resolve().parent / "static" / "genomic.gff",
]


def _load_gene_list_and_db():
    """Load gene list only; DB path for on-demand reads. Returns (gene_list, db_path) or (None, None)."""
    if not GENE_LIST_PATH.exists() or not GENE_REGIONS_DB_PATH.exists():
        return None, None
    try:
        with open(GENE_LIST_PATH, encoding="utf-8") as f:
            gene_list = json.load(f)
        return gene_list, str(GENE_REGIONS_DB_PATH)
    except Exception as e:
        logging.warning("[API] Failed to load gene list / DB: %s", e)
        return None, None


def _resolve_gene_name(gene_list: list[str], query: str) -> str:
    """Resolve user query to canonical gene_id (same logic as get_gene_methylation_from_cached)."""
    q = query.strip()
    if not q:
        raise KeyError("gene_id is required")
    # Exact
    if q in gene_list:
        return q
    for k in gene_list:
        if k == q or k.replace("gene-", "") == q or k.endswith("_" + q):
            return k
    # Case-insensitive substring
    lower = q.lower()
    for k in gene_list:
        if lower in k.lower() or k.replace("gene-", "").lower() == lower:
            return k
    raise KeyError(f"Gene not found: {query}")


def _load_one_gene_regions(db_path: str, gene_id: str):
    """Load one gene's regs dict from SQLite. Raises KeyError if not found."""
    with sqlite3.connect(db_path) as conn:
        row = conn.execute("SELECT data FROM genes WHERE gene_id = ?", (gene_id,)).fetchone()
    if not row:
        raise KeyError(f"Gene not found: {gene_id}")
    return pickle.loads(row[0])


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
                    gff_path, promoter_up=2000, downstream_down=2000
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


def _gene_regions_ready():
    """True if gene data is available (lazy or legacy)."""
    if _state.get("gene_list") is not None and _state.get("gene_regions_db_path"):
        return True
    return _state.get("gene_regions") is not None


def _genes_count():
    if _state.get("gene_list") is not None:
        return len(_state["gene_list"])
    gr = _state.get("gene_regions")
    return len(gr) if gr else 0


# In-memory store: job-based BED (for Render: upload returns fast, process in background)
# _jobs[job_id] = { "status": "processing"|"ready", "bed_df": None|DataFrame, "bed_path": str, "created_at": float }
_jobs: dict[str, dict] = {}
_state = {
    "current_job_id": None,  # latest ready job_id; /api/gene uses this job's bed_df
    "gene_regions": None,
    "gene_list": None,
    "gene_regions_db_path": None,
}


def _get_current_bed_df():
    """Return bed_df for the current job if ready, else None."""
    jid = _state.get("current_job_id")
    if not jid:
        return None
    job = _jobs.get(jid)
    if not job or job.get("status") != "ready" or job.get("bed_df") is None:
        return None
    return job["bed_df"]


def _evict_other_jobs(keep_job_id: str):
    """Keep only one job with bed_df in memory (for 512MB limit). Drop bed_df from others."""
    for jid, job in list(_jobs.items()):
        if jid != keep_job_id and job.get("bed_df") is not None:
            job["bed_df"] = None
            logging.info("[API] Evicted job %s from memory (keep %s)", jid, keep_job_id)


async def _process_job(job_id: str, bed_path: Path):
    """Load BED in a thread (CPU-bound); then set job ready and evict previous job."""
    try:
        loop = asyncio.get_event_loop()
        bed_df = await loop.run_in_executor(None, load_bed_dataframe, bed_path)
        _jobs[job_id]["bed_df"] = bed_df
        _jobs[job_id]["status"] = "ready"
        _jobs[job_id]["bed_rows"] = len(bed_df)
        _state["current_job_id"] = job_id
        _evict_other_jobs(job_id)
        logging.info("[API] Job %s ready: %d rows", job_id, len(bed_df))
    except Exception as e:
        logging.exception("[API] Job %s failed: %s", job_id, e)
        _jobs[job_id]["status"] = "failed"
        _jobs[job_id]["error"] = str(e)


@app.on_event("startup")
def startup_load_gene_regions():
    UPLOAD_DIR.mkdir(parents=True, exist_ok=True)
    # Prefer lazy-load format (gene_list + DB) for low memory (e.g. Render free tier)
    gene_list, db_path = _load_gene_list_and_db()
    if gene_list is not None and db_path:
        _state["gene_list"] = gene_list
        _state["gene_regions_db_path"] = db_path
        logging.info("[API] Lazy-load mode: gene list %d genes, DB at %s", len(gene_list), db_path)
        return
    # Fallback: full in-memory (legacy pkl)
    _state["gene_regions"] = _load_bundled_gene_regions()
    if _state["gene_regions"] is not None:
        logging.info("[API] Loaded gene regions (legacy): %d genes", len(_state["gene_regions"]))
    else:
        logging.warning(
            "[API] No gene regions. Place genomic.gff in project root or data/ and restart; "
            "or run: python scripts/build_gene_regions.py"
        )


@app.get("/health")
async def health():
    """Health check for load balancers and deployment. Returns 200 when the app is up."""
    return {"status": "ok"}


@app.get("/api/status")
async def api_status():
    """Return app state: whether gene regions and BED are loaded. Useful for debugging and UI."""
    bed_df = _get_current_bed_df()
    return {
        "gene_regions_loaded": _gene_regions_ready(),
        "genes_count": _genes_count(),
        "bed_loaded": bed_df is not None,
        "bed_rows": len(bed_df) if bed_df is not None else 0,
        "current_job_id": _state.get("current_job_id"),
    }


@app.post("/api/upload")
async def upload_files(
    bed: UploadFile = File(..., description="modkit BED (e.g. *_m_chr.bed)"),
):
    """Upload BED: save file and return job_id immediately. BED is loaded in background (avoids timeout on Render)."""
    if not _gene_regions_ready():
        raise HTTPException(
            status_code=503,
            detail="Gene regions not loaded. Place genomic.gff in project root or data/ and restart the app.",
        )
    job_id = str(uuid.uuid4())
    bed_path = UPLOAD_DIR / f"{job_id}.bed"
    try:
        content = await bed.read()
        bed_path.write_bytes(content)
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Failed to save BED: {e}")

    _jobs[job_id] = {
        "status": "processing",
        "bed_df": None,
        "bed_path": str(bed_path),
        "created_at": time.time(),
    }
    asyncio.create_task(_process_job(job_id, bed_path))
    logging.info("[API] Upload job %s created; processing in background", job_id)
    return {"job_id": job_id, "status": "processing"}


@app.get("/api/jobs/{job_id}")
async def get_job_status(job_id: str):
    """Poll job status after upload. Returns processing | ready | failed."""
    if job_id not in _jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    job = _jobs[job_id]
    out = {"job_id": job_id, "status": job["status"]}
    if job["status"] == "ready":
        out["bed_rows"] = job.get("bed_rows", 0)
    if job["status"] == "failed":
        out["error"] = job.get("error", "Unknown error")
    return out


@app.get("/api/genes")
async def list_genes(q: str | None = Query(None, description="Filter genes by substring (case-insensitive)")):
    """Return list of gene identifiers from bundled gene regions. Optional query param 'q' filters by substring."""
    if not _gene_regions_ready():
        raise HTTPException(
            status_code=503,
            detail="Gene regions not loaded. Place genomic.gff in project root or data/ and restart the app.",
        )
    if _state.get("gene_list") is not None:
        genes = _state["gene_list"]  # already sorted from build
    else:
        genes = sorted(_state["gene_regions"].keys())
    if q and q.strip():
        lower = q.strip().lower()
        genes = [g for g in genes if lower in g.lower()]
    return {"genes": genes}


@app.get("/api/gene/{gene_id:path}")
async def get_gene_methylation(gene_id: str):
    """Return methylation sites and region boundaries for the given gene (for interactive plot)."""
    logging.info("[API] Gene plot request: %s", gene_id)
    bed_df = _get_current_bed_df()
    if not _gene_regions_ready():
        raise HTTPException(status_code=503, detail="Gene regions not loaded. Place genomic.gff in project root or data/ and restart.")
    if bed_df is None:
        raise HTTPException(
            status_code=400,
            detail="No BED ready. Upload a BED file and wait until job status is 'ready' (poll /api/jobs/{job_id}).",
        )

    if not gene_id or not gene_id.strip():
        raise HTTPException(status_code=400, detail="gene_id is required.")

    try:
        if _state.get("gene_list") is not None and _state.get("gene_regions_db_path"):
            canonical = _resolve_gene_name(_state["gene_list"], gene_id.strip())
            regs = _load_one_gene_regions(_state["gene_regions_db_path"], canonical)
            gene_regions_one = {canonical: regs}
            data = get_gene_methylation_from_cached(gene_regions_one, bed_df, canonical)
        else:
            data = get_gene_methylation_from_cached(
                _state["gene_regions"], bed_df, gene_id.strip()
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
