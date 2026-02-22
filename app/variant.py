"""
Variant search API: VCF list from config, POST search runs auto_variant logic.
No AI; deterministic search by gene or location.
"""
from __future__ import annotations

import os
from pathlib import Path

from fastapi import APIRouter, HTTPException
from pydantic import BaseModel, Field

router = APIRouter()

# Project root: resolve relative VCF_DIR from here so it works regardless of process cwd
_APP_ROOT = Path(__file__).resolve().parent.parent

# When source is "paths", list is from env (cached). When "dir", we re-scan the folder on each request.
_VCF_PATHS: list[tuple[str, str]] = []  # (label, absolute path)
_VCF_SOURCE: str = ""  # "paths" | "dir"
_VCF_DIR_PATH: Path | None = None  # resolved dir when source is "dir"


def _get_vcf_paths_env() -> str:
    return os.environ.get("VCF_PATHS", "").strip()


def _get_vcf_dir_env() -> str:
    return os.environ.get("VCF_DIR", "").strip()


def _is_vcf_file(name: str) -> bool:
    """Only .vcf and .vcf.gz; ignore other files in the folder."""
    return name.endswith(".vcf.gz") or (name.endswith(".vcf") and not name.endswith(".vcf.gz"))


def _discover_vcf_in_dir(dir_path: Path) -> list[tuple[str, str]]:
    """Scan directory for .vcf and .vcf.gz only; return list of (label, absolute path)."""
    out = []
    if not dir_path.is_dir():
        return out
    for f in sorted(dir_path.iterdir()):
        if f.is_file() and _is_vcf_file(f.name):
            out.append((f.name, str(f.resolve())))
    return out


# Default folders to try when VCF_DIR and VCF_PATHS are not set (first that has VCFs wins)
_DEFAULT_VCF_DIRS = ("data", "/data")


def init_app() -> None:
    """Determine VCF source (VCF_PATHS or VCF_DIR / defaults). Idempotent."""
    global _VCF_PATHS, _VCF_SOURCE, _VCF_DIR_PATH
    raw_paths = _get_vcf_paths_env()
    if raw_paths:
        out = []
        for part in raw_paths.split(","):
            path_str = part.strip()
            if not path_str:
                continue
            p = Path(path_str)
            if not p.is_absolute():
                p = (_APP_ROOT / path_str).resolve()
            else:
                p = p.resolve()
            label = p.name or path_str
            out.append((label, str(p)))
        _VCF_PATHS = out
        _VCF_SOURCE = "paths"
        _VCF_DIR_PATH = None
        return
    raw_dir = _get_vcf_dir_env()
    if raw_dir:
        dir_path = Path(raw_dir)
        if not dir_path.is_absolute():
            dir_path = (_APP_ROOT / raw_dir).resolve()
        else:
            dir_path = dir_path.resolve()
        _VCF_SOURCE = "dir"
        _VCF_DIR_PATH = dir_path
        _VCF_PATHS = _discover_vcf_in_dir(dir_path)  # initial list; get_vcf_list will re-scan
        return
    # No env set: try default dirs
    for raw_dir in _DEFAULT_VCF_DIRS:
        dir_path = Path(raw_dir)
        if not dir_path.is_absolute():
            dir_path = (_APP_ROOT / raw_dir).resolve()
        else:
            dir_path = dir_path.resolve()
        found = _discover_vcf_in_dir(dir_path)
        if found:
            _VCF_SOURCE = "dir"
            _VCF_DIR_PATH = dir_path
            _VCF_PATHS = found
            return
    _VCF_SOURCE = ""
    _VCF_DIR_PATH = None
    _VCF_PATHS = []


def get_vcf_list() -> list[dict]:
    """Return list of { path, label } for dropdown. Re-scans folder on each call when using a directory."""
    init_app()
    if _VCF_SOURCE == "dir" and _VCF_DIR_PATH is not None:
        # Re-read folder so newly added .vcf / .vcf.gz files appear
        items = _discover_vcf_in_dir(_VCF_DIR_PATH)
        return [{"path": path, "label": label} for label, path in items]
    return [{"path": path, "label": label} for label, path in _VCF_PATHS]


def get_vcf_path_for_request(vcf_path: str) -> str | None:
    """Return resolved absolute path if vcf_path is in allowed list, else None. Re-scans dir when using a directory."""
    init_app()
    if _VCF_SOURCE == "dir" and _VCF_DIR_PATH is not None:
        allowed = {path for _, path in _discover_vcf_in_dir(_VCF_DIR_PATH)}
    else:
        allowed = {path for _, path in _VCF_PATHS}
    candidate = str(Path(vcf_path).resolve())
    return candidate if candidate in allowed else None


# --- Pydantic schemas ---


class SearchRequest(BaseModel):
    """Request body for variant search."""

    query: str = Field(..., min_length=1, description="Gene symbol (e.g. FAH) or location (e.g. 80472479 or chr15:80472479)")
    vcf_path: str = Field(..., min_length=1, description="VCF path from the configured list")


class VcfItem(BaseModel):
    """Single VCF entry for the dropdown."""

    path: str = Field(..., description="Resolved absolute path to the VCF file")
    label: str = Field(..., description="Display label (usually filename)")


class VcfListResponse(BaseModel):
    """Response for GET /api/vcf-list."""

    vcf_list: list[VcfItem] = Field(..., description="List of VCFs available for search")


class SearchResponse(BaseModel):
    """Response for POST /api/search."""

    message: str = Field(..., description="7-step variant message (SOP)")
    query: str = Field(..., description="Query that was searched")
    vcf_path: str = Field(..., description="VCF path that was used")


# --- API ---


@router.get("/api/debug/vcf-config")
def api_debug_vcf_config():
    """Debug: VCF env, app root, and resulting vcf list (fresh scan when using dir). No auth."""
    init_app()
    current = get_vcf_list()
    return {
        "vcf_dir_env": _get_vcf_dir_env() or "(not set)",
        "vcf_paths_env": _get_vcf_paths_env() or "(empty)",
        "app_root": str(_APP_ROOT),
        "vcf_source": _VCF_SOURCE or "(none)",
        "vcf_list_count": len(current),
        "vcf_list_labels": [item["label"] for item in current],
    }


@router.get("/api/vcf-list", response_model=VcfListResponse)
def api_vcf_list() -> VcfListResponse:
    """Return list of VCFs (.vcf and .vcf.gz only). Re-scans folder on each call when using a directory."""
    items = [VcfItem(path=item["path"], label=item["label"]) for item in get_vcf_list()]
    return VcfListResponse(vcf_list=items)


@router.post("/api/search", response_model=SearchResponse)
def api_search(body: SearchRequest) -> SearchResponse:
    """Run variant search on the selected VCF; return 7-step message."""
    from app.auto_variant import run_variant_search

    resolved = get_vcf_path_for_request(body.vcf_path)
    if not resolved:
        raise HTTPException(
            status_code=400,
            detail="Selected VCF path is not in the configured list (VCF_PATHS).",
        )
    query = body.query.strip()
    try:
        message = run_variant_search(resolved, query, sample_name=None)
        return SearchResponse(message=message, query=query, vcf_path=resolved)
    except FileNotFoundError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
