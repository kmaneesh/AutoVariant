"""
AutoVariant FastAPI app: web (static), variant search API. No AI.
Routes:
  GET  /              – Index page (web)
  GET  /health        – Health check (web)
  GET  /static/*      – Static assets
  GET  /api/vcf-list  – VCF list from config (variant)
  POST /api/search    – Variant search: query + vcf_path → 7-step message (variant)
"""
from contextlib import asynccontextmanager
from pathlib import Path

from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from app import variant, web


@asynccontextmanager
async def lifespan(app: FastAPI):
    """Load VCF list from config at startup."""
    variant.init_app()
    yield


app = FastAPI(
    title="AutoVariant",
    description="Variant search by gene or location; 7-step SOP message. No AI.",
    lifespan=lifespan,
)

APP_DIR = Path(__file__).resolve().parent
STATIC_DIR = APP_DIR / "static"
TEMPLATES_DIR = APP_DIR / "templates"

web.STATIC_DIR = STATIC_DIR if STATIC_DIR.exists() else None
web.TEMPLATES_DIR = TEMPLATES_DIR if TEMPLATES_DIR.exists() else None

app.include_router(web.router)
app.include_router(variant.router)

if STATIC_DIR.exists():
    app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
