"""
MyVariant.info: fetch full variant annotation and filter by source (gnomAD, ExAC, etc.).
Uses hg19 genomic HGVS, e.g. chr15:g.80472479C>T. No API key required.
"""
from __future__ import annotations

import urllib.parse
from typing import Any

import httpx

from lookup_types import PopulationFreq

MYVARIANT_BASE = "https://myvariant.info/v1/variant"


def hg19_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build MyVariant hg19 genomic HGVS (chrN:g.POSref>alt)."""
    c = str(chrom).strip()
    if not c.lower().startswith("chr"):
        c = "chr" + c
    return f"{c}:g.{pos}{ref}>{alt}"


def fetch_variant(chrom: str, pos: int, ref: str, alt: str) -> dict[str, Any] | None:
    """
    Fetch full variant annotation from MyVariant.info (all fields).
    Returns the raw JSON object or None on error.
    """
    vid = hg19_variant_id(chrom, pos, ref, alt)
    url = MYVARIANT_BASE + "/" + urllib.parse.quote(vid, safe="")
    try:
        with httpx.Client(timeout=20.0) as client:
            r = client.get(url)
            r.raise_for_status()
            return r.json()
    except Exception:
        return None


def get_gnomad(data: dict[str, Any] | None) -> dict[str, Any]:
    """Extract all gnomAD-related fields from a MyVariant variant object."""
    if not data or not isinstance(data, dict):
        return {}
    return {
        k: v
        for k, v in data.items()
        if k.startswith("gnomad")
    }


def get_exac(data: dict[str, Any] | None) -> dict[str, Any]:
    """Extract ExAC field from a MyVariant variant object."""
    if not data or not isinstance(data, dict):
        return {}
    exac = data.get("exac")
    return {"exac": exac} if exac is not None else {}


def get_clinvar(data: dict[str, Any] | None) -> dict[str, Any]:
    """Extract ClinVar-related fields from a MyVariant variant object."""
    if not data or not isinstance(data, dict):
        return {}
    out = {}
    if "clinvar" in data:
        out["clinvar"] = data["clinvar"]
    return out


def get_dbsnp(data: dict[str, Any] | None) -> dict[str, Any]:
    """Extract dbSNP field from a MyVariant variant object."""
    if not data or not isinstance(data, dict):
        return {}
    out = {}
    if "dbsnp" in data:
        out["dbsnp"] = data["dbsnp"]
    return out


def get_cadd(data: dict[str, Any] | None) -> dict[str, Any]:
    """Extract CADD field from a MyVariant variant object."""
    if not data or not isinstance(data, dict):
        return {}
    out = {}
    if "cadd" in data:
        out["cadd"] = data["cadd"]
    return out


def get_dbnsfp(data: dict[str, Any] | None) -> dict[str, Any]:
    """Extract dbNSFP field from a MyVariant variant object."""
    if not data or not isinstance(data, dict):
        return {}
    out = {}
    if "dbnsfp" in data:
        out["dbnsfp"] = data["dbnsfp"]
    return out


def get_all_sources(data: dict[str, Any] | None) -> dict[str, Any]:
    """
    Return variant object with only annotation sources (no _id, _version, _score).
    Useful to see which sources are present.
    """
    if not data or not isinstance(data, dict):
        return {}
    skip = {"_id", "_version", "_score"}
    return {k: v for k, v in data.items() if k not in skip and not k.startswith("_")}


def _float_af(obj: Any) -> float | None:
    if obj is None:
        return None
    if isinstance(obj, (int, float)):
        return float(obj)
    if isinstance(obj, dict) and "af" in obj:
        return _float_af(obj["af"])
    return None


def lookup_gnomad_exac(chrom: str, pos: int, ref: str, alt: str) -> PopulationFreq | None:
    """
    Fetch variant from MyVariant and return PopulationFreq (gnomAD max AF, ExAC).
    Convenience for report Step 4; uses fetch_variant + get_gnomad/get_exac.
    """
    data = fetch_variant(chrom, pos, ref, alt)
    if not data:
        return None
    out = PopulationFreq()
    gnomad = get_gnomad(data)
    gnomad_af_max = None
    for key in ("gnomad_exome", "gnomad_genome"):
        g = gnomad.get(key)
        if not isinstance(g, dict):
            continue
        out.raw_gnomad = g
        af = _float_af(g.get("af")) or g.get("allele_freq")
        if af is not None:
            try:
                af = float(af)
                if gnomad_af_max is None or af > gnomad_af_max:
                    gnomad_af_max = af
            except (TypeError, ValueError):
                pass
        if out.gnomad_ac is None and g.get("ac") is not None:
            try:
                out.gnomad_ac = int(g["ac"])
            except (TypeError, ValueError):
                pass
        if out.gnomad_an is None and g.get("an") is not None:
            try:
                out.gnomad_an = int(g["an"])
            except (TypeError, ValueError):
                pass
    out.gnomad_af_max = gnomad_af_max

    exac = get_exac(data).get("exac")
    if isinstance(exac, dict):
        out.raw_exac = exac
        af = _float_af(exac.get("af"))
        if af is not None:
            out.exac_af = float(af)
        if exac.get("ac") is not None:
            try:
                out.exac_ac = int(exac["ac"])
            except (TypeError, ValueError):
                pass
        if exac.get("an") is not None:
            try:
                out.exac_an = int(exac["an"])
            except (TypeError, ValueError):
                pass
    return out
