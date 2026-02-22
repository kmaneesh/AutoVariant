"""
Variant lookups: MyVariant.info (gnomAD, ExAC, etc.), VarSome, Franklin Genoox.
Re-exports from myvariant, varsome, franklin, lookup_types.
"""
from __future__ import annotations

from lookup_types import LookupResult, PopulationFreq
from myvariant import fetch_variant, get_gnomad, get_exac, lookup_gnomad_exac
from varsome import lookup_varsome
from franklin import lookup_franklin

__all__ = [
    "fetch_variant",
    "get_gnomad",
    "get_exac",
    "lookup_gnomad_exac",
    "lookup_varsome",
    "lookup_franklin",
    "LookupResult",
    "PopulationFreq",
]


def run_lookups(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    varsome_api_key: str | None = None,
) -> LookupResult:
    """Run all available lookups for one variant (chr, pos, ref, alt)."""
    result = LookupResult()
    result.gnomad_exac = lookup_gnomad_exac(chrom, pos, ref, alt)
    result.varsome = lookup_varsome(chrom, pos, ref, alt, api_key=varsome_api_key)
    result.franklin = lookup_franklin(chrom, pos, ref, alt)
    return result
