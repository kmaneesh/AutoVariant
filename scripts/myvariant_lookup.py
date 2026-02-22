#!/usr/bin/env python3
"""
Standalone MyVariant.info lookup: fetch all variant annotations and optionally
filter by source (gnomAD, ExAC, ClinVar, dbSNP, CADD, dbNSFP).

Usage:
  python scripts/myvariant_lookup.py "chr15:80472479 C T"
  python scripts/myvariant_lookup.py 15 80472479 C T
  python scripts/myvariant_lookup.py chr15:80472479-C-T --exac --gnomad
  python scripts/myvariant_lookup.py chr15:80472479-C-T -o variant.json

Variant can be:
  - "chr15:80472479 C T" or "15 80472479 C T"
  - "chr15:80472479-C-T" or "15-80472479-C-T"
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

# Allow running from repo root: scripts/ on path so lookup modules are found
_scripts_dir = Path(__file__).resolve().parent
if str(_scripts_dir) not in sys.path:
    sys.path.insert(0, str(_scripts_dir))
from myvariant import (
    fetch_variant,
    get_gnomad,
    get_exac,
    get_clinvar,
    get_dbsnp,
    get_cadd,
    get_dbnsfp,
    get_all_sources,
    hg19_variant_id,
)


def parse_variant_arg(s: str) -> tuple[str, int, str, str] | None:
    """Parse 'chr15:80472479 C T' or '15-80472479-C-T' -> (chrom, pos, ref, alt)."""
    s = s.strip()
    # chr15:80472479-C-T or 15-80472479-C-T
    m = re.match(r"^(chr)?(\d+|X|Y|MT|M)[:\-](\d+)[\-\s]+([ACGT]+)[\-\s]+([ACGT]+)$", s, re.I)
    if m:
        chrom = (m.group(1) or "") + m.group(2)
        if not chrom.lower().startswith("chr"):
            chrom = "chr" + chrom
        return (chrom, int(m.group(3)), m.group(4).upper(), m.group(5).upper())
    # chr15:80472479 C T (space-separated)
    parts = s.replace(":", " ").replace("-", " ").split()
    if len(parts) >= 4:
        c, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
        if not c.upper().startswith("CHR"):
            c = "chr" + c
        try:
            return (c, int(pos_s), ref.upper(), alt.upper())
        except ValueError:
            pass
    return None


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Fetch full variant annotation from MyVariant.info and optionally filter by source."
    )
    ap.add_argument(
        "variant",
        nargs="?",
        help="Variant: 'chr15:80472479 C T' or '15 80472479 C T' or 'chr15:80472479-C-T'",
    )
    ap.add_argument(
        "pos",
        nargs="?",
        type=int,
        help="Position (if variant given as chrom only).",
    )
    ap.add_argument("ref", nargs="?", help="Reference allele (if chrom pos ref alt given separately).")
    ap.add_argument("alt", nargs="?", help="Alternate allele.")
    ap.add_argument("--gnomad", action="store_true", help="Print only gnomAD-related fields.")
    ap.add_argument("--exac", action="store_true", help="Print only ExAC field.")
    ap.add_argument("--clinvar", action="store_true", help="Print only ClinVar field.")
    ap.add_argument("--dbsnp", action="store_true", help="Print only dbSNP field.")
    ap.add_argument("--cadd", action="store_true", help="Print only CADD field.")
    ap.add_argument("--dbnsfp", action="store_true", help="Print only dbNSFP field.")
    ap.add_argument("--sources", action="store_true", help="Print only source keys present (no _id/_version).")
    ap.add_argument("-o", "--output", metavar="FILE", help="Write JSON to file instead of stdout.")
    ap.add_argument("--indent", type=int, default=2, metavar="N", help="JSON indent (default 2).")
    args = ap.parse_args()

    # Build chrom, pos, ref, alt
    if args.variant and args.pos is not None and args.ref and args.alt:
        c = args.variant.strip()
        if not c.lower().startswith("chr"):
            c = "chr" + c
        chrom, pos, ref, alt = c, args.pos, args.ref.upper(), args.alt.upper()
    elif args.variant:
        parsed = parse_variant_arg(args.variant)
        if not parsed:
            print("Could not parse variant. Use e.g. 'chr15:80472479 C T' or '15-80472479-C-T'", file=sys.stderr)
            return 1
        chrom, pos, ref, alt = parsed
    else:
        ap.print_help()
        return 0

    data = fetch_variant(chrom, pos, ref, alt)
    if not data:
        print(f"Variant not found or error: {hg19_variant_id(chrom, pos, ref, alt)}", file=sys.stderr)
        return 1

    filters = []
    if args.gnomad:
        filters.append(("gnomAD", get_gnomad))
    if args.exac:
        filters.append(("ExAC", get_exac))
    if args.clinvar:
        filters.append(("ClinVar", get_clinvar))
    if args.dbsnp:
        filters.append(("dbSNP", get_dbsnp))
    if args.cadd:
        filters.append(("CADD", get_cadd))
    if args.dbnsfp:
        filters.append(("dbNSFP", get_dbnsfp))
    if args.sources:
        filters.append(("sources", get_all_sources))

    if filters:
        out = {}
        for name, fn in filters:
            out[name] = fn(data)
        payload = out
    else:
        payload = data

    indent = args.indent if args.indent >= 0 else None
    text = json.dumps(payload, indent=indent, default=str)

    if args.output:
        Path(args.output).write_text(text, encoding="utf-8")
        print(f"Wrote {args.output}", file=sys.stderr)
    else:
        print(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
