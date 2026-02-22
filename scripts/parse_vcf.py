#!/usr/bin/env python3
"""Parse VCF row by row using cyvcf2; show INFO, FORMAT, Exomiser and FUNC.

Uses only cyvcf2 for parsing (no custom VCF/line parsing).
Exomiser and FUNC are extracted as (header_tuple, list of value tuples) so you can
convert to dict when required: dict(zip(header, row)).

VCF columns:
  - FORMAT: per-sample field names (GT, GQ, DP, ...).
  - Last column: sample format values (from cyvcf2).
"""

import ast
import re
import sys
from pathlib import Path


def get_all_info_descriptions(vcf):
    """
    Build INFO key -> description from the VCF header once.
    Single-use: call once per file, reuse the returned dict for all rows.
    """
    out = {}
    raw = getattr(vcf, "raw_header", None) or ""
    # ##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">
    for m in re.finditer(r"##INFO=<ID=([^,]+),[^>]*Description=\"([^\"]*)\"", raw):
        key, desc = m.group(1), m.group(2)
        out[key] = desc.strip('"') if desc else ""
    return out

from cyvcf2 import VCF


# Exomiser: pipe-separated fields per VCF ##INFO Exomiser description
EXOMISER_HEADER = (
    "RANK",
    "ID",
    "GENE_SYMBOL",
    "ENTREZ_GENE_ID",
    "MOI",
    "P-VALUE",
    "EXOMISER_GENE_COMBINED_SCORE",
    "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE",
    "EXOMISER_VARIANT_SCORE",
    "CONTRIBUTING_VARIANT",
    "WHITELIST_VARIANT",
    "FUNCTIONAL_CLASS",
    "HGVS",
    "EXOMISER_ACMG_CLASSIFICATION",
    "EXOMISER_ACMG_EVIDENCE",
    "EXOMISER_ACMG_DISEASE_ID",
    "EXOMISER_ACMG_DISEASE_NAME",
)

# FUNC: canonical key order for tuple rows (convert to dict with dict(zip(FUNC_HEADER, row)))
FUNC_HEADER = (
    "origPos",
    "origRef",
    "normalizedRef",
    "gene",
    "normalizedPos",
    "normalizedAlt",
    "gt",
    "coding",
    "transcript",
    "function",
    "protein",
    "location",
    "origAlt",
    "exon",
    "codon",
    "polyphen",
    "sift",
    "grantham",
    "CLNREVSTAT1",
    "CLNACC1",
    "CLNSIG1",
    "CLNID1",
)


def info_to_dict(record):
    """Build dict from record INFO using cyvcf2 (no custom parsing)."""
    return dict(record.INFO)


def _format_value_for_display(record, key, idx):
    """One sample's value for a FORMAT key using cyvcf2; for display."""
    val = record.format(key)
    if val is None or (hasattr(val, "__len__") and len(val) == 0):
        return None
    v0 = val[idx]
    if key == "GT" and getattr(record, "genotypes", None) is not None:
        g = record.genotypes[idx]
        return "/".join(str(x) for x in g[:2]) if len(g) >= 2 else str(v0)
    if hasattr(v0, "item"):
        try:
            return v0.item()
        except ValueError:
            return list(v0) if hasattr(v0, "__len__") else v0
    if hasattr(v0, "__len__") and not isinstance(v0, (str, bytes)):
        return [x.item() if hasattr(x, "item") else x for x in v0]
    return str(v0) if isinstance(v0, (bytes,)) else v0


def extract_exomiser(record):
    """
    Extract Exomiser from INFO as (header_tuple, list of value tuples).
    Convert to dict when required: dict(zip(EXOMISER_HEADER, row)).
    Returns (EXOMISER_HEADER, [tuple, ...]); empty list if no Exomiser.
    """
    raw = record.INFO.get("Exomiser")
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return (EXOMISER_HEADER, [])

    raw = str(raw).strip()
    if not raw:
        return (EXOMISER_HEADER, [])

    rows = []
    # One or more blocks like {a|b|c|...}; strip braces and split by |
    for block in re.split(r"\}\s*;\s*\{", raw):
        block = block.strip(" {}")
        if not block:
            continue
        vals = block.split("|")
        # Pad or trim to header length
        n = len(EXOMISER_HEADER)
        if len(vals) < n:
            vals = list(vals) + [None] * (n - len(vals))
        else:
            vals = vals[:n]
        rows.append(tuple(vals))
    return (EXOMISER_HEADER, rows)


def extract_func(record):
    """
    Extract FUNC from INFO as (header_tuple, list of value tuples).
    FUNC is stored as a Python literal list of dicts; we parse and convert to rows.
    Convert to dict when required: dict(zip(FUNC_HEADER, row)).
    Returns (FUNC_HEADER, [tuple, ...]); empty list if no FUNC or parse error.
    """
    raw = record.INFO.get("FUNC")
    if raw is None:
        return (FUNC_HEADER, [])

    if isinstance(raw, (list, dict)):
        items = raw if isinstance(raw, list) else [raw]
    else:
        try:
            items = ast.literal_eval(raw)
        except (ValueError, SyntaxError):
            return (FUNC_HEADER, [])
        if not isinstance(items, list):
            items = [items] if items is not None else []

    rows = []
    for d in items:
        if not isinstance(d, dict):
            continue
        row = tuple(d.get(k) for k in FUNC_HEADER)
        rows.append(row)
    return (FUNC_HEADER, rows)


def print_row(
    record,
    info_dict,
    row_num,
    *,
    info_descriptions=None,
    print_description=True,
    sample_index=-1,
    exomiser_header=None,
    exomiser_rows=None,
    func_header=None,
    func_rows=None,
):
    """Print one VCF row. FORMAT and per-sample values from cyvcf2 (record.FORMAT, record.format).
    Exomiser/FUNC shown as (header, rows); to dict: dict(zip(header, row)).
    info_descriptions: optional dict of INFO key -> header description (for the # comment)."""
    print(f"\n--- Row {row_num} ---")
    print(f"CHROM={record.CHROM} POS={record.POS} REF={record.REF} ALT={record.ALT} QUAL={record.QUAL} FILTER={record.FILTER}")

    # FORMAT and per-sample values from cyvcf2 directly
    if record.FORMAT:
        keys = list(record.FORMAT)
        print("\nFORMAT (field names):")
        print("  " + ":".join(keys))
        raw = record.format(keys[0])
        n_s = len(raw) if raw is not None else 0
        if n_s > 0:
            idx = (n_s + sample_index) if sample_index < 0 else sample_index
            idx = max(0, min(idx, n_s - 1))
            print("\nPer-sample (format values):")
            for k in keys:
                v = _format_value_for_display(record, k, idx)
                print(f"  {k}: {v}")

    print("\nINFO (as dict):")
    desc_map = info_descriptions if info_descriptions is not None else {}
    for k, v in sorted(info_dict.items()):
        if k in ("Exomiser", "FUNC"):
            line = f"  {k}: (see below as header + rows)"
        else:
            line = f"  {k}: {v}"
        if print_description:
            desc = desc_map.get(k, "")
            if desc:
                line = f"{line}  # {desc}"
        print(line)

    if exomiser_header is not None and exomiser_rows is not None:
        print("\nExomiser (header, rows) — to dict: dict(zip(header, row)):")
        print("  header:", exomiser_header)
        if not exomiser_rows:
            print("  (no entries)")
        for j, row in enumerate(exomiser_rows):
            print(f"  row[{j}]: {row}")
            if j == 0 and row:
                print(f"  -> dict(zip(header, row)): {dict(zip(exomiser_header, row))}")

    if func_header is not None and func_rows is not None:
        print("\nFUNC (header, rows) — to dict: dict(zip(header, row)):")
        print("  header:", func_header)
        if not func_rows:
            print("  (no entries)")
        for j, row in enumerate(func_rows):
            print(f"  row[{j}]: {row}")
            if j == 0 and row:
                print(f"  -> dict(zip(header, row)): {dict(zip(func_header, row))}")


def main():
    vcf_path = sys.argv[1] if len(sys.argv) > 1 else "WES166_26032846_BSC_IN.vcf"
    path = Path(vcf_path)
    if not path.exists():
        print(f"File not found: {path}", file=sys.stderr)
        sys.exit(1)

    vcf = VCF(str(path))
    # Build INFO descriptions once from header (single-use), then reuse for every row
    info_descriptions = get_all_info_descriptions(vcf)

    for i, record in enumerate(vcf):
        info_dict = info_to_dict(record)
        exomiser_header, exomiser_rows = extract_exomiser(record)
        func_header, func_rows = extract_func(record)

        print_row(
            record,
            info_dict,
            i + 1,
            info_descriptions=info_descriptions,
            print_description=True,
            sample_index=-1,
            exomiser_header=exomiser_header,
            exomiser_rows=exomiser_rows,
            func_header=func_header,
            func_rows=func_rows,
        )

        try:
            input("\n[Press Enter for next row, Ctrl+C to stop] ")
        except (KeyboardInterrupt, EOFError):
            print("\nStopped.")
            break

    vcf.close()


if __name__ == "__main__":
    main()
