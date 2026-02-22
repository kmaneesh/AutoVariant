#!/usr/bin/env python3
"""
AutoVariant: search a VCF by gene symbol or location and output
a final variant message following the WES Variant Details SOP (7 steps).
Step 2 uses amino acid residue properties (AutoVariant.mdc) for wording.
Step 4: optional gnomAD/ExAC via MyVariant.info (httpx). Self-contained; no lookup package.
"""
import argparse
import ast
import re
import sys
import urllib.parse
from pathlib import Path

try:
    from cyvcf2 import VCF
    CYVCF2_AVAILABLE = True
except ImportError:
    CYVCF2_AVAILABLE = False

try:
    from vcfpy import Reader
    VCFPY_AVAILABLE = True
except ImportError:
    VCFPY_AVAILABLE = False

try:
    import httpx
    HTTPX_AVAILABLE = True
except ImportError:
    HTTPX_AVAILABLE = False


def _step4_population_message(chrom: str, pos: int, ref: str, alt: str) -> str | None:
    """Fetch MyVariant.info for variant; return Step 4 sentence (gnomAD/ExAC) or None."""
    if not HTTPX_AVAILABLE:
        return None
    c = str(chrom).strip()
    if not c.lower().startswith("chr"):
        c = "chr" + c
    vid = f"{c}:g.{pos}{ref}>{alt}"
    url = "https://myvariant.info/v1/variant/" + urllib.parse.quote(vid, safe="")
    try:
        with httpx.Client(timeout=20.0) as client:
            r = client.get(url)
            r.raise_for_status()
            data = r.json()
    except Exception:
        return None

    def _af(obj):
        if obj is None:
            return None
        if isinstance(obj, (int, float)):
            return float(obj)
        if isinstance(obj, dict) and "af" in obj:
            return _af(obj["af"])
        return None

    gnomad_af_max = None
    for key in ("gnomad_exome", "gnomad_genome"):
        g = data.get(key)
        if not isinstance(g, dict):
            continue
        af = _af(g.get("af")) or g.get("allele_freq")
        if af is not None:
            try:
                af = float(af)
                if gnomad_af_max is None or af > gnomad_af_max:
                    gnomad_af_max = af
            except (TypeError, ValueError):
                pass
    exac_af = None
    exac = data.get("exac")
    if isinstance(exac, dict):
        exac_af = _af(exac.get("af"))

    def _pct(af_val):
        if af_val is None:
            return None
        pct = af_val * 100
        return f"{pct:.2f}%" if pct >= 0.01 else f"{pct:.4f}%"

    if gnomad_af_max is not None and exac_af is not None:
        return (
            f"This variant has minor allele frequency of {_pct(gnomad_af_max)} "
            f"and {_pct(exac_af)} in gnomAD (max) and ExAC database respectively."
        )
    if gnomad_af_max is not None:
        return (
            f"This variant has minor allele frequency of {_pct(gnomad_af_max)} in gnomAD (max)."
        )
    if exac_af is not None:
        return f"This variant has minor allele frequency of {_pct(exac_af)} in ExAC database."
    return None


# Exomiser pipe-separated field indices (from VCF header)
EXOMISER_KEYS = [
    "RANK", "ID", "GENE_SYMBOL", "ENTREZ_GENE_ID", "MOI", "P-VALUE",
    "EXOMISER_GENE_COMBINED_SCORE", "EXOMISER_GENE_PHENO_SCORE",
    "EXOMISER_GENE_VARIANT_SCORE", "EXOMISER_VARIANT_SCORE",
    "CONTRIBUTING_VARIANT", "WHITELIST_VARIANT", "FUNCTIONAL_CLASS", "HGVS",
    "EXOMISER_ACMG_CLASSIFICATION", "EXOMISER_ACMG_EVIDENCE",
    "EXOMISER_ACMG_DISEASE_ID", "EXOMISER_ACMG_DISEASE_NAME",
]
EXOMISER_HEADER = tuple(EXOMISER_KEYS)

# FUNC: canonical key order for tuple rows (dict(zip(FUNC_HEADER, row)))
FUNC_HEADER = (
    "origPos", "origRef", "normalizedRef", "gene", "normalizedPos", "normalizedAlt",
    "gt", "coding", "transcript", "function", "protein", "location", "origAlt",
    "exon", "codon", "polyphen", "sift", "grantham",
    "CLNREVSTAT1", "CLNACC1", "CLNSIG1", "CLNID1",
)


def _format_value_for_display(record, key: str, idx: int):
    """One sample's value for a FORMAT key using cyvcf2."""
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
    record: cyvcf2 Variant. To dict: dict(zip(EXOMISER_HEADER, row)).
    """
    raw = record.INFO.get("Exomiser")
    if raw is None or (isinstance(raw, str) and not raw.strip()):
        return (EXOMISER_HEADER, [])
    raw = str(raw).strip()
    if not raw:
        return (EXOMISER_HEADER, [])
    rows = []
    for block in re.split(r"\}\s*;\s*\{", raw):
        block = block.strip(" {}")
        if not block:
            continue
        vals = block.split("|")
        n = len(EXOMISER_HEADER)
        vals = list(vals) + [None] * (n - len(vals)) if len(vals) < n else vals[:n]
        rows.append(tuple(vals))
    return (EXOMISER_HEADER, rows)


def extract_func(record):
    """
    Extract FUNC from INFO as (header_tuple, list of value tuples).
    record: cyvcf2 Variant. To dict: dict(zip(FUNC_HEADER, row)).
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
    rows = [tuple(d.get(k) for k in FUNC_HEADER) for d in items if isinstance(d, dict)]
    return (FUNC_HEADER, rows)


# Amino acid residue properties for Step 2 wording (from AutoVariant.mdc / biochemistry refs)
# 3-letter code -> (full name, property phrase for "which is ...")
AMINO_ACID_PROPERTIES = {
    "Ala": ("alanine", "nonpolar, hydrophobic"),
    "Arg": ("arginine", "basic, positively charged"),
    "Asn": ("asparagine", "polar uncharged (amide)"),
    "Asp": ("aspartic acid", "acidic, negatively charged"),
    "Cys": ("cysteine", "polar uncharged (thiol, can form disulfide)"),
    "Gln": ("glutamine", "polar uncharged (amide)"),
    "Glu": ("glutamic acid", "acidic, negatively charged"),
    "Gly": ("glycine", "small, neutral (no side chain)"),
    "His": ("histidine", "basic, positively charged (imidazole)"),
    "Ile": ("isoleucine", "nonpolar, hydrophobic"),
    "Leu": ("leucine", "nonpolar, hydrophobic"),
    "Lys": ("lysine", "basic, positively charged"),
    "Met": ("methionine", "nonpolar, hydrophobic (sulfur)"),
    "Phe": ("phenylalanine", "nonpolar, hydrophobic (aromatic)"),
    "Pro": ("proline", "nonpolar, cyclic (rigid)"),
    "Ser": ("serine", "polar uncharged (hydroxyl)"),
    "Thr": ("threonine", "polar uncharged (hydroxyl)"),
    "Trp": ("tryptophan", "nonpolar, hydrophobic (aromatic)"),
    "Tyr": ("tyrosine", "polar uncharged (aromatic hydroxyl)"),
    "Val": ("valine", "nonpolar, hydrophobic"),
}


def _parse_protein_substitution(protein: str) -> tuple[str, int, str] | None:
    """Parse p.X###Y or p.(X###Y) to (ref_3letter, position, alt_3letter). Returns None if not a simple substitution."""
    if not protein or "?" in protein or "fs" in protein.lower() or "*" in protein:
        return None
    protein = protein.strip().strip("()")
    if not protein.startswith("p."):
        return None
    rest = protein[2:].strip()
    # Match 3-letter + digits + 3-letter (e.g. Thr325Met, Gly1098Asp)
    m = re.match(r"^([A-Za-z]{3})(\d+)([A-Za-z]{3})$", rest)
    if m:
        ref_aa = m.group(1)
        pos = int(m.group(2))
        alt_aa = m.group(3)
        if ref_aa in AMINO_ACID_PROPERTIES and alt_aa in AMINO_ACID_PROPERTIES:
            return (ref_aa, pos, alt_aa)
    return None


def _format_amino_acid_change(protein: str) -> str:
    """Format Step 2 sentence using residue properties when possible."""
    parsed = _parse_protein_substitution(protein)
    if parsed:
        ref_aa, pos, alt_aa = parsed
        ref_name, ref_prop = AMINO_ACID_PROPERTIES[ref_aa]
        alt_name, alt_prop = AMINO_ACID_PROPERTIES[alt_aa]
        return (
            f"This variant ({protein}) replaces {ref_name}, which is {ref_prop}, "
            f"with {alt_name}, which is {alt_prop}, at the {pos}th amino acid position."
        )
    return (
        f"This variant ({protein}) is reported at the protein level. "
        "Refer to transcript for residue properties (e.g. neutral vs acidic/polar)."
    )


def parse_func(func_str: str) -> list[dict]:
    """Parse INFO FUNC string (Ion Reporter style: list of dicts as string)."""
    if not func_str or func_str == ".":
        return []
    # FUNC is often stored as [{'gene':'FAH', ...}] with single quotes
    try:
        return ast.literal_eval(func_str)
    except (ValueError, SyntaxError):
        return []


def parse_exomiser(exo_str: str) -> list[dict]:
    """Parse INFO Exomiser pipe-separated value. Returns list of dicts per allele."""
    if not exo_str or exo_str == ".":
        return []
    result = []
    # Exomiser can be {1|15-80472479-C-T_AR|FAH|...}; strip braces and split by };
    raw = exo_str.strip("{}")
    for block in raw.split("};"):
        block = block.strip().strip("{")
        if not block:
            continue
        parts = block.split("|")
        row = {}
        for i, key in enumerate(EXOMISER_KEYS):
            if i < len(parts):
                row[key] = parts[i].strip('"')
        result.append(row)
    return result


def parse_gt(gt_str: str) -> tuple[str, str]:
    """Return (zygosity_label, gt_raw)."""
    if not gt_str or gt_str in (".", "./.", ".|."):
        return "unknown", gt_str or "."
    gt = gt_str.replace("|", "/")
    if gt in ("0/0", "0|0"):
        return "homozygous reference", gt
    if gt in ("1/1", "1|1"):
        return "homozygous", gt
    if gt in ("0/1", "1/0", "0|1", "1|0"):
        return "heterozygous", gt
    if "/" in gt and "1" in gt and "0" not in gt:
        return "homozygous (alternate)", gt
    return "heterozygous (multi-alt)", gt


def _parse_info(info_str: str) -> dict:
    """Parse INFO string; split only on ; that start a new KEY= (so values can contain ;)."""
    if not info_str or info_str == ".":
        return {}
    # Split on semicolon followed by a valid INFO key (KEY=)
    parts = re.split(r";(?=[A-Za-z_][A-Za-z0-9_]*=)", info_str)
    info = {}
    for part in parts:
        if "=" not in part:
            continue
        key, _, value = part.partition("=")
        if key and value:
            info[key] = value
    return info


def _parse_format_sample(format_keys: list[str], sample_str: str) -> dict:
    """Parse FORMAT and sample column into sample dict."""
    if not format_keys or not sample_str:
        return {}
    values = sample_str.split(":")
    return dict(zip(format_keys, values)) if len(values) >= len(format_keys) else {}


def _normalize_chrom(c: str) -> str:
    """Normalize chromosome to chrN form for comparison (lowercase)."""
    if not c:
        return c
    c = str(c).strip()
    if c.upper().startswith("CHR"):
        return ("chr" + c[3:]) if c[0] != "c" else c
    return "chr" + c


def _match_record_dict(record: dict, search: str) -> bool:
    """Return True if record (dict with chrom, pos, info) matches search."""
    search = search.strip().upper()
    chrom, pos = record.get("chrom", ""), record.get("pos", 0)
    info = record.get("info", {})
    if ":" in search:
        chrom_part, pos_part = search.split(":", 1)
        chrom_part = _normalize_chrom(chrom_part.strip())
        chrom_norm = _normalize_chrom(chrom)
        try:
            pos_int = int(pos_part.strip())
            if chrom_norm == chrom_part and pos == pos_int:
                return True
        except ValueError:
            pass
    else:
        try:
            if pos == int(search):
                return True
        except (ValueError, TypeError):
            pass
    func_str = info.get("FUNC")
    if func_str:
        for f in parse_func(str(func_str)):
            if isinstance(f, dict) and f.get("gene", "").upper() == search:
                return True
    exo_str = info.get("Exomiser")
    if exo_str:
        for exo in parse_exomiser(str(exo_str)):
            if exo.get("GENE_SYMBOL", "").upper() == search:
                return True
    return False


def _vcfpy_record_to_dict(record) -> dict:
    """Convert a vcfpy Record to our standard record dict."""
    alts = []
    if record.ALT:
        for a in record.ALT:
            if hasattr(a, "value"):
                alts.append(str(a.value))
            else:
                alts.append(str(a))
    samples = {}
    if record.calls:
        for call in record.calls:
            data = dict(call.data) if hasattr(call, "data") else {}
            # GT may be a tuple (0, 1); normalize to "0/1"
            gt = data.get("GT")
            if gt is not None and not isinstance(gt, str) and hasattr(gt, "__len__") and len(gt) >= 2:
                data = dict(data)
                data["GT"] = f"{gt[0]}/{gt[1]}" if gt[0] is not None and gt[1] is not None else "."
            samples[call.sample] = data
    return {
        "chrom": record.CHROM,
        "pos": record.POS,
        "ref": record.REF,
        "alts": alts,
        "info": dict(record.INFO),
        "samples": samples,
        "format": list(record.FORMAT) if record.FORMAT else [],
    }


def _read_vcf_vcfpy(vcf_path: Path) -> tuple[list[str], list[dict]]:
    """Read VCF with vcfpy; return (sample_names, list of record dicts). Raises on parse/encoding error."""
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)  # vcfpy FieldInfoNotFound for non-standard INFO
        with open(vcf_path, encoding="utf-8", errors="replace") as stream:
            reader = Reader.from_stream(stream)
            sample_names = list(reader.header.samples.names)
            records = [_vcfpy_record_to_dict(rec) for rec in reader]
            reader.close()
    return sample_names, records


def _cyvcf2_record_to_dict(record, sample_names: list[str]) -> dict:
    """Convert a cyvcf2 Variant to our standard record dict (same shape as _vcfpy_record_to_dict)."""
    alts = [str(a) for a in (record.ALT or [])]
    info = dict(record.INFO)
    format_keys = list(record.FORMAT) if record.FORMAT else []
    samples = {}
    if sample_names and format_keys:
        raw = record.format(format_keys[0])
        n_s = len(raw) if raw is not None else 0
        for si, sname in enumerate(sample_names):
            if si >= n_s:
                break
            samples[sname] = {
                k: _format_value_for_display(record, k, si) for k in format_keys
            }
    return {
        "chrom": record.CHROM,
        "pos": record.POS,
        "ref": record.REF,
        "alts": alts,
        "info": info,
        "samples": samples,
        "format": format_keys,
    }


def _read_vcf_cyvcf2(vcf_path: Path) -> tuple[list[str], list[dict]]:
    """Read VCF with cyvcf2; return (sample_names, list of record dicts)."""
    vcf = VCF(str(vcf_path))
    sample_names = list(vcf.samples)
    records = [_cyvcf2_record_to_dict(rec, sample_names) for rec in vcf]
    vcf.close()
    return sample_names, records


def _read_vcf_records(vcf_path: Path) -> tuple[list[str], list[dict]]:
    """Read VCF file with built-in text parser (tolerates Ion Reporter FUNC / semicolons in INFO)."""
    sample_names = []
    records = []
    with open(vcf_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n\r")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM") or line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 10:
                    sample_names = parts[9:]
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            chrom, pos_s, ref, alt = parts[0], parts[1], parts[3], parts[4]
            try:
                pos = int(pos_s)
            except ValueError:
                continue
            info = _parse_info(parts[7])
            format_keys = parts[8].split(":") if len(parts) > 8 else []
            samples = {}
            for i, sname in enumerate(sample_names):
                if 9 + i < len(parts):
                    samples[sname] = _parse_format_sample(format_keys, parts[9 + i])
            records.append({
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alts": [a for a in alt.split(",") if a and a != "."],
                "info": info,
                "samples": samples,
                "format": format_keys,
            })
    return sample_names, records


def get_sample_dp(record, sample_name: str) -> str:
    """Get sample-specific DP from FORMAT; fallback to INFO DP. record is dict or object."""
    if hasattr(record, "samples"):
        s = record.samples.get(sample_name) if sample_name else None
        if s is not None and hasattr(record, "format") and "DP" in (record.format or []):
            try:
                dp = s.get("DP") if isinstance(s, dict) else (s["DP"] if "DP" in s else None)
                if dp is not None:
                    return str(dp)
            except (IndexError, KeyError, TypeError):
                pass
    info = record.info if hasattr(record, "info") else record.get("info", {})
    info_dp = info.get("DP") or info.get("FDP")
    if info_dp is not None:
        if isinstance(info_dp, (list, tuple)):
            info_dp = info_dp[0] if info_dp else None
        return str(info_dp) if info_dp is not None else ""
    return ""


def match_record(record, search: str) -> bool:
    """Return True if record matches search. record can be dict or object with .chrom, .pos, .info."""
    if isinstance(record, dict):
        return _match_record_dict(record, search)
    search = search.strip().upper()
    if ":" in search:
        chrom_part, pos_part = search.split(":", 1)
        chrom_part = chrom_part.strip()
        if not chrom_part.upper().startswith("CHR"):
            chrom_part = "chr" + chrom_part
        try:
            pos = int(pos_part.strip())
            if record.chrom == chrom_part and record.pos == pos:
                return True
        except ValueError:
            pass
    else:
        try:
            if record.pos == int(search):
                return True
        except (ValueError, TypeError):
            pass
    func_str = record.info.get("FUNC") if hasattr(record, "info") else None
    if func_str:
        for f in parse_func(str(func_str)):
            if isinstance(f, dict) and f.get("gene", "").upper() == search:
                return True
    exo_str = record.info.get("Exomiser") if hasattr(record, "info") else None
    if exo_str:
        for exo in parse_exomiser(str(exo_str)):
            if exo.get("GENE_SYMBOL", "").upper() == search:
                return True
    return False


def _record_info(record) -> dict:
    """Normalize record to .info-style access."""
    if isinstance(record, dict):
        return record.get("info", {})
    return getattr(record, "info", {})


def _record_samples(record) -> dict:
    if isinstance(record, dict):
        return record.get("samples", {})
    return getattr(record, "samples", {})


def _record_format(record) -> list:
    if isinstance(record, dict):
        return record.get("format", []) or []
    return list(getattr(record, "format", [])) or []


def build_report(record, sample_name: str | None = None) -> str:
    """Build the 7-step final variant message from one VCF record (dict or object)."""
    info = _record_info(record)
    func_str = info.get("FUNC") or ""
    funcs = parse_func(str(func_str))
    exo_str = info.get("Exomiser") or ""
    exo_list = parse_exomiser(str(exo_str))
    exo = exo_list[0] if exo_list else {}

    func = funcs[0] if (funcs and isinstance(funcs[0], dict)) else {}
    if not func and exo:
        hgvs = exo.get("HGVS", "")
        parts = hgvs.split(":") if hgvs else []
        coding = next((p for p in parts if p.strip().startswith("c.")), parts[-1] if parts else "")
        protein = next((p for p in parts if p.strip().startswith("p.")), "")
        func = {"gene": exo.get("GENE_SYMBOL", ""), "coding": coding, "protein": protein.strip("()"), "function": exo.get("FUNCTIONAL_CLASS", ""), "location": "exonic"}

    gene = func.get("gene") or exo.get("GENE_SYMBOL") or "Unknown"
    coding = func.get("coding") or ""
    protein = func.get("protein", "").strip("()") or ""
    function = func.get("function") or exo.get("FUNCTIONAL_CLASS") or "variant"
    location = func.get("location", "exonic")
    exon = func.get("exon", "")

    samples = _record_samples(record)
    samples_list = list(samples) if samples else []
    sname = sample_name or (samples_list[0] if samples_list else None)
    gt_raw = "."
    if sname and sname in samples:
        sdata = samples[sname]
        gt_raw = sdata.get("GT") if isinstance(sdata, dict) else "."
        if gt_raw and hasattr(gt_raw, "__len__") and not isinstance(gt_raw, str) and len(gt_raw) >= 2:
            gt_raw = f"{gt_raw[0]}/{gt_raw[1]}" if gt_raw[0] is not None and gt_raw[1] is not None else "."
        else:
            gt_raw = str(gt_raw) if gt_raw else "."
    zygosity_label, _ = parse_gt(gt_raw)

    chrom = record.get("chrom") if isinstance(record, dict) else getattr(record, "chrom", "")
    pos = record.get("pos") if isinstance(record, dict) else getattr(record, "pos", 0)
    ref = record.get("ref") if isinstance(record, dict) else getattr(record, "ref", "")
    alts = record.get("alts") if isinstance(record, dict) else (list(getattr(record, "alts", [])) or [])
    locus = f"{chrom}:{pos}"
    alt = ",".join(alts) if alts else "."
    depth = get_sample_dp(record, sname) if sname else (info.get("DP") or info.get("FDP") or "")
    if isinstance(depth, (list, tuple)) and depth:
        depth = depth[0]
    depth_str = f"{depth}x" if depth else ""

    # In-silico and ClinVar from FUNC/Exomiser
    polyphen = func.get("polyphen", "")
    sift = func.get("sift", "")
    clnacc = func.get("CLNACC1", "")
    clnid = func.get("CLNID1", "")
    acmg = exo.get("EXOMISER_ACMG_CLASSIFICATION", "")
    disease = (exo.get("EXOMISER_ACMG_DISEASE_NAME") or "").replace("_", " ")
    # Sanitize when value was split by strict INFO parsing (e.g. vcfpy)
    for bad in ['"', "'", "}]", "}", "']"]:
        disease = disease.split(bad)[0].strip()
    disease = disease.strip()

    # Human-readable variant effect
    effect_map = {
        "missense": "missense",
        "missense_variant": "missense",
        "frameshiftDeletion": "frameshift deletion",
        "frameshift_truncation": "frameshift truncation",
        "splice_region_variant": "splice region",
        "splice_donor_variant": "splice donor",
        "synonymous": "synonymous",
        "synonymous_variant": "synonymous",
        "intronic": "intronic",
        "exonic": "exonic",
    }
    effect_text = effect_map.get(function, function.replace("_", " ") if function else "variant")

    lines = []

    # 1. Variant sentence
    line1 = f"A {zygosity_label} {effect_text} variant {coding} in {gene} gene ({locus}; Depth:{depth_str}) was detected."
    lines.append("1.")
    lines.append(line1)
    lines.append("")

    # 2. Amino acid change (use residue properties from AutoVariant.mdc when available)
    if protein and protein != "p?" and "?" not in protein:
        lines.append("2.")
        lines.append(_format_amino_acid_change(protein))
    else:
        lines.append("2.")
        lines.append("No protein change reported or variant is non-coding.")
    lines.append("")

    # 3. Gene/protein impact
    lines.append("3.")
    if disease:
        lines.append(f"This variant is associated with {disease}. Consult ClinVar, Franklin Genoox, VarSome and literature (PubMed / population databases) for detailed evidence.")
    else:
        lines.append("Consult ClinVar, Franklin Genoox, VarSome and literature (PubMed / population databases) for how this variant affects the gene and protein.")
    lines.append("")

    # 4. Population frequencies (MyVariant.info when available)
    lines.append("4.")
    step4_msg = None
    if alts:
        try:
            step4_msg = _step4_population_message(chrom, pos, ref, alts[0])
        except Exception:
            pass
    if step4_msg:
        lines.append(step4_msg)
    else:
        lines.append("Population frequencies (gnomAD max, ExAC) are not available in this VCF. Check Franklin Genoox variant assessment or state: absent from major databases.")
    lines.append("")

    # 5. ClinVar
    lines.append("5.")
    if clnacc or clnid:
        lines.append(f"The variant has been reported in ClinVar (Variation ID: {clnacc}; {clnid}).")
    else:
        lines.append("The variant has not been reported in ClinVar database.")
    lines.append("")

    # 6. In-silico
    lines.append("6.")
    if polyphen or sift:
        pp = f"PolyPhen-2: {polyphen}" if polyphen else ""
        sf = f"SIFT: {sift}" if sift else ""
        preds = "; ".join(x for x in [pp, sf] if x)
        if (polyphen and float(polyphen) >= 0.85) or (sift and float(sift) <= 0.05):
            lines.append(f"In-silico prediction tools suggest a deleterious effect ({preds}).")
        else:
            lines.append(f"In-silico prediction tools: {preds or 'SIFT, REVEL, Polyphen-2, CADD, MVP — check Franklin Genoox for full assessment'}.")
    else:
        lines.append("In-silico data (SIFT, REVEL, MT, Polyphen-2, CADD, MVP) not present in this VCF; check Franklin Genoox variant assessment.")
    lines.append("")

    # 7. ACMG
    lines.append("7.")
    acmg_upper = (acmg or "").upper()
    if "PATHOGENIC" in acmg_upper and "LIKELY" not in acmg_upper:
        lines.append("Based on the aforementioned evidence, this variant is classified as pathogenic according to the ACMG AMP guidelines.")
    elif "LIKELY_PATHOGENIC" in acmg_upper:
        lines.append("Based on the aforementioned evidence, this variant is classified as likely pathogenic according to the ACMG AMP guidelines.")
    elif "UNCERTAIN" in acmg_upper or "VUS" in acmg_upper or "SIGNIFICANCE" in acmg_upper:
        lines.append("Based on the aforementioned evidence, this variant is classified as variant of uncertain significance according to the ACMG AMP guidelines.")
    elif "BENIGN" in acmg_upper:
        lines.append("Based on the aforementioned evidence, this variant is classified as benign / likely benign according to the ACMG AMP guidelines.")
    else:
        lines.append(f"ACMG classification from VCF: {acmg or 'Not available'}.")
    lines.append("")

    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description="AutoVariant: search VCF by gene or location and print final variant message (SOP 7 steps).")
    parser.add_argument("search", help="Gene symbol (e.g. FAH) or location (e.g. chr15:80472479 or 80472479)")
    parser.add_argument("--vcf", "-v", default="WES166_26032846_BSC.vcf", help="Path to VCF file")
    parser.add_argument("--sample", "-s", default=None, help="Sample name (default: first sample)")
    args = parser.parse_args()

    vcf_path = Path(args.vcf)
    if not vcf_path.is_file():
        print(f"Error: VCF file not found: {vcf_path}", file=sys.stderr)
        return 1

    sample_names, records = [], []
    # Prefer cyvcf2; fall back to text parser then vcfpy
    if CYVCF2_AVAILABLE:
        try:
            sample_names, records = _read_vcf_cyvcf2(vcf_path)
        except Exception:
            pass
    if not records:
        try:
            sample_names, records = _read_vcf_records(vcf_path)
        except Exception:
            pass
    if not records and VCFPY_AVAILABLE:
        try:
            sample_names, records = _read_vcf_vcfpy(vcf_path)
        except Exception:
            pass
    if not records:
        print("Error: no VCF records could be read.", file=sys.stderr)
        return 1

    found = None
    for record in records:
        if match_record(record, args.search):
            found = record
            break

    if found is None:
        print(f"No variant found for search: {args.search}", file=sys.stderr)
        return 1

    report = build_report(found, args.sample or (sample_names[0] if sample_names else None))
    print(report)
    return 0


if __name__ == "__main__":
    sys.exit(main())
