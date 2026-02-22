# AutoVariant — Project Plan

## Objective

Generate **variant information reports** by pulling data from a VCF file. The program takes a **single input (search)**—either a **gene symbol** or a **genomic location**—and produces a final narrative message that follows the WES Variant Details SOP.

---

## Inputs & Data Sources

| Item | Description |
|------|-------------|
| **User input** | One search term: **gene name** (e.g. `COL4A5`, `CLCN6`) or **location** (e.g. `chrX:107898607`, `chr1:11888529` or range). |
| **VCF file** | Primary data source. Read using **vembrane** (Python/VCF tooling). |
| **Sample VCF** | `WES166_26032846_BSC.vcf` — Ion Reporter–style VCF (hg19) with INFO fields `FUNC`, `Exomiser`, `DP`, etc. |

---

## Output: Final Message Format

The final message must follow the **SOP for WES Variant Details** described in `AutoVariant-example.md`:

1. **Step 1 — Variant sentence**  
   Zygosity, variant effect, coding change, gene name, locus, and coverage.  
   *e.g. “A hemizygous missense variant c.3293G>A in COL4A5 gene (chrX:107898607; Depth:20x) was detected.”*

2. **Step 2 — Amino acid change**  
   Protein change and residue properties (e.g. neutral vs acidic/polar).  
   *e.g. “This variant (p.Gly1098Asp) replaces glycine, which is neutral and non-polar, with aspartic acid, which is acidic and polar, at 1098th amino acid position.”*

3. **Step 3 — Gene/protein impact & evidence**  
   How the variant affects the gene/protein; evidence or literature (PubMed, population DBs) if available.

4. **Step 4 — Population frequencies**  
   gnomAD (max) and ExAC (or statement if absent from major DBs).

5. **Step 5 — ClinVar**  
   Variation ID if reported, or “not reported in ClinVar”.

6. **Step 6 — In-silico predictions**  
   SIFT, REVEL, MT, Polyphen-2, CADD, MVP — deleterious vs uncertain.

7. **Step 7 — ACMG classification**  
   Pathogenic / Likely pathogenic / Variant of uncertain significance (per Exomiser/ACMG-AMP).

---

## VCF Data Mapping (from `WES166_26032846_BSC.vcf`)

Relevant fields for building the report:

| SOP need | VCF source |
|----------|------------|
| Gene name | `INFO.FUNC` (list of dicts) → `gene`; or `INFO.Exomiser` (pipe-separated) → GENE_SYMBOL |
| Locus | `CHROM`, `POS` → e.g. `chr1:11888529` |
| Depth / coverage | `INFO.DP` or FORMAT `DP` |
| Zygosity | FORMAT `GT` (0/0, 0/1, 1/1, 1/2, etc.) → hemizygous/heterozygous/homozygous |
| Variant effect | `FUNC[].function` (e.g. missense, frameshiftDeletion, splice_region_variant); or Exomiser FUNCTIONAL_CLASS |
| Coding change | `FUNC[].coding` (e.g. c.975delA) |
| Protein change | `FUNC[].protein` (e.g. p.Lys325AsnfsTer15) |
| Transcript | `FUNC[].transcript` |
| In-silico (partial) | `FUNC[].polyphen`, `FUNC[].sift`, `FUNC[].grantham` when present |
| ACMG classification | `INFO.Exomiser` → EXOMISER_ACMG_CLASSIFICATION (e.g. LIKELY_PATHOGENIC, UNCERTAIN_SIGNIFICANCE) |
| ClinVar (when present) | `FUNC[]` keys: `CLNACC1` (accession), `CLNID1` (rs/ID), `CLNSIG1`, `CLNREVSTAT1` |

Population frequencies (gnomAD, ExAC) and detailed literature are **not** in this VCF; they can be placeholders or future external API/DB lookups.

---

## High-Level Architecture

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────────┐
│  User search    │────▶│  VCF read        │────▶│  Match by gene or    │
│  (gene/location)│     │  (vembrane)      │     │  location (filter)   │
└─────────────────┘     └──────────────────┘     └──────────┬──────────┘
                                                             │
                                                             ▼
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────────┐
│  Final message  │◀────│  Report builder   │◀────│  Extract fields from │
│  (SOP steps 1–7)│     │  (SOP template)   │     │  INFO/FUNC/Exomiser   │
└─────────────────┘     └──────────────────┘     └─────────────────────┘
```

---

## Implementation Phases

### Phase 1: VCF read and search

- Use **vembrane** (or its backend: pysam/cyvcf2) to read the VCF.
- Implement **search**:
  - **By gene:** filter records where any annotation (e.g. `FUNC` or `Exomiser`) matches the given gene symbol.
  - **By location:** filter by `CHROM` and `POS` (exact or range, e.g. `chr1:11888529` or `chr1:11888000-11889000`).
- Support the sample file `WES166_26032846_BSC.vcf` (Ion Reporter INFO structure).

### Phase 2: Parsing INFO (FUNC, Exomiser)

- Parse **INFO.FUNC**: list of dicts with `gene`, `coding`, `protein`, `function`, `location`, `polyphen`, `sift`, `grantham`, `exon`, and optional `CLN*` keys.
- Parse **INFO.Exomiser**: pipe-separated segments; extract GENE_SYMBOL, HGVS, EXOMISER_ACMG_CLASSIFICATION, EXOMISER_ACMG_EVIDENCE, EXOMISER_ACMG_DISEASE_*.
- Derive **zygosity** from FORMAT `GT` (and sample sex if needed for “hemizygous” on chrX).
- Map **variant effect** and **ACMG** to the wording expected in the SOP (Step 1 and Step 7).

### Phase 3: Report builder (final message)

- Implement a **report builder** that outputs text for each of the 7 SOP steps.
- **Steps 1–2:** Use only VCF-derived data (variant sentence + amino acid change).
- **Steps 3–6:** Use VCF where available (e.g. ClinVar IDs and in-silico from FUNC); add placeholders or “not available” for missing data (e.g. population frequencies, literature).
- **Step 7:** Use Exomiser ACMG classification from VCF to choose Pathogenic / Likely pathogenic / VUS wording.

### Phase 4 (optional): External data and polish

- Add optional lookups for population frequencies (gnomAD, ExAC) and ClinVar Variation ID if not in VCF.
- Add optional literature/PubMed references.
- Improve amino acid property text (e.g. neutral/polar vs acidic) from protein change when possible.

---

## Deliverables

| Deliverable | Description |
|-------------|-------------|
| **CLI or script** | Single entry point: input = search (gene or location); output = final message (e.g. to stdout or file). |
| **VCF reader** | VCF opened and iterated via vembrane (or agreed backend); support for Ion Reporter–style INFO. |
| **Search** | Filter variants by gene symbol or by genomic location. |
| **Report generator** | Assemble the 7-step narrative from extracted fields and write the final message. |
| **Tests** | At least one run on `WES166_26032846_BSC.vcf` for a known gene/location (e.g. CLCN6, chr1:11888529) to validate end-to-end output. |

---

## Reference Files

- **Output format (SOP):** `AutoVariant-example.md`
- **Sample VCF:** `WES166_26032846_BSC.vcf`
- **Conventions / tooling:** AutoHPO repo (project structure, patterns); `.cursor/rules/` (venv, vembrane, vcf).

---

## Summary

| Aspect | Choice |
|--------|--------|
| **Input** | One search: **gene** or **location** |
| **Data** | VCF read with **vembrane** |
| **Output** | Final message following **AutoVariant-example.md** (7 steps) |
| **Sample** | `WES166_26032846_BSC.vcf` for gene/location-based extraction |
