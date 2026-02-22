# Variant at chr15:80472479 (FAH) — INFO and FORMAT fields

Single-nucleotide variant **C→T** at position 80472479 on chromosome 15 (hg19). This is the FAH missense variant **c.974C>T**, **p.Thr325Met**, associated with Tyrosinemia type I.

---

## Fixed columns (per VCF spec)

| Column | Value | Meaning |
|--------|--------|--------|
| CHROM | chr15 | Chromosome |
| POS | 80472479 | Position (1-based, hg19) |
| ID | . | No dbSNP/ID in this column (rs770713168 is in FUNC/ClinVar) |
| REF | C | Reference allele |
| ALT | T | Alternate allele |
| QUAL | 1603.11 | Phred-scaled quality of the variant call |
| FILTER | PASS | Variant passed all filters |

---

## INFO fields (variant-level)

| INFO key | Value | What it tells you |
|----------|--------|-------------------|
| **AF** | 0.988024 | **Allele frequency** in this sample: ~98.8% of reads show T (alternate). Consistent with homozygous alternate (1/1). |
| **AO** | 162 | **Alternate allele observations**: 162 reads support T. |
| **DP** | 164 | **Total read depth** at the locus (164 reads). |
| **Exomiser** | See [Exomiser field (expanded)](#exomiser-field-expanded) below. | Phenotype-driven variant prioritization; one pipe-separated block per compatible mode of inheritance (MOI). |
| **FAO** | 165 | **Flow Evaluator** alternate observations (Ion Torrent flow-based count). |
| **FDP** | 167 | Flow Evaluator total depth at locus. |
| **FDVR** | -1 | Flow disruption of alt vs reference (negative = less disruption). |
| **FR** | . | Filter reason (empty = not filtered). |
| **FRO** | 2 | Flow Evaluator **reference** observations (only 2 ref reads). |
| **FSAF** | 123 | Alternate observations on **forward** strand (flow). |
| **FSAR** | 42 | Alternate observations on **reverse** strand (flow). |
| **FSRF** | 2 | Reference observations on forward strand. |
| **FSRR** | 0 | Reference observations on reverse strand. |
| **FUNC** | [{'gene':'FAH','coding':'c.974C>T',…}] | **Functional annotation** (Ion Reporter): gene **FAH**, **c.974C>T**, **p.Thr325Met**, transcript NM_000137.4, **missense**, exon 12, **PolyPhen 1.0** (probably damaging), **SIFT 0.0** (deleterious), **Grantham 81**, **ClinVar** accession 552826, **rs770713168**, clinical significance **Conflicting interpretations of pathogenicity**, review status **criteria provided, conflicting interpretations**. |
| **FWDB** | -0.00446849 | Forward-strand bias in prediction (small = minimal bias). |
| **FXX** | 0 | Flow Evaluator failed-read ratio. |
| **GCM** | 0 | Not a GC-motif artifact. |
| **HRUN** | 1 | Homopolymer run length for alt allele. |
| **HS_ONLY** | 0 | Not hotspot-only. |
| **LEN** | 1 | Allele length change (1 bp substitution). |
| **MLLD** | 169.709 | **Mean log-likelihood delta** per read (support for variant). |
| **OALT** | T | Original alternate base (same as ALT). |
| **OID** | . | Original hotspot ID (none). |
| **OMAPALT** | T | Mapping of original alt to ALT. |
| **OPOS** | 80472479 | Original position. |
| **OREF** | C | Original reference base. |
| **PB** | 0.5 | **Position bias**: ref vs alt read position (0.5 = neutral). |
| **PBP** | 1 | P-value for position bias (1 = no evidence of bias). |
| **PPD** | 0 | No prefix padding on alt. |
| **QD** | 38.3979 | **Quality by depth** (QUAL/FDP-like): high confidence. |
| **RBI** | 0.0695339 | Distance of bias parameters from zero. |
| **REFB** | -0.162362 | Reference-hypothesis bias. |
| **REVB** | 0.0693902 | Reverse-strand bias. |
| **RO** | 2 | **Reference allele observations**: only 2 ref reads (rest are T). |
| **SAF** | 123 | Alternate observations on **forward** strand (standard). |
| **SAR** | 39 | Alternate observations on **reverse** strand. |
| **SPD** | 0 | No suffix padding. |
| **SRF** | 2 | Reference observations on forward strand. |
| **SRR** | 0 | Reference observations on reverse strand. |
| **SSEN** | 0 | Strand-specific error (negative strand). |
| **SSEP** | 0 | Strand-specific error (positive strand). |
| **SSSB** | -0.00767438 | **Strand-specific strand bias** for allele (small = minimal). |
| **STB** | 0.504032 | **Strand bias** (variant vs ref); ~0.5 = balanced. |
| **STBP** | 0.434 | P-value for strand bias. |
| **TYPE** | snp | Variant type: **single-nucleotide polymorphism**. |
| **VARB** | 9.49604E-4 | Variant-hypothesis bias. |

---

## Exomiser field (expanded)

The VCF header defines the **Exomiser** INFO field as:

> *"A pipe-separated set of values for the proband allele(s) from the record with one per compatible MOI following the format: {RANK|ID|GENE_SYMBOL|ENTREZ_GENE_ID|MOI|P-VALUE|EXOMISER_GENE_COMBINED_SCORE|EXOMISER_GENE_PHENO_SCORE|EXOMISER_GENE_VARIANT_SCORE|EXOMISER_VARIANT_SCORE|CONTRIBUTING_VARIANT|WHITELIST_VARIANT|FUNCTIONAL_CLASS|HGVS|EXOMISER_ACMG_CLASSIFICATION|EXOMISER_ACMG_EVIDENCE|EXOMISER_ACMG_DISEASE_ID|EXOMISER_ACMG_DISEASE_NAME}"*

For this variant there is **one** block (one compatible MOI: AR). Raw value:

`{1|15-80472479-C-T_AR|FAH|2184|AR|0.0001|0.9824|0.7831|0.9985|0.9985|1|0|missense_variant|FAH:ENST00000261755.5:c.974C>T:p.(Thr325Met)|LIKELY_PATHOGENIC|PM2_Supporting,PP3_Strong,PP4_Moderate,BP1|OMIM:276700|"Tyrosinemia,_type_I"}`

### Field-by-field interpretation (per header order)

| # | Field name (from VCF header) | Value | What it tells you |
|---|------------------------------|--------|-------------------|
| 1 | **RANK** | 1 | **Prioritization rank**: this variant is **#1** in the Exomiser result for this sample under the given MOI — the top candidate given the phenotype. |
| 2 | **ID** | 15-80472479-C-T_AR | **Variant identifier** for this record: chrom 15, position 80472479, ref C, alt T, with MOI suffix **_AR** (autosomal recessive). |
| 3 | **GENE_SYMBOL** | FAH | **Gene symbol**: **FAH** (fumarylacetoacetate hydrolase). |
| 4 | **ENTREZ_GENE_ID** | 2184 | **Entrez Gene ID** for FAH (NCBI). |
| 5 | **MOI** | AR | **Mode of inheritance** considered: **AR** = autosomal recessive. Exomiser scored this variant under the AR model. |
| 6 | **P-VALUE** | 0.0001 | **P-value** for the gene–phenotype association under this MOI; very low (0.0001) — strong statistical support that this gene fits the phenotype. |
| 7 | **EXOMISER_GENE_COMBINED_SCORE** | 0.9824 | **Gene combined score** (0–1): combines how well the gene matches the phenotype and variant data; 0.98 = very strong gene-level support. |
| 8 | **EXOMISER_GENE_PHENO_SCORE** | 0.7831 | **Gene phenotype score** (0–1): how well the **gene** matches the patient phenotype; 0.78 = good match. |
| 9 | **EXOMISER_GENE_VARIANT_SCORE** | 0.9985 | **Gene variant score** (0–1): support from **variant data** (frequency, predicted impact, etc.) for this gene; 0.9985 = very high. |
| 10 | **EXOMISER_VARIANT_SCORE** | 0.9985 | **Variant score** (0–1): overall support for **this specific variant**; 0.9985 = very strong variant-level support. |
| 11 | **CONTRIBUTING_VARIANT** | 1 | **Contributing variant flag**: 1 = this variant **contributes** to the gene’s prioritization (used in the analysis). |
| 12 | **WHITELIST_VARIANT** | 0 | **Whitelist flag**: 0 = not on a user whitelist; prioritization is from Exomiser’s own scoring. |
| 13 | **FUNCTIONAL_CLASS** | missense_variant | **Functional class**: **missense_variant** — single amino acid substitution (Sequence Ontology / VEP style). |
| 14 | **HGVS** | FAH:ENST00000261755.5:c.974C>T:p.(Thr325Met) | **HGVS** (Human Genome Variation Society): gene, transcript (ENST00000261755.5), **c.974C>T** (coding), **p.(Thr325Met)** (protein). Unambiguous nomenclature for the variant. |
| 15 | **EXOMISER_ACMG_CLASSIFICATION** | LIKELY_PATHOGENIC | **ACMG classification** from Exomiser: **LIKELY_PATHOGENIC** — variant is classified as likely pathogenic per ACMG/AMP guidelines. |
| 16 | **EXOMISER_ACMG_EVIDENCE** | PM2_Supporting,PP3_Strong,PP4_Moderate,BP1 | **ACMG evidence codes** applied: **PM2** (Supporting) = absent/infrequent in population; **PP3** (Strong) = multiple lines of computational support for deleterious; **PP4** (Moderate) = phenotype highly specific for single-gene disorder; **BP1** (Benign supporting) = gene associated only with dominant disease (here, BP1 may reduce weight in some workflows). Net result supports likely pathogenic. |
| 17 | **EXOMISER_ACMG_DISEASE_ID** | OMIM:276700 | **Disease identifier**: **OMIM:276700** — Tyrosinemia type I (FAH-related) in OMIM. |
| 18 | **EXOMISER_ACMG_DISEASE_NAME** | "Tyrosinemia,_type_I" | **Disease name**: **Tyrosinemia, type I** (underscores in VCF represent spaces). The phenotype/disease used for prioritization and ACMG context. |

### Summary of what the Exomiser field holds

- **Prioritization**: Rank 1, under AR; gene and variant scores ~0.98–0.99; p-value 0.0001 — Exomiser considers this the top hit.
- **Variant identity**: FAH, c.974C>T, p.Thr325Met (missense), with full HGVS and transcript.
- **ACMG**: LIKELY_PATHOGENIC with evidence PM2, PP3, PP4, BP1.
- **Disease**: OMIM:276700, Tyrosinemia type I — the disorder this variant is evaluated against.

---

## FORMAT and sample column (WES166_26032846_BSC)

Format string: `GT:AF:AO:DP:FAO:FDP:FRO:FSAF:FSAR:FSRF:FSRR:GQ:RO:SAF:SAR:SRF:SRR`

| FORMAT | Value | What it tells you |
|--------|--------|-------------------|
| **GT** | 1/1 | **Genotype**: homozygous alternate (T/T). Both chromosomes carry the variant. |
| **AF** | 0.988024 | Sample **allele frequency** for alt: ~98.8% of reads are T. |
| **AO** | 162 | Alternate allele observation count (162 reads). |
| **DP** | 164 | **Depth** at this site for this sample (164 reads). |
| **FAO** | 165 | Flow Evaluator alternate count. |
| **FDP** | 167 | Flow Evaluator depth. |
| **FRO** | 2 | Flow Evaluator reference count. |
| **FSAF** | 123 | Flow: alt on forward strand. |
| **FSAR** | 42 | Flow: alt on reverse strand. |
| **FSRF** | 2 | Flow: ref on forward. |
| **FSRR** | 0 | Flow: ref on reverse. |
| **GQ** | 53 | **Genotype quality** (Phred): ~99.95% confidence in 1/1. |
| **RO** | 2 | Reference allele observations. |
| **SAF** | 123 | Alt on forward strand. |
| **SAR** | 39 | Alt on reverse strand. |
| **SRF** | 2 | Ref on forward. |
| **SRR** | 0 | Ref on reverse. |

---

## Summary in plain language

- **Variant**: C→T at chr15:80472479 (hg19), in **FAH**, **c.974C>T**, **p.Thr325Met** (missense).
- **Genotype**: **Homozygous** (1/1): both copies carry T; ~162 alt reads vs 2 ref, depth 164.
- **Quality**: High QUAL, PASS, good QD and GQ — call is technically reliable.
- **Interpretation**: Exomiser ranks it #1, **LIKELY_PATHOGENIC**, with PM2/PP3/PP4/BP1; associated with **Tyrosinemia type I** (OMIM:276700). ClinVar lists **conflicting interpretations**; rs770713168.
- **In silico**: PolyPhen 1.0 (probably damaging), SIFT 0.0 (deleterious), Grantham 81 — all support deleterious missense.
- **Strand/position**: No strong strand or position bias; flow and standard counts agree.
