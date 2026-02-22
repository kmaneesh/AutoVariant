"""Verify we can read both .vcf and .vcf.gz with cyvcf2 (cyvcf2 supports both natively)."""
import gzip
import shutil
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent

# Minimal valid VCF that cyvcf2 can parse (simple INFO, no semicolons in values)
# One variant at chr1 1000, gene FAH in FUNC for run_variant_search
MINIMAL_VCF = """##fileformat=VCFv4.1
##contig=<ID=chr1,length=249250621>
##INFO=<ID=FUNC,Number=.,Type=String,Description="Functional">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	1000	.	A	G	30	PASS	FUNC=[{"gene":"FAH","protein":"p.Thr325Met"}]	GT:DP	0/1:50
"""


@pytest.fixture
def vcf_and_gz(tmp_path):
    """Write minimal .vcf and .vcf.gz; return (path_vcf, path_gz)."""
    vcf = tmp_path / "test.vcf"
    vcf.write_text(MINIMAL_VCF, encoding="utf-8")
    gz = tmp_path / "test.vcf.gz"
    with open(vcf, "rb") as f_in:
        with gzip.open(gz, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    return vcf, gz


def test_read_vcf_with_cyvcf2(vcf_and_gz):
    """_read_vcf_cyvcf2 reads plain .vcf and returns records."""
    from scripts.auto_variant import _read_vcf_cyvcf2

    vcf_path, _ = vcf_and_gz
    sample_names, records = _read_vcf_cyvcf2(vcf_path)
    assert len(sample_names) >= 1
    assert len(records) == 1
    assert records[0]["chrom"] == "chr1"
    assert records[0]["pos"] == 1000


def test_read_vcf_gz_with_cyvcf2(vcf_and_gz):
    """_read_vcf_cyvcf2 reads .vcf.gz and returns the same records as .vcf (cyvcf2 supports both)."""
    from scripts.auto_variant import _read_vcf_cyvcf2

    vcf_path, gz_path = vcf_and_gz
    _, records_vcf = _read_vcf_cyvcf2(vcf_path)
    _, records_gz = _read_vcf_cyvcf2(gz_path)
    assert len(records_vcf) == len(records_gz)
    assert records_vcf[0]["chrom"] == records_gz[0]["chrom"]
    assert records_vcf[0]["pos"] == records_gz[0]["pos"]


def test_run_variant_search_works_with_vcf_and_vcf_gz(vcf_and_gz):
    """run_variant_search works for both .vcf and .vcf.gz (same query, same report)."""
    from scripts.auto_variant import run_variant_search

    vcf_path, gz_path = vcf_and_gz
    report_vcf = run_variant_search(vcf_path, "FAH")
    report_gz = run_variant_search(gz_path, "FAH")
    assert "FAH" in report_vcf and "FAH" in report_gz
    assert "1." in report_vcf and "1." in report_gz
