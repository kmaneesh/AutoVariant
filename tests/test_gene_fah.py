"""Tests for AutoVariant report: search by gene FAH."""
import subprocess
import sys
from pathlib import Path

import pytest

# Project root
ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "scripts" / "auto_variant.py"
VCF = ROOT / "WES166_26032846_BSC.vcf"


def _run_report(search: str, vcf: Path = VCF) -> tuple[str, str, int]:
    """Run auto_variant.py; return (stdout, stderr, returncode)."""
    cmd = [sys.executable, str(SCRIPT), search, "--vcf", str(vcf)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    return r.stdout, r.stderr, r.returncode


@pytest.mark.skipif(not VCF.exists(), reason="Sample VCF not found")
def test_search_by_gene_fah_finds_variant():
    """Searching by gene FAH returns exit 0 and a report."""
    stdout, stderr, code = _run_report("FAH")
    assert code == 0, f"Expected exit 0, got {code}. stderr: {stderr}"
    assert "FAH" in stdout
    assert "1." in stdout and "7." in stdout


@pytest.mark.skipif(not VCF.exists(), reason="Sample VCF not found")
def test_fah_report_contains_expected_content():
    """FAH variant report contains locus, coding change, protein, zygosity, depth."""
    stdout, stderr, code = _run_report("FAH")
    assert code == 0
    # Locus chr15:80472479
    assert "chr15:80472479" in stdout or "80472479" in stdout
    # Coding and protein (from FUNC)
    assert "c.974C>T" in stdout
    assert "Thr325Met" in stdout or "p.Thr325Met" in stdout
    # Zygosity: this sample is 1/1 homozygous
    assert "homozygous" in stdout.lower()
    # Depth
    assert "164" in stdout or "Depth:" in stdout
    # ACMG
    assert "likely pathogenic" in stdout.lower() or "pathogenic" in stdout.lower() or "ACMG" in stdout


@pytest.mark.skipif(not VCF.exists(), reason="Sample VCF not found")
def test_fah_final_message_has_seven_steps():
    """Final variant message includes all 7 SOP steps."""
    stdout, stderr, code = _run_report("FAH")
    assert code == 0
    for i in range(1, 8):
        assert f"{i}." in stdout, f"Missing step {i} in output"
