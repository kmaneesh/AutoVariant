"""Tests for AutoVariant report: search by location chr15:80472479 (FAH variant)."""
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
VCF = ROOT / "data" / "WES166_26032846_BSC_IN.vcf"


def _run_report(search: str, vcf: Path = VCF) -> tuple[str, str, int]:
    cmd = [sys.executable, "-m", "app.auto_variant", search, "--vcf", str(vcf)]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=str(ROOT))
    return r.stdout, r.stderr, r.returncode


@pytest.mark.skipif(not VCF.exists(), reason="Sample VCF not found")
def test_search_by_location_80472479_finds_variant():
    """Searching by chr15:80472479 returns exit 0 and a report."""
    stdout, stderr, code = _run_report("chr15:80472479")
    assert code == 0, f"Expected exit 0, got {code}. stderr: {stderr}"
    assert "80472479" in stdout
    assert "1." in stdout


@pytest.mark.skipif(not VCF.exists(), reason="Sample VCF not found")
def test_search_by_position_only_80472479_finds_same_variant():
    """Searching by position 80472479 (no chrom) finds the variant at chr15:80472479."""
    stdout, stderr, code = _run_report("80472479")
    assert code == 0
    assert "80472479" in stdout
    assert "FAH" in stdout
    assert "c.974C>T" in stdout


@pytest.mark.skipif(not VCF.exists(), reason="Sample VCF not found")
def test_location_80472479_report_is_fah():
    """Variant at chr15:80472479 is the FAH missense c.974C>T, p.Thr325Met."""
    stdout, stderr, code = _run_report("chr15:80472479")
    assert code == 0
    assert "FAH" in stdout
    assert "c.974C>T" in stdout
    assert "Thr325Met" in stdout or "p.Thr325Met" in stdout
    assert "missense" in stdout.lower()
    assert "164" in stdout  # depth
    assert "homozygous" in stdout.lower()
