"""Shared types for lookup results."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class PopulationFreq:
    """Population frequency summary (gnomAD max, ExAC, etc.)."""
    gnomad_af_max: float | None = None
    gnomad_ac: int | None = None
    gnomad_an: int | None = None
    exac_af: float | None = None
    exac_ac: int | None = None
    exac_an: int | None = None
    raw_gnomad: dict | None = None
    raw_exac: dict | None = None

    def summary(self) -> str:
        parts = []
        if self.gnomad_af_max is not None:
            parts.append(f"gnomAD max AF: {self.gnomad_af_max:.2e}")
        if self.exac_af is not None:
            parts.append(f"ExAC AF: {self.exac_af:.2e}")
        if not parts:
            return "No population frequency data from gnomAD/ExAC."
        return "; ".join(parts)

    def _af_percent(self, af: float) -> str:
        """Format allele frequency as percentage (e.g. 0.0001 -> 0.01%)."""
        pct = af * 100
        if pct >= 0.01:
            return f"{pct:.2f}%"
        if pct >= 0.0001:
            return f"{pct:.4f}%"
        return f"{pct:.4f}%" if pct != 0 else "0%"

    def step4_message(self) -> str:
        """Sentence for report Step 4: minor allele frequency in gnomAD (max) and ExAC."""
        g = self.gnomad_af_max is not None
        e = self.exac_af is not None
        if g and e:
            return (
                f"This variant has minor allele frequency of {self._af_percent(self.gnomad_af_max)} "
                f"and {self._af_percent(self.exac_af)} in gnomAD (max) and ExAC database respectively."
            )
        if g:
            return (
                f"This variant has minor allele frequency of {self._af_percent(self.gnomad_af_max)} "
                "in gnomAD (max)."
            )
        if e:
            return (
                f"This variant has minor allele frequency of {self._af_percent(self.exac_af)} "
                "in ExAC database."
            )
        return "No population frequency data from gnomAD/ExAC."


@dataclass
class LookupResult:
    """Aggregated result from all lookup sources."""
    gnomad_exac: PopulationFreq | None = None
    varsome: dict | str | None = None  # JSON or error message
    franklin: dict | str | None = None  # JSON or error message
