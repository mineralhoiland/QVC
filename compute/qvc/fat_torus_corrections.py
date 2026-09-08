"""Fat-torus regime: locked F(x)=Σ K₀(nx) and useful expansions.

Regime
------
x = 2π r / R. Preprint (r,R)=(8,50) nm ⇒ x ≈ 1.005 → **fat torus**
(neither thin-ring PFA nor thick-tube limit).

Locked spectral weight
----------------------
    F(x) = Σ_{n=1}^∞ K₀(n x)

Asymptotics (derived, numerically verified in suite):
  • Large-x (thick / strong decay): F(x) ∼ √(π/(2x)) e^{-x}  (n=1 dominates)
  • Small-x (thin torus): F(x) ∼ (1/x)(γ_E + ln(2/x))/2 + … (harmonic-like)
  • Fat torus x∼O(1): no closed elementary reduction — use direct sum.

PFA is not a lock here; use F-sum + Tier-C form only.
"""

from __future__ import annotations

import math
from typing import Any

from scipy.special import k0

from qvc.params import PARAMS
from qvc.vacuum_freeze import bessel_mode_sum, vacuum_self_energy_corrected


def F_large_x_leading(x: float) -> float:
    """Leading n=1 Macdonald asymptote: √(π/(2x)) e^{-x}."""
    if x <= 0:
        return float("nan")
    return math.sqrt(math.pi / (2.0 * x)) * math.exp(-x)


def F_small_x_proxy(x: float, n_cut: int = 2000) -> float:
    """Thin-torus numerical proxy (truncated sum); not a closed form."""
    return bessel_mode_sum(x, nmax=n_cut)


def fat_torus_F_diagnostics(x: float | None = None) -> dict[str, Any]:
    """Compare exact F to large-x asymptote at fat-torus x."""
    p = PARAMS
    if x is None:
        x = 2.0 * math.pi * p.r_preprint_um / p.R_preprint_um
    F = bessel_mode_sum(x)
    F1 = float(k0(x))
    F_asym = F_large_x_leading(x)
    return {
        "x": x,
        "regime": "fat_torus",
        "F_exact": F,
        "F_n1_only": F1,
        "F_large_x_asym": F_asym,
        "rel_err_n1_vs_exact": abs(F1 - F) / F if F else float("nan"),
        "rel_err_asym_vs_n1": abs(F_asym - F1) / F1 if F1 else float("nan"),
        "notes": (
            "At x∼1, n=1 carries most weight but higher n still O(10%). "
            "Do not replace F by PFA or Jacobi θ₃."
        ),
    }


def tierC_from_geometry(
    r_um: float,
    R_um: float,
    hbar_cs: float | None = None,
) -> dict[str, float]:
    """E_vac = (1/4π²)(ℏc_s/r) F(2π r/R) — locked Tier-C evaluation."""
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    x = 2.0 * math.pi * r_um / R_um
    F = bessel_mode_sum(x)
    E = vacuum_self_energy_corrected(hbar_cs, r_um, F)
    return {"x": x, "F": F, "E_tierC_meV": E, "hbar_cs": hbar_cs, "r_um": r_um, "R_um": R_um}


def scan_fat_torus_aspect(
    *,
    R_um: float | None = None,
    hbar_cs: float | None = None,
    x_values: list[float] | None = None,
) -> list[dict[str, float]]:
    """Sweep aspect ratio at fixed major radius — fat-torus correction curve."""
    R_um = PARAMS.R_preprint_um if R_um is None else R_um
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    if x_values is None:
        x_values = [0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
    rows = []
    for x in x_values:
        r = x * R_um / (2.0 * math.pi)
        row = tierC_from_geometry(r, R_um, hbar_cs)
        diag = fat_torus_F_diagnostics(x)
        rows.append(
            {
                **row,
                "F_n1_only": diag["F_n1_only"],
                "F_large_x_asym": diag["F_large_x_asym"],
            }
        )
    return rows
