"""Vacuum freeze computational suite: fat-torus analysis + inverse freeze.

Derives and numerically evaluates:
  • F(x)=Σ K₀(nx) convergence
  • Tier-C E_vac vs old Casimir on fat-torus / moiré geometries
  • PFA vs full Bessel (fat-torus correction factor)
  • Inverse freeze: which ℏc_s (or r) would hit a target meV under Tier-C
  • SI ↔ package unit cross-check

Does NOT promote a lock to frozen — that remains an explicit decision.
"""

from __future__ import annotations

import math
from typing import Any

from qvc.params import PARAMS
from qvc.vacuum_freeze import (
    FrozenVacuumSpec,
    assert_preprint_0p244_not_reproduced,
    bessel_mode_sum,
    casimir_preprint_formula,
    evaluate_vacuum,
    jacobi_theta3_alternating,
    vacuum_self_energy_corrected,
)


def F_partial_sums(x: float, nmax: int = 80) -> list[dict[str, float]]:
    """Convergence table for F(x)=Σ_{n=1}^{N} K₀(nx)."""
    from scipy.special import k0

    s = 0.0
    rows = []
    for n in range(1, nmax + 1):
        s += float(k0(n * x))
        if n in {1, 2, 5, 10, 20, 50, 80} or n == nmax:
            rows.append({"n_max": float(n), "F_partial": s})
    return rows


def pfa_energy_estimate(hbar_cs: float, r: float, R: float) -> float:
    """Crude parallel-plate PFA scale (overestimates fat torus).

    E_PFA ∼ (π²/720) (ℏ c_s) (2π R) / r³   [dimensional proxy]
    Used only for ratio diagnostics, not as a lock.
    """
    return (math.pi**2 / 720.0) * hbar_cs * (2.0 * math.pi * R) / (r**3)


def fat_torus_report(
    r_um: float | None = None,
    R_um: float | None = None,
    hbar_cs: float | None = None,
) -> dict[str, Any]:
    """Full fat-torus vacuum report on preprint (r,R)=(8,50) nm by default."""
    p = PARAMS
    r = p.r_preprint_um if r_um is None else r_um
    R = p.R_preprint_um if R_um is None else R_um
    hbar_cs = p.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    x = 2.0 * math.pi * r / R
    F = bessel_mode_sum(x)
    E_tierC = vacuum_self_energy_corrected(hbar_cs, r, F)
    E_old = casimir_preprint_formula(hbar_cs, r, R, F)
    E_pfa = pfa_energy_estimate(hbar_cs, r, R)
    # Paper notes E_old / E_PFA ≈ 0.47 for fat torus with F≈0.596
    ratio_old_pfa = E_old / E_pfa if E_pfa else float("nan")
    ratio_tierC_pfa = E_tierC / E_pfa if E_pfa else float("nan")
    audit = assert_preprint_0p244_not_reproduced(r=r, R=R, hbar_cs=hbar_cs)
    return {
        "geometry": "preprint_fat_torus" if r_um is None else "custom",
        "r_nm": r * 1e3,
        "R_nm": R * 1e3,
        "x": x,
        "regime": "fat_torus" if 0.5 < x < 2.0 else ("thin_torus" if x < 0.5 else "thick_tube"),
        "F": F,
        "F_partial": F_partial_sums(x),
        "E_tierC_meV": E_tierC,
        "E_old_meV": E_old,
        "E_pfa_proxy_meV": E_pfa,
        "E_old_over_E_pfa": ratio_old_pfa,
        "E_tierC_over_E_pfa": ratio_tierC_pfa,
        "claimed_0p244_audit": audit,
        "jacobi_theta3_at_moire_alpha": jacobi_theta3_alternating(8.0 / 181.0),
        "notes": (
            "Tier-C form locked; meV provisional. "
            "PFA proxy is diagnostic only. Jacobi θ₃ ≠ F."
        ),
    }


def inverse_freeze_hbar_cs(
    target_meV: float,
    r_um: float,
    R_um: float,
    *,
    formula: str = "tierC",
) -> dict[str, Any]:
    """Solve for ℏc_s that yields target energy under locked F(x).

    Tier-C: E = (1/4π²)(ℏc_s/r) F  ⇒  ℏc_s = E · 4π² · r / F
    Old:    E = (1/4π)(ℏc_s/r)(R/r) F ⇒ ℏc_s = E · 4π · r² / (R F)
    """
    x = 2.0 * math.pi * r_um / R_um
    F = bessel_mode_sum(x)
    if formula == "tierC":
        hbar_cs = target_meV * (4.0 * math.pi**2) * r_um / F
        E_check = vacuum_self_energy_corrected(hbar_cs, r_um, F)
    elif formula == "old":
        hbar_cs = target_meV * (4.0 * math.pi) * (r_um**2) / (R_um * F)
        E_check = casimir_preprint_formula(hbar_cs, r_um, R_um, F)
    else:
        raise ValueError("formula must be 'tierC' or 'old'")
    return {
        "target_meV": target_meV,
        "formula": formula,
        "r_nm": r_um * 1e3,
        "R_nm": R_um * 1e3,
        "x": x,
        "F": F,
        "hbar_cs_required_meV_um": hbar_cs,
        "E_reproduced_meV": E_check,
        "package_hbar_cs": PARAMS.hbar_cs_meV_um,
        "ratio_vs_package": hbar_cs / PARAMS.hbar_cs_meV_um,
        "status": "inverse_diagnostic_not_a_freeze",
    }


def inverse_freeze_r_for_target(
    target_meV: float,
    R_um: float,
    hbar_cs: float | None = None,
    *,
    r_lo_um: float = 1e-5,
    r_hi_um: float = 0.05,
) -> dict[str, Any]:
    """Binary search r such that Tier-C E_vac(r,R)=target (fixed R, ℏc_s)."""
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs

    def E(r: float) -> float:
        x = 2.0 * math.pi * r / R_um
        F = bessel_mode_sum(x)
        return vacuum_self_energy_corrected(hbar_cs, r, F)

    lo, hi = r_lo_um, r_hi_um
    # E decreases with r typically; ensure bracket
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if E(mid) > target_meV:
            lo = mid
        else:
            hi = mid
    r_star = 0.5 * (lo + hi)
    return {
        "target_meV": target_meV,
        "R_nm": R_um * 1e3,
        "hbar_cs": hbar_cs,
        "r_star_nm": r_star * 1e3,
        "E_at_r_star": E(r_star),
        "x_star": 2.0 * math.pi * r_star / R_um,
        "status": "inverse_diagnostic_not_a_freeze",
    }


def geometry_comparison_table(hbar_cs: float | None = None) -> list[dict[str, Any]]:
    """Side-by-side vacuum evaluation for freeze candidates."""
    p = PARAMS
    hbar_cs = p.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    rows = []
    for label, r, R in [
        ("preprint_fat_torus", p.r_preprint_um, p.R_preprint_um),
        ("corrections_moire", p.r_corr_um, p.R_corr_um),
    ]:
        spec = FrozenVacuumSpec(
            formula_id="corrected_self_energy_2d",
            r=r,
            R=R,
            hbar_cs=hbar_cs,
            geometry_label=label,
            status="provisional",
        )
        ev = evaluate_vacuum(spec)
        rows.append(
            {
                "geometry": label,
                "r_nm": r * 1e3,
                "R_nm": R * 1e3,
                "x": ev["x"],
                "F": ev["F"],
                "E_tierC_meV": ev["E_corrected_self_energy_meV"],
                "E_old_meV": ev["E_old_casimir_meV"],
                "old_matches_0p244": ev["old_reproduces_0.244"],
            }
        )
    return rows


def si_cross_check(
    r_nm: float = 0.335,
    R_nm: float = 13.4,
    cs_m_s: float = 2.1e4,
) -> dict[str, Any]:
    """SI unit evaluation of Tier-C form (sensitive to c_s calibration)."""
    hbar = 1.054571817e-34  # J·s
    r_m = r_nm * 1e-9
    R_m = R_nm * 1e-9
    x = 2.0 * math.pi * r_m / R_m
    F = bessel_mode_sum(x)
    E_J = (1.0 / (4.0 * math.pi**2)) * (hbar * cs_m_s / r_m) * F
    E_meV = E_J / 1.602176634e-22
    return {
        "r_nm": r_nm,
        "R_nm": R_nm,
        "cs_m_s": cs_m_s,
        "x": x,
        "F": F,
        "E_meV_SI": E_meV,
        "note": "Order-of-magnitude cross-check; c_s must be platform-calibrated",
    }


def run_vacuum_freeze_suite() -> dict[str, Any]:
    """Aggregate freeze diagnostics for reporting."""
    from qvc.fat_torus_corrections import fat_torus_F_diagnostics, scan_fat_torus_aspect
    from qvc.self_consistency import bridge_scan_row, fold_in_E_vac_bound

    p = PARAMS
    fat = fat_torus_report()
    inv_hbar = inverse_freeze_hbar_cs(
        0.244, p.r_preprint_um, p.R_preprint_um, formula="tierC"
    )
    inv_r = inverse_freeze_r_for_target(0.244, p.R_preprint_um)
    geom = geometry_comparison_table()
    bridge_scan = [
        bridge_scan_row(row["E_tierC_meV"], row["geometry"]) for row in geom
    ]
    return {
        "fat_torus": fat,
        "fat_torus_F_diagnostics": fat_torus_F_diagnostics(),
        "fat_torus_aspect_scan": scan_fat_torus_aspect(),
        "geometry_table": geom,
        "fold_in_bounds": fold_in_E_vac_bound(),
        "bridge_scan": bridge_scan,
        "inverse_hbar_cs_for_0p244_tierC_fat": inv_hbar,
        "inverse_r_for_0p244_tierC_fat": inv_r,
        "si_corrections": si_cross_check(),
        "freeze_open_decisions": [
            "geometry: preprint_fat_torus vs corrections_moire",
            "formula: corrected_self_energy_2d (recommended)",
            "hbar_cs calibration vs SI/platform",
            "n_modes entering A_eff=√n A_single (vs bare bridge)",
            "bridge choice: package Tier-C already past fold on fat torus",
            "promote status provisional → frozen only after 1–5 fixed",
        ],
    }
