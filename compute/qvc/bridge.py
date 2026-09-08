"""Derived α ↔ E_vac ↔ A_ZPF bridge (frozen-cavity gap equation).

PROVEN ALGEBRA (from the locked gap equation)
--------------------------------------------
    Δ = Δ₀ + ℏ c_s A [ ln(Δ a_H / (2 ℏ c_s)) + γ_E ],   Δ₀ = 2J.

Frozen cavity (preprint geometry): a_H = r = const, not ℏ c_s/Δ.
Then with δ = Δ/Δ₀,

    ℓ(r) = ln(Δ₀ r / (2 ℏ c_s)) + γ_E,
    α    = (ℏ c_s A) / Δ₀,
    δ    = 1 + α (ln δ + ℓ).

Geometric amplitude from the locked Tier-C form:

    A_geom = (1/4π²) (1/r) F(2π r/R) = E_vac / (ℏ c_s).

Unit participation (λ=1) therefore gives the unique dimensional map

    α = E_vac / Δ₀.

√n mode quadrature is NOT inserted here: F(x)=Σ K₀(nx) already sums
toroidal windings. Independent-channel √n is a separate optional factor.

ℓ_Corrections ≈ −2.791 is not an independent IR fit — it is ℓ(r) for the
preprint fat torus to O(0.2%).

Status
------
Algebra: PROVEN.
Geometry + ℏc_s: FROZEN_WORKING (fat torus, package 0.07 meV·µm).
meV E_vac: COMPUTED from frozen triple (replaces folklore 0.244).
Δ*: fold theorem — may be absent at λ=1 (supercritical).
"""

from __future__ import annotations

import math
from typing import Any

from qvc.gap_fold import (
    ELL_CORRECTIONS,
    fold_condition,
    solve_fold_branches,
    tree_level_ell,
)
from qvc.params import PARAMS
from qvc.vacuum_freeze import (
    FrozenVacuumSpec,
    bessel_mode_sum,
    casimir_preprint_formula,
    vacuum_self_energy_corrected,
)


def ell_frozen_cavity(
    r_um: float | None = None,
    *,
    Delta0: float | None = None,
    hbar_cs: float | None = None,
) -> float:
    """ℓ = ln(Δ₀ r / (2 ℏ c_s)) + γ_E  for frozen healing radius r."""
    p = PARAMS
    r_um = p.r_preprint_um if r_um is None else r_um
    Delta0 = p.Delta0_meV if Delta0 is None else Delta0
    hbar_cs = p.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    arg = Delta0 * r_um / (2.0 * hbar_cs)
    if arg <= 0:
        raise ValueError("log argument must be positive")
    return math.log(arg) + p.gamma_E


def A_geom_from_F(r_um: float, F: float) -> float:
    """A_geom = (1/4π²) (1/r) F   [µm⁻¹]."""
    return (1.0 / (4.0 * math.pi**2)) * (1.0 / r_um) * F


def alpha_from_A(A_um_inv: float, *, Delta0: float | None = None, hbar_cs: float | None = None) -> float:
    """α = ℏ c_s A / Δ₀  (proven rewrite of the gap equation)."""
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    return (hbar_cs * A_um_inv) / Delta0


def alpha_from_E_vac(E_vac_meV: float, *, Delta0: float | None = None, lam: float = 1.0) -> float:
    """α = λ E_vac / Δ₀  with participation λ (λ=1 is unit coupling)."""
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    return lam * E_vac_meV / Delta0


def alpha_for_target_delta(delta: float, ell: float) -> float:
    """Invert fold map: α = (δ − 1) / (ln δ + ℓ)."""
    if delta <= 0 or delta >= 1:
        raise ValueError("target δ must lie in (0, 1)")
    return (delta - 1.0) / (math.log(delta) + ell)


def derive_bridge(
    r_um: float | None = None,
    R_um: float | None = None,
    hbar_cs: float | None = None,
    *,
    lam: float = 1.0,
) -> dict[str, Any]:
    """Full derived bridge + fold evaluation for one frozen cavity."""
    p = PARAMS
    r_um = p.r_preprint_um if r_um is None else r_um
    R_um = p.R_preprint_um if R_um is None else R_um
    hbar_cs = p.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    Delta0 = p.Delta0_meV

    x = 2.0 * math.pi * r_um / R_um
    F = bessel_mode_sum(x)
    E_vac = vacuum_self_energy_corrected(hbar_cs, r_um, F)
    E_old = casimir_preprint_formula(hbar_cs, r_um, R_um, F)
    A = A_geom_from_F(r_um, F)
    # Consistency: E_vac = ℏ c_s A
    E_from_A = hbar_cs * A
    ell = ell_frozen_cavity(r_um, Delta0=Delta0, hbar_cs=hbar_cs)
    delta_c, alpha_c = fold_condition(ell)
    alpha = alpha_from_E_vac(E_vac, Delta0=Delta0, lam=lam)
    alpha_from_amp = alpha_from_A(A, Delta0=Delta0, hbar_cs=hbar_cs) * lam
    br = solve_fold_branches(alpha, ell)

    Delta_star = None if br["phys"] is None else br["phys"] * Delta0
    gap_pct = None if Delta_star is None else 100.0 * (1.0 - br["phys"])

    # Inverse: ℏc_s that sits exactly at the fold (λ=1)
    hbar_cs_fold = alpha_c * Delta0 * (4.0 * math.pi**2) * r_um / F
    # Participation that sits at 0.90 α_c
    lam_sub = (0.90 * alpha_c / (E_vac / Delta0)) if E_vac else float("nan")
    # α that would give the withdrawn 35% reduction (δ=0.389/0.6)
    delta_folklore = 0.389 / 0.600
    alpha_35 = alpha_for_target_delta(delta_folklore, ell)
    E_for_35 = alpha_35 * Delta0

    return {
        "status_algebra": "PROVEN",
        "status_geometry": "FROZEN_WORKING_fat_torus",
        "status_meV": "COMPUTED_FROM_FROZEN_TRIPLE",
        "status_Delta_star": "fold_theorem",
        "freeze_option": "A_star_frozen_cavity_fat_torus",
        "r_nm": r_um * 1e3,
        "R_nm": R_um * 1e3,
        "hbar_cs_meV_um": hbar_cs,
        "Delta0_meV": Delta0,
        "x": x,
        "F": F,
        "F_preprint_claim": 0.596,
        "E_vac_meV": E_vac,
        "E_vac_replaces_0p244": 0.244,
        "E_old_meV": E_old,
        "A_geom_um_inv": A,
        "A_times_R": A * R_um,
        "A_times_xi_preprint_claim": 0.296,
        "E_vac_equals_hbar_cs_A": abs(E_vac - E_from_A) < 1e-12,
        "ell_frozen": ell,
        "ell_corrections": ELL_CORRECTIONS,
        "ell_tree": tree_level_ell(),
        "ell_rel_err_vs_corrections": abs(ell - ELL_CORRECTIONS) / abs(ELL_CORRECTIONS),
        "participation_lambda": lam,
        "alpha": alpha,
        "alpha_from_A_check": alpha_from_amp,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "alpha_over_ac": alpha / alpha_c if alpha_c else float("nan"),
        "folded_out": br["phys"] is None,
        "delta_star": br["phys"],
        "delta_unst": br["unst"],
        "Delta_star_meV": Delta_star,
        "gap_reduction_pct": gap_pct,
        "hbar_cs_at_fold_meV_um": hbar_cs_fold,
        "lambda_for_0p90_ac": lam_sub,
        "alpha_for_withdrawn_35pct": alpha_35,
        "E_vac_for_withdrawn_35pct_meV": E_for_35,
        "bridge": "alpha = lambda * E_vac / Delta0,  A = E_vac / (hbar_cs),  ell = ln(Delta0 r / (2 hbar_cs)) + gamma_E",
        "notes": (
            "0.244 meV is withdrawn. Unit-participation fat-torus Tier-C is "
            "supercritical (α>α_c) at package ℏc_s. Running healing is not "
            "the Corrections-ℓ theory."
        ),
    }


def frozen_fat_torus_spec() -> FrozenVacuumSpec:
    """Working freeze: fat torus + Tier-C + package ℏc_s."""
    p = PARAMS
    return FrozenVacuumSpec(
        formula_id="corrected_self_energy_2d",
        r=p.r_preprint_um,
        R=p.R_preprint_um,
        hbar_cs=p.hbar_cs_meV_um,
        F_definition="bessel_K0_sum",
        geometry_label="preprint_fat_torus",
        status="frozen_working",
        notes=(
            "FROZEN WORKING STANDARD: fat-torus cavity (8,50) nm, Tier-C form, "
            "derived bridge α=E_vac/Δ₀, frozen-cavity ℓ. meV computed; "
            "unit-participation Δ* is fold-out (theorem). Replaces 0.244 folklore."
        ),
    )
