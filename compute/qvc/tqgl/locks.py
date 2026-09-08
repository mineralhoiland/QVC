"""Assemble frozen Option A* numbers from existing QVCCompute APIs.

No folklore 0.244 / 0.389 / 0.032 is treated as a physical result.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any

from qvc.bridge import derive_bridge
from qvc.gap_fold import fold_condition, solve_fold_branches
from qvc.gap_reduction import lambda_fold_critical, static_at_lambda
from qvc.lens_chern_simons import fractional_hall_charge
from qvc.params import PARAMS

# ℏ in the (meV, ps) unit system used by the GPE / DCE integrators.
HBAR_MEV_PS = 0.6582119569

# Withdrawn preprint Casimir number — used only as a drive-ratio denominator.
E_CAS_WITHDRAWN_MEV = 0.244

# Withdrawn folklore (must not appear as suite results).
DELTA_QVC_FOLKLORE_MEV = 0.389
DELTA_1PS_FOLKLORE_MEV = 0.032


@dataclass(frozen=True)
class FrozenTQGL:
    """Immutable Option A* lock consumed by the GPE / DCE suite."""

    r_um: float
    R_um: float
    r_nm: float
    R_nm: float
    hbar_cs_meV_um: float
    Delta0_meV: float
    N_layers: int
    Q_H: float
    x: float
    F: float
    E_vac_meV: float
    A_geom_um_inv: float
    ell: float
    alpha_c: float
    delta_c: float
    Delta_c_meV: float
    max_static_reduction_pct: float
    lambda_fold: float
    lambda_star: float
    alpha_star: float
    alpha_star_over_ac: float
    Delta_star_meV: float
    delta_star: float
    gap_reduction_star_pct: float
    A_c_fold_um_inv: float
    A_star_um_inv: float
    drive_scale_vs_old_0244: float
    hbar_meV_ps: float
    kB_meV_per_K: float

    def as_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_frozen_tqgl(*, lambda_star_frac: float = 0.90) -> FrozenTQGL:
    """Pull vacuum, fold, and λ★ from locked modules (not hardcoded folklore)."""
    p = PARAMS
    b = derive_bridge(lam=1.0)
    E_vac = float(b["E_vac_meV"])
    ell = float(b["ell_frozen"])
    Delta0 = float(b["Delta0_meV"])
    alpha_c = float(b["alpha_c"])
    delta_c = float(b["delta_c"])
    lam_fold = lambda_fold_critical(E_vac, alpha_c, Delta0)
    lam_star = lambda_star_frac * lam_fold
    star = static_at_lambda(lam_star, E_vac, ell, Delta0)
    if star["folded_out"] or star["Delta_star_meV"] is None:
        raise RuntimeError("λ★ must sit on the physical fold branch")

    A_geom = float(b["A_geom_um_inv"])
    A_c = alpha_c * Delta0 / p.hbar_cs_meV_um
    A_star = lam_star * A_geom

    return FrozenTQGL(
        r_um=p.r_preprint_um,
        R_um=p.R_preprint_um,
        r_nm=p.r_preprint_um * 1e3,
        R_nm=p.R_preprint_um * 1e3,
        hbar_cs_meV_um=p.hbar_cs_meV_um,
        Delta0_meV=Delta0,
        N_layers=p.N_layers,
        Q_H=fractional_hall_charge(p.N_layers),
        x=float(b["x"]),
        F=float(b["F"]),
        E_vac_meV=E_vac,
        A_geom_um_inv=A_geom,
        ell=ell,
        alpha_c=alpha_c,
        delta_c=delta_c,
        Delta_c_meV=delta_c * Delta0,
        max_static_reduction_pct=100.0 * (1.0 - delta_c),
        lambda_fold=lam_fold,
        lambda_star=lam_star,
        alpha_star=float(star["alpha"]),
        alpha_star_over_ac=float(star["alpha_over_ac"]),
        Delta_star_meV=float(star["Delta_star_meV"]),
        delta_star=float(star["delta_star"]),
        gap_reduction_star_pct=float(star["gap_reduction_pct"]),
        A_c_fold_um_inv=A_c,
        A_star_um_inv=A_star,
        drive_scale_vs_old_0244=E_vac / E_CAS_WITHDRAWN_MEV,
        hbar_meV_ps=HBAR_MEV_PS,
        kB_meV_per_K=p.kB_meV_per_K,
    )


def fold_delta_of_A(A_um_inv: float, lock: FrozenTQGL) -> dict[str, Any]:
    """Static fold branches at amplitude A (α = ℏ c_s A / Δ₀)."""
    alpha = lock.hbar_cs_meV_um * A_um_inv / lock.Delta0_meV
    br = solve_fold_branches(alpha, lock.ell)
    phys = br["phys"]
    return {
        "A_um_inv": A_um_inv,
        "alpha": alpha,
        "alpha_over_ac": alpha / lock.alpha_c if lock.alpha_c else float("nan"),
        "folded_out": phys is None,
        "delta_phys": phys,
        "Delta_phys_meV": None if phys is None else phys * lock.Delta0_meV,
        "delta_unst": br["unst"],
        "past_fold": bool(phys is None or (phys is not None and phys < lock.delta_c)),
    }


def fold_condition_numbers(lock: FrozenTQGL) -> tuple[float, float]:
    """Re-export (δ_c, α_c) from the locked ℓ."""
    return fold_condition(lock.ell)
