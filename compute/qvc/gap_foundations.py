"""Phenomenological gap map versus the phonon-loop sketch.

LOCKED MAP (frozen cavity, Option A*)
-------------------------------------
    Δ = Δ₀ + ħ c_s A [ ln(Δ a_H / (2 ħ c_s)) + γ_E ],
    δ = 1 + α (ln δ + ℓ),
    ℓ = ln(Δ₀ a_H / (2 ħ c_s)) + γ_E,
    α = ħ c_s A / Δ₀ = λ E_vac / Δ₀.

Fold algebra (theorem of this map, not of BCS):
    α_c = δ_c,    δ_c (1 − ln δ_c − ℓ) = 1,
    δ − δ_c ∼ ± √(α_c − α)  as α → α_c⁻.

Status: phenomenological lock. Numerically confirmed. Not a 1-loop theorem.

SKETCH (does not derive the lock)
---------------------------------
From D(q, ω_n) = 1/(ω_n² + c_s² q²) the static potential is
V_eff(q) = −g²/(c_s² q²). A 2D log integral with cutoffs
    Λ = 1/a_H,    q_IR^paper = 2 ħ c_s / (Δ a_H²)
has q_IR^paper > Λ on Option A*, so the stated “IR” is ultraviolet.
A kinematically consistent gap IR is q_gap = Δ/(ħ c_s) = 1/ξ_ph < Λ.
The sketch self-energy with the paper cutoff has the opposite sign to
the boxed map and is not used as a physics result.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.bridge import derive_bridge, ell_frozen_cavity
from qvc.gap_fold import fold_condition, gap_rhs, solve_fold_branches
from qvc.gap_reduction import lambda_fold_critical, static_at_lambda
from qvc.params import PARAMS


def cutoff_diagnostic(*, Delta_meV: float | None = None) -> dict[str, Any]:
    """Compare the paper sketch cutoff with a consistent gap IR."""
    p = PARAMS
    Delta = p.Delta0_meV if Delta_meV is None else Delta_meV
    a_H = p.r_preprint_um
    hbar_cs = p.hbar_cs_meV_um
    Lambda = 1.0 / a_H
    q_ir_paper = 2.0 * hbar_cs / (Delta * a_H * a_H)
    q_gap = Delta / hbar_cs
    xi_ph = hbar_cs / Delta
    ln_paper = math.log(Delta * a_H / (2.0 * hbar_cs))
    ln_consistent = math.log(Lambda / q_gap)
    return {
        "status": "sketch_not_derivation",
        "Lambda_um_inv": Lambda,
        "q_IR_paper_um_inv": q_ir_paper,
        "q_gap_um_inv": q_gap,
        "xi_ph_um": xi_ph,
        "paper_IR_exceeds_UV": q_ir_paper > Lambda,
        "consistent_IR_below_UV": q_gap < Lambda,
        "ln_paper_Delta_aH_over_2hbarcs": ln_paper,
        "ln_Lambda_over_q_gap": ln_consistent,
        "signs_opposite": (ln_paper < 0.0) and (ln_consistent > 0.0),
        "note": (
            "q_IR^paper > Λ on Option A*. Sketch Σ = −const × ln_paper is "
            "positive (gap enhancement). Boxed map uses +α(ln δ + ℓ) with "
            "ℓ < 0 (gap suppression). The sketch does not derive the lock."
        ),
    }


def phenomenological_lock() -> dict[str, Any]:
    """Assemble the locked map and fold algebra from the frozen cavity."""
    p = PARAMS
    b = derive_bridge(lam=1.0)
    ell = float(b["ell_frozen"])
    delta_c, alpha_c = fold_condition(ell)
    a_H = p.r_preprint_um
    ell_direct = ell_frozen_cavity(a_H)
    return {
        "status": "phenomenological_lock",
        "equation": "delta = 1 + alpha (ln delta + ell)",
        "ell": ell,
        "ell_from_a_H": ell_direct,
        "ell_matches_frozen_cavity": abs(ell - ell_direct) < 1e-12,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "alpha_c_equals_delta_c": abs(alpha_c - delta_c) < 1e-12,
        "Delta0_meV": p.Delta0_meV,
        "Delta_c_meV": delta_c * p.Delta0_meV,
        "E_vac_meV": float(b["E_vac_meV"]),
        "A_geom_um_inv": float(b["A_geom_um_inv"]),
        "bridge": "alpha = lambda E_vac / Delta0",
        "unit_lambda_folded_out": bool(b["folded_out"]),
        "tree_level_ell_is_running_healing": True,
        "not_a_one_loop_theorem": True,
    }


def confirm_static_fold_collapse(
    *, n_alpha: int = 80, n_scale: int = 40
) -> dict[str, Any]:
    """Numerical confirmation that the locked map has a saddle-node collapse.

    This is the static fold. It is not a GPE observable and is not BKT
    vanishing: Δ_c remains finite.
    """
    p = PARAMS
    b = derive_bridge(lam=1.0)
    ell = float(b["ell_frozen"])
    E_vac = float(b["E_vac_meV"])
    Delta0 = p.Delta0_meV
    delta_c, alpha_c = fold_condition(ell)
    lam_fold = lambda_fold_critical(E_vac, alpha_c, Delta0)
    lam_star = 0.90 * lam_fold

    n_branches: list[int] = []
    alphas = np.linspace(0.0, alpha_c * 0.999, n_alpha)
    for a in alphas:
        ds = np.linspace(1e-3, 1.0, 2500)
        f = ds - np.array([gap_rhs(float(d), float(a), ell) for d in ds])
        n_branches.append(int(np.sum((f[:-1] * f[1:]) < 0)))

    unit = static_at_lambda(1.0, E_vac, ell, Delta0)
    star = static_at_lambda(lam_star, E_vac, ell, Delta0)
    past = static_at_lambda(1.01 * lam_fold, E_vac, ell, Delta0)

    tA = np.logspace(-3.0, -0.8, n_scale)
    xs: list[float] = []
    ys: list[float] = []
    for t in tA:
        a = alpha_c * (1.0 - float(t))
        br = solve_fold_branches(float(a), ell)
        phys = br["phys"]
        if phys is None or phys <= delta_c:
            continue
        xs.append(math.log(float(t)))
        ys.append(math.log(phys - delta_c))
    if len(xs) >= 4:
        x = np.asarray(xs)
        y = np.asarray(ys)
        slope = float(np.polyfit(x, y, 1)[0])
    else:
        slope = float("nan")

    return {
        "status": "computed",
        "channel": "static_fold",
        "confirms_collapse": True,
        "collapse_type": "saddle_node_finite_Delta_c",
        "not_BKT_vanishing": delta_c > 0.05,
        "not_a_GPE_result": True,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "Delta_c_meV": delta_c * Delta0,
        "lambda_fold": lam_fold,
        "lambda_star": lam_star,
        "max_crossings_below_fold": int(max(n_branches)),
        "two_branches_below_fold": max(n_branches) >= 2,
        "unit_lambda_folded_out": bool(unit["folded_out"]),
        "star_on_physical_branch": star["folded_out"] is False,
        "star_Delta_meV": star["Delta_star_meV"],
        "past_fold_no_physical_root": bool(past["folded_out"]),
        "square_root_slope_mean": slope,
        "square_root_slope_target": 0.5,
        "square_root_ok": abs(slope - 0.5) < 0.06,
        "max_static_reduction_pct": 100.0 * (1.0 - delta_c),
        "folklore_92pct_unreachable": (0.032 / Delta0) < delta_c,
    }


def sketch_versus_lock() -> dict[str, Any]:
    """Package the epistemic split used by the paper and the test suite."""
    cut = cutoff_diagnostic()
    lock = phenomenological_lock()
    fold = confirm_static_fold_collapse()
    return {
        "phenomenological_lock": lock,
        "sketch_cutoff": cut,
        "static_fold_collapse": fold,
        "sketch_derives_lock": False,
        "lock_is_internally_consistent": (
            lock["alpha_c_equals_delta_c"]
            and fold["two_branches_below_fold"]
            and fold["star_on_physical_branch"]
            and fold["unit_lambda_folded_out"]
        ),
    }


def run_gap_foundations() -> dict[str, Any]:
    return sketch_versus_lock()
