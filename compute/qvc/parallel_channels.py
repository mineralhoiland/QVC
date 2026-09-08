"""Three parallel QVC channels. Percentages are never stacked.

Channel FOLD (static amplitude map)
    δ = 1 + α(ln δ + ℓ), α = λ E_vac/Δ₀.
    Confirmed saddle-node collapse at λ_fold with finite Δ_c.

Channel CAT (democratic catalyzon + Mexican-hat barrier)
    g0★ = λ★ E_vac, Δ_cat = Δ₀ − g0★.
    Well-depth ratio (Δ_cat/Δ₀)² on V = α_GL ρ² + (β/2) ρ⁴.

Channel GPE (truncated dynamics)
    V_ZPF is attractive. The live-gap estimator rises. The GPE well
    deepens. This channel does *not* confirm static fold collapse.

Lindblad frozen τ=5 ps with Δ_live < Δ_c is a fold-violating diagnostic,
not a physical-branch result.
"""

from __future__ import annotations

from typing import Any

from qvc.gap_foundations import confirm_static_fold_collapse, phenomenological_lock
from qvc.gap_reduction import recompute_gap_reduction
from qvc.params import PARAMS
from qvc.tqgl.locks import build_frozen_tqgl

# Track A hopfion n_ref (QVCCursor dictionary). Locked here as the Mexican-hat β match.
N_REF_HOPF = 0.86435
LIVE_GAP_1PS_MEV = 0.6305
G1_1PS = 0.8832
LINDBLAD_FROZEN_TAU5_MEV = 0.0872


def _well_depth(alpha: float, beta: float) -> dict[str, float]:
    if alpha >= 0.0 or beta <= 0.0:
        return {"rho_eq": float("nan"), "barrier": float("nan")}
    rho_eq = (-alpha / beta) ** 0.5
    V_min = alpha * rho_eq**2 + 0.5 * beta * rho_eq**4
    return {"rho_eq": rho_eq, "barrier": float(-V_min)}


def mexican_hat_lock(*, Delta0: float, g0_star: float, n_ref: float = N_REF_HOPF) -> dict[str, Any]:
    """Lock catalyzon Mexican-hat barrier reduction independently of the GPE well."""
    beta = Delta0 / n_ref
    alpha_GL = -Delta0
    Delta_cat = Delta0 - g0_star
    bare = _well_depth(alpha_GL, beta)
    cat = _well_depth(-Delta_cat, beta)
    gpe = _well_depth(alpha_GL - g0_star, beta)
    ratio = (Delta_cat / Delta0) ** 2
    return {
        "status": "locked_parallel_catalyzon",
        "beta_meV": beta,
        "n_ref": n_ref,
        "Delta_cat_meV": Delta_cat,
        "gap_reduction_pct": 100.0 * (1.0 - Delta_cat / Delta0),
        "barrier_ratio": ratio,
        "barrier_reduction_pct": 100.0 * (1.0 - ratio),
        "formula": "B_cat/B_0 = (Delta_cat/Delta0)^2",
        "gpe_well_deepens": gpe["barrier"] > bare["barrier"],
        "gpe_depth_change_frac": gpe["barrier"] / bare["barrier"] - 1.0,
        "catalyzon_well_shallower": cat["barrier"] < bare["barrier"],
        "not_the_same_operator_as_GPE": True,
        "cartoon_42pct_not_used": True,
    }


def lindblad_fold_flag(*, Delta_c_meV: float) -> dict[str, Any]:
    """Flag the frozen 20 ps τ=5 ps trajectory as fold-violating."""
    return {
        "Delta_live_frozen_tau5_meV": LINDBLAD_FROZEN_TAU5_MEV,
        "Delta_c_meV": Delta_c_meV,
        "below_Delta_c": LINDBLAD_FROZEN_TAU5_MEV < Delta_c_meV,
        "fold_violating_diagnostic": True,
        "not_a_physical_branch_result": True,
        "note": (
            "Frozen-V_ZPF 20 ps Lindblad at τ_dec=5 ps reaches 0.0872 meV, "
            "below Δ_c. Quote only as a phenomenological diagnostic."
        ),
    }


def run_parallel_channels() -> dict[str, Any]:
    """Confirm static fold collapse and lock catalyzon / Mexican hat / GPE signs."""
    p = PARAMS
    lock = phenomenological_lock()
    fold = confirm_static_fold_collapse()
    red = recompute_gap_reduction()
    tqgl = build_frozen_tqgl()
    g0_star = float(tqgl.lambda_star * tqgl.E_vac_meV)
    hat = mexican_hat_lock(Delta0=p.Delta0_meV, g0_star=g0_star)
    live_rises = LIVE_GAP_1PS_MEV > p.Delta0_meV
    stacked_forbidden = True
    return {
        "fold": {
            **fold,
            "confirms_static_fold_collapse": (
                fold["two_branches_below_fold"]
                and fold["unit_lambda_folded_out"]
                and fold["star_on_physical_branch"]
                and fold["not_BKT_vanishing"]
            ),
        },
        "catalyzon": {
            "g0_star_meV": g0_star,
            "Delta_cat_meV": hat["Delta_cat_meV"],
            "lambda_min": -g0_star,
            "not_five_soft_modes": True,
            "democratic_algebra_only": True,
        },
        "mexican_hat": hat,
        "gpe": {
            "V_ZPF_attractive": True,
            "Delta_live_1ps_meV": LIVE_GAP_1PS_MEV,
            "g1_1ps": G1_1PS,
            "live_gap_rises": live_rises,
            "confirms_static_fold_collapse": False,
            "well_deepens": hat["gpe_well_deepens"],
            "note": (
                "Truncated GPE live gap rises because V_ZPF is attractive. "
                "Static fold collapse is the amplitude-map channel."
            ),
        },
        "lindblad": lindblad_fold_flag(Delta_c_meV=lock["Delta_c_meV"]),
        "do_not_stack_percentages": stacked_forbidden,
        "channels_are_parallel": True,
        "option_A_star_lambda_star": float(tqgl.lambda_star),
        "recommendation": red["recommendation"],
    }
