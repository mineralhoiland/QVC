"""Mexican-hat potential on the GPE ρ axis.

V(ρ) = α_GL ρ² + (β/2) ρ⁴ − λ★ E_vac ρ²

α_GL = −Δ₀, β = Δ₀ / n_ref from hopfion_initial.
The cartoon (α0=−0.5, α_QVC=−0.35, 42% drop) is not used.

Two evaluations share the same ρ axis and the same β:

1. GPE-sign vacuum shift −λ★ E_vac ρ² (attractive V_ZPF). The well deepens.
2. Parallel catalyzon dressing |α| → Δ_cat = Δ₀ − g0★ (shallower hat).
   Barrier reduction 1 − (Δ_cat/Δ₀)² replaces the cartoon 42%.
   Do not stack this percentage with fold or live-gap channels.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc_l1.locks import CARTOON_BARRIER_PCT, L1Locks


def V_hat(rho: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """V(ρ) = α ρ² + (β/2) ρ⁴."""
    r2 = np.asarray(rho, dtype=float) ** 2
    return alpha * r2 + 0.5 * beta * r2 * r2


def well_depth(alpha: float, beta: float) -> dict[str, float]:
    """Height from the Mexican-hat minimum up to ρ=0. Requires α<0, β>0."""
    if alpha >= 0.0 or beta <= 0.0:
        return {
            "rho_eq": float("nan"),
            "V_min": float("nan"),
            "barrier": float("nan"),
        }
    rho_eq = float(np.sqrt(-alpha / beta))
    V_min = float(alpha * rho_eq**2 + 0.5 * beta * rho_eq**4)
    return {
        "rho_eq": rho_eq,
        "V_min": V_min,
        "barrier": float(-V_min),  # V(0)=0
    }


def mexican_hat(locks: L1Locks, *, beta_meV: float, n_ref: float) -> dict[str, Any]:
    """Recompute the hat on the GPE ρ axis. Replace cartoon 42%."""
    alpha_GL = -locks.Delta0_meV
    g0 = float(locks.lambda_star * locks.E_vac_meV)
    Delta_cat = locks.Delta0_meV - g0

    # 1. As written: V = α_GL ρ² + (β/2) ρ⁴ − λ★ E_vac ρ²
    alpha_star_gpe = alpha_GL - g0  # more negative
    bare = well_depth(alpha_GL, beta_meV)
    gpe_shift = well_depth(alpha_star_gpe, beta_meV)

    # 2. Parallel catalyzon dressing of the Landau mass (shallower hat)
    alpha_cat = -Delta_cat
    cat = well_depth(alpha_cat, beta_meV)
    barrier_reduction_cat = 1.0 - cat["barrier"] / bare["barrier"]

    rho_max = 1.6 * max(bare["rho_eq"], cat["rho_eq"], gpe_shift["rho_eq"])
    rho = np.linspace(0.0, rho_max, 401)
    V_bare = V_hat(rho, alpha_GL, beta_meV)
    V_gpe = V_hat(rho, alpha_star_gpe, beta_meV)
    V_cat = V_hat(rho, alpha_cat, beta_meV)

    return {
        "formula": "V(ρ) = α_GL ρ² + (β/2) ρ⁴ − λ★ E_vac ρ²",
        "alpha_GL_meV": alpha_GL,
        "beta_meV": float(beta_meV),
        "n_ref": float(n_ref),
        "g0_star_meV": g0,
        "lambda_star": locks.lambda_star,
        "rho": rho,
        "V_bare": V_bare,
        "V_gpe_sign": V_gpe,
        "V_catalyzon_dressed": V_cat,
        "bare": bare,
        "gpe_sign_shift": {
            **gpe_shift,
            "alpha_eff_meV": alpha_star_gpe,
            "well_deepens": gpe_shift["barrier"] > bare["barrier"],
            "depth_change_frac": gpe_shift["barrier"] / bare["barrier"] - 1.0,
            "status": "computed",
            "note": (
                "The −λ★ E_vac ρ² term is the GPE V_ZPF sign (attractive). "
                "The well deepens; this is not the cartoon 42% reduction."
            ),
        },
        "catalyzon_dressed": {
            **cat,
            "alpha_eff_meV": alpha_cat,
            "Delta_cat_meV": Delta_cat,
            "barrier_reduction_frac": barrier_reduction_cat,
            "barrier_reduction_pct": 100.0 * barrier_reduction_cat,
            "formula": "1 − (Δ_cat/Δ₀)² with Δ_cat = Δ₀ − g0★",
            "status": "computed",
            "replaces_cartoon_42pct": True,
            "cartoon_42pct_not_used": CARTOON_BARRIER_PCT,
            "parallel_channel_only": True,
            "note": (
                "Landau mass tracks the catalyzon-dressed gap. Quadratic well "
                "depth scales as α². This percentage is not stacked with the "
                "fold or live-gap channels."
            ),
        },
        "cartoon_not_used": {
            "alpha0": -0.5,
            "alpha_QVC": -0.35,
            "barrier_pct": CARTOON_BARRIER_PCT,
            "used": False,
        },
    }
