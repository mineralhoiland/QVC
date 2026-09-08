"""Gap-map branch solvers and fold critical table — Tier C lock.

CLOSED-FORM LOCK
----------------
δ = 1 + α (ln δ + ℓ)
At fold: α_c = δ_c and δ_c (1 - ln δ_c - ℓ) = 1
Near fold: δ - δ_c ∼ ± √(α_c - α)

Corrections-canon IR: ℓ ≈ -2.791 → δ_c ≈ 0.182
Tree-level: ℓ_tree = ln(1/2) + γ_E ≈ -0.11593
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.params import PARAMS

# Corrections-canon IR coefficient (from fold inversion at δ_c ≈ 0.182)
ELL_CORRECTIONS = -2.791


def gap_rhs(delta: float, alpha: float, ell: float) -> float:
    """Right-hand side of δ = 1 + α (ln δ + ℓ)."""
    if delta <= 0:
        return float("nan")
    return 1.0 + alpha * (math.log(delta) + ell)


def fold_residual(delta: float, alpha: float, ell: float) -> float:
    """δ - [1 + α(ln δ + ℓ)]."""
    return delta - gap_rhs(delta, alpha, ell)


def fold_condition(ell: float) -> tuple[float, float]:
    """Solve δ_c (1 - ln δ_c - ℓ) = 1; return (δ_c, α_c) with α_c = δ_c."""
    lo, hi = 1e-12, 1.0

    def f(d: float) -> float:
        return d * (1.0 - math.log(d) - ell) - 1.0

    for _ in range(120):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    delta_c = 0.5 * (lo + hi)
    return delta_c, delta_c


def fold_square_root_scaling(
    A: float, A_c: float, Delta_c: float, sign: float = 1.0
) -> float:
    """Δ ≈ Δ_c ± √(A_c - A) (leading saddle-node form; A ≤ A_c)."""
    if A > A_c:
        return float("nan")
    return Delta_c + sign * math.sqrt(max(A_c - A, 0.0))


def tree_level_ell(gamma_E: float | None = None) -> float:
    """ℓ_tree = ln(1/2) + γ_E."""
    gE = PARAMS.gamma_E if gamma_E is None else gamma_E
    return math.log(0.5) + gE


def ell_from_delta_c(delta_c: float) -> float:
    """Invert fold condition: ℓ = 1 - ln δ_c - 1/δ_c."""
    if delta_c <= 0:
        raise ValueError("delta_c must be positive")
    return 1.0 - math.log(delta_c) - 1.0 / delta_c


def solve_fold_branches(
    alpha: float, ell: float, *, n_scan: int = 8000
) -> dict[str, float | None]:
    """Find physical (larger) and unstable (smaller) roots of the fold map."""
    delta_c, alpha_c = fold_condition(ell)
    if alpha <= 0:
        return {"phys": 1.0, "unst": None, "alpha_c": alpha_c, "delta_c": delta_c}
    if alpha >= alpha_c * (1.0 - 1e-12):
        return {"phys": None, "unst": None, "alpha_c": alpha_c, "delta_c": delta_c}

    roots: list[float] = []
    d_lo, d_hi = 1e-6, 1.0
    prev_d = d_lo
    prev_f = fold_residual(prev_d, alpha, ell)
    for i in range(1, n_scan + 1):
        d = d_lo + (d_hi - d_lo) * i / n_scan
        f = fold_residual(d, alpha, ell)
        if prev_f == 0.0:
            roots.append(prev_d)
        elif prev_f * f < 0.0:
            t = -prev_f / (f - prev_f)
            roots.append(prev_d + t * (d - prev_d))
        prev_d, prev_f = d, f
    roots = sorted(set(round(r, 12) for r in roots), reverse=True)
    return {
        "phys": roots[0] if roots else None,
        "unst": roots[1] if len(roots) > 1 else None,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
    }


def fold_critical_table(
    ell: float = ELL_CORRECTIONS,
    Delta0_meV: float | None = None,
) -> dict[str, Any]:
    """Canonical fold numbers for appendix / freeze suite."""
    Delta0 = PARAMS.Delta0_meV if Delta0_meV is None else Delta0_meV
    delta_c, alpha_c = fold_condition(ell)
    return {
        "ell": ell,
        "ell_tree": tree_level_ell(),
        "delta_c": delta_c,
        "alpha_c": alpha_c,
        "Delta0_meV": Delta0,
        "Delta_c_meV": delta_c * Delta0,
        "scaling": "delta - delta_c ~ ± sqrt(alpha_c - alpha)",
        "bkt_status": "contingent_on_vortex_RG",
        "status": "locked",
    }


def fold_curve(
    ell: float = ELL_CORRECTIONS, n: int = 100
) -> dict[str, Any]:
    """Sample physical/unstable branches vs α/α_c for plotting/tables."""
    delta_c, alpha_c = fold_condition(ell)
    alphas: list[float] = []
    phys: list[float | None] = []
    unst: list[float | None] = []
    for i in range(n):
        a = alpha_c * 0.995 * i / max(n - 1, 1)
        br = solve_fold_branches(a, ell)
        alphas.append(a / alpha_c)
        phys.append(br["phys"])
        unst.append(br["unst"])
    return {
        "alpha_over_ac": alphas,
        "delta_phys": phys,
        "delta_unst": unst,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "ell": ell,
    }
