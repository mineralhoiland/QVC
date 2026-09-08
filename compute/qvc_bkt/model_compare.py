"""C5: model comparison √t_A vs exp(−C/√t_A) on fold and vortex-implied Δ.

Over ≳1.5 decades in t_A. Fold-only C5 is a statement about the locked
amplitude map; it is not a BKT disproof. Vortex-implied Δ is a C4 construct.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.gap_fold import fold_condition, solve_fold_branches  # noqa: E402

from qvc_bkt.locks import BKTLocks, build_bkt_locks


def _mask_decades(
    tA: np.ndarray,
    y: np.ndarray,
    *,
    t_lo: float = 1e-3,
    t_hi: float = 10 ** (-0.5),
) -> np.ndarray:
    decades = math.log10(t_hi) - math.log10(t_lo)
    if decades < 1.5:
        raise ValueError("t_A window must span ≳1.5 decades")
    return (
        np.isfinite(tA)
        & np.isfinite(y)
        & (tA > t_lo)
        & (tA < t_hi)
        & (y > 0.0)
    )


def fit_sqrt_branch(tA: np.ndarray, delta: np.ndarray, delta_c: float) -> dict[str, float]:
    """ln(δ−δ_c) = c0 + s ln t_A. Fold theorem wants s=1/2."""
    y = delta - delta_c
    mask = _mask_decades(tA, y)
    t = tA[mask]
    z = np.log(np.maximum(y[mask], 1e-16))
    n = int(t.size)
    if n < 8:
        return {"n": float(n), "slope_free": float("nan"), "ll": float("nan")}
    X = np.column_stack([np.ones_like(t), np.log(t)])
    coef, *_ = np.linalg.lstsq(X, z, rcond=None)
    resid = z - X @ coef
    s2 = float(np.mean(resid**2))
    ll = -0.5 * n * math.log(2.0 * math.pi * max(s2, 1e-30)) - 0.5 * n
    # Forced slope 1/2
    Xh = np.column_stack([np.ones_like(t), 0.5 * np.log(t)])
    c_h, *_ = np.linalg.lstsq(Xh, z, rcond=None)
    r_h = z - Xh @ c_h
    s2h = float(np.mean(r_h**2))
    ll_h = -0.5 * n * math.log(2.0 * math.pi * max(s2h, 1e-30)) - 0.5 * n
    t_span = float(np.log10(t.max()) - np.log10(t.min()))
    return {
        "n": float(n),
        "slope_free": float(coef[1]),
        "intercept": float(coef[0]),
        "ll": float(ll),
        "ll_forced_half": float(ll_h),
        "mse": s2,
        "tA_decades": t_span,
        "status_expected": "slope_near_0.5_on_fold",
    }


def fit_essential(tA: np.ndarray, delta: np.ndarray) -> dict[str, float]:
    """ln δ = c0 − C / √t_A  (continuous vanishing; uses ln δ not ln(δ−δ_c))."""
    mask = _mask_decades(tA, delta)
    t = tA[mask]
    z = np.log(np.maximum(delta[mask], 1e-16))
    n = int(t.size)
    if n < 8:
        return {"n": float(n), "C": float("nan"), "ll": float("nan")}
    X = np.column_stack([np.ones_like(t), -1.0 / np.sqrt(t)])
    coef, *_ = np.linalg.lstsq(X, z, rcond=None)
    resid = z - X @ coef
    s2 = float(np.mean(resid**2))
    ll = -0.5 * n * math.log(2.0 * math.pi * max(s2, 1e-30)) - 0.5 * n
    t_span = float(np.log10(t.max()) - np.log10(t.min()))
    return {
        "n": float(n),
        "C": float(coef[1]),
        "intercept": float(coef[0]),
        "ll": float(ll),
        "mse": s2,
        "tA_decades": t_span,
    }


def compare_models(
    tA: np.ndarray,
    delta: np.ndarray,
    delta_c: float,
    *,
    label: str,
) -> dict[str, Any]:
    sqrt_fit = fit_sqrt_branch(tA, delta, delta_c)
    bkt_fit = fit_essential(tA, delta)
    dll = float(sqrt_fit.get("ll", float("nan")) - bkt_fit.get("ll", float("nan")))
    return {
        "label": label,
        "sqrt": sqrt_fit,
        "essential": bkt_fit,
        "delta_ll_sqrt_minus_essential": dll,
        "preferred": (
            "sqrt_fold"
            if np.isfinite(dll) and dll > 0
            else "essential_BKT_shape"
            if np.isfinite(dll)
            else "undecided"
        ),
        "tA_window": "[1e-3, 10^{-0.5}]  (≳2 decades)",
        "note": (
            "Gaussian iid residual likelihood; not a field-theory evidence ratio. "
            "Fold-only preference for √t_A is expected from the locked map and "
            "is not a BKT disproof."
        ),
    }


def fold_branch_curve(locks: BKTLocks, n: int = 240) -> dict[str, np.ndarray]:
    delta_c, alpha_c = fold_condition(locks.ell)
    tA = np.logspace(-3.0, 0.0, n)
    alpha = alpha_c * (1.0 - tA)
    d_phys = np.full_like(tA, np.nan)
    for i, a in enumerate(alpha):
        if a <= 0.0:
            d_phys[i] = 1.0
            continue
        br = solve_fold_branches(float(a), locks.ell)
        if br["phys"] is not None:
            d_phys[i] = float(br["phys"])
    return {
        "t_A": tA,
        "delta": d_phys,
        "Delta_meV": d_phys * locks.Delta0_meV,
        "delta_c": np.full_like(tA, delta_c),
    }


def run_c5(
    c4: dict[str, Any],
    locks: BKTLocks | None = None,
) -> dict[str, Any]:
    locks = locks or build_bkt_locks()
    fold = fold_branch_curve(locks)
    cmp_fold = compare_models(fold["t_A"], fold["delta"], locks.delta_c, label="locked_fold_delta")

    vortex_cmp = None
    vortex_note = "no finite vortex-implied Δ on the primary T* scan"
    prim = c4["primary_Tstar"]
    if prim["unbinding_occurs"]:
        rows = prim["rows"]
        tA = np.array([r["t_A"] for r in rows], dtype=float)
        # Use Δ_vort/Δ₀ as δ_vort; essential fit uses ln δ so small Δ is ok.
        Dv = np.array([r["Delta_vortex_meV"] for r in rows], dtype=float)
        dv = Dv / locks.Delta0_meV
        finite = np.array([r["finite_unbinding"] for r in rows], dtype=bool) & np.isfinite(dv)
        if np.count_nonzero(finite) >= 8:
            vortex_cmp = compare_models(
                tA[finite], dv[finite], 0.0, label="vortex_implied_Delta_from_C4"
            )
            # For essential vs sqrt, delta_c=0 is the H-VRG1 continuous vanishing.
            vortex_note = "compared on C4 vortex-implied Δ (unbound points)"
        else:
            vortex_note = "unbinding flagged but too few finite Δ_vort points for OLS"

    # Analytic local √ form for plotting
    tA = fold["t_A"]
    d = fold["delta"]
    finite = np.isfinite(d)
    i0 = int(np.where(finite)[0][0])
    B = (d[i0] - locks.delta_c) / math.sqrt(max(tA[i0], 1e-16))
    d_sqrt = locks.delta_c + B * np.sqrt(tA)
    C_bkt = 0.35
    d_bkt = np.exp(-C_bkt / np.sqrt(np.maximum(tA, 1e-12)))

    return {
        "task": "C5",
        "fold_comparison": cmp_fold,
        "vortex_comparison": vortex_cmp,
        "vortex_note": vortex_note,
        "curves": {
            "t_A": tA,
            "delta_phys": d,
            "delta_sqrt_local": d_sqrt,
            "delta_bkt_toy": d_bkt,
            "B_sqrt": B,
            "C_bkt_toy": C_bkt,
        },
        "disclaimer": (
            "Fold-only C5 preferring √t_A is the locked amplitude theorem and "
            "is not a BKT disproof. A BKT shape would have to win on a "
            "vortex-implied Δ constructed in C4, under stated C1 closures."
        ),
    }
