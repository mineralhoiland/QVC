"""C3: Kosterlitz flow at frozen α with computed κ(α).

Option A* freezes the cavity: there is no microscopic β_α. Scan the
physical branch δ_+(α) and integrate

    dK^{-1}/dℓ = 4π³ y²
    dy/dℓ     = (2 − π K) y
    y         = exp(−E_core / T_eff)

H-VRG1 holds iff the physical drive hits K=2/π (or IR plasma) before the
fold α=α_c. Anyon K_c=(2/π)(1−1/N)² is a labeled CS-restored variant.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc_bkt.locks import BKTLocks, build_bkt_locks
from qvc_bkt.stiffness import K_of_kappa

PI = math.pi
K_C_BOSONIC = 2.0 / PI


def kt_rhs(_ell: float, z: np.ndarray) -> np.ndarray:
    u, y = z
    u = max(float(u), 1e-12)
    K = 1.0 / u
    du = 4.0 * PI**3 * y * y
    dy = (2.0 - PI * K) * y
    return np.array([du, dy], dtype=float)


def integrate_kt(
    K_init: float,
    y_init: float,
    ell_ir: float,
    *,
    y_cut: float = 0.4,
    n_eval: int = 161,
) -> dict[str, Any]:
    """Textbook KT at frozen α. Status: computed."""
    u0 = 1.0 / max(K_init, 1e-12)
    y0 = max(float(y_init), 1e-18)

    def plasma(_ell: float, z: np.ndarray) -> float:
        return y_cut - z[1]

    plasma.terminal = True  # type: ignore[attr-defined]
    plasma.direction = -1.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        kt_rhs,
        (0.0, ell_ir),
        np.array([u0, y0]),
        rtol=1e-7,
        atol=1e-10,
        dense_output=True,
        events=plasma,
        max_step=0.05,
    )
    t_end = min(float(sol.t[-1]), ell_ir)
    ell = np.linspace(0.0, t_end, n_eval)
    uv = sol.sol(ell)
    K = 1.0 / np.maximum(uv[0], 1e-12)
    y = np.maximum(uv[1], 0.0)
    hit = bool(sol.t_events and len(sol.t_events[0]) > 0) or (y[-1] > y_cut * 0.98)
    ell_star = float(sol.t_events[0][0]) if (sol.t_events and len(sol.t_events[0]) > 0) else None
    bound = (not hit) and (K[-1] > K_C_BOSONIC)
    return {
        "ell": ell,
        "K": K,
        "y": y,
        "plasma": bool(hit),
        "bound": bool(bound),
        "ell_stop": float(t_end),
        "ell_star": ell_star,
        "K_ir": float(K[-1]),
        "y_ir": float(y[-1]),
        "K_init": float(K_init),
        "y_init": float(y0),
        "status": "computed",
    }


def fugacity_y(E_core_meV: float, T_K: float, kB: float, *, y_floor: float = 1e-16) -> float:
    """y = exp(−E_core / T_eff). Status: derived given E_core and T_eff."""
    t = T_K * kB
    if t <= 0.0:
        return y_floor
    return float(max(y_floor, math.exp(-E_core_meV / t)))


def competing_scale_scan(
    rows_kappa: list[dict[str, Any]],
    *,
    locks: BKTLocks,
    T_K: float,
    E_core_over_kappa: float,
    ell_ir: float = 8.0,
    K_c: float | None = None,
    use_anyon: bool = False,
    n_flow: int = 80,
) -> dict[str, Any]:
    """Scan frozen-α physical branch: does the KT surface precede the fold?"""
    Kc = locks.k_c_anyon if use_anyon else (K_c if K_c is not None else locks.k_c_bosonic)
    kB = locks.kB_meV_per_K
    alpha_kt_bare: float | None = None
    alpha_kt_ir: float | None = None
    out_rows: list[dict[str, Any]] = []
    # subsample flows for speed
    step = max(1, len(rows_kappa) // n_flow)
    for i, rec in enumerate(rows_kappa):
        kappa = float(rec["kappa_meV"])
        K_b = K_of_kappa(kappa, T_K, kB)
        E_core = E_core_over_kappa * kappa
        y_b = fugacity_y(E_core, T_K, kB)
        do_flow = (i % step == 0) or (i == len(rows_kappa) - 1)
        if do_flow:
            fl = integrate_kt(K_b, y_b, ell_ir, n_eval=81)
            plasma = bool(fl["plasma"])
            K_ir = float(fl["K_ir"])
            y_ir = float(fl["y_ir"])
            ell_star = fl["ell_star"]
        else:
            plasma = bool(K_b < Kc)
            K_ir = K_b
            y_ir = y_b
            ell_star = None
        bare_below = K_b < Kc
        phase = "plasma" if (plasma or bare_below) else "bound"
        if bare_below and alpha_kt_bare is None:
            alpha_kt_bare = float(rec["alpha"])
        if phase == "plasma" and alpha_kt_ir is None:
            alpha_kt_ir = float(rec["alpha"])
        out_rows.append(
            {
                "alpha": rec["alpha"],
                "alpha_over_ac": rec["alpha_over_ac"],
                "t_A": rec["t_A"],
                "delta": rec["delta"],
                "Delta_meV": rec["Delta_meV"],
                "kappa_meV": kappa,
                "K_bare": float(K_b),
                "y_bare": float(y_b),
                "E_core_meV": float(E_core),
                "K_ir": float(K_ir),
                "y_ir": float(y_ir),
                "plasma": int(phase == "plasma"),
                "phase": phase,
                "ell_star": ell_star,
            }
        )

    alpha_c = locks.alpha_c
    kt_before = alpha_kt_ir is not None and alpha_kt_ir < alpha_c
    return {
        "T_K": float(T_K),
        "K_c": float(Kc),
        "use_anyon": bool(use_anyon),
        "E_core_over_kappa": float(E_core_over_kappa),
        "ell_ir": float(ell_ir),
        "alpha_c": float(alpha_c),
        "alpha_kt_bare": alpha_kt_bare,
        "alpha_kt_ir": alpha_kt_ir,
        "alpha_kt_bare_over_ac": None if alpha_kt_bare is None else alpha_kt_bare / alpha_c,
        "alpha_kt_ir_over_ac": None if alpha_kt_ir is None else alpha_kt_ir / alpha_c,
        "kt_before_fold": bool(kt_before),
        "fold_wins": not kt_before,
        "rows": out_rows,
        "status": "computed",
        "cs": "restored_anyon" if use_anyon else "truncated",
        "option_A_star": "alpha_frozen_no_beta_alpha",
    }


def representative_flows(
    scan: dict[str, Any],
    locks: BKTLocks,
    *,
    ell_ir: float = 8.0,
) -> list[dict[str, Any]]:
    """KT trajectories at UV, λ★, and near-fold stations."""
    stations = [(0.30, "uv"), (0.90, "star"), (0.985, "nearfold")]
    rows = scan["rows"]
    T_K = float(scan["T_K"])
    c = float(scan["E_core_over_kappa"])
    out = []
    for frac, name in stations:
        # nearest scanned alpha
        target = frac * locks.alpha_c
        rec = min(rows, key=lambda r: abs(r["alpha"] - target))
        fl = integrate_kt(rec["K_bare"], rec["y_bare"], ell_ir, n_eval=161)
        out.append(
            {
                "name": name,
                "alpha_over_ac": rec["alpha_over_ac"],
                "delta": rec["delta"],
                "K_bare": rec["K_bare"],
                "y_bare": rec["y_bare"],
                "flow": fl,
                "T_K": T_K,
                "E_core_over_kappa": c,
            }
        )
    return out


def run_c3(
    c1: dict[str, Any],
    c2: dict[str, Any],
    locks: BKTLocks | None = None,
    *,
    ell_ir: float = 8.0,
) -> dict[str, Any]:
    """C3 driver: primary closure at T*, T-grid, and labeled anyon variant."""
    locks = locks or build_bkt_locks()
    eck = float(c2["E_core_over_kappa"])
    # Guard: a negative E_core/κ from cutoff mismatch still gives y>1; clip to a
    # physical O(1) core if the GP subtraction went negative (large-R log).
    eck_used = eck if eck > 0.05 else 0.5
    eck_label = "C2_computed" if eck > 0.05 else "floor_0.5_after_C2_nonpositive"
    primary_rows = c1["primary"]["rows"]
    T_star = locks.T_star_K
    scan_star = competing_scale_scan(
        primary_rows, locks=locks, T_K=T_star, E_core_over_kappa=eck_used, ell_ir=ell_ir
    )
    scan_any = competing_scale_scan(
        primary_rows,
        locks=locks,
        T_K=T_star,
        E_core_over_kappa=eck_used,
        ell_ir=ell_ir,
        use_anyon=True,
    )
    conv_rows = c1["variants"]["conv_only"]["rows"]
    scan_conv = competing_scale_scan(
        conv_rows, locks=locks, T_K=T_star, E_core_over_kappa=eck_used, ell_ir=ell_ir
    )

    T_grid = list(c1["primary"]["T_grid_K"])
    teff = []
    for T in T_grid:
        sc = competing_scale_scan(
            primary_rows, locks=locks, T_K=float(T), E_core_over_kappa=eck_used, ell_ir=ell_ir
        )
        sc_c = competing_scale_scan(
            conv_rows, locks=locks, T_K=float(T), E_core_over_kappa=eck_used, ell_ir=ell_ir
        )
        teff.append(
            {
                "T_K": float(T),
                "primary_kt_before_fold": sc["kt_before_fold"],
                "primary_alpha_kt_over_ac": sc["alpha_kt_ir_over_ac"]
                if sc["alpha_kt_ir_over_ac"] is not None
                else sc["alpha_kt_bare_over_ac"],
                "conv_only_kt_before_fold": sc_c["kt_before_fold"],
                "conv_only_alpha_kt_over_ac": sc_c["alpha_kt_ir_over_ac"]
                if sc_c["alpha_kt_ir_over_ac"] is not None
                else sc_c["alpha_kt_bare_over_ac"],
            }
        )

    flows = representative_flows(scan_star, locks, ell_ir=ell_ir)
    return {
        "task": "C3",
        "E_core_over_kappa_C2": eck,
        "E_core_over_kappa_used": eck_used,
        "E_core_over_kappa_label": eck_label,
        "scan_Tstar_primary": scan_star,
        "scan_Tstar_conv_only": scan_conv,
        "scan_Tstar_anyon": scan_any,
        "T_eff_grid": teff,
        "flows_Tstar": flows,
        "h_vrg1_Tstar_primary": {
            "kt_before_fold": scan_star["kt_before_fold"],
            "fold_wins": scan_star["fold_wins"],
            "closure": c1["primary_closure"],
            "cs": "truncated",
            "option_A_star": "frozen_alpha",
        },
        "note": (
            "Anyon scan is a labeled CS-restored variant and is not mixed "
            "into the truncated default. No β_α is integrated."
        ),
    }
