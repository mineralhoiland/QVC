r"""Diagnostic vortex RG for the QVC fold vs contingent BKT (H-VRG1).

Status
------
This module is a *diagnostic closure*, not a microscopic derivation of KT
flow from the TQGL functional.

- Fold numbers are **locked** (``qvc.gap_fold``).
- Two-sector split and competing-scale algebra are **derived here**.
- Integrals and finite-size illustrations are **computed**.
- H-VRG1 promotion and anyon-plasma shift are **conjecture**.
- Chern--Simons is **truncated** from the XY--KT sector (same truncation
  as the live GPE), unless the anyon correction is switched on explicitly.

Kosterlitz normalisation
------------------------
We use \(K=\kappa/T_{\mathrm{eff}}\) so that the textbook flow
\(\partial_\ell K^{-1}=4\pi^3 y^2\), \(\partial_\ell y=(2-\pi K)y\)
has critical point \(K_c=2/\pi\).  Then \(\pi\kappa/T_{\mathrm{eff}}=\pi K\)
equals 2 on the universal jump.  Do not mix this \(K\) with the
José--Kadanoff convention \(\tilde K=\pi\kappa/T\) (\(\tilde K_c=2\)).

Diagnostic stiffness
--------------------
On the *physical* fold branch, postulate \(\kappa(\alpha)\propto\Delta(\alpha)^2
=\Delta_0^2\delta_+(\alpha)^2\).  Bare \(K(\alpha)=K_0[\delta_+(\alpha)]^2\)
with \(K_0=(2/\pi)(T_{\mathrm{BKT}}^{\mathrm{bare}}/T_{\mathrm{eff}})\),
i.e. the undriven stiffness is assigned a BKT temperature in the TEVC
window as an *input*, not a computed \(\kappa\).
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import solve_ivp

from qvc.gap_fold import (
    ELL_CORRECTIONS,
    fold_condition,
    solve_fold_branches,
)
from qvc.params import PARAMS

PI = math.pi
K_C_BOSONIC = 2.0 / PI  # 2/π in Kosterlitz normalisation
N_LAYERS = 5
THETA_STAT = PI / N_LAYERS  # π/N; conjecture if CS restored
# Anyon-plasma shift (Wen–Zee / statistical angle); conjecture.
K_C_ANYON = K_C_BOSONIC * (1.0 - THETA_STAT / PI) ** 2

# Locked Option A* (must match macros.tex)
ALPHA_C_LOCK = 0.1818
DELTA_C_MEV = 0.109
DELTA0_MEV = 0.600
LAMBDA_STAR = 0.7628
DELTA_STAR_MEV = 0.2324
T_STAR_K = 0.1187  # Schwinger diagnostic; not T_BKT
T_BKT_BARE_DEFAULT = 1.0  # TEVC window input, not derived from κ


@dataclass(frozen=True)
class VortexRGLocks:
    """Immutable diagnostic parameters (fold locked; rest labelled)."""

    ell_gap: float = ELL_CORRECTIONS
    alpha_c: float = ALPHA_C_LOCK
    delta_c: float = ALPHA_C_LOCK
    delta0_meV: float = DELTA0_MEV
    delta_c_meV: float = DELTA_C_MEV
    lambda_star: float = LAMBDA_STAR
    delta_star_meV: float = DELTA_STAR_MEV
    T_star_K: float = T_STAR_K
    T_BKT_bare_K: float = T_BKT_BARE_DEFAULT
    T_eff_K: float = T_STAR_K
    y_uv: float = 1e-4  # numerical floor; physical y=exp(-c K)
    core_coeff: float = 1.0  # E_core / κ ~ O(1)
    ell_ir: float = 8.0  # ln(L/a_c) ~ a few decades
    n_layers: int = N_LAYERS
    theta_stat: float = THETA_STAT
    k_c_bosonic: float = K_C_BOSONIC
    k_c_anyon: float = K_C_ANYON
    status: str = "diagnostic_closure"


LOCKS = VortexRGLocks()


def k_c_bosonic() -> float:
    return K_C_BOSONIC


def k_c_anyon(n: int = N_LAYERS) -> float:
    """Conjectural anyon-plasma shift: K_c=(2/π)(1-θ_stat/π)^2, θ_stat=π/n."""
    theta = PI / n
    return K_C_BOSONIC * (1.0 - theta / PI) ** 2


def K0_from_temperatures(T_BKT_bare: float, T_eff: float) -> float:
    """Bare K at α=0: K_0=(2/π)(T_BKT^bare / T_eff).  Diagnostic input."""
    if T_eff <= 0.0:
        return float("inf")
    return K_C_BOSONIC * (T_BKT_bare / T_eff)


def physical_delta(alpha: float, ell: float = ELL_CORRECTIONS) -> float | None:
    """Physical (upper) fold-branch δ_+(α).  Locked solver."""
    br = solve_fold_branches(float(alpha), ell)
    phys = br["phys"]
    return float(phys) if phys is not None else None


def d_delta_d_alpha(delta: float, alpha: float, ell: float = ELL_CORRECTIONS) -> float:
    """Implicit derivative of δ=1+α(ln δ+ℓ).  Diverges at the fold."""
    if delta <= 0.0:
        return float("nan")
    denom = 1.0 - alpha / delta
    if abs(denom) < 1e-14:
        return float("inf")
    return (math.log(delta) + ell) / denom


def stiffness_K(delta: float, K0: float, delta0: float = 1.0) -> float:
    """Diagnostic: K(α)=K_0 [δ_+(α)/δ_0]^2 with δ_0=1 at α=0."""
    return K0 * (delta / delta0) ** 2


def fugacity_y(delta: float, K0: float, y_uv: float, core_coeff: float) -> float:
    """Textbook fugacity y=exp(-E_core/T) with E_core=c κ, hence y=exp(-c K).

    A floor y_uv is applied only so the ODE is non-singular at large K.
    Do *not* pin y(α=0) to a T-independent O(10^{-2}): that inflates y near
    the fold by exp(K_0) and forces a spurious plasma.
    """
    K = stiffness_K(delta, K0)
    return float(max(y_uv, math.exp(-core_coeff * K)))


def kt_rhs(ell: float, z: np.ndarray) -> np.ndarray:
    """Textbook KT vector field on (u,y)=(K^{-1},y).  α frozen."""
    u, y = z
    u = max(u, 1e-12)
    K = 1.0 / u
    du = 4.0 * PI**3 * y * y
    dy = (2.0 - PI * K) * y
    return np.array([du, dy], dtype=float)


def beta_alpha_diagnostic(
    K: float,
    y: float,
    alpha: float,
    delta: float,
    K0: float,
    ell_gap: float = ELL_CORRECTIONS,
) -> float:
    """Chain-rule diagnostic β_α, not a microscopic beta function.

    Identify running K with K_id(α)=K_0 δ_+(α)^2 on the physical branch.
    KT gives dK/dℓ=-4π³ y² K².  Then
        β_α = (dK/dℓ) / (dK_id/dα).
    Both numerator and dK_id/dα are negative on the physical branch, so
    β_α>0: vortex suppression is read as drive toward the fold.
    Frozen-cavity Option A* has microscopic β_α=0; this closure is extra.
    """
    if delta is None or delta <= 0.0 or K0 <= 0.0:
        return 0.0
    ddel = d_delta_d_alpha(delta, alpha, ell_gap)
    dK_dalpha = 2.0 * K0 * delta * ddel
    if not math.isfinite(dK_dalpha) or abs(dK_dalpha) < 1e-18:
        return 0.0
    dK_dl = -4.0 * PI**3 * (y**2) * (K**2)
    beta = dK_dl / dK_dalpha
    # α cannot run past the fold in this diagnostic.
    if alpha >= 0.999 * fold_condition(ell_gap)[1]:
        return 0.0
    return float(max(beta, 0.0))


def coupled_rhs(ell: float, z: np.ndarray, K0: float, ell_gap: float) -> np.ndarray:
    """(u, y, α) with diagnostic β_α."""
    u, y, alpha = z
    u = max(float(u), 1e-12)
    y = max(float(y), 0.0)
    K = 1.0 / u
    du = 4.0 * PI**3 * y * y
    dy = (2.0 - PI * K) * y
    delta = physical_delta(alpha, ell_gap)
    if delta is None:
        ba = 0.0
    else:
        ba = beta_alpha_diagnostic(K, y, alpha, delta, K0, ell_gap)
    return np.array([du, dy, ba], dtype=float)


def _plasma_event_factory(y_cut: float = 0.4):
    def event(ell: float, z: np.ndarray) -> float:
        return y_cut - z[1]

    event.terminal = True  # type: ignore[attr-defined]
    event.direction = -1.0  # type: ignore[attr-defined]
    return event


def integrate_kt_fixed_alpha(
    K_init: float,
    y_init: float,
    ell_ir: float,
    *,
    y_cut: float = 0.4,
    n_eval: int = 241,
) -> dict[str, Any]:
    """Integrate textbook KT at frozen α.  Plasma if y hits y_cut."""
    u0 = 1.0 / max(K_init, 1e-12)
    y0 = max(y_init, 1e-16)
    ev = _plasma_event_factory(y_cut)
    sol = solve_ivp(
        kt_rhs,
        (0.0, ell_ir),
        np.array([u0, y0]),
        rtol=1e-7,
        atol=1e-9,
        dense_output=True,
        events=ev,
        max_step=0.05,
    )
    ell = np.linspace(0.0, min(sol.t[-1], ell_ir), n_eval)
    uv = sol.sol(ell)
    K = 1.0 / np.maximum(uv[0], 1e-12)
    y = uv[1]
    plasma = bool(sol.t_events and len(sol.t_events[0]) > 0) or (
        y[-1] > y_cut * 0.98
    )
    bound = (not plasma) and (K[-1] > K_C_BOSONIC) and (y[-1] < y0 * 2.0 or y[-1] < 0.05)
    return {
        "ell": ell,
        "K": K,
        "y": np.maximum(y, 0.0),
        "plasma": plasma,
        "bound": bool(bound),
        "ell_stop": float(sol.t[-1]),
        "K_ir": float(K[-1]),
        "y_ir": float(y[-1]),
    }


def integrate_coupled(
    K_init: float,
    y_init: float,
    alpha_init: float,
    K0: float,
    ell_ir: float,
    ell_gap: float = ELL_CORRECTIONS,
    *,
    y_cut: float = 0.4,
    n_eval: int = 241,
) -> dict[str, Any]:
    """Coupled (K,y,α) diagnostic.  Stops on plasma or fold."""
    u0 = 1.0 / max(K_init, 1e-12)
    y0 = max(y_init, 1e-16)
    delta_c, alpha_c = fold_condition(ell_gap)

    already_plasma = (y0 >= y_cut) or (K_init < K_C_BOSONIC and y0 > 0.05)
    if already_plasma:
        ell = np.linspace(0.0, min(0.05, ell_ir), 8)
        K = np.full_like(ell, K_init)
        y = np.full_like(ell, y0)
        alpha = np.full_like(ell, alpha_init)
        return {
            "ell": ell,
            "K": K,
            "y": y,
            "alpha": alpha,
            "hit_plasma": True,
            "hit_fold": False,
            "winner": "KT",
            "alpha_ir": float(alpha_init),
            "K_ir": float(K_init),
            "y_ir": float(y0),
            "alpha_c": float(alpha_c),
            "delta_c": float(delta_c),
            "already_plasma": True,
        }

    def rhs(ell: float, z: np.ndarray) -> np.ndarray:
        return coupled_rhs(ell, z, K0, ell_gap)

    def plasma(ell: float, z: np.ndarray) -> float:
        return y_cut - z[1]

    plasma.terminal = True  # type: ignore[attr-defined]
    plasma.direction = -1.0  # type: ignore[attr-defined]

    def fold_hit(ell: float, z: np.ndarray) -> float:
        return 0.999 * alpha_c - z[2]

    fold_hit.terminal = True  # type: ignore[attr-defined]
    fold_hit.direction = -1.0  # type: ignore[attr-defined]

    sol = solve_ivp(
        rhs,
        (0.0, ell_ir),
        np.array([u0, max(y_init, 1e-16), alpha_init]),
        rtol=1e-6,
        atol=1e-9,
        dense_output=True,
        events=(plasma, fold_hit),
        max_step=0.05,
    )
    t_end = min(float(sol.t[-1]), ell_ir)
    ell = np.linspace(0.0, t_end, n_eval)
    uv = sol.sol(ell)
    K = 1.0 / np.maximum(uv[0], 1e-12)
    y = np.maximum(uv[1], 0.0)
    alpha = uv[2]
    hit_plasma = bool(sol.t_events and len(sol.t_events[0]) > 0)
    hit_fold = bool(sol.t_events and len(sol.t_events[1]) > 0)
    winner = "undecided"
    if hit_plasma and not hit_fold:
        winner = "KT"
    elif hit_fold and not hit_plasma:
        winner = "fold"
    elif hit_plasma and hit_fold:
        t_p = sol.t_events[0][0]
        t_f = sol.t_events[1][0]
        winner = "KT" if t_p < t_f else "fold"
    elif alpha[-1] >= 0.995 * alpha_c:
        winner = "fold"
    elif y[-1] > y_cut * 0.9:
        winner = "KT"
    else:
        # Reached ℓ_IR still bound and below fold.
        winner = "fold_first_at_IR" if K[-1] > K_C_BOSONIC else "KT"
    return {
        "ell": ell,
        "K": K,
        "y": y,
        "alpha": alpha,
        "hit_plasma": hit_plasma,
        "hit_fold": hit_fold,
        "winner": winner,
        "alpha_ir": float(alpha[-1]),
        "K_ir": float(K[-1]),
        "y_ir": float(y[-1]),
        "alpha_c": float(alpha_c),
        "delta_c": float(delta_c),
    }


def competing_scale_scan(
    locks: VortexRGLocks = LOCKS,
    n_alpha: int = 80,
    use_anyon: bool = False,
) -> dict[str, Any]:
    """Scan the physical branch: does IR KT plasma occur before α_c?"""
    ell_gap = locks.ell_gap
    delta_c, alpha_c = fold_condition(ell_gap)
    K0 = K0_from_temperatures(locks.T_BKT_bare_K, locks.T_eff_K)
    Kc = locks.k_c_anyon if use_anyon else locks.k_c_bosonic
    alphas = np.linspace(1e-4, 0.995 * alpha_c, n_alpha)
    rows: list[dict[str, float | int | str]] = []
    alpha_kt: float | None = None
    for a in alphas:
        d = physical_delta(float(a), ell_gap)
        if d is None:
            continue
        K_b = stiffness_K(d, K0)
        y_b = fugacity_y(d, K0, locks.y_uv, locks.core_coeff)
        flow = integrate_kt_fixed_alpha(K_b, y_b, locks.ell_ir)
        # Bare (pre-RG) crossing of K_c:
        bare_below = K_b < Kc
        ir_plasma = bool(flow["plasma"])
        phase = "bound"
        if ir_plasma or bare_below:
            phase = "plasma"
        tA = (alpha_c - a) / alpha_c
        if phase == "plasma" and alpha_kt is None:
            alpha_kt = float(a)
        rows.append(
            {
                "alpha": float(a),
                "alpha_over_ac": float(a / alpha_c),
                "t_A": float(tA),
                "delta": float(d),
                "Delta_meV": float(d * locks.delta0_meV),
                "K_bare": float(K_b),
                "y_bare": float(y_b),
                "K_ir": float(flow["K_ir"]),
                "y_ir": float(flow["y_ir"]),
                "plasma": int(ir_plasma or bare_below),
                "phase": phase,
            }
        )
    fold_wins = alpha_kt is None
    # Also the bare-K crossing (no RG), for the criterion theorem.
    alpha_bare_kt: float | None = None
    delta_bare_kt: float | None = None
    for row in rows:
        if row["K_bare"] < Kc:  # type: ignore[operator]
            alpha_bare_kt = float(row["alpha"])  # type: ignore[arg-type]
            delta_bare_kt = float(row["delta"])  # type: ignore[arg-type]
            break
    return {
        "K0": K0,
        "Kc": Kc,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "alpha_kt_ir": alpha_kt,
        "alpha_kt_bare": alpha_bare_kt,
        "delta_kt_bare": delta_bare_kt,
        "fold_wins": fold_wins,
        "kt_before_fold": (not fold_wins)
        and (alpha_kt is not None)
        and (alpha_kt < alpha_c),
        "use_anyon": use_anyon,
        "T_eff_K": locks.T_eff_K,
        "T_BKT_bare_K": locks.T_BKT_bare_K,
        "rows": rows,
        "status": "computed_diagnostic",
    }


def T_eff_crossover_for_fold_win(
    T_BKT_bare: float = T_BKT_BARE_DEFAULT,
    delta_c: float | None = None,
) -> float:
    """T_eff below which even K(α_c)=K_0 δ_c² stays above 2/π (bare).

    K_fold = (2/π)(T_BKT/T_eff) δ_c² > 2/π  ⇒  T_eff < T_BKT δ_c².
    """
    dc = delta_c if delta_c is not None else fold_condition(ELL_CORRECTIONS)[0]
    return T_BKT_bare * (dc**2)


def eta_spinwave(K: float) -> float:
    """η=1/(2πK).  Equals 1/4 at K=2/π."""
    if K <= 0.0:
        return float("inf")
    return 1.0 / (2.0 * PI * K)


def g_r_spinwave(r: np.ndarray, eta: float, a: float = 1.0) -> np.ndarray:
    """g(r)~(a/r)^η, Villain spin-wave (no vortices)."""
    rr = np.maximum(np.asarray(r, dtype=float), a)
    return (a / rr) ** eta


def g_r_plasma(r: np.ndarray, xi: float, a: float = 1.0) -> np.ndarray:
    """Exponential plasma correlator (illustration)."""
    rr = np.maximum(np.asarray(r, dtype=float), a)
    return np.exp(-rr / xi) * (a / rr) ** 0.25


def eta_eff_finite_size(g_half: float, L: float, a: float = 1.0) -> float:
    """η_eff(L)=-ln g(L/2)/ln(L/(2a))."""
    if g_half <= 0.0 or L <= 2.0 * a:
        return float("nan")
    return float(-math.log(g_half) / math.log(L / (2.0 * a)))


def xy_metropolis_gL(
    L: int,
    K_lat: float,
    *,
    n_therm: int = 400,
    n_meas: int = 800,
    seed: int = 0,
) -> dict[str, float]:
    """Modest 2D XY Metropolis; program/computed, not experiment.

    Hamiltonian -K_lat ∑_<ij> cos(θ_i-θ_j).  Reports g(L/2) along a lattice
    axis and η_eff.  Equilibration is intentionally modest.
    """
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0.0, 2.0 * PI, size=(L, L))
    delta = 0.5 * PI

    def energy_change(i: int, j: int, dth: float) -> float:
        tn = theta[i, j] + dth
        s = 0.0
        for ni, nj in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
            s += math.cos(tn - theta[ni % L, nj % L]) - math.cos(
                theta[i, j] - theta[ni % L, nj % L]
            )
        return -K_lat * s

    def sweep() -> None:
        nonlocal delta
        acc = 0
        n = L * L
        for _ in range(n):
            i = int(rng.integers(0, L))
            j = int(rng.integers(0, L))
            dth = float(rng.uniform(-delta, delta))
            dE = energy_change(i, j, dth)
            if dE <= 0.0 or rng.random() < math.exp(-dE):
                theta[i, j] = (theta[i, j] + dth) % (2.0 * PI)
                acc += 1
        rate = acc / n
        if rate < 0.3:
            delta *= 0.9
        elif rate > 0.5:
            delta *= 1.1
        delta = min(max(delta, 0.05), PI)

    for _ in range(n_therm):
        sweep()

    g_acc = 0.0
    n_acc = 0
    half = L // 2
    for _ in range(n_meas):
        sweep()
        # g(r=L/2) along x, averaged over origins
        dth = theta[:, :] - np.roll(theta, shift=half, axis=1)
        g_acc += float(np.mean(np.cos(dth)))
        n_acc += 1
    g_half = g_acc / max(n_acc, 1)
    eta = eta_eff_finite_size(max(g_half, 1e-16), float(L))
    return {
        "L": float(L),
        "K_lat": float(K_lat),
        "g_half": float(g_half),
        "eta_eff": float(eta),
    }


def fold_vs_essential_curves(
    n: int = 200,
    C_bkt: float = 0.35,
    ell: float = ELL_CORRECTIONS,
) -> dict[str, np.ndarray]:
    """Diagnostic Δ(t_A): locked √ vs contingent essential singularity.

    Essential-singularity curve is a *target shape*, not a fit to data,
    and is not a QVC result.
    """
    delta_c, alpha_c = fold_condition(ell)
    tA = np.logspace(-3.0, 0.0, n)
    alpha = alpha_c * (1.0 - tA)
    d_phys = np.full_like(tA, np.nan)
    for i, a in enumerate(alpha):
        d = physical_delta(float(a), ell)
        if d is not None:
            d_phys[i] = d
    # Local √ form: δ-δ_c = B √t_A.  Amplitude from the smallest t_A with a root.
    finite = np.where(np.isfinite(d_phys))[0]
    if finite.size:
        i0 = int(finite[0])
        B = (d_phys[i0] - delta_c) / math.sqrt(max(tA[i0], 1e-16))
    else:
        B = 0.5
    d_sqrt = delta_c + B * np.sqrt(tA)
    d_bkt = np.exp(-C_bkt / np.sqrt(np.maximum(tA, 1e-12)))
    # Scale BKT toy so it is visible, not so it matches Δ_c (it vanishes).
    return {
        "t_A": tA,
        "delta_phys": d_phys,
        "delta_sqrt": d_sqrt,
        "delta_bkt": d_bkt,
        "delta_c": np.full_like(tA, delta_c),
        "B_sqrt": np.array([B]),
    }


def log_likelihood_scaling(
    tA: np.ndarray,
    delta: np.ndarray,
    delta_c: float,
) -> dict[str, float]:
    """C5 diagnostic: ½ ln t_A (fold) vs 1/√t_A (essential) on ln(δ-δ_c).

    Uses only points with δ>δ_c and t_A in (10^{-3}, 10^{-0.5}).
    """
    mask = (
        np.isfinite(delta)
        & np.isfinite(tA)
        & (delta > delta_c + 1e-6)
        & (tA > 1e-3)
        & (tA < 10 ** (-0.5))
    )
    t = tA[mask]
    y = np.log(np.maximum(delta[mask] - delta_c, 1e-16))
    if t.size < 8:
        return {"n": float(t.size), "ll_sqrt": float("nan"), "ll_bkt": float("nan")}
    # Fold: ln(δ-δ_c)=c0 + 0.5 ln t_A
    X1 = np.column_stack([np.ones_like(t), 0.5 * np.log(t)])
    # Force slope 1 on the second column by OLS intercept-only after subtracting
    coef1, *_ = np.linalg.lstsq(X1, y, rcond=None)
    r1 = y - X1 @ coef1
    # Free slope as a check (should be ~0.5 for a fold).
    Xf = np.column_stack([np.ones_like(t), np.log(t)])
    coef_free, *_ = np.linalg.lstsq(Xf, y, rcond=None)
    # Essential: ln δ = c0 - C/√t  (δ→0, so use ln δ not ln(δ-δ_c))
    yb = np.log(np.maximum(delta[mask], 1e-16))
    X2 = np.column_stack([np.ones_like(t), -1.0 / np.sqrt(t)])
    coef2, *_ = np.linalg.lstsq(X2, yb, rcond=None)
    r2 = yb - X2 @ coef2
    n = float(t.size)
    s1 = float(np.mean(r1**2))
    s2 = float(np.mean(r2**2))
    ll1 = -0.5 * n * math.log(2.0 * PI * max(s1, 1e-30)) - 0.5 * n
    ll2 = -0.5 * n * math.log(2.0 * PI * max(s2, 1e-30)) - 0.5 * n
    return {
        "n": n,
        "ll_sqrt": ll1,
        "ll_bkt": ll2,
        "delta_ll": ll1 - ll2,
        "slope_free": float(coef_free[1]),
        "C_bkt_fit": float(coef2[1]),
        "status": "computed_on_fold_branch",
    }


def write_dat(path: Path, header: str, columns: dict[str, np.ndarray]) -> None:
    """Write a PGFPlots table with one ``#`` header line (column names).

    Matches ``fold_phys.dat``: pgfplots uses ``skip first n=1`` because ``#``
    is not a TeX comment in that parser.
    """
    keys = list(columns.keys())
    arr = np.column_stack([np.asarray(columns[k], dtype=float) for k in keys])
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, arr, header="  ".join(keys), comments="# ")


def run_diagnostic(
    data_dir: Path,
    *,
    T_eff_K: float | None = None,
    T_BKT_bare_K: float | None = None,
    do_mc: bool = True,
) -> dict[str, Any]:
    """Run the full diagnostic suite and write PGFPlots tables."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    locks = VortexRGLocks(
        T_eff_K=LOCKS.T_star_K if T_eff_K is None else T_eff_K,
        T_BKT_bare_K=LOCKS.T_BKT_bare_K if T_BKT_bare_K is None else T_BKT_bare_K,
    )
    delta_c, alpha_c = fold_condition(locks.ell_gap)
    K0 = K0_from_temperatures(locks.T_BKT_bare_K, locks.T_eff_K)
    T_fold_win = T_eff_crossover_for_fold_win(locks.T_BKT_bare_K, delta_c)

    scan = competing_scale_scan(locks, n_alpha=100, use_anyon=False)
    scan_any = competing_scale_scan(locks, n_alpha=100, use_anyon=True)

    # Representative KT flows at three stations on the physical branch.
    stations = []
    for frac, name in ((0.30, "uv"), (0.90, "star"), (0.985, "nearfold")):
        a = frac * alpha_c
        d = physical_delta(a)
        assert d is not None
        Kb = stiffness_K(d, K0)
        yb = fugacity_y(d, K0, locks.y_uv, locks.core_coeff)
        fl = integrate_kt_fixed_alpha(Kb, yb, locks.ell_ir, n_eval=201)
        stations.append((name, a, d, fl))
        write_dat(
            data_dir / f"vortex_rg_flow_{name}.dat",
            f"fixed-alpha KT diagnostic  alpha/ac={frac:.3f}  status=computed",
            {"ell": fl["ell"], "K": fl["K"], "y": fl["y"], "Kinv": 1.0 / fl["K"]},
        )

    # Coupled flow starting at λ* working point (α/α_c=0.90).
    a_star = 0.90 * alpha_c
    d_star = physical_delta(a_star)
    assert d_star is not None
    coupled = integrate_coupled(
        stiffness_K(d_star, K0),
        fugacity_y(d_star, K0, locks.y_uv, locks.core_coeff),
        a_star,
        K0,
        locks.ell_ir,
        locks.ell_gap,
        n_eval=241,
    )
    write_dat(
        data_dir / "vortex_rg_coupled.dat",
        "coupled (K,y,alpha) diagnostic closure from lambda*  status=computed",
        {
            "ell": coupled["ell"],
            "K": coupled["K"],
            "y": coupled["y"],
            "alpha": coupled["alpha"],
            "alpha_over_ac": coupled["alpha"] / alpha_c,
        },
    )

    # Also coupled from the UV (small α) to see the global toy.
    a_uv = 0.20 * alpha_c
    d_uv = physical_delta(a_uv)
    assert d_uv is not None
    coupled_uv = integrate_coupled(
        stiffness_K(d_uv, K0),
        fugacity_y(d_uv, K0, locks.y_uv, locks.core_coeff),
        a_uv,
        K0,
        locks.ell_ir,
        locks.ell_gap,
        n_eval=241,
    )
    write_dat(
        data_dir / "vortex_rg_coupled_uv.dat",
        "coupled diagnostic from UV  status=computed",
        {
            "ell": coupled_uv["ell"],
            "K": coupled_uv["K"],
            "y": coupled_uv["y"],
            "alpha": coupled_uv["alpha"],
            "alpha_over_ac": coupled_uv["alpha"] / alpha_c,
        },
    )

    rows = scan["rows"]
    write_dat(
        data_dir / "competing_scales.dat",
        "physical-branch competing-scale scan  status=computed_diagnostic",
        {
            "alpha_over_ac": np.array([r["alpha_over_ac"] for r in rows]),
            "t_A": np.array([r["t_A"] for r in rows]),
            "delta": np.array([r["delta"] for r in rows]),
            "K_bare": np.array([r["K_bare"] for r in rows]),
            "y_bare": np.array([r["y_bare"] for r in rows]),
            "K_ir": np.array([r["K_ir"] for r in rows]),
            "y_ir": np.array([r["y_ir"] for r in rows]),
            "plasma": np.array([r["plasma"] for r in rows], dtype=float),
        },
    )
    write_dat(
        data_dir / "competing_scales_anyon.dat",
        "same scan with conjectural anyon K_c  status=conjecture",
        {
            "alpha_over_ac": np.array([r["alpha_over_ac"] for r in scan_any["rows"]]),
            "K_bare": np.array([r["K_bare"] for r in scan_any["rows"]]),
            "plasma": np.array([r["plasma"] for r in scan_any["rows"]], dtype=float),
        },
    )

    curves = fold_vs_essential_curves()
    write_dat(
        data_dir / "delta_vs_tA.dat",
        "fold branch vs contingent essential singularity  diagnostic only",
        {
            "t_A": curves["t_A"],
            "delta_phys": curves["delta_phys"],
            "delta_sqrt": curves["delta_sqrt"],
            "delta_bkt": curves["delta_bkt"],
            "delta_c": curves["delta_c"],
        },
    )
    ll = log_likelihood_scaling(
        curves["t_A"], curves["delta_phys"], float(curves["delta_c"][0])
    )

    # g(r): spin-wave at several K along the branch + plasma illustration.
    r = np.linspace(1.0, 64.0, 128)
    # Working point, KT, fold-adjacent:
    d_work = physical_delta(0.90 * alpha_c)
    d_nf = physical_delta(0.985 * alpha_c)
    K_work = stiffness_K(d_work, K0)  # type: ignore[arg-type]
    K_nf = stiffness_K(d_nf, K0)  # type: ignore[arg-type]
    eta_w = eta_spinwave(K_work)
    eta_c = 0.25
    eta_nf = eta_spinwave(max(K_nf, 1e-3))
    write_dat(
        data_dir / "g_r_scaling.dat",
        "Villain spin-wave g(r); plasma toy.  status=program_computed",
        {
            "r": r,
            "g_work": g_r_spinwave(r, eta_w),
            "g_kt": g_r_spinwave(r, eta_c),
            "g_nearfold_sw": g_r_spinwave(r, min(eta_nf, 2.0)),
            "g_plasma": g_r_plasma(r, xi=6.0),
        },
    )

    # η_eff vs α (spin-wave finite-size formula) for L=16,32,64.
    Ls = (16, 32, 64)
    eta_cols: dict[str, np.ndarray] = {
        "alpha_over_ac": np.array([r["alpha_over_ac"] for r in rows]),
        "K_bare": np.array([r["K_bare"] for r in rows]),
    }
    for L in Ls:
        etas = []
        for rrow in rows:
            K = float(rrow["K_bare"])
            eta = eta_spinwave(max(K, 1e-6))
            g_half = g_r_spinwave(np.array([L / 2.0]), eta)[0]
            # If plasma flagged, mix in exponential suppression.
            if int(rrow["plasma"]) == 1:
                g_half *= math.exp(-(L / 2.0) / 8.0)
            etas.append(eta_eff_finite_size(max(g_half, 1e-18), float(L)))
        eta_cols[f"eta_L{L}"] = np.array(etas)
    write_dat(
        data_dir / "eta_eff.dat",
        "eta_eff(L) along physical branch  spin-wave+plasma mix  program",
        eta_cols,
    )

    mc_rows: list[dict[str, float]] = []
    if do_mc:
        # L=16 only: L=32 was not equilibrated at these sweep counts.
        # Status: program/computed illustration, not experiment (C7 open).
        for Klat, seed in ((1.40, 1), (1.10, 2), (0.90, 3), (0.70, 4), (0.50, 5)):
            mc_rows.append(
                xy_metropolis_gL(16, Klat, n_therm=800, n_meas=1200, seed=seed)
            )
        if mc_rows:
            write_dat(
                data_dir / "g_r_mc.dat",
                "modest XY Metropolis  status=program_computed_not_experiment",
                {
                    "L": np.array([m["L"] for m in mc_rows]),
                    "K_lat": np.array([m["K_lat"] for m in mc_rows]),
                    "g_half": np.array([m["g_half"] for m in mc_rows]),
                    "eta_eff": np.array([m["eta_eff"] for m in mc_rows]),
                },
            )

    # T_eff scan: who wins?
    T_grid = np.array([0.02, 0.033, 0.05, 0.1187, 0.25, 0.50, 1.00])
    winners = []
    alpha_kts = []
    for T in T_grid:
        lk = VortexRGLocks(T_eff_K=float(T), T_BKT_bare_K=locks.T_BKT_bare_K)
        sc = competing_scale_scan(lk, n_alpha=60, use_anyon=False)
        winners.append(0.0 if sc["fold_wins"] else 1.0)
        # Prefer the IR-flow crossing; fall back to bare K=K_c.
        akt = sc["alpha_kt_ir"] if sc["alpha_kt_ir"] is not None else sc["alpha_kt_bare"]
        alpha_kts.append(float(akt) / alpha_c if akt is not None else math.nan)
    write_dat(
        data_dir / "Teff_scan.dat",
        "competing-scale winner vs T_eff  1=KT_before_fold  0=fold_wins",
        {
            "T_eff": T_grid,
            "kt_before_fold": np.array(winners),
            "alpha_kt_over_ac": np.array(alpha_kts),
        },
    )

    summary = {
        "status": "computed_diagnostic_not_a_BKT_theorem",
        "locks": asdict(locks),
        "K0": K0,
        "K_c_bosonic": K_C_BOSONIC,
        "K_c_anyon": K_C_ANYON,
        "theta_stat": THETA_STAT,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "T_eff_fold_win_K": T_fold_win,
        "path_scan_fold_wins": scan["fold_wins"],
        "path_scan_kt_before_fold": scan["kt_before_fold"],
        "alpha_kt_bare_over_ac": (
            None
            if scan["alpha_kt_bare"] is None
            else scan["alpha_kt_bare"] / alpha_c
        ),
        "delta_kt_bare": scan["delta_kt_bare"],
        "anyon_kt_before_fold": scan_any["kt_before_fold"],
        "alpha_kt_anyon_over_ac": (
            None
            if scan_any["alpha_kt_bare"] is None
            else scan_any["alpha_kt_bare"] / alpha_c
        ),
        "coupled_from_star_winner": coupled["winner"],
        "coupled_from_uv_winner": coupled_uv["winner"],
        "coupled_star_alpha_ir_over_ac": coupled["alpha_ir"] / alpha_c,
        "loglik": ll,
        "mc_points": len(mc_rows),
        "note": (
            "Failure of H-VRG1 does not invalidate the fold. "
            "A KT-before-fold toy hit is not a BKT theorem."
        ),
    }
    # Human-readable sidecar
    lines = [
        "QVC vortex-RG diagnostic summary",
        f"status: {summary['status']}",
        f"K0={K0:.4f}  K_c(bosonic)={K_C_BOSONIC:.4f}  K_c(anyon)={K_C_ANYON:.4f}",
        f"T_eff={locks.T_eff_K:.4f} K  T_BKT_bare={locks.T_BKT_bare_K:.4f} K",
        f"T_eff below which bare fold wins: {T_fold_win:.4f} K",
        f"path scan KT before fold (bosonic): {scan['kt_before_fold']}",
        f"alpha_kt_bare/alpha_c={summary['alpha_kt_bare_over_ac']}",
        f"delta_kt_bare={scan['delta_kt_bare']}",
        f"anyon KT before fold: {scan_any['kt_before_fold']}",
        f"coupled from lambda* winner: {coupled['winner']}",
        f"coupled from UV winner: {coupled_uv['winner']}",
        f"C5 delta_ll (fold minus BKT)={ll.get('delta_ll')}",
        f"C5 free slope on ln(delta-delta_c) vs ln t_A={ll.get('slope_free')}",
        summary["note"],
    ]
    (data_dir / "vortex_rg_summary.txt").write_text("\n".join(lines) + "\n")
    return summary


if __name__ == "__main__":
    import argparse
    import json

    p = argparse.ArgumentParser(description="QVC vortex-RG diagnostic")
    p.add_argument(
        "--data-dir",
        type=Path,
        default=Path("tex/QVC_arXiv_v4/data"),
    )
    p.add_argument("--no-mc", action="store_true")
    args = p.parse_args()
    out = run_diagnostic(args.data_dir, do_mc=not args.no_mc)
    print(json.dumps({k: v for k, v in out.items() if k != "locks"}, indent=2, default=str))
