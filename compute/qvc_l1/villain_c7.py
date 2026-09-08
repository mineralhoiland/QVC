"""Villain reduction (algebra) and C7 vectorized XY g(r).

Polar substitution Ψ = ρ e^{iθ} in F_bulk with ρ frozen at the saddle
produces an XY stiffness κ. The Villain approximation replaces

    exp(K cos φ)  →  Σ_{n∈Z} exp( −(K_V/2) (φ − 2π n)² )

with K_V ∼ K at large K. Poisson summation dualizes integer link
variables to a height field; the curl is the integer vortex gas

    H/T = 2π K Σ_{i<j} m_i m_j ln(|r_i−r_j|/a_c) + (E_core/T) Σ m_i²

with Kosterlitz K = κ/T_eff and K_c = 2/π. Defects are 2D phase vortices,
not hopfions and not CS fluxons.

C7 measures g(r) = ⟨e^{i[θ(r)−θ(0)]}⟩ on an L×L XY torus by checkerboard
Metropolis. Spin-wave theory gives η_sw = 1/(2π K); at the jump K=2/π
one has η=1/4 in the thermodynamic limit.

Status: derived as effective theory; computed as finite-L MC.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_bkt.kt_flow import K_C_BOSONIC


def villain_reduction_card(*, K: float, T_eff: float, kappa: float) -> dict[str, Any]:
    """Algebra of the Villain dualization. No Monte Carlo here."""
    return {
        "status": "derived_effective_ansatz_for_TEVC",
        "polar_split": "Ψ = ρ e^{iθ}",
        "Villain_kernel": "exp(K cos φ) → Σ_n exp(−(K_V/2)(φ−2π n)²)",
        "Poisson_dual": "integer link currents → height; m = curl n",
        "Coulomb_gas": "H/T = 2π K Σ_{i<j} m_i m_j ln(|r_i−r_j|/a) + E_core/T Σ m²",
        "K": K,
        "K_definition": "K = κ / T_eff (Kosterlitz, not José–Kadanoff πκ/T)",
        "K_c": K_C_BOSONIC,
        "eta_spinwave": 1.0 / (2.0 * math.pi * K) if K > 0 else float("inf"),
        "eta_KT_jump": 0.25,
        "kappa_meV": kappa,
        "T_eff_meV": T_eff,
        "universal_jump": "κ(T_BKT^-) = 2 T / π",
        "defect": "2D phase vortex, not hopfion, not CS fluxon",
        "note": (
            "Villain dualization is textbook once the phase stiffness exists. "
            "The TEVC identification of κ is C1, still partly ansatz."
        ),
    }


def _checkerboard_sweep(th: np.ndarray, K: float, rng: np.random.Generator, delta: float) -> None:
    """Simultaneous Metropolis on one sublattice (neighbors are the other)."""
    L = th.shape[0]
    i = np.arange(L)[:, None]
    j = np.arange(L)[None, :]
    for parity in (0, 1):
        mask = ((i + j) % 2) == parity
        dth = rng.uniform(-delta, delta, size=th.shape)
        new = np.where(mask, th + dth, th)

        def bonds(center: np.ndarray) -> np.ndarray:
            return (
                np.cos(center - np.roll(th, 1, axis=0))
                + np.cos(center - np.roll(th, -1, axis=0))
                + np.cos(center - np.roll(th, 1, axis=1))
                + np.cos(center - np.roll(th, -1, axis=1))
            )

        de = -K * (bonds(new) - bonds(th))
        accept = ((de <= 0.0) | (rng.random(th.shape) < np.exp(-np.clip(de, -50.0, 50.0)))) & mask
        th[:] = np.where(accept, new, th)


def _g_of_r(th: np.ndarray) -> np.ndarray:
    L = th.shape[0]
    g = np.empty(L // 2 + 1, dtype=float)
    for r in range(L // 2 + 1):
        g[r] = float(
            0.5
            * (
                np.mean(np.cos(th - np.roll(th, r, axis=1)))
                + np.mean(np.cos(th - np.roll(th, r, axis=0)))
            )
        )
    return g


def _eta_log_fit(g: np.ndarray, rmin: int = 2) -> float:
    """η from g(r) ∼ r^{−η} on r ∈ [rmin, L/4]."""
    Lh = g.size - 1
    rmax = max(rmin + 1, Lh // 2)
    r = np.arange(g.size, dtype=float)
    sl = slice(rmin, rmax + 1)
    gp = np.clip(g[sl], 1e-12, None)
    rr = np.clip(r[sl], 1.0, None)
    x = np.log(rr)
    y = np.log(gp)
    varx = float(np.var(x))
    if varx < 1e-18:
        return float("nan")
    slope = float(np.cov(x, y, bias=True)[0, 1] / varx)
    return float(-slope)


def xy_metropolis_g(
    *,
    L: int = 24,
    K: float = 1.0,
    n_therm: int = 600,
    n_samp: int = 1000,
    seed: int = 0,
    delta: float = 1.2,
) -> dict[str, Any]:
    """Vectorized checkerboard Metropolis XY on a torus."""
    rng = np.random.default_rng(seed)
    th = rng.uniform(0.0, 2.0 * math.pi, size=(L, L))
    for _ in range(n_therm):
        _checkerboard_sweep(th, K, rng, delta)
    acc = np.zeros(L // 2 + 1)
    e_acc = 0.0
    for _ in range(n_samp):
        _checkerboard_sweep(th, K, rng, delta)
        acc += _g_of_r(th)
        e_acc += float(
            np.mean(1.0 - np.cos(th - np.roll(th, 1, axis=0)))
            + np.mean(1.0 - np.cos(th - np.roll(th, 1, axis=1)))
        )
    g = acc / n_samp
    r = np.arange(L // 2 + 1)
    g_half = max(float(g[L // 2]), 1e-12)
    eta_half = -math.log(g_half) / math.log(L / 2.0)
    eta_fit = _eta_log_fit(g)
    eta_sw = 1.0 / (2.0 * math.pi * K)
    return {
        "status": "computed_checkerboard_xy",
        "L": L,
        "K": K,
        "n_therm": n_therm,
        "n_samp": n_samp,
        "g_r": g.tolist(),
        "r": r.tolist(),
        "eta_eff_from_g_Lhalf": eta_half,
        "eta_log_fit": eta_fit,
        "eta_spinwave": eta_sw,
        "eta_KT_target": 0.25,
        "energy_per_bond": e_acc / n_samp,
        "equilibrated_villain": False,
        "note": (
            "Checkerboard Metropolis XY, not Villain dual. Finite L. "
            "η→1/4 at K=2/π is the thermodynamic target."
        ),
    }


def eta_scan(
    *,
    L: int = 24,
    Ks: tuple[float, ...] | None = None,
    n_therm: int = 500,
    n_samp: int = 800,
) -> dict[str, Any]:
    if Ks is None:
        Ks = (2.4, 1.6, 1.0, K_C_BOSONIC, 0.45)
    rows = []
    for i, K in enumerate(Ks):
        mc = xy_metropolis_g(
            L=L, K=float(K), n_therm=n_therm, n_samp=n_samp, seed=10 + i
        )
        rows.append(
            {
                "K": float(K),
                "eta_log_fit": mc["eta_log_fit"],
                "eta_Lhalf": mc["eta_eff_from_g_Lhalf"],
                "eta_spinwave": mc["eta_spinwave"],
                "g_r": mc["g_r"],
                "r": mc["r"],
            }
        )
    return {"L": L, "rows": rows, "Ks": [float(k) for k in Ks]}


def run_villain_c7(
    *,
    kappa_meV: float,
    T_eff_meV: float,
    L: int = 24,
) -> dict[str, Any]:
    K = kappa_meV / max(T_eff_meV, 1e-12)
    card = villain_reduction_card(K=K, T_eff=T_eff_meV, kappa=kappa_meV)
    K_uv_mc = min(max(K, 1.2), 3.0)
    # Physical K at T* is typically ≫ 2, so UV MC uses a capped XY coupling.
    uv = xy_metropolis_g(L=L, K=K_uv_mc, n_therm=500, n_samp=800, seed=1)
    kt = xy_metropolis_g(L=L, K=K_C_BOSONIC, n_therm=600, n_samp=1000, seed=2)
    scan = eta_scan(L=L, n_therm=400, n_samp=600)
    return {
        "villain_card": card,
        "mc_uv": uv,
        "mc_near_KT": kt,
        "eta_scan": scan,
        "K_physical": K,
        "K_uv_capped_for_MC": K_uv_mc,
        "K_c": K_C_BOSONIC,
    }
