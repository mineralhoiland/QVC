"""Coupled GPE + DCE: live gap from ψ feeds DCE; A(t) feeds V_ZPF.

Two timescales (documented, not faked as one run):
  1. Split-step GPE on ~1 ps (dt ~ 2 fs) with DCE RK4 at the GPE step.
     Δ_live(ψ) enters dA/dt; A(t)/A(0) scales V_ZPF.
  2. Longer DCE window (~80 ps) slaved to the *static fold* Δ(A).
     Justified because Γ_cav^{-1} ~ tens of ps ≫ 1 ps GPE window.

CS / Berry / Skyrme are not in either EOM.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc.tqgl.dce import DCEConfig, make_dce, run_dce_fold_slaved
from qvc.tqgl.gpe import (
    GPEGrid,
    hopfion_initial,
    kinetic_coeff_mev_um2,
    kinetic_propagator,
    l2_norm,
    make_grid,
    nl_potential,
    phonon_potential,
    split_step_advance,
    zpf_potential,
)
from qvc.tqgl.live_gap import coherence_g1, delta_live_from_psi, mean_density_weighted, mu_nl_from_psi
from qvc.tqgl.locks import FrozenTQGL, fold_delta_of_A


def _fold_delta_meV(A: float, lock: FrozenTQGL, table: tuple[np.ndarray, np.ndarray] | None = None) -> float | None:
    if table is not None:
        A_g, D_g = table
        if A >= lock.A_c_fold_um_inv * (1.0 - 1e-12):
            return None
        return float(np.interp(A, A_g, D_g))
    info = fold_delta_of_A(A, lock)
    return info["Delta_phys_meV"]


def _fold_table(lock: FrozenTQGL, n: int = 80) -> tuple[np.ndarray, np.ndarray]:
    A_c = lock.A_c_fold_um_inv
    A = np.linspace(0.0, A_c * 0.999, n)
    D = np.empty_like(A)
    for i, a in enumerate(A):
        info = fold_delta_of_A(float(a), lock)
        D[i] = np.nan if info["Delta_phys_meV"] is None else float(info["Delta_phys_meV"])
    return A, D


def run_coupled_gpe_dce(
    lock: FrozenTQGL,
    *,
    n_grid: int = 128,
    L_um: float = 0.200,
    dt_ps: float = 0.002,
    t_end_ps: float = 1.0,
    n0: float = 1.0,
    seed: int = 42,
    pump_mode: str = "resonant",
    dce_cfg: DCEConfig | None = None,
    snapshot_every: int | None = None,
) -> dict[str, Any]:
    """1 ps GPE with live-gap DCE. Returns timeseries + final |ψ|².

    ``snapshot_every`` (steps) additionally records (t, |ψ|², arg ψ) frames
    for the workbench replay; physics is unchanged.
    """
    grid: GPEGrid = make_grid(lock, n_grid=n_grid, L_um=L_um)
    ke = kinetic_coeff_mev_um2(lock)
    prop_K = kinetic_propagator(grid, dt_ps, ke, lock.hbar_meV_ps)
    psi = hopfion_initial(grid, lock, n0=n0, seed=seed)
    n_ref = mean_density_weighted(psi, grid.weight)
    if n_ref <= 0.0:
        raise RuntimeError("n_ref vanished on hopfion IC")

    alpha = -lock.Delta0_meV
    beta = lock.Delta0_meV / n0
    g_ph = 0.10 * lock.Delta0_meV
    omega_ph = 2.0 * lock.Delta0_meV / lock.hbar_meV_ps

    cfg = DCEConfig(pump_mode=pump_mode) if dce_cfg is None else dce_cfg
    dce = make_dce(lock, cfg)
    # 1 ps field starts at Δ_live(0)=Δ₀, not the static Δ★. Resonant pump
    # tracks the IC catalyzon; detuning keeps the +n_γ Γ offset from that.
    w_ic = 2.0 * lock.Delta0_meV / lock.hbar_meV_ps
    dce.Gamma_cav = w_ic / cfg.Q_moire
    if cfg.pump_mode == "resonant":
        dce.omega_pump_half = w_ic
    else:
        dce.omega_pump_half = w_ic + cfg.n_gamma_detune * dce.Gamma_cav
    A_ref = dce.A0

    n_steps = int(round(t_end_ps / dt_ps))
    t = np.arange(n_steps + 1, dtype=float) * dt_ps
    Delta_live = np.empty(n_steps + 1)
    Delta_fold = np.empty(n_steps + 1)
    g1 = np.empty(n_steps + 1)
    A = np.empty(n_steps + 1)
    L = np.empty(n_steps + 1)
    mu = np.empty(n_steps + 1)
    folded_live = np.zeros(n_steps + 1, dtype=bool)
    nrm0 = l2_norm(psi, grid.dx)
    fold_tab = _fold_table(lock)

    def record(i: int) -> None:
        A[i] = dce.A
        Delta_live[i] = delta_live_from_psi(
            psi, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
        )
        fd = _fold_delta_meV(dce.A, lock, fold_tab)
        Delta_fold[i] = np.nan if fd is None else fd
        g1[i] = coherence_g1(psi)
        L[i] = dce.lorentzian_now(Delta_live[i])
        mu[i] = mu_nl_from_psi(psi, weight=grid.weight, alpha=alpha, beta=beta)
        folded_live[i] = bool(Delta_live[i] < lock.Delta_c_meV)

    record(0)
    psi0_density = np.abs(psi) ** 2
    snapshots: list[dict[str, Any]] = []
    if snapshot_every:
        snapshots.append({"t_ps": 0.0, "density": psi0_density.copy(), "phase": np.angle(psi)})

    for step in range(n_steps):
        t_now = t[step]
        V = (
            nl_potential(psi, alpha, beta)
            + zpf_potential(dce.A, A_ref, grid.gauss_ring, lock.E_vac_meV)
            + phonon_potential(t_now, grid.gauss_ring, g_ph=g_ph, omega_ps=omega_ph)
        )
        psi = split_step_advance(
            psi, prop_K=prop_K, V=V, dt_ps=dt_ps, hbar=lock.hbar_meV_ps
        )
        d_live = delta_live_from_psi(
            psi, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
        )
        dce.step(dt_ps, d_live)
        record(step + 1)
        if snapshot_every and ((step + 1) % snapshot_every == 0 or step + 1 == n_steps):
            snapshots.append({"t_ps": float(t[step + 1]), "density": np.abs(psi) ** 2, "phase": np.angle(psi)})

    nrm1 = l2_norm(psi, grid.dx)
    live_below_c = bool(np.any(folded_live))
    return {
        "pump_mode": pump_mode,
        "n_grid": n_grid,
        "L_um": L_um,
        "dt_ps": dt_ps,
        "t_end_ps": t_end_ps,
        "n_ref": n_ref,
        "n0": n0,
        "alpha_nl": alpha,
        "beta_nl": beta,
        "g_ph": g_ph,
        "omega_ph_ps_inv": omega_ph,
        "ke_coef_meV_um2": ke,
        "A0_um_inv": A_ref,
        "t_ps": t,
        "Delta_live_meV": Delta_live,
        "Delta_fold_of_A_meV": Delta_fold,
        "g1": g1,
        "A_um_inv": A,
        "lorentzian": L,
        "mu_nl_meV": mu,
        "live_below_delta_c": live_below_c,
        "any_live_below_delta_c": live_below_c,
        "Delta_live_final_meV": float(Delta_live[-1]),
        "Delta_live_min_meV": float(np.min(Delta_live)),
        "g1_final": float(g1[-1]),
        "A_final_um_inv": float(A[-1]),
        "A_final_over_Ac_fold": float(A[-1] / lock.A_c_fold_um_inv),
        "norm_rel_drift": abs(nrm1 - nrm0) / nrm0 if nrm0 else float("nan"),
        "density_final": np.abs(psi) ** 2,
        "density_initial": psi0_density,
        "snapshots": snapshots,
        "x_um": grid.x,
        "y_um": grid.y,
        "truncation": (
            "EOM is kinetic + cubic-quintic + V_ZPF(A(t)) + V_phonon. "
            "Chern–Simons, Berry, and Skyrme/Hopf terms are not discretised."
        ),
        "live_gap_uses_psi": True,
        "live_gap_has_no_time_argument": True,
        "V_ZPF_uses_running_A": True,
        "drive_is_E_vac_not_0p244": True,
    }


def run_two_timescale_suite(
    lock: FrozenTQGL,
    *,
    n_grid: int = 128,
    t_gpe_ps: float = 1.0,
    dt_gpe_ps: float = 0.002,
    n_berry: int = 2000,
) -> dict[str, Any]:
    """GPE+DCE (resonant) plus fold-slaved DCE (resonant and detuning)."""
    from qvc.tqgl.diagnostics import berry_goe_gue, schwinger_twin_channel

    coupled = run_coupled_gpe_dce(
        lock, n_grid=n_grid, t_end_ps=t_gpe_ps, dt_ps=dt_gpe_ps, pump_mode="resonant"
    )
    dce_res = run_dce_fold_slaved(lock, DCEConfig(pump_mode="resonant"))
    dce_det = run_dce_fold_slaved(lock, DCEConfig(pump_mode="detuning"))
    sch = schwinger_twin_channel(lock)
    berry = berry_goe_gue(n_samples=n_berry, g_bar=abs(lock.lambda_star * lock.E_vac_meV))
    return {
        "lock": lock.as_dict(),
        "coupled_gpe_dce": coupled,
        "dce_fold_resonant": dce_res,
        "dce_fold_detuning": dce_det,
        "schwinger": sch,
        "berry": berry,
        "two_timescale": {
            "gpe_window_ps": t_gpe_ps,
            "dce_window_ps": dce_res["t_ps"][-1] if len(dce_res["t_ps"]) else None,
            "justification": (
                "GPE resolves fs–ps condensate response. DCE cavity time is "
                "Γ_cav^{-1} ~ Q ℏ / (2Δ★) ~ tens of ps, so the 80 ps window "
                "slaves Δ to the static fold branch Δ(A) rather than advancing "
                "the 128² GPE for 4e4 extra steps."
            ),
        },
    }
