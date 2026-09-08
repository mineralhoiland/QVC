"""B2: number-damping and dephasing on the truncated split-step GPE.

The EOM remains the truncated 2D GPE (kinetic + cubic-quintic + V_ZPF + V_ph).
Chern–Simons, Berry, and Skyrme/Hopf terms are not discretised.

Lindblad channels (phenomenological, applied after each unitary split-step):

- Number damping:  ψ ← ψ exp(−γ_n Δt / 2)     [no-jump loss; N(t)∼e^{−γ_n t}]
- Dephasing:       ψ ← ψ exp(i √(γ_φ Δt) ξ)   [local phase diffusion; ξ∼N(0,1)]

γ = 1/τ_dec. Scan τ_dec ∈ {5, 20, 50} ps versus the coherent (γ=0) baseline.

Live gap uses ψ: Δ_live = Δ0 √(n_w / n_ref). There is no time schedule and
no restored 0.600→0.389→0.032 / 92% folklore.

1 ps window keeps coupled DCE so the coherent run is comparable to the
published Δ_live≈0.6305 meV, g^{(1)}≈0.8832 (n_grid=128). A 20 ps window
is also run (V_ZPF frozen at A(0), no DCE) because γ τ ≲ 1 ps is a small
fraction of τ_dec.
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.tqgl.dce import DCEConfig, make_dce  # noqa: E402
from qvc.tqgl.gpe import (  # noqa: E402
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
from qvc.tqgl.live_gap import (  # noqa: E402
    coherence_g1,
    delta_live_from_psi,
    mean_density_weighted,
    mu_nl_from_psi,
)
from qvc.tqgl.locks import FrozenTQGL  # noqa: E402
from qvc_l1.locks import (  # noqa: E402
    COHERENT_DELTA_LIVE_1PS_MEV,
    COHERENT_G1_1PS,
    L1Locks,
    build_l1_locks,
)

TAU_DEC_PS = (5.0, 20.0, 50.0)


def apply_lindblad(
    psi: np.ndarray,
    *,
    dt_ps: float,
    gamma_n: float,
    gamma_phi: float,
    rng: np.random.Generator | None,
) -> np.ndarray:
    """Number damping and/or local dephasing. Does not restore a gap schedule."""
    out = psi
    if gamma_n > 0.0:
        out = out * np.exp(-0.5 * gamma_n * dt_ps)
    if gamma_phi > 0.0:
        if rng is None:
            raise ValueError("dephasing requires an RNG")
        xi = rng.standard_normal(out.shape)
        out = out * np.exp(1j * np.sqrt(gamma_phi * dt_ps) * xi)
    return out


def run_gpe_lindblad(
    lock: FrozenTQGL,
    *,
    n_grid: int = 128,
    L_um: float = 0.200,
    dt_ps: float = 0.002,
    t_end_ps: float = 1.0,
    n0: float = 1.0,
    seed: int = 42,
    gamma_n: float = 0.0,
    gamma_phi: float = 0.0,
    couple_dce: bool = True,
    pump_mode: str = "resonant",
    record_stride: int = 1,
    dephasing_seed: int | None = None,
) -> dict[str, Any]:
    """Truncated GPE + optional DCE + Lindblad. Live gap from ψ."""
    grid = make_grid(lock, n_grid=n_grid, L_um=L_um)
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

    dce = None
    A_ref = 1.0
    A_frozen = 1.0
    if couple_dce:
        cfg = DCEConfig(pump_mode=pump_mode)
        dce = make_dce(lock, cfg)
        w_ic = 2.0 * lock.Delta0_meV / lock.hbar_meV_ps
        dce.Gamma_cav = w_ic / cfg.Q_moire
        if cfg.pump_mode == "resonant":
            dce.omega_pump_half = w_ic
        else:
            dce.omega_pump_half = w_ic + cfg.n_gamma_detune * dce.Gamma_cav
        A_ref = dce.A0
        A_frozen = dce.A0

    rng = None
    if gamma_phi > 0.0:
        rng = np.random.default_rng(seed if dephasing_seed is None else dephasing_seed)

    n_steps = int(round(t_end_ps / dt_ps))
    rec_idx = np.arange(0, n_steps + 1, record_stride, dtype=int)
    if rec_idx[-1] != n_steps:
        rec_idx = np.append(rec_idx, n_steps)
    n_rec = rec_idx.size
    t = np.empty(n_rec)
    Delta_live = np.empty(n_rec)
    g1 = np.empty(n_rec)
    nrm = np.empty(n_rec)
    A_arr = np.empty(n_rec)

    nrm0 = l2_norm(psi, grid.dx)
    rec_set = {int(i): k for k, i in enumerate(rec_idx)}

    def record(step: int) -> None:
        k = rec_set[step]
        t[k] = step * dt_ps
        Delta_live[k] = delta_live_from_psi(
            psi, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
        )
        g1[k] = coherence_g1(psi)
        nrm[k] = l2_norm(psi, grid.dx)
        A_arr[k] = dce.A if dce is not None else A_frozen

    record(0)

    for step in range(n_steps):
        t_now = step * dt_ps
        A_now = dce.A if dce is not None else A_frozen
        V = (
            nl_potential(psi, alpha, beta)
            + zpf_potential(A_now, A_ref, grid.gauss_ring, lock.E_vac_meV)
            + phonon_potential(t_now, grid.gauss_ring, g_ph=g_ph, omega_ps=omega_ph)
        )
        psi = split_step_advance(
            psi, prop_K=prop_K, V=V, dt_ps=dt_ps, hbar=lock.hbar_meV_ps
        )
        if gamma_n > 0.0 or gamma_phi > 0.0:
            psi = apply_lindblad(
                psi, dt_ps=dt_ps, gamma_n=gamma_n, gamma_phi=gamma_phi, rng=rng
            )
        if dce is not None:
            d_live = delta_live_from_psi(
                psi, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
            )
            dce.step(dt_ps, d_live)
        if (step + 1) in rec_set:
            record(step + 1)

    nrm1 = l2_norm(psi, grid.dx)
    i_1ps = int(np.argmin(np.abs(t - 1.0))) if t[-1] >= 1.0 - 0.5 * dt_ps else -1
    return {
        "n_grid": n_grid,
        "L_um": L_um,
        "dt_ps": dt_ps,
        "t_end_ps": t_end_ps,
        "couple_dce": couple_dce,
        "gamma_n_ps_inv": gamma_n,
        "gamma_phi_ps_inv": gamma_phi,
        "tau_n_ps": (1.0 / gamma_n) if gamma_n > 0.0 else None,
        "tau_phi_ps": (1.0 / gamma_phi) if gamma_phi > 0.0 else None,
        "n_ref": n_ref,
        "t_ps": t,
        "Delta_live_meV": Delta_live,
        "g1": g1,
        "norm": nrm,
        "A_um_inv": A_arr,
        "Delta_live_final_meV": float(Delta_live[-1]),
        "g1_final": float(g1[-1]),
        "Delta_live_1ps_meV": float(Delta_live[i_1ps]) if i_1ps >= 0 else None,
        "g1_1ps": float(g1[i_1ps]) if i_1ps >= 0 else None,
        "norm_rel_change": abs(nrm1 - nrm0) / nrm0 if nrm0 else float("nan"),
        "norm_ratio_final": nrm1 / nrm0 if nrm0 else float("nan"),
        "live_gap_uses_psi": True,
        "live_gap_has_no_time_argument": True,
        "no_prescribed_Delta_t": True,
        "truncation": (
            "EOM is kinetic + cubic-quintic + V_ZPF + V_phonon + Lindblad. "
            "Chern–Simons, Berry, and Skyrme/Hopf terms are not discretised."
        ),
        "channel": (
            "coherent"
            if gamma_n == 0.0 and gamma_phi == 0.0
            else (
                "number+dephasing"
                if gamma_n > 0.0 and gamma_phi > 0.0
                else ("number" if gamma_n > 0.0 else "dephasing")
            )
        ),
    }


def _at_time(run: dict[str, Any], t_ps: float) -> dict[str, float | None]:
    t = np.asarray(run["t_ps"])
    if t.size == 0 or t[-1] < t_ps - 1e-9:
        return {"t_ps": t_ps, "Delta_live_meV": None, "g1": None}
    i = int(np.argmin(np.abs(t - t_ps)))
    return {
        "t_ps": float(t[i]),
        "Delta_live_meV": float(run["Delta_live_meV"][i]),
        "g1": float(run["g1"][i]),
    }


def _row(
    *,
    label: str,
    tau_dec_ps: float | None,
    channel: str,
    run: dict[str, Any],
    t_report_ps: float,
) -> dict[str, Any]:
    pt = _at_time(run, t_report_ps)
    end = _at_time(run, float(run["t_end_ps"]))
    return {
        "label": label,
        "channel": channel,
        "tau_dec_ps": tau_dec_ps,
        "t_report_ps": t_report_ps,
        "Delta_live_meV": pt["Delta_live_meV"],
        "g1": pt["g1"],
        "Delta_live_end_meV": end["Delta_live_meV"],
        "g1_end": end["g1"],
        "n_grid": run["n_grid"],
        "couple_dce": run["couple_dce"],
        "norm_ratio_final": run["norm_ratio_final"],
    }


def run_lindblad_scan(
    locks: L1Locks | None = None,
    *,
    tau_dec_ps: tuple[float, ...] = TAU_DEC_PS,
    n_grid_1ps: int = 128,
    n_grid_20ps: int = 96,
    dt_ps: float = 0.002,
    t_short_ps: float = 1.0,
    t_long_ps: float = 20.0,
    include_channel_split: bool = True,
    seed: int = 42,
) -> dict[str, Any]:
    """Coherent baseline + τ_dec scan. Combined channels; optional split."""
    if locks is None:
        locks = build_l1_locks()
    lock = locks.tqgl

    def _short(**kw: Any) -> dict[str, Any]:
        return run_gpe_lindblad(
            lock,
            n_grid=n_grid_1ps,
            dt_ps=dt_ps,
            t_end_ps=t_short_ps,
            couple_dce=True,
            seed=seed,
            **kw,
        )

    def _long(**kw: Any) -> dict[str, Any]:
        stride = max(1, int(round(0.05 / dt_ps)))
        return run_gpe_lindblad(
            lock,
            n_grid=n_grid_20ps,
            dt_ps=dt_ps,
            t_end_ps=t_long_ps,
            couple_dce=False,
            seed=seed,
            record_stride=stride,
            **kw,
        )

    runs_1ps: dict[str, dict[str, Any]] = {}
    runs_20ps: dict[str, dict[str, Any]] = {}

    jobs_1: dict[str, tuple] = {"coherent": (0.0, 0.0, seed)}
    jobs_20: dict[str, tuple] = {"coherent": (0.0, 0.0, seed)}
    for tau in tau_dec_ps:
        g = 1.0 / tau
        jobs_1[f"tau_{tau:g}"] = (g, g, seed + int(10 * tau))
        jobs_20[f"tau_{tau:g}"] = (g, g, seed + int(10 * tau))
    if include_channel_split and tau_dec_ps:
        tau_s = float(tau_dec_ps[0])
        g = 1.0 / tau_s
        jobs_1[f"number_tau_{tau_s:g}"] = (g, 0.0, seed)
        jobs_1[f"dephasing_tau_{tau_s:g}"] = (0.0, g, seed + 7)

    def _launch_short(item: tuple[str, tuple]) -> tuple[str, dict[str, Any]]:
        name, (gn, gp, sd) = item
        return name, _short(gamma_n=gn, gamma_phi=gp, dephasing_seed=sd)

    def _launch_long(item: tuple[str, tuple]) -> tuple[str, dict[str, Any]]:
        name, (gn, gp, sd) = item
        return name, _long(gamma_n=gn, gamma_phi=gp, dephasing_seed=sd)

    n_workers = min(4, max(1, len(jobs_1) + len(jobs_20)))
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        short_f = [pool.submit(_launch_short, item) for item in jobs_1.items()]
        long_f = [pool.submit(_launch_long, item) for item in jobs_20.items()]
        for fut in short_f:
            name, run = fut.result()
            runs_1ps[name] = run
        for fut in long_f:
            name, run = fut.result()
            runs_20ps[name] = run

    table_1ps = [
        _row(
            label="coherent",
            tau_dec_ps=None,
            channel="coherent",
            run=runs_1ps["coherent"],
            t_report_ps=t_short_ps,
        )
    ]
    table_20ps = [
        _row(
            label="coherent",
            tau_dec_ps=None,
            channel="coherent",
            run=runs_20ps["coherent"],
            t_report_ps=t_long_ps,
        )
    ]
    for tau in tau_dec_ps:
        key = f"tau_{tau:g}"
        table_1ps.append(
            _row(
                label=key,
                tau_dec_ps=tau,
                channel="number+dephasing",
                run=runs_1ps[key],
                t_report_ps=t_short_ps,
            )
        )
        table_20ps.append(
            _row(
                label=key,
                tau_dec_ps=tau,
                channel="number+dephasing",
                run=runs_20ps[key],
                t_report_ps=t_long_ps,
            )
        )

    split_rows: list[dict[str, Any]] = []
    if include_channel_split and tau_dec_ps:
        tau_s = float(tau_dec_ps[0])
        split_rows = [
            _row(
                label=f"number_tau_{tau_s:g}",
                tau_dec_ps=tau_s,
                channel="number",
                run=runs_1ps[f"number_tau_{tau_s:g}"],
                t_report_ps=t_short_ps,
            ),
            _row(
                label=f"dephasing_tau_{tau_s:g}",
                tau_dec_ps=tau_s,
                channel="dephasing",
                run=runs_1ps[f"dephasing_tau_{tau_s:g}"],
                t_report_ps=t_short_ps,
            ),
        ]

    coh_1 = table_1ps[0]
    return {
        "status": "PROGRAM_COMPUTED_truncated_GPE_Lindblad",
        "truncation": (
            "CS / Berry / Skyrme not in the EOM. No prescribed Δ(t). "
            "Do not quote 92% or 0.032/0.389 as results."
        ),
        "published_coherent_baseline": {
            "Delta_live_1ps_meV": COHERENT_DELTA_LIVE_1PS_MEV,
            "g1_1ps": COHERENT_G1_1PS,
            "note": "n_grid=128 coupled GPE+DCE from output/tqgl_v4; comparison target, not a schedule.",
        },
        "this_run_coherent_1ps": {
            "Delta_live_meV": coh_1["Delta_live_meV"],
            "g1": coh_1["g1"],
            "n_grid": n_grid_1ps,
        },
        "tau_dec_ps": list(tau_dec_ps),
        "table_1ps": table_1ps,
        "table_20ps": table_20ps,
        "channel_split_1ps": split_rows,
        "runs_1ps": runs_1ps,
        "runs_20ps": runs_20ps,
        "windows": {
            "short_ps": t_short_ps,
            "long_ps": t_long_ps,
            "short_couple_dce": True,
            "long_couple_dce": False,
            "long_note": (
                "20 ps window freezes V_ZPF at A(0) (no DCE) so the scan isolates γ. "
                "1 ps window keeps coupled DCE to match the coherent baseline."
            ),
        },
        "computed_vs_ansatz": {
            "computed": [
                "Δ_live(t) from ring-weighted |ψ|²",
                "g^{(1)}(t) = |⟨ψ⟩|² / ⟨|ψ|²⟩",
                "split-step truncated GPE with Lindblad after each step",
            ],
            "phenomenological": [
                "γ_n = γ_φ = 1/τ_dec (equal rates in the combined scan)",
                "dephasing as local Wiener phase kicks (one unraveling)",
            ],
            "not": [
                "restored 92% gap collapse",
                "prescribed Δ(t) schedule 0.600→0.389→0.032",
                "CS / Berry / Skyrme in the EOM",
            ],
        },
    }


def run_b2(
    locks: L1Locks | None = None,
    **kwargs: Any,
) -> dict[str, Any]:
    return run_lindblad_scan(locks, **kwargs)
