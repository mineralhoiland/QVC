"""Honest TQGL v4 suite tests. Fast grids; no folklore-as-result."""

from __future__ import annotations

import inspect
import math

import numpy as np

from qvc.tqgl.dce import DCEConfig, make_dce, run_dce_fold_slaved
from qvc.tqgl.gpe import (
    hopfion_initial,
    kinetic_coeff_mev_um2,
    kinetic_propagator,
    l2_norm,
    make_grid,
    nl_potential,
    split_step_advance,
)
from qvc.tqgl.live_gap import coherence_g1, delta_live_from_psi, mean_density_weighted
from qvc.tqgl.locks import (
    DELTA_1PS_FOLKLORE_MEV,
    DELTA_QVC_FOLKLORE_MEV,
    E_CAS_WITHDRAWN_MEV,
    build_frozen_tqgl,
)


def test_locks_option_a_star():
    lock = build_frozen_tqgl(lambda_star_frac=0.90)
    assert abs(lock.r_nm - 8.0) < 1e-9
    assert abs(lock.R_nm - 50.0) < 1e-9
    assert abs(lock.E_vac_meV - 0.1287) < 5e-4
    assert abs(lock.F - 0.58068) < 5e-4
    assert abs(lock.alpha_c - 0.1818) < 5e-4
    assert abs(lock.Delta_c_meV - 0.109) < 5e-3
    assert abs(lock.Delta0_meV - 0.600) < 1e-12
    assert abs(lock.lambda_fold - 0.8475) < 5e-3
    assert abs(lock.lambda_star - 0.7628) < 5e-3
    assert abs(lock.Delta_star_meV - 0.232) < 5e-3
    assert lock.E_vac_meV != E_CAS_WITHDRAWN_MEV
    assert lock.drive_scale_vs_old_0244 == lock.E_vac_meV / E_CAS_WITHDRAWN_MEV


def test_live_gap_uses_psi_no_time_arg():
    sig = inspect.signature(delta_live_from_psi)
    assert "t" not in sig.parameters
    assert "time" not in sig.parameters
    assert "psi" in sig.parameters
    src = inspect.getsource(delta_live_from_psi)
    # Body must not interpolate the withdrawn schedule (docstring may mention it).
    body = src.split('"""', 2)[-1]
    assert "return" in body
    assert "0.389" not in body
    assert "0.032" not in body

    lock = build_frozen_tqgl()
    grid = make_grid(lock, n_grid=32, L_um=0.200)
    psi = hopfion_initial(grid, lock, seed=1)
    n_ref = mean_density_weighted(psi, grid.weight)
    d0 = delta_live_from_psi(psi, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV)
    assert abs(d0 - lock.Delta0_meV) < 1e-9
    d_half = delta_live_from_psi(
        0.5 * psi, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
    )
    assert abs(d_half - 0.5 * lock.Delta0_meV) < 1e-9
    d_other = delta_live_from_psi(
        psi * 1.2, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
    )
    assert d_other != d0


def test_gpe_norm_and_live_gap_from_field():
    from qvc.tqgl.coupled import run_coupled_gpe_dce

    lock = build_frozen_tqgl()
    out = run_coupled_gpe_dce(
        lock, n_grid=32, t_end_ps=0.02, dt_ps=0.002, pump_mode="resonant"
    )
    assert out["live_gap_uses_psi"]
    assert out["live_gap_has_no_time_argument"]
    assert out["norm_rel_drift"] < 1e-8
    assert abs(out["Delta_live_meV"][0] - lock.Delta0_meV) < 1e-9
    grid = make_grid(lock, n_grid=32)
    psi = hopfion_initial(grid, lock)
    ke = kinetic_coeff_mev_um2(lock)
    prop = kinetic_propagator(grid, 0.002, ke)
    V = nl_potential(psi, -lock.Delta0_meV, lock.Delta0_meV)
    n0 = l2_norm(psi, grid.dx)
    psi2 = split_step_advance(psi, prop_K=prop, V=V, dt_ps=0.002)
    assert abs(l2_norm(psi2, grid.dx) - n0) / n0 < 1e-10
    assert 0.0 <= coherence_g1(psi) <= 1.0


def test_dce_does_not_backcompute_0p52():
    lock = build_frozen_tqgl()
    out = run_dce_fold_slaved(
        lock, DCEConfig(pump_mode="resonant", t_end_ps=4.0, dt_ps=0.05)
    )
    assert out["A_c_is_fold_not_BKT"]
    assert out["A_c_not_backcomputed_from_0p52"]
    assert out["paper_ratio_0p52_not_imposed"]
    assert not out["Lorentzian_identically_zero"]
    assert abs(out["A_c_fold_um_inv"] - out["A_star_um_inv"] / 0.52) > 0.05
    assert abs(out["A_c_fold_um_inv"] - lock.A_c_fold_um_inv) < 1e-12
    integ = make_dce(lock, DCEConfig(pump_mode="detuning"))
    assert integ.Lam_loss > 0.0


def test_no_folklore_as_suite_results():
    from qvc.tqgl.coupled import run_coupled_gpe_dce

    lock = build_frozen_tqgl()
    out = run_coupled_gpe_dce(lock, n_grid=32, t_end_ps=0.016, dt_ps=0.002)
    folklore = {E_CAS_WITHDRAWN_MEV, DELTA_QVC_FOLKLORE_MEV, DELTA_1PS_FOLKLORE_MEV}
    for key in ("Delta_live_final_meV", "Delta_live_min_meV", "A_final_um_inv"):
        val = float(out[key])
        for f in folklore:
            assert abs(val - f) > 1e-4
    assert abs(out["Delta_live_final_meV"] - 0.032) > 1e-3
    assert abs(out["Delta_live_final_meV"] - 0.389) > 1e-3


def test_schwinger_berry_not_democratic_g():
    from qvc.tqgl.diagnostics import berry_goe_gue, schwinger_twin_channel

    lock = build_frozen_tqgl()
    sch = schwinger_twin_channel(lock)
    assert sch["lambda_min_is_minus_g0"]
    assert sch["E_pair_not_0p172_unless_rederived"]
    assert abs(sch["E_pair_meV"] - 0.172) > 1e-4
    assert abs(sch["g0_democratic_meV"] - lock.lambda_star * lock.E_vac_meV) < 1e-12
    berry = berry_goe_gue(n_samples=64, seed=0)
    assert berry["separate_from_democratic_G"]
    assert math.isfinite(berry["R_mean"])
    assert berry["R_mean"] > 0.0
