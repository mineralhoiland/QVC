"""Tests for L1 G_ij (B1) and Lindblad GPE (B2).

Fast grids. Folklore 0.244/0.389/0.032/92% must not appear as physics results.
G_ij must not be hand-inserted as g0(J−I).
"""

from __future__ import annotations

import inspect
from pathlib import Path

import numpy as np
import pytest

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.coupling_matrix import democratic_spectrum  # noqa: E402
from qvc_l1.gij_ring import (  # noqa: E402
    run_ring_fourier_G,
    run_ws_gaussian_G,
)
from qvc_l1.lindblad import apply_lindblad, run_gpe_lindblad  # noqa: E402
from qvc_l1.locks import (  # noqa: E402
    FOLKLORE_DELTA_1PS_MEV,
    FOLKLORE_DELTA_QVC_MEV,
    FOLKLORE_E_CAS_MEV,
    FOLKLORE_PCT,
    build_l1_locks,
)

FOLKLORE = (FOLKLORE_E_CAS_MEV, FOLKLORE_DELTA_QVC_MEV, FOLKLORE_DELTA_1PS_MEV, FOLKLORE_PCT)


@pytest.fixture(scope="module")
def locks():
    return build_l1_locks()


def test_option_a_star_locks(locks):
    assert abs(locks.r_nm - 8.0) < 1e-12
    assert abs(locks.R_nm - 50.0) < 1e-12
    assert abs(locks.hbar_cs_meV_um - 0.07) < 1e-12
    assert abs(locks.Delta0_meV - 0.600) < 1e-12
    assert abs(locks.E_vac_meV - 0.1287) < 5e-4
    assert abs(locks.lambda_star - 0.7628) < 5e-4
    assert abs(locks.g0_star_meV - locks.lambda_star * locks.E_vac_meV) < 1e-12
    assert abs(locks.g0_star_meV - 0.0982) < 5e-4
    assert locks.N_layers == 5
    assert locks.cs_status == "truncated"
    assert abs(locks.Q_H - 1.0 / 5.0) < 1e-12


def test_democratic_algebra_not_minus_4g0(locks):
    spec = democratic_spectrum(locks.N_layers, locks.g0_star_meV)
    assert abs(spec["lambda_max"] - 4.0 * locks.g0_star_meV) < 1e-12
    assert abs(spec["lambda_min"] + locks.g0_star_meV) < 1e-12
    assert spec["lambda_min_mult"] == 4
    assert abs(locks.dem_lambda_max_over_g0 - 4.0) < 1e-12
    assert abs(locks.dem_lambda_min_over_g0 + 1.0) < 1e-12


def test_ws_finer_grid_computed_not_JI(locks):
    ws = run_ws_gaussian_G(locks, n_grid=32)
    G = np.asarray(ws["G"])
    assert G.shape == (5, 5)
    assert ws["acceptance"]["hermitian"]
    assert ws["acceptance"]["diag_zero"]
    assert ws["acceptance"]["lambda_min_lt_0"]
    assert ws["g0_eff_is_not_g0_star"]
    # Kernel number is not the BRIDGE meV identification.
    assert abs(ws["g0_eff_kernel"] - locks.g0_star_meV) > 1e-4
    # Must not be exactly g0(J−I) at g0★ (would mean hand insertion).
    j_minus_i = np.ones((5, 5)) - np.eye(5)
    assert not np.allclose(G, locks.g0_star_meV * j_minus_i)
    src = inspect.getsource(run_ws_gaussian_G)
    assert "democratic_coupling_matrix" not in src
    assert "run_tqgl_ws_G" in src


def test_ring_fourier_hermitian_diag_zero(locks):
    ring = run_ring_fourier_G(locks, n_grid=48)
    G = np.asarray(ring["G"])
    assert G.shape == (5, 5)
    assert ring["acceptance"]["hermitian"]
    assert ring["acceptance"]["diag_zero"]
    assert ring["acceptance"]["lambda_min_lt_0"]
    assert ring["harmonics"] == [0, 1, 2, 3, 4]
    assert ring["g0_eff_is_not_g0_star"]
    assert "not_torsion" in ring["status"]
    assert "not_TQGL_theorem" in ring["status"]
    src = inspect.getsource(run_ring_fourier_G)
    assert "ones((N, N)) - np.eye" not in src
    assert "ones((n, n)) - np.eye" not in src
    assert ring["gpe_quadratic"]["status"].startswith("diagnostic")
    assert np.allclose(np.diag(ring["gpe_quadratic"]["G"]), 0.0, atol=1e-10)


def test_number_damping_reduces_norm(locks):
    coh = run_gpe_lindblad(
        locks.tqgl, n_grid=32, t_end_ps=0.04, dt_ps=0.002, couple_dce=False
    )
    damp = run_gpe_lindblad(
        locks.tqgl,
        n_grid=32,
        t_end_ps=0.04,
        dt_ps=0.002,
        couple_dce=False,
        gamma_n=5.0,
        gamma_phi=0.0,
    )
    assert coh["norm_ratio_final"] == pytest.approx(1.0, rel=1e-6, abs=1e-6)
    assert damp["norm_ratio_final"] < 0.95
    assert damp["Delta_live_final_meV"] < coh["Delta_live_final_meV"]
    assert abs(damp["Delta_live_meV"][0] - locks.Delta0_meV) < 1e-9


def test_dephasing_reduces_g1(locks):
    coh = run_gpe_lindblad(
        locks.tqgl, n_grid=32, t_end_ps=0.08, dt_ps=0.002, couple_dce=False, seed=1
    )
    deph = run_gpe_lindblad(
        locks.tqgl,
        n_grid=32,
        t_end_ps=0.08,
        dt_ps=0.002,
        couple_dce=False,
        gamma_n=0.0,
        gamma_phi=8.0,
        seed=1,
        dephasing_seed=3,
    )
    assert deph["g1_final"] < coh["g1_final"] - 0.02
    assert deph["live_gap_uses_psi"]
    assert deph["no_prescribed_Delta_t"]


def test_live_gap_no_time_argument():
    from qvc.tqgl.live_gap import delta_live_from_psi

    sig = inspect.signature(delta_live_from_psi)
    assert "t" not in sig.parameters
    assert "time" not in sig.parameters


def test_no_folklore_as_results(locks):
    out = run_gpe_lindblad(
        locks.tqgl, n_grid=32, t_end_ps=0.02, dt_ps=0.002, couple_dce=True
    )
    for val in (out["Delta_live_final_meV"], out["g1_final"], out["Delta_live_meV"][0]):
        for f in FOLKLORE[:3]:
            assert abs(float(val) - f) > 1e-3
    assert abs(out["Delta_live_final_meV"] - 0.032) > 1e-3
    assert abs(out["Delta_live_final_meV"] - 0.389) > 1e-3


def test_apply_lindblad_channels_independent():
    rng = np.random.default_rng(0)
    psi = np.ones((4, 4), dtype=np.complex128)
    n0 = np.sum(np.abs(psi) ** 2)
    lost = apply_lindblad(psi, dt_ps=0.1, gamma_n=2.0, gamma_phi=0.0, rng=None)
    assert np.sum(np.abs(lost) ** 2) < n0
    deph = apply_lindblad(psi, dt_ps=0.1, gamma_n=0.0, gamma_phi=4.0, rng=rng)
    assert np.sum(np.abs(deph) ** 2) == pytest.approx(n0, rel=1e-12)
    assert not np.allclose(deph, psi)


def test_gij_source_does_not_insert_JI():
    from qvc_l1 import gij_ring

    src = Path(gij_ring.__file__).read_text()
    assert "g0 * (np.ones" not in src
    assert "g0 * (np.ones((N, N)) - np.eye" not in src
    assert "democratic_coupling_matrix(" in src  # comparison only
    assert "overlap_G_complex" in src
    assert "run_tqgl_ws_G" in src
    assert "hopfion_initial" in src
