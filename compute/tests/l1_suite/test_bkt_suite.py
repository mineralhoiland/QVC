"""Tests for the QVC / TEVC BKT C1–C6 suite.

Fast grids. Folklore 0.244/0.389/0.032/92% must not appear as physics results.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc_bkt.locks import (  # noqa: E402
    FOLKLORE_DELTA_1PS_MEV,
    FOLKLORE_DELTA_QVC_MEV,
    FOLKLORE_E_CAS_MEV,
    build_bkt_locks,
)
from qvc_bkt.stiffness import (  # noqa: E402
    K_of_kappa,
    kappa_GL_mev_um2,
    kappa_conv_meV,
    n0_healing_disk_um2,
    run_c1,
    twist_stiffness_uniform,
)
from qvc_bkt.vortex_core import run_c2  # noqa: E402
from qvc_bkt.kt_flow import K_C_BOSONIC, fugacity_y, integrate_kt  # noqa: E402
from qvc_bkt.model_compare import fit_sqrt_branch, fold_branch_curve  # noqa: E402
from qvc_bkt.phonon_damping import ahns_criterion  # noqa: E402
from qvc_bkt.suite import run_all  # noqa: E402


FOLKLORE = (FOLKLORE_E_CAS_MEV, FOLKLORE_DELTA_QVC_MEV, FOLKLORE_DELTA_1PS_MEV, 92.0)


@pytest.fixture(scope="module")
def locks():
    return build_bkt_locks()


def test_option_a_star_locks(locks):
    assert abs(locks.r_um - 0.008) < 1e-12
    assert abs(locks.R_um - 0.050) < 1e-12
    assert abs(locks.hbar_cs_meV_um - 0.07) < 1e-12
    assert abs(locks.F - 0.58068) < 5e-4
    assert abs(locks.E_vac_meV - 0.1287) < 5e-4
    assert abs(locks.alpha_c - 0.1818) < 5e-4
    assert abs(locks.delta_c - locks.alpha_c) < 1e-9
    assert abs(locks.Delta_c_meV - 0.109) < 5e-3
    assert abs(locks.Delta0_meV - 0.600) < 1e-12
    assert abs(locks.lambda_star - 0.7628) < 5e-4
    assert abs(locks.Delta_star_meV - 0.2324) < 5e-4
    assert abs(locks.T_star_K * 1e3 - 118.7) < 0.5
    assert abs(locks.E_pair_meV - 0.1963) < 5e-4


def test_democratic_lambda_min_and_QH(locks):
    assert abs(locks.lambda_min_meV + locks.g0_meV) < 1e-12
    assert abs(locks.Q_H - 1.0 / locks.N_layers) < 1e-12
    assert locks.N_layers == 5
    assert locks.cs_status == "truncated"


def test_kosterlitz_convention():
    assert abs(K_C_BOSONIC - 2.0 / math.pi) < 1e-15


def test_healing_disk_cancels_aH(locks):
    dens = n0_healing_disk_um2(locks.r_um)
    k_gl = kappa_GL_mev_um2(locks.tqgl)
    k_uv = kappa_conv_meV(kappa_GL=k_gl, n0_um2=dens["n0_um2"], rho2=1.0)
    assert abs(k_uv - 2.0 * locks.Delta0_meV / math.pi) < 1e-9


def test_twist_matches_two_sector(locks):
    k_gl = kappa_GL_mev_um2(locks.tqgl)
    n0 = n0_healing_disk_um2(locks.r_um)["n0_um2"]
    tw = twist_stiffness_uniform(kappa_GL=k_gl, n0_um2=n0, rho2=1.0)
    assert tw["rel_err"] < 0.05


def test_c1_no_T_BKT_bare_input(locks):
    c1 = run_c1(locks, n_alpha=24)
    dumped = json.dumps(c1["primary"]["closures"])
    assert "T_BKT_bare" in dumped  # mentioned as NOT used
    assert "NOT an input" in dumped
    k_fold = c1["primary"]["K_fold_Tstar"]
    fake = (2.0 / math.pi) * (1.0 / locks.T_star_K) * (locks.delta_c**2)
    assert abs(k_fold - fake) > 0.2


def test_c2_vortex_core_oom(locks):
    c2 = run_c2(locks, n_grid=240, R_over_xi=24.0)
    assert c2["defect"] == "2d_phase_vortex_winding_1"
    assert c2["radial"]["not_hopfion"]
    ratio = c2["E_core_over_kappa"]
    assert math.isfinite(ratio)
    assert abs(ratio) < 20.0
    assert c2["radial"]["a_half_um"] > 0.0
    assert c2["radial"]["xi_um"] > 0.0


def test_kt_flow_bound_vs_plasma():
    bound = integrate_kt(2.0, 1e-4, 8.0)
    assert bound["bound"]
    assert not bound["plasma"]
    plasma = integrate_kt(0.3, 0.05, 8.0)
    assert plasma["plasma"]


def test_fugacity_from_core_energy(locks):
    y = fugacity_y(0.1, locks.T_star_K, locks.kB_meV_per_K)
    t = locks.kT_star_meV
    assert abs(y - math.exp(-0.1 / t)) < 1e-12


def test_c3_frozen_alpha_no_beta(locks):
    c1 = run_c1(locks, n_alpha=24)
    c2 = {
        "E_core_over_kappa": 0.5,
        "a_c_um": locks.r_um,
        "a_c_label": "test",
    }
    from qvc_bkt.kt_flow import run_c3

    c3 = run_c3(c1, c2, locks, ell_ir=4.0)
    assert c3["scan_Tstar_primary"]["option_A_star"] == "alpha_frozen_no_beta_alpha"
    assert c3["scan_Tstar_primary"]["cs"] == "truncated"
    assert c3["scan_Tstar_anyon"]["cs"] == "restored_anyon"
    assert c3["scan_Tstar_anyon"]["K_c"] == pytest.approx(
        (2.0 / math.pi) * (1.0 - 1.0 / 5.0) ** 2, rel=1e-12
    )


def test_c5_fold_slope_near_half(locks):
    cur = fold_branch_curve(locks, n=120)
    fit = fit_sqrt_branch(cur["t_A"], cur["delta"], locks.delta_c)
    assert fit["tA_decades"] >= 1.5
    assert abs(fit["slope_free"] - 0.5) < 0.08


def test_ahns_algebra():
    under = ahns_criterion(0.01, 1.0)
    assert not under["overdamped"]
    assert under["verdict"] == "sharp_jump_survives"
    over = ahns_criterion(10.0, 0.1)
    assert over["overdamped"]


def test_full_suite_no_folklore(tmp_path: Path):
    bundle = run_all(out_dir=tmp_path, n_alpha=24, make_figures=False)
    kn = bundle["key_numbers"]
    for bad in FOLKLORE:
        for v in kn.values():
            if isinstance(v, float) and abs(v - float(bad)) < 1e-9:
                raise AssertionError(f"folklore {bad} leaked into key_numbers")
    assert kn["kt_before_fold_Tstar_primary"] in (True, False)
    assert (tmp_path / "summary.json").is_file()
    assert (tmp_path / "summary.md").is_file()
    text = (tmp_path / "summary.md").read_text()
    assert "0.244" not in text
    assert "0.389" not in text
    assert "not a BKT theorem" in text.lower() or "not** a BKT theorem" in text
    assert "kt_before_fold" in json.dumps(bundle["h_vrg1"])
    assert bundle["computed_vs_assumed"]["C1"]["not_a_theorem"] is True
    assert bundle["C3"]["scan_Tstar_primary"]["option_A_star"] == "alpha_frozen_no_beta_alpha"


def test_K_equals_kappa_over_T(locks):
    K = K_of_kappa(0.2, locks.T_star_K, locks.kB_meV_per_K)
    assert abs(K - 0.2 / locks.kT_star_meV) < 1e-12
