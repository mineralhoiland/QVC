"""Light tests for Tier C closed-form locks."""

from __future__ import annotations

import math

import numpy as np

from qvc.coupling_matrix import (
    catalyzon_dark_basis,
    democratic_coupling_matrix,
    democratic_spectrum,
)
from qvc.gap_fold import fold_condition, tree_level_ell
from qvc.lens_chern_simons import chern_simons_flat_u1, linking_form_self
from qvc.retardation import eta_ret_from_healing
from qvc.vacuum_freeze import assert_preprint_0p244_not_reproduced


def test_democratic_spectrum():
    N, g0 = 5, 1.0
    eigs = np.sort(np.linalg.eigvalsh(democratic_coupling_matrix(N, g0)))
    spec = democratic_spectrum(N, g0)
    assert np.isclose(eigs[-1], spec["lambda_max"])
    assert np.allclose(eigs[:-1], spec["lambda_min"])


def test_helmert_orthonormal():
    B = catalyzon_dark_basis(5)
    assert np.allclose(B @ B.T, np.eye(4), atol=1e-12)
    assert np.allclose(B.sum(axis=1), 0.0, atol=1e-12)


def test_cs_linking():
    assert chern_simons_flat_u1(1, 5) == 0.2
    assert linking_form_self(5) == 0.2


def test_tree_level_ell():
    assert abs(tree_level_ell() - (math.log(0.5) + 0.5772156649015329)) < 1e-12


def test_fold_alpha_equals_delta():
    dc, ac = fold_condition(-2.79)
    assert abs(dc - ac) < 1e-12


def test_eta_ret_healing():
    assert abs(eta_ret_from_healing(0.6, 0.07) - 2.0) < 1e-12


def test_preprint_0p244_fails():
    audit = assert_preprint_0p244_not_reproduced()
    assert audit["claim_reproduced"] is False


def test_bridge_ell_matches_corrections():
    from qvc.bridge import derive_bridge
    from qvc.gap_fold import ELL_CORRECTIONS

    b = derive_bridge()
    assert b["E_vac_equals_hbar_cs_A"]
    assert abs(b["ell_frozen"] - ELL_CORRECTIONS) / abs(ELL_CORRECTIONS) < 0.01
    assert abs(b["alpha"] - b["E_vac_meV"] / b["Delta0_meV"]) < 1e-12
    assert b["folded_out"] is True


def test_tqgl_ws_acceptance():
    from qvc.gij_tqgl_ws import run_tqgl_ws_G

    g = run_tqgl_ws_G(n_grid=28)
    acc = g["acceptance"]
    assert acc["hermitian"]
    assert acc["diag_zero"]
    assert acc["lambda_min_lt_0"]
    assert acc["near_democratic"]


def test_gap_reduction_lambda_fold():
    from qvc.gap_reduction import recompute_gap_reduction

    r = recompute_gap_reduction()
    f = r["fold_channel"]
    assert f["lambda_fold"] < 1.0
    assert f["delta_1ps_below_delta_c"] is True
    assert f["max_static_reduction_pct"] < 92.0
    assert r["recommendation"]["lambda_star"] < f["lambda_fold"]
    assert r["recommendation"]["fold"]["folded_out"] is False
    # unit λ remains supercritical
    assert r["fold_channel"]["unit_lambda"]["folded_out"] is True


def test_lambda_fold_bound_and_partition():
    from qvc.gap_reduction import evaluate_partition, run_gap_reduction_suite

    s = run_gap_reduction_suite()
    assert s["lambda_alpha_max"] < 1.0
    assert s["g0_working_exceeds_Evac"] is True
    unit = s["scenarios"]["unit_alpha_foldout"]
    assert unit["static"]["folded_out"] is True
    eq = s["scenarios"]["equal_split"]
    assert eq["static"]["folded_out"] is False
    assert abs(eq["catalyzon"]["g0_meV"] - 0.5 * s["bridge"]["E_vac_meV"]) < 1e-12
    sat = s["scenarios"]["fold_saturated_leftover_G"]
    assert sat["static"]["folded_out"] is False
    row = evaluate_partition(0.4, lam_G=0.4)
    assert abs(row["lambda_tot"] - 0.8) < 1e-12
