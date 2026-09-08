"""Tests for the Track A L1 parameter dictionary.

Fail if ξ=181 nm is mixed into a_H=8 nm, if η_ret=3.1 is treated as the unique
lock, or if g0=0.21 is emitted as derived.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.retardation import eta_ret  # noqa: E402

from qvc_l1.acoustic import acoustic_matching  # noqa: E402
from qvc_l1.catalyzon import catalyzon_card, g0_from_lambda  # noqa: E402
from qvc_l1.gl_coefficients import gl_coefficients, kappa_GL_mev_um2  # noqa: E402
from qvc_l1.hall import E2_OVER_H_S, hall_convention  # noqa: E402
from qvc_l1.kinematics import locked_kinematics  # noqa: E402
from qvc_l1.lengths import length_dictionary  # noqa: E402
from qvc_l1.locks import (  # noqa: E402
    CARTOON_BARRIER_PCT,
    ETA_RET_ALTERNATE,
    FOLKLORE_DELTA_1PS_MEV,
    FOLKLORE_DELTA_QVC_MEV,
    FOLKLORE_E_CAS_MEV,
    HOLOGRAPHIC_GAMMA_QM_MEV_NM2,
    RETIRED_G0_MEV,
    XI_ALTERNATE_UM,
    build_l1_locks,
)
from qvc_l1.mexican_hat import mexican_hat  # noqa: E402
from qvc_l1.dictionary import build_dictionary, run_all  # noqa: E402

FOLKLORE = (FOLKLORE_E_CAS_MEV, FOLKLORE_DELTA_QVC_MEV, FOLKLORE_DELTA_1PS_MEV, 92.0)


@pytest.fixture(scope="module")
def locks():
    return build_l1_locks()


@pytest.fixture(scope="module")
def bundle():
    return build_dictionary()


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
    assert abs(locks.lambda_fold - 0.8475) < 5e-4
    assert abs(locks.lambda_star - 0.7628) < 5e-4
    assert abs(locks.Delta_star_meV - 0.2324) < 5e-4
    assert abs(locks.T_star_K * 1e3 - 118.7) < 0.5
    assert abs(locks.E_pair_meV - 0.1963) < 5e-4


def test_a_H_is_r_not_xi_181(locks):
    """FAIL if someone mixes ξ=181 nm into a_H=8 nm."""
    a_H = locks.a_H_um
    assert abs(a_H - 0.008) < 1e-15
    assert abs(a_H - XI_ALTERNATE_UM) > 0.1
    assert abs(1e3 * a_H - 8.0) < 1e-12
    assert abs(1e3 * a_H - 181.0) > 100.0
    lengths = length_dictionary(locks)
    assert lengths["a_H_is_not_xi_181"]
    assert lengths["lengths"]["a_H"]["is_lock"]
    assert not lengths["lengths"]["xi_181"]["is_lock"]


def test_three_lengths_are_distinct(locks):
    lengths = length_dictionary(locks)
    a_H = lengths["lengths"]["a_H"]["value_um"]
    xi_ph = lengths["lengths"]["xi_ph"]["value_um"]
    xi_181 = lengths["lengths"]["xi_181"]["value_um"]
    assert abs(a_H - xi_ph) > 0.05
    assert abs(a_H - xi_181) > 0.05
    assert abs(xi_ph - xi_181) > 0.02
    assert abs(xi_ph - locks.hbar_cs_meV_um / locks.Delta0_meV) < 1e-15


def test_eta_ret_identity_and_not_unique_3p1(locks):
    """FAIL if η_ret=3.1 is treated as the unique lock."""
    lengths = length_dictionary(locks)
    eta_aH = lengths["eta_ret_a_H"]
    eta_ph = lengths["eta_ret_xi_ph"]
    eta_181 = lengths["eta_ret_xi_181"]
    assert abs(eta_ph - 2.0) < 1e-12
    assert abs(eta_aH - 0.13714285714285712) < 1e-9
    assert abs(eta_aH - 2.0) > 1.0
    assert abs(eta_aH - ETA_RET_ALTERNATE) > 2.0
    assert lengths["unique_lock_is_not_3p1"]
    assert not lengths["lengths"]["xi_181"]["eta_ret_boxed_3p1_is_unique"]
    # 3.1 is the alternate evaluation, not a lock
    assert abs(eta_181 - ETA_RET_ALTERNATE) < 0.02
    card = lengths["lengths"]
    assert card["a_H"]["status"] == "locked"
    assert card["xi_ph"]["status"] == "derived"
    assert card["xi_181"]["status"] == "alternate"


def test_eta_ret_formula_matches_qvccompute(locks):
    a_H = locks.a_H_um
    expected = eta_ret(locks.Delta0_meV, a_H, locks.hbar_cs_meV_um)
    got = length_dictionary(locks)["eta_ret_a_H"]
    assert abs(got - expected) < 1e-15


def test_g0_star_is_lambda_Evac_not_0p21(locks):
    """FAIL if g0=0.21 is treated as derived."""
    cat = catalyzon_card(locks)
    g0 = cat["g0_star_meV"]
    assert abs(g0 - locks.lambda_star * locks.E_vac_meV) < 1e-15
    assert abs(g0 - 0.0982) < 5e-4
    assert abs(g0 - RETIRED_G0_MEV) > 0.05
    assert cat["at_lambda_star"]["retired_g0_0p21_not_used"]
    assert cat["lambda_min_equals_minus_g0"]
    assert abs(cat["Delta_cat_star_meV"] - (locks.Delta0_meV - g0)) < 1e-15
    assert abs(cat["at_lambda_star"]["reduction_pct"] - 16.4) < 0.2
    # fold point: g0 = Δ_c ≈ α_c Δ₀
    assert cat["g0_fold_equals_Delta_c"]
    assert cat["reduction_fold_equals_alpha_c"]
    assert abs(g0_from_lambda(locks, 1.0) - locks.E_vac_meV) < 1e-15


def test_kappa_GL_from_a_H_not_gamma_QM(locks):
    k_gl = kappa_GL_mev_um2(locks)
    assert abs(k_gl - locks.Delta0_meV * locks.a_H_um**2) < 1e-18
    assert abs(k_gl - 3.84e-5) < 1e-12
    # holographic γ_QM is meV·nm²; converting 3.73 meV·nm² → meV·µm² is 3.73e-6
    gamma_as_um2 = HOLOGRAPHIC_GAMMA_QM_MEV_NM2 * 1e-6
    assert abs(k_gl - gamma_as_um2) > 1e-5
    gl = gl_coefficients(locks)
    assert gl["kappa_GL"]["not_gamma_QM"]
    assert gl["kappa_S"]["not_holographic_gamma_QM"]
    assert gl["cartoon_not_used"]["used"] is False


def test_beta_from_hopfion_n_ref_not_cartoon(locks):
    gl = gl_coefficients(locks)
    n_ref = gl["n_ref"]["n_ref"]
    beta = gl["beta"]["beta_meV"]
    assert n_ref > 0.5
    assert n_ref < 1.0  # core depletion on the ring
    assert abs(beta - locks.Delta0_meV / n_ref) < 1e-12
    assert abs(gl["alpha_GL"]["value_meV"] + locks.Delta0_meV) < 1e-15
    hat = mexican_hat(locks, beta_meV=beta, n_ref=n_ref)
    assert hat["cartoon_not_used"]["used"] is False
    assert abs(hat["cartoon_not_used"]["barrier_pct"] - CARTOON_BARRIER_PCT) < 1e-12
    red = hat["catalyzon_dressed"]["barrier_reduction_pct"]
    assert abs(red - CARTOON_BARRIER_PCT) > 5.0
    assert hat["catalyzon_dressed"]["replaces_cartoon_42pct"]
    # quadratic identity 1 - (Δ_cat/Δ0)²
    g0 = locks.lambda_star * locks.E_vac_meV
    expected = 100.0 * (1.0 - ((locks.Delta0_meV - g0) / locks.Delta0_meV) ** 2)
    assert abs(red - expected) < 1e-6


def test_hall_one_formula(locks):
    hall = hall_convention(locks)
    assert hall["k_CS"] == 2 * locks.N_layers
    assert abs(hall["Q_H"] - 1.0 / locks.N_layers) < 1e-15
    assert hall["pairing_identity"]
    assert abs(hall["delta_sigma_xy_S"] - E2_OVER_H_S / hall["k_CS"]) < 1e-18
    assert abs(hall["delta_sigma_xy_S"] - hall["Q_H"] * E2_OVER_H_S / 2.0) < 1e-18
    assert abs(hall["sigma_xy_parent_S"] - locks.N_layers * E2_OVER_H_S) < 1e-18
    assert hall["not_in_GPE_EOM"]
    assert hall["pairing_is_lock_product"]


def test_R_WS_matching_not_withdrawn_0p91(locks):
    ac = acoustic_matching(locks)
    assert abs(ac["R_WS_um"] - 2.18) < 0.05
    assert abs(ac["R_WS_um"] - 0.91) > 1.0
    assert ac["not_a_duality"]
    assert ac["Berry_R_mean_is_ratio"]
    # withdrawn 285 mK reproduces the old ~0.91 µm
    assert abs(ac["R_WS_withdrawn_285mK_um"] - 0.91) < 0.02


def test_parameter_card_rejects_folklore_and_insertions(bundle):
    card = bundle["parameter_card"]
    g0 = card["g0_star_meV"]
    assert abs(g0 - RETIRED_G0_MEV) > 0.05
    assert abs(card["eta_ret_a_H"] - ETA_RET_ALTERNATE) > 2.0
    assert abs(card["a_H_um"] - XI_ALTERNATE_UM) > 0.1
    dumped = json.dumps(card["folklore_not_results"])
    for val in FOLKLORE[:3]:
        assert str(val) in dumped  # recorded as folklore, not as a result key
    assert "E_Cas_meV" not in card
    assert "Delta_QVC_meV" not in card
    assert card["kappa_GL_meV_um2"] == pytest.approx(3.84e-5, rel=1e-12)
    ledger = bundle["status_ledger"]
    retired = " ".join(ledger["retired"])
    assert "0.21" in retired
    assert "3.1" in " ".join(ledger["alternate_not_lock"])


def test_run_all_writes_outputs(tmp_path):
    bundle = run_all(out_dir=tmp_path, make_figures=True, mirror=False)
    assert (tmp_path / "summary.json").is_file()
    assert (tmp_path / "summary.md").is_file()
    assert (tmp_path / "fig_lengths.png").is_file()
    assert (tmp_path / "fig_mexican_hat.png").is_file()
    assert (tmp_path / "fig_spec_G.png").is_file()
    payload = json.loads((tmp_path / "summary.json").read_text())
    text = (tmp_path / "summary.md").read_text()
    assert "not the lock" in text
    assert "g_0^\\star" in text or "g0" in json.dumps(payload)
    assert abs(payload["parameter_card"]["g0_star_meV"] - 0.21) > 0.05
    assert payload["lengths"]["unique_lock_is_not_3p1"]
    assert "Hubbard" in " ".join(bundle["status_ledger"]["out_of_scope"])
    assert bundle["figures"]


def test_locked_kinematics_two_velocities_and_vk(locks):
    kin = locked_kinematics(locks.tqgl, locks)
    assert kin["velocities_not_identified"]
    assert abs(kin["m_star_over_m_e"] - 1.0) < 0.02
    assert abs(kin["hbar_omega_ph_aH_meV"] - 8.75) < 1e-6
    assert abs(kin["drive_Delta0"]["nu_THz"] - 0.290) < 0.002
    assert abs(kin["drive_Deltastar"]["nu_THz"] - 0.112) < 0.002
    assert abs(kin["drive_Deltastar"]["q_over_K_M"] - 0.102) < 0.005
    assert kin["drive_Deltastar"]["q_over_K_M"] < 0.15
    assert abs(kin["C_VK_min_meV"] - 0.430) < 0.005
    assert kin["not_E_cas_0244"]
    # withdrawn 290 THz and VK 0.816 are not the lock
    assert kin["drive_Delta0"]["nu_THz"] < 1.0
    assert kin["C_VK_min_meV"] < 0.6

