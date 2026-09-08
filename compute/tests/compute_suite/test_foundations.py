"""Foundations: phenomenological gap lock, parallel channels, F_moiré, Schwinger."""

from __future__ import annotations

from qvc.cavity_moire import apply_F_moire_to_A_geom
from qvc.gap_foundations import run_gap_foundations
from qvc.parallel_channels import run_parallel_channels
from qvc.schwinger_rate import boxed_schwinger_rate


def test_sketch_does_not_derive_lock():
    g = run_gap_foundations()
    assert g["sketch_derives_lock"] is False
    assert g["sketch_cutoff"]["paper_IR_exceeds_UV"]
    assert g["phenomenological_lock"]["alpha_c_equals_delta_c"]
    assert g["lock_is_internally_consistent"]


def test_static_fold_collapse_confirmed():
    fold = run_gap_foundations()["static_fold_collapse"]
    assert fold["confirms_collapse"]
    assert fold["not_BKT_vanishing"]
    assert fold["not_a_GPE_result"]
    assert fold["two_branches_below_fold"]
    assert fold["unit_lambda_folded_out"]
    assert fold["star_on_physical_branch"]
    assert fold["folklore_92pct_unreachable"]
    assert abs(fold["square_root_slope_mean"] - 0.5) < 0.06


def test_parallel_channels_mexican_hat_and_gpe_signs():
    ch = run_parallel_channels()
    assert ch["fold"]["confirms_static_fold_collapse"]
    assert ch["channels_are_parallel"]
    assert ch["do_not_stack_percentages"]
    hat = ch["mexican_hat"]
    assert hat["catalyzon_well_shallower"]
    assert hat["gpe_well_deepens"]
    assert abs(hat["barrier_reduction_pct"] - 30.0) < 0.5
    assert ch["gpe"]["live_gap_rises"]
    assert ch["gpe"]["confirms_static_fold_collapse"] is False
    assert ch["lindblad"]["fold_violating_diagnostic"]


def test_f_moire_on_A_geom_rechecks_lambda_fold():
    cav = apply_F_moire_to_A_geom()
    assert cav["withdrawn_A_ZPF_0p018_not_used"]
    assert cav["withdrawn_A_cav_0p0201_not_used"]
    assert abs(cav["F_moire"] - 1.119) < 5e-3
    assert abs(cav["lambda_fold_F"] - cav["lambda_fold"] / cav["F_moire"]) < 1e-12
    assert cav["physical_branch_exists_for_lambda_le_lambda_fold_F"]
    assert cav["retuned_star_F_folded_out"] is False
    assert cav["A_geom_um_inv"] > 1.0


def test_boxed_schwinger_matches_T_star():
    s = boxed_schwinger_rate()
    assert s["ratio_is_alpha_star"]
    assert s["T_star_matches_alpha"]
    assert s["F_moire_not_in_exponent"]
    assert abs(s["T_star_K"] * 1e3 - 118.7) < 0.5
    assert s["Gamma_sch_Hz"] > 100.0
    assert s["Gamma_sch_Hz"] < 2000.0
