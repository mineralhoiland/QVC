"""S9 acceptance tests for the TEVC / QVC–TDVT workbench.

Every viewport must agree with the Option A* locks (docs/PHASE_TWO_PLAN.md):
E_vac = 0.1287 meV, Δ* = 0.2324 meV, λ_min = −g₀ (×N−1), δσ_xy = e²/(2Nh),
Q_H = 1/N. A visualization that disagrees with the locks is a bug; a viewport
that reproduces a withdrawn folklore number (0.244 meV, 92 %, 285 mK, …) is a bug.

Heavy pieces (GPE replay, Panel app) run on small grids and are marked ``slow``.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

from qvc.workbench.locks import (
    DEPRECATED,
    DeprecatedNumberError,
    Layer,
    assert_not_deprecated,
    check_against_phase_two,
    get_locks,
    stamp,
)

L = get_locks()


# --------------------------------------------------------------------------- S0
def test_locks_match_phase_two_table():
    r = check_against_phase_two()
    failed = {k: v for k, v in r["rows"].items() if not v["pass"]}
    assert r["all_pass"], failed


def test_lock_numbers():
    assert math.isclose(L.E_vac_meV, 0.1287, rel_tol=2e-3)
    assert math.isclose(L.Delta_star_meV, 0.2324, rel_tol=2e-3)
    assert math.isclose(L.lambda_min_meV, -L.g0_star_meV, rel_tol=1e-9)
    assert L.lambda_min_multiplicity == L.N - 1
    assert L.Q_H == 1.0 / L.N
    assert math.isclose(L.hall_step_e2_over_h, 1.0 / (2 * L.N))
    assert math.isclose(L.hall_step_siemens, 3.874e-6, rel_tol=1e-3)


def test_deprecated_numbers_are_rejected():
    for rec in DEPRECATED.values():
        with pytest.raises(DeprecatedNumberError):
            assert_not_deprecated(float(rec["value"]), unit=rec["unit"])
    # the frozen numbers themselves must pass the guard
    assert_not_deprecated(L.E_vac_meV, unit="meV")
    assert_not_deprecated(L.Delta_star_meV, unit="meV")
    assert_not_deprecated(L.tqgl.gap_reduction_star_pct, unit="%")


def test_layer_stamp_prefix():
    assert stamp(Layer.L1, "x").startswith("[L1]")
    assert stamp(Layer.L2, "x").startswith("[L2]")
    assert stamp(Layer.L3, "x").startswith("[L3]")
    for _, _, layer in L.table():
        assert layer == "L1"  # all locks are laboratory-layer numbers


# --------------------------------------------------------------------------- S1
def test_tevc_geometry_and_exporters(tmp_path):
    from qvc.viz.structure import build_tevc_geometry, export_tevc, hopfion_bcc_cif

    g = build_tevc_geometry(1.10, L.N, fragment_radius_nm=1.0)
    assert g.n_atoms > 1000
    assert math.isclose(g.lambda_M_A / 10.0, 12.81, rel_tol=2e-2)  # a/(2 sin θ/2) at 1.10°
    assert len(g.layer_positions) == L.N
    r = export_tevc(tmp_path, g)
    assert set(r["files"]) >= {"cif", "poscar", "lammps", "fragment_cif"}
    for f in r["files"].values():
        assert Path(f).exists(), f
    assert "TEVC_AA_fragment" in (tmp_path / "TEVC_fragment.cif").read_text()
    hopfion_bcc_cif(tmp_path / "bcc.cif")
    txt = (tmp_path / "bcc.cif").read_text()
    assert "I m -3 m" in txt and "229" in txt
    assert "50.0" in txt  # a = ξ = 50 nm written as 50 display-Å
    assert "not a chemical species" in txt  # 'Ho' is a display proxy


# --------------------------------------------------------------------------- S2
@pytest.fixture(scope="module")
def fields():
    from qvc.viz.fields3d import compute_static_fields

    return compute_static_fields(n_grid=40, R_nm=L.R_nm, r_nm=L.r_nm, E_vac_meV=L.E_vac_meV, Q_H=L.Q_H)


def test_hopfion_unit_norm_and_vacuum_boundary(fields):
    n = fields["n"]
    assert np.allclose(np.linalg.norm(n, axis=-1), 1.0, atol=1e-8)
    assert fields["acceptance"]["boundary_is_vacuum"]


def test_whitehead_charge_converges_to_one():
    from qvc.viz.fields3d import whitehead_convergence

    c = whitehead_convergence(grids=(32, 48, 64), core_nm=L.R_nm)
    Qs = [r["Q"] for r in c["rows"]]
    assert all(np.diff(Qs) > 0), Qs  # monotone in resolution
    assert abs(c["richardson_Q"] - 1.0) < 0.05, c


def test_old_cube_header_normalisation_is_rejected(fields):
    # The QVCCursor hopfion_charge_density.cube header claimed Q≈3.81 with a 1/4π² factor.
    assert abs(fields["Q_if_4pi2_normalisation"] - 1.0) > 0.5
    assert fields["acceptance"]["old_4pi2_header_rejected"]


def test_cs_flux_is_lens_space_lock(fields):
    assert abs(fields["Q_H_flux_2d"] - L.Q_H) < 0.02
    assert abs(fields["Q_H_from_cut"] - L.Q_H) < 0.02


def test_faddeev_niemi_vk_bound(fields):
    vk = fields["VK"]
    assert vk["satisfied"]
    assert 1.0 < vk["ratio"] < 10.0  # unrelaxed stereographic ansatz sits above the relaxed 1.22


def test_casimir_envelope_equals_E_vac_lock(fields):
    cas = fields["casimir_envelope_meV"]
    # envelope −E_vac·exp(−ρ²/2r²): the grid samples the ring within (dx/r)² ≈ 1e-4 of the lock
    assert math.isclose(float(np.min(cas)), -L.E_vac_meV, rel_tol=1e-3)
    assert float(np.max(cas)) <= 0.0


def test_static_field_exporters(tmp_path, fields):
    from qvc.viz.fields3d import export_static_fields

    files = export_static_fields(fields, tmp_path)
    assert any(str(p).endswith(".cube") for p in files.values())
    assert any(str(p).endswith(".vti") for p in files.values())
    head = open(files["cube_rho_H"]).read(600)
    assert "16pi^2" in head  # header states the 1/16π² normalisation (old 1/4π² header rejected)
    assert f"CS-flux Q_H={fields['Q_H_flux_2d']:.3f}" in head  # computed flux, within 0.02 of the 1/N lock


# --------------------------------------------------------------------------- S3
def test_democratic_spectrum_panel():
    from qvc.workbench import catalyzon

    d = catalyzon.democratic_panel()
    assert d["N"] == L.N
    assert math.isclose(d["lambda_max_meV"], (L.N - 1) * L.g0_star_meV, rel_tol=1e-9)
    assert math.isclose(d["lambda_min_meV"], -L.g0_star_meV, rel_tol=1e-9)
    assert d["lambda_min_multiplicity"] == L.N - 1
    assert d["dark_subspace_residual"] < 1e-10
    assert math.isclose(d["Delta_cat_meV"], L.Delta_cat_meV)


def test_ws_overlap_is_compared_not_forced():
    from qvc.workbench import catalyzon

    ws = catalyzon.ws_overlap_comparison(n_grid=24, xi_over_lambda_values=(0.5, 3.85))
    rows = ws["rows"]
    assert len(rows) == 2 and ws["layer"] == "L1"
    for r in rows:
        # G^(H) is symmetric with N eigenvalues; λ_max/g0_eff ≈ N−1 (democratic-like), λ_min/g0_eff ≈ −1
        assert len(r["eigenvalues"]) == L.N
        assert r["lambda_max_over_g0eff"] > 0 > r["lambda_min_over_g0eff"]
    # short-range kernel deviates from J − I; the comparison is reported, not forced
    assert rows[0]["rel_frobenius_to_democratic"] > rows[1]["rel_frobenius_to_democratic"] or not rows[0]["near_democratic"]
    assert "not forced" in ws["caption"]


# --------------------------------------------------------------------------- S4
def test_fold_map_and_operating_point():
    from qvc.workbench import barriers

    fm = barriers.fold_map()
    assert math.isclose(fm["alpha_c"], 0.1818, rel_tol=2e-3)
    assert math.isclose(fm["delta_c"], fm["alpha_c"], rel_tol=1e-9)
    op = fm["operating_point"]
    assert op["alpha_over_ac"] < 1.0  # λ* sits inside the fold
    assert math.isclose(op["delta"] * fm["Delta0_meV"], L.Delta_star_meV, rel_tol=1e-6)
    assert fm["unit_lambda"]["alpha_over_ac"] > 1.0  # λ = 1 is fold-out


def test_mexican_hat_tracks_are_parallel_channels():
    from qvc.workbench import barriers

    mh = barriers.mexican_hat_tracks()
    t = mh["tracks"]
    d_bare = t["bare_Delta0"]["well_depth_meV"]
    d_cat = t["catalyzon_Delta_cat"]["well_depth_meV"]
    d_fold = t["static_fold_Delta_star"]["well_depth_meV"]
    assert d_bare > d_cat > d_fold > 0
    assert math.isclose(mh["B_cat_over_B0"], mh["B_cat_over_B0_expected"], rel_tol=1e-6)


def test_u_eff_reshaping_at_fixed_energy():
    from qvc.workbench import barriers

    A_c = L.tqgl.A_c_fold_um_inv
    u0 = barriers.u_eff_of_A(0.0)
    u1 = barriers.u_eff_of_A(0.9 * A_c)
    assert not u1["folded_out"] and math.isclose(u1["A_over_Ac"], 0.9)
    # fixed total energy: the pinned landscape has the same minimum as the bare one
    assert math.isclose(float(np.min(u1["V_meV_fixed_energy"])), float(np.min(u0["V_meV"])), rel_tol=1e-9)
    assert 0 < u1["well_depth_ratio"] < 1  # A_ZPF shallows the well (shape only)
    assert u1["Delta_track_meV"] < L.Delta0_meV
    u2 = barriers.u_eff_of_A(1.1 * A_c)
    assert u2["folded_out"] and u2["Delta_track_meV"] is None


# --------------------------------------------------------------------------- S6
def test_hall_staircase_matches_lock():
    from qvc.workbench import topology

    h = topology.hall_panel()
    assert h["lock_step_matches"]
    c = h["by_N"][L.N]
    assert math.isclose(c["step_e2h"], 1.0 / (2 * L.N))
    assert math.isclose(c["base_e2h"], float(L.N))
    assert math.isclose(c["step_S"], L.hall_step_siemens, rel_tol=1e-6)


def test_preimage_linking_is_one(fields):
    from qvc.workbench import topology

    pp = topology.preimage_panel(fields["n"], fields["grid"], n_fibres=3)
    assert abs(pp["mean_abs_linking"] - 1.0) < 0.25
    assert pp["layer"] == "L2"
    assert pp["Q_H_lock"] == L.Q_H and pp["expected_linking"] == 1.0  # kept distinct


# --------------------------------------------------------------------------- S7 (L2)
def test_analogue_panel_is_diagnostic_only():
    from qvc.workbench import analogue

    am = analogue.acoustic_metric_slice()
    assert am["layer"] == "L2"
    assert am["horizon_inside_core"]  # pure fractional vortex: no sonic horizon
    amd = analogue.acoustic_metric_slice(drain_v0_over_cs=1.2)
    assert not amd["horizon_inside_core"]
    assert not math.isclose(amd["T_H_mK"], L.T_star_mK, rel_tol=0.1)  # T_H is not T*
    ein = analogue.h_ein1_overlap()
    assert not ein["is_graviton_like"]
    tp = analogue.torsion_proxy(np.ones((4, 4, 4)), h_tor1_toggle=True)
    assert tp["feeds_Gij"] is False and "parked" in tp["note"].lower()


# --------------------------------------------------------------------------- S5 / S8 (slow)
@pytest.mark.slow
def test_gpe_replay_small():
    from qvc.workbench import dynamics

    res = dynamics.run_replay(n_grid=32, t_end_ps=0.1, n_frames=4)
    assert len(res["snapshots"]) >= 4
    assert abs(res["norm_rel_drift"]) < 1e-6
    rep = dynamics.dynamics_report(res)
    assert math.isclose(rep["Delta_star_static_fold_meV"], L.Delta_star_meV)
    fr = dynamics.replay_frames_3d(res, n_phi=24, n_chi=8)
    assert fr[0]["X_um"].shape == (24, 8)


@pytest.mark.slow
def test_panel_app_builds():
    pytest.importorskip("panel")
    from qvc.workbench import app

    tpl = app.build_app(fields_grid=32, gpe_grid=32, gpe_t_ps=0.1, gpe_frames=4)
    assert tpl is not None
    html = tpl.main.__repr__()
    assert "Column" in html or "Row" in html
