"""QVCCursor tests for foundations / parallel-channel locks imported from QVCCompute."""

from __future__ import annotations

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.cavity_moire import apply_F_moire_to_A_geom  # noqa: E402
from qvc.gap_foundations import run_gap_foundations  # noqa: E402
from qvc.parallel_channels import run_parallel_channels  # noqa: E402
from qvc.schwinger_rate import boxed_schwinger_rate  # noqa: E402
from qvc_l1.locks import build_l1_locks  # noqa: E402
from qvc_l1.mexican_hat import mexican_hat  # noqa: E402
from qvc_l1.gl_coefficients import gl_coefficients  # noqa: E402


def test_cursor_fold_collapse_matches_compute():
    g = run_gap_foundations()
    locks = build_l1_locks()
    fold = g["static_fold_collapse"]
    assert fold["confirms_collapse"]
    assert abs(fold["lambda_star"] - locks.lambda_star) < 1e-9
    assert abs(fold["star_Delta_meV"] - locks.Delta_star_meV) < 1e-9


def test_cursor_mexican_hat_parallel_to_fold():
    locks = build_l1_locks()
    gl = gl_coefficients(locks)
    n_ref = float(gl["n_ref"]["n_ref"])
    beta = float(gl["beta"]["beta_meV"])
    hat = mexican_hat(locks, beta_meV=beta, n_ref=n_ref)
    ch = run_parallel_channels()
    assert hat["catalyzon_dressed"]["parallel_channel_only"]
    assert hat["gpe_sign_shift"]["well_deepens"]
    assert abs(hat["catalyzon_dressed"]["barrier_reduction_pct"] - ch["mexican_hat"]["barrier_reduction_pct"]) < 0.6


def test_cursor_cavity_and_schwinger():
    cav = apply_F_moire_to_A_geom()
    s = boxed_schwinger_rate()
    assert cav["not_the_option_A_star_static_lock"]
    assert s["ratio_is_alpha_star"]
    locks = build_l1_locks()
    assert abs(s["T_star_K"] - locks.T_star_K) < 1e-9
