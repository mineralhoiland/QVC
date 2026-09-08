"""Tests for L1 closure suite: catalyzon 1–5, CS, C1, Track B, Villain/C7, Hopf."""

from __future__ import annotations

import math

import numpy as np

from qvc_l1.c1_metric import qi_wu_zhang_metric_kappa, run_c1_metric
from qvc_l1.cs_eom import background_a, k_c_anyon
from qvc_l1.hopf_integral import hopfion_n_hat, s3_whitehead_hopf, stereographic_hopf_n
from qvc_l1.linearize import galerkin_basis, tqgl_energy
from qvc_l1.locks import FOLKLORE_PCT, build_l1_locks
from qvc_l1.track_b import run_track_b
from qvc_l1.villain_c7 import villain_reduction_card
from qvc.tqgl.gpe import hopfion_initial, make_grid, zpf_potential


def test_energy_finite():
    locks = build_l1_locks()
    grid = make_grid(locks.tqgl, n_grid=32, L_um=0.200)
    psi = hopfion_initial(grid, locks.tqgl, noise_amp=0.0, seed=0)
    V = zpf_potential(1.0, 1.0, grid.gauss_ring, locks.E_vac_meV)
    F = tqgl_energy(psi, grid, locks.tqgl, V_zpf=V)
    assert np.isfinite(F)


def test_galerkin_orthonormal():
    locks = build_l1_locks()
    grid = make_grid(locks.tqgl, n_grid=32, L_um=0.200)
    b = galerkin_basis(grid, locks.tqgl)
    Phi, dA = b["Phi"], b["dA"]
    G = np.zeros((Phi.shape[0], Phi.shape[0]), dtype=complex)
    for i in range(Phi.shape[0]):
        for j in range(Phi.shape[0]):
            G[i, j] = np.sum(np.conj(Phi[i]) * Phi[j]) * dA
    assert np.allclose(G, np.eye(Phi.shape[0]), atol=1e-6)


def test_anyon_kc_lower():
    assert k_c_anyon(5) < 2.0 / np.pi


def test_track_b_not_holographic():
    tb = run_track_b()
    assert abs(tb["t_over_U"] - 1.0) < 1e-9
    assert tb["status"] == "matched_not_hubbard"
    assert tb["not_BM_Schrieffer_Wolff"] is True
    assert tb["t_lambda_over_U"] < 1.0


def test_qwz_chern_order_one():
    toy = qi_wu_zhang_metric_kappa(Delta_meV=0.6, lambda_M_um=0.013, n_k=24)
    assert abs(toy["chern_numerical"] - 1.0) < 0.35


def test_c1_branch_kappa_drops_toward_fold():
    c1 = run_c1_metric(n_alpha=12)
    br = c1["branch_frozen_toy"]
    assert br["kappa_fold_meV"] < br["kappa_uv_meV"]
    assert br["K_always_above_Kc_Tstar"] is True
    assert br["KT_hit_before_fold_at_Tstar"] is False


def test_hopfion_unit_length():
    ax = np.linspace(-2.0, 2.0, 8)
    Z, Y, X = np.meshgrid(ax, ax, ax, indexing="ij")
    nx, ny, nz = hopfion_n_hat(X, Y, Z, R=1.0, r=0.35)
    n2 = nx * nx + ny * ny + nz * nz
    assert np.allclose(n2, 1.0, atol=1e-5)
    sx, sy, sz = stereographic_hopf_n(X, Y, Z)
    s2 = sx * sx + sy * sy + sz * sz
    assert np.allclose(s2, 1.0, atol=1e-5)


def test_s3_whitehead_unit_hopf():
    s3w = s3_whitehead_hopf(n_chi=48)
    assert abs(s3w["Q_Hopf_map"] - 1.0) < 5e-4


def test_cs_analytic_flux_attachment():
    locks = build_l1_locks()
    grid = make_grid(locks.tqgl, n_grid=32, L_um=0.200)
    ax, ay = background_a(grid, locks.Q_H)
    rho = np.sqrt(grid.X**2 + grid.Y**2)
    mask = rho > 2.0 * grid.dx
    vx = locks.Q_H * (-np.sin(grid.theta) / np.maximum(rho, 1e-16)) - ax
    vy = locks.Q_H * (np.cos(grid.theta) / np.maximum(rho, 1e-16)) - ay
    rms = float(np.sqrt(np.mean(vx[mask] ** 2 + vy[mask] ** 2)))
    assert rms < 1e-12


def test_villain_kc_algebra():
    card = villain_reduction_card(K=2.0 / math.pi, T_eff=1.0, kappa=2.0 / math.pi)
    assert abs(card["eta_spinwave"] - 0.25) < 1e-12


def test_no_folklore_92():
    assert FOLKLORE_PCT == 92.0
