"""Chern–Simons constraint as a covariant kinetic EOM.

Level k = 2N attaches flux Φ = 2π Q_H = 2π/N. The background

    a_θ = Q_H / ρ ,   a = a_θ ê_θ

is the mean-field constraint, not a dynamical 2+1D gauge field. The
truncated GPE uses −κ ∇². The CS-restored kinetic operator is

    H_kin = κ (−i∇ − a)²
          = κ [ −∇² + 2 i a·∇ + i (∇·a) + |a|² ] .

A Strang split advances half a covariant kinetic RK4 step, the local
potential, then the other kinetic half-step, so a is inside the stepper.
On the hopfion, ∇θ = Q_H ê_θ/ρ, so the superfluid velocity ∇θ − a is
cancelled by flux attachment; rms |j_cov| ≪ rms |j_trunc| is the
diagnostic. Hall conductance remains CS algebra σ_xy ∼ 1/k unless a
transported edge is simulated (it is not).

Status: computed background-CS covariant GPE. Not a dynamical CS theorem.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.lens_chern_simons import (  # noqa: E402
    anyonic_exchange_angle,
    fractional_hall_charge,
)
from qvc.tqgl.gpe import (  # noqa: E402
    hopfion_initial,
    kinetic_coeff_mev_um2,
    kinetic_propagator,
    make_grid,
    nl_potential,
    split_step_advance,
    zpf_potential,
)
from qvc.tqgl.live_gap import coherence_g1, delta_live_from_psi  # noqa: E402
from qvc.tqgl.locks import HBAR_MEV_PS  # noqa: E402

from qvc_l1.locks import L1Locks, build_l1_locks


def background_a(grid, Q_H: float) -> tuple[np.ndarray, np.ndarray]:
    """a_θ = Q_H / ρ  ⇒  a_x = −a_θ sinθ, a_y = a_θ cosθ. Regularized at origin."""
    rho = np.sqrt(grid.X**2 + grid.Y**2)
    a_th = Q_H / np.maximum(rho, 0.5 * grid.dx)
    ax = -a_th * np.sin(grid.theta)
    ay = a_th * np.cos(grid.theta)
    return ax, ay


def _spectral_grad(psi: np.ndarray, grid) -> tuple[np.ndarray, np.ndarray]:
    kx = 2.0 * np.pi * np.fft.fftfreq(grid.N, d=grid.dx)
    KX, KY = np.meshgrid(kx, kx, indexing="xy")
    pk = np.fft.fft2(psi)
    return np.fft.ifft2(1j * KX * pk), np.fft.ifft2(1j * KY * pk)


def covariant_kinetic_apply(
    psi: np.ndarray, ax: np.ndarray, ay: np.ndarray, grid
) -> np.ndarray:
    """Apply (−i∇ − a)² to ψ. Spectral derivatives, periodic box."""
    dx_psi, dy_psi = _spectral_grad(psi, grid)
    lap = np.fft.ifft2(-grid.K2 * np.fft.fft2(psi))
    dax_x, _ = _spectral_grad(ax.astype(np.complex128), grid)
    _, day_y = _spectral_grad(ay.astype(np.complex128), grid)
    div_a = dax_x.real + day_y.real
    a2 = ax * ax + ay * ay
    return -lap + 2j * (ax * dx_psi + ay * dy_psi) + 1j * div_a * psi + a2 * psi


def covariant_k_energy(
    psi: np.ndarray, grid, ax: np.ndarray, ay: np.ndarray, kappa: float
) -> float:
    dx_psi, dy_psi = _spectral_grad(psi, grid)
    Dx = dx_psi - 1j * ax * psi
    Dy = dy_psi - 1j * ay * psi
    dA = grid.dx * grid.dx
    return kappa * float(np.sum(np.abs(Dx) ** 2 + np.abs(Dy) ** 2)) * dA


def current_density(
    psi: np.ndarray, grid, ax: np.ndarray, ay: np.ndarray, kappa: float
) -> tuple[np.ndarray, np.ndarray]:
    """j = 2 κ Im(ψ* (∇ − i a) ψ)."""
    dx_psi, dy_psi = _spectral_grad(psi, grid)
    Dx = dx_psi - 1j * ax * psi
    Dy = dy_psi - 1j * ay * psi
    jx = 2.0 * kappa * np.imag(np.conj(psi) * Dx)
    jy = 2.0 * kappa * np.imag(np.conj(psi) * Dy)
    return jx.real, jy.real


def flux_enclosed(ax: np.ndarray, ay: np.ndarray, grid) -> float:
    """∫ (∂x ay − ∂y ax) d²r."""
    dx_ay = np.gradient(ay, grid.dx, axis=1)
    dy_ax = np.gradient(ax, grid.dx, axis=0)
    return float(np.sum(dx_ay - dy_ax) * grid.dx * grid.dx)


def k_c_anyon(N: int) -> float:
    """Labeled CS-restored KT criterion. Conjecture unless derived from this EOM."""
    return (2.0 / math.pi) * (1.0 - 1.0 / N) ** 2


def _rk4_kinetic(
    psi: np.ndarray,
    ax: np.ndarray,
    ay: np.ndarray,
    grid,
    kappa: float,
    dt_ps: float,
    hbar: float,
) -> np.ndarray:
    """Unitary-leaning RK4 for iℏ ∂t ψ = κ (−i∇−a)² ψ."""
    c = -1j * kappa / hbar

    def f(p: np.ndarray) -> np.ndarray:
        return c * covariant_kinetic_apply(p, ax, ay, grid)

    k1 = f(psi)
    k2 = f(psi + 0.5 * dt_ps * k1)
    k3 = f(psi + 0.5 * dt_ps * k2)
    k4 = f(psi + dt_ps * k3)
    return psi + (dt_ps / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def covariant_split_step(
    psi: np.ndarray,
    *,
    ax: np.ndarray,
    ay: np.ndarray,
    grid,
    kappa: float,
    V: np.ndarray,
    dt_ps: float,
    hbar: float = HBAR_MEV_PS,
) -> np.ndarray:
    """Strang: kinetic/2 (RK4), potential multiply, kinetic/2."""
    psi = _rk4_kinetic(psi, ax, ay, grid, kappa, 0.5 * dt_ps, hbar)
    psi = psi * np.exp(-1j * V * dt_ps / hbar)
    return _rk4_kinetic(psi, ax, ay, grid, kappa, 0.5 * dt_ps, hbar)


def run_cs_eom(
    locks: L1Locks | None = None,
    *,
    n_grid: int = 64,
    t_end_ps: float = 0.25,
    dt_ps: float = 0.002,
) -> dict[str, Any]:
    """Covariant CS-background GPE versus truncated FFT GPE."""
    locks = locks or build_l1_locks()
    lock = locks.tqgl
    grid = make_grid(lock, n_grid=n_grid, L_um=0.200)
    kappa = kinetic_coeff_mev_um2(lock)
    ax, ay = background_a(grid, lock.Q_H)
    zax, zay = np.zeros_like(ax), np.zeros_like(ay)
    flux = flux_enclosed(ax, ay, grid)
    psi0 = hopfion_initial(grid, lock, noise_amp=0.0, seed=0)
    n_ref = float(np.sum(np.abs(psi0) ** 2 * grid.weight) / np.sum(grid.weight))
    dA = grid.dx * grid.dx
    n0_mass = float(np.sum(np.abs(psi0) ** 2) * dA)

    jx_tr, jy_tr = current_density(psi0, grid, zax, zay, kappa)
    jx_cs, jy_cs = current_density(psi0, grid, ax, ay, kappa)
    rms_tr = math.sqrt(float(np.mean(jx_tr**2 + jy_tr**2)))
    rms_cs = math.sqrt(float(np.mean(jx_cs**2 + jy_cs**2)))
    rho = np.sqrt(grid.X**2 + grid.Y**2)
    mask = rho > 2.0 * grid.dx
    vx = lock.Q_H * (-np.sin(grid.theta) / np.maximum(rho, 1e-16)) - ax
    vy = lock.Q_H * (np.cos(grid.theta) / np.maximum(rho, 1e-16)) - ay
    attach_rms = math.sqrt(float(np.mean((vx[mask] ** 2 + vy[mask] ** 2))))

    prop = kinetic_propagator(grid, dt_ps, kappa)
    n_steps = int(round(t_end_ps / dt_ps))
    alpha_nl = -lock.Delta0_meV
    beta_nl = lock.Delta0_meV
    Vz = zpf_potential(1.0, 1.0, grid.gauss_ring, lock.E_vac_meV)

    psi_tr = psi0.copy()
    psi_cs = psi0.copy()
    for _ in range(n_steps):
        Vtr = nl_potential(psi_tr, alpha_nl, beta_nl) + Vz
        psi_tr = split_step_advance(
            psi_tr, prop_K=prop, V=Vtr, dt_ps=dt_ps, hbar=HBAR_MEV_PS
        )
        Vcs = nl_potential(psi_cs, alpha_nl, beta_nl) + Vz
        psi_cs = covariant_split_step(
            psi_cs, ax=ax, ay=ay, grid=grid, kappa=kappa, V=Vcs, dt_ps=dt_ps
        )

    mass_cs = float(np.sum(np.abs(psi_cs) ** 2) * dA)
    N = locks.N_layers
    d_live_tr = delta_live_from_psi(
        psi_tr, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
    )
    d_live_cs = delta_live_from_psi(
        psi_cs, weight=grid.weight, n_ref=n_ref, Delta0_meV=lock.Delta0_meV
    )
    flux_target = 2.0 * math.pi * lock.Q_H
    return {
        "status": "computed_covariant_stepper",
        "a_in_stepper": True,
        "k_cs": 2 * N,
        "Q_H": fractional_hall_charge(N),
        "theta_stat": anyonic_exchange_angle(N),
        "K_c_anyon": k_c_anyon(N),
        "K_c_bosonic": 2.0 / math.pi,
        "flux_enclosed": flux,
        "flux_target_2pi_QH": flux_target,
        "flux_rel_err": abs(flux - flux_target) / (flux_target + 1e-18),
        "E_kinetic_truncated_meV": covariant_k_energy(psi_tr, grid, zax, zay, kappa),
        "E_kinetic_covariant_on_truncated_meV": covariant_k_energy(
            psi_tr, grid, ax, ay, kappa
        ),
        "E_kinetic_covariant_evolved_meV": covariant_k_energy(
            psi_cs, grid, ax, ay, kappa
        ),
        "j_rms_truncated_IC": rms_tr,
        "j_rms_covariant_IC": rms_cs,
        "j_cancellation_ratio": rms_cs / (rms_tr + 1e-18),
        "analytic_attachment_rms": attach_rms,
        "psi_current_polluted_by_QH_branch_cut": True,
        "mass_rel_drift_cs": abs(mass_cs - n0_mass) / (n0_mass + 1e-18),
        "Delta_live_truncated_meV": d_live_tr,
        "Delta_live_covariant_meV": d_live_cs,
        "g1_truncated": coherence_g1(psi_tr),
        "g1_covariant": coherence_g1(psi_cs),
        "t_end_ps": t_end_ps,
        "n_grid": n_grid,
        "Hall_status": "algebra_not_transported_edge",
        "sigma_xy_over_e2_h": 1.0 / (2 * N),
        "note": (
            "Covariant Strang+RK4 puts a inside the kinetic EOM. "
            "Flux-attachment cancellation is the IC current diagnostic. "
            "Hall conductance is still the CS algebra 1/k, not an edge transport. "
            "Anyon K_c is labeled, not derived from this stepper."
        ),
        "xy_kt_still_truncated_default": True,
    }
