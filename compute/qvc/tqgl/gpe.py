"""Split-step Fourier 2D GPE (truncated TQGL; no CS / Berry / Skyrme in EOM).

Units: energy meV, length µm, time ps.

    i ℏ ∂ψ/∂t = [ −(ℏ²/2m*) ∇² + α |ψ|² + β |ψ|⁴ + V_ZPF(r,t) + V_ph(r,t) ] ψ

Chern–Simons, Berry connection, and Skyrme/Hopf terms are *not* discretised.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from qvc.tqgl.live_gap import ring_weight
from qvc.tqgl.locks import FrozenTQGL, HBAR_MEV_PS


@dataclass
class GPEGrid:
    """Uniform 2D periodic box and ring geometry."""

    N: int
    L_um: float
    R_um: float
    a_H_um: float
    x: np.ndarray
    y: np.ndarray
    X: np.ndarray
    Y: np.ndarray
    dx: float
    K2: np.ndarray
    gauss_ring: np.ndarray
    weight: np.ndarray
    theta: np.ndarray


def make_grid(
    lock: FrozenTQGL,
    *,
    n_grid: int = 128,
    L_um: float = 0.200,
) -> GPEGrid:
    """Square box of side ``L_um`` (default 200 nm) with ring at major radius R."""
    dx = L_um / n_grid
    x = np.linspace(-0.5 * L_um, 0.5 * L_um - dx, n_grid)
    y = x.copy()
    X, Y = np.meshgrid(x, y, indexing="xy")
    kx = 2.0 * np.pi * np.fft.fftfreq(n_grid, d=dx)
    KX, KY = np.meshgrid(kx, kx, indexing="xy")
    rho = np.sqrt(X * X + Y * Y)
    gauss = np.exp(-((rho - lock.R_um) ** 2) / (2.0 * lock.r_um * lock.r_um))
    return GPEGrid(
        N=n_grid,
        L_um=L_um,
        R_um=lock.R_um,
        a_H_um=lock.r_um,
        x=x,
        y=y,
        X=X,
        Y=Y,
        dx=dx,
        K2=KX * KX + KY * KY,
        gauss_ring=gauss,
        weight=ring_weight(X, Y, lock.R_um, lock.r_um),
        theta=np.arctan2(Y, X),
    )


def kinetic_coeff_mev_um2(lock: FrozenTQGL) -> float:
    """ℏ²/(2m*) from healing a_H: ξ = a_H = ℏ/√(2m*|α|) with |α|=Δ₀."""
    return lock.Delta0_meV * (lock.r_um ** 2)


def kinetic_propagator(
    grid: GPEGrid, dt_ps: float, ke_coef: float, hbar: float = HBAR_MEV_PS
) -> np.ndarray:
    """Half-step kinetic factor exp(−i KE k² dt / 2ℏ)."""
    return np.exp(-1j * ke_coef * grid.K2 * (0.5 * dt_ps) / hbar)


def hopfion_initial(
    grid: GPEGrid,
    lock: FrozenTQGL,
    *,
    n0: float = 1.0,
    core_depletion: float = 0.10,
    noise_amp: float = 0.03,
    seed: int = 42,
) -> np.ndarray:
    """Ring hopfion: core depletion, phase Q_H θ, small complex noise."""
    amp = np.sqrt(n0) * (1.0 - core_depletion * grid.gauss_ring)
    phase = np.exp(1j * lock.Q_H * grid.theta)
    rng = np.random.default_rng(seed)
    noise = noise_amp * (
        rng.standard_normal(grid.gauss_ring.shape)
        + 1j * rng.standard_normal(grid.gauss_ring.shape)
    )
    return (amp * phase + noise).astype(np.complex128)


def nl_potential(psi: np.ndarray, alpha: float, beta: float) -> np.ndarray:
    """α n + β n²  (cubic-quintic potential multiplying ψ)."""
    n = np.abs(psi) ** 2
    return alpha * n + beta * n * n


def phonon_potential(
    t_ps: float,
    gauss_ring: np.ndarray,
    *,
    g_ph: float,
    omega_ps: float,
) -> np.ndarray:
    return g_ph * np.sin(omega_ps * t_ps) * gauss_ring


def zpf_potential(
    A: float,
    A_ref: float,
    gauss_ring: np.ndarray,
    E_vac_meV: float,
) -> np.ndarray:
    """V_ZPF = −E_vac (A/A_ref) × ring Gaussian. Drive is locked E_vac, not 0.244."""
    scale = 0.0 if A_ref <= 0.0 else A / A_ref
    return -E_vac_meV * scale * gauss_ring


def split_step_advance(
    psi: np.ndarray,
    *,
    prop_K: np.ndarray,
    V: np.ndarray,
    dt_ps: float,
    hbar: float = HBAR_MEV_PS,
) -> np.ndarray:
    """One Strang split-step: kinetic/2, potential, kinetic/2. Real V ⇒ unitary."""
    psi_k = np.fft.fft2(psi)
    psi_k *= prop_K
    psi = np.fft.ifft2(psi_k)
    psi *= np.exp(-1j * V * dt_ps / hbar)
    psi_k = np.fft.fft2(psi)
    psi_k *= prop_K
    return np.fft.ifft2(psi_k)


def l2_norm(psi: np.ndarray, dx: float) -> float:
    return float(np.sqrt(np.sum(np.abs(psi) ** 2) * dx * dx))
