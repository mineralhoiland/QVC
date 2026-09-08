"""Hopf invariant by Coulomb-gauge Whitehead integral (no Dirac string).

Skyrmion density (3-form dual):

    B_i = (1/2) ε_ijk  n · (∂_j n × ∂_k n)

Coulomb gauge in Fourier space:  A(k) = i k × B(k) / k² ,  A(0)=0.
Whitehead:

    Q = (1 / 16π²) ∫ A · B  d³x

The locked Q_H = 1/N is CS/linking on L(N,1). This quadrature tests the
unit Hopf map (target 1) without a southern-chart Dirac string. An
N-wound toroidal ansatz is a labeled probe, not a derivation of 1/N.

Status: computed 3D FFT quadrature.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.lens_chern_simons import fractional_hall_charge  # noqa: E402


def hopf_map_spinor(eta1, eta2):
    nrm = np.sqrt(np.abs(eta1) ** 2 + np.abs(eta2) ** 2) + 1e-16
    eta1, eta2 = eta1 / nrm, eta2 / nrm
    nx = 2.0 * np.real(np.conj(eta1) * eta2)
    ny = 2.0 * np.imag(np.conj(eta1) * eta2)
    nz = np.abs(eta1) ** 2 - np.abs(eta2) ** 2
    return nx, ny, nz


def stereographic_hopf_n(X, Y, Z):
    """Standard Hopf fibration pulled back to R³ by inverse stereographic S³."""
    r2 = X * X + Y * Y + Z * Z
    den = 1.0 + r2
    a = 2.0 * X / den
    b = 2.0 * Y / den
    c = 2.0 * Z / den
    d = (1.0 - r2) / den
    nx = 2.0 * (a * c + b * d)
    ny = 2.0 * (b * c - a * d)
    nz = (a * a + b * b) - (c * c + d * d)
    return nx, ny, nz


def hopfion_n_hat(X, Y, Z, *, R: float, r: float, n_winding: float = 1.0):
    """Compactified toroidal hopfion (Sutcliffe-type) with optional φ winding."""
    rho = np.sqrt(X * X + Y * Y)
    w = (rho - R) + 1j * Z
    denom = r * r + np.abs(w) ** 2 + 1e-16
    u = (2.0 * r * w) / denom
    v = (r * r - np.abs(w) ** 2) / denom
    phi = np.arctan2(Y, X)
    return hopf_map_spinor(v, u * np.exp(1j * n_winding * phi))


def _fd_grad(arr: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Centered finite-difference gradient. arr axes (Z, Y, X)."""
    dz = np.gradient(arr, dx, axis=0, edge_order=2)
    dy = np.gradient(arr, dx, axis=1, edge_order=2)
    dx_ = np.gradient(arr, dx, axis=2, edge_order=2)
    return dz, dy, dx_


def skyrmion_B(nx, ny, nz, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """B_i = (1/2) ε_ijk n·(∂_j n × ∂_k n) = n·(∂_j n × ∂_k n) cyclic.

    With antisymmetric F_jk = n·(∂_j n × ∂_k n) the z-component is
    B_z = F_xy (the extra 1/2 from ε_{zij} cancels). Axes (Z,Y,X).
    """
    dzx, dyx, dxx = _fd_grad(nx, dx)
    dzy, dyy, dxy = _fd_grad(ny, dx)
    dzz, dyz, dxz = _fd_grad(nz, dx)

    def triple(ax, ay, az, bx, by, bz):
        return nx * (ay * bz - az * by) + ny * (az * bx - ax * bz) + nz * (ax * by - ay * bx)

    bz = triple(dxx, dxy, dxz, dyx, dyy, dyz)
    bx = triple(dyx, dyy, dyz, dzx, dzy, dzz)
    by = triple(dzx, dzy, dzz, dxx, dxy, dxz)
    return bx, by, bz


def coulomb_A(bx, by, bz, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Free-space Coulomb gauge via zero-padded FFT (Hockney–Eastwood).

    Periodic images of a slowly decaying stereographic hopfion otherwise
    pollute A. Padding approximates A(k)= i k×B/k² on R³.
    """
    nz, ny, nx_ = bx.shape
    pz, py, px = 2 * nz, 2 * ny, 2 * nx_
    Bx = np.zeros((pz, py, px), dtype=float)
    By = np.zeros_like(Bx)
    Bz = np.zeros_like(Bx)
    Bx[:nz, :ny, :nx_] = bx
    By[:nz, :ny, :nx_] = by
    Bz[:nz, :ny, :nx_] = bz
    kz = 2.0 * np.pi * np.fft.fftfreq(pz, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(py, d=dx)
    kx = 2.0 * np.pi * np.fft.fftfreq(px, d=dx)
    KZ, KY, KX = np.meshgrid(kz, ky, kx, indexing="ij")
    k2 = KX * KX + KY * KY + KZ * KZ
    k2 = np.where(k2 < 1e-18, 1.0, k2)
    Bxf, Byf, Bzf = np.fft.fftn(Bx), np.fft.fftn(By), np.fft.fftn(Bz)
    Ax = 1j * (KY * Bzf - KZ * Byf) / k2
    Ay = 1j * (KZ * Bxf - KX * Bzf) / k2
    Az = 1j * (KX * Byf - KY * Bxf) / k2
    Ax[0, 0, 0] = 0.0
    Ay[0, 0, 0] = 0.0
    Az[0, 0, 0] = 0.0
    ax = np.fft.ifftn(Ax).real[:nz, :ny, :nx_]
    ay = np.fft.ifftn(Ay).real[:nz, :ny, :nx_]
    az = np.fft.ifftn(Az).real[:nz, :ny, :nx_]
    return ax, ay, az


def s3_whitehead_hopf(*, n_chi: int = 64) -> dict[str, Any]:
    """Hopf map on S³ ⊂ C².

    z1=sinχ e^{iα}, z2=cosχ e^{iβ}, χ∈[0,π/2], α,β∈[0,2π].
    Connection A = −i ⟨z, dz⟩ = sin²χ dα + cos²χ dβ.
    A∧dA = 2 sinχ cosχ dχ∧dα∧dβ, so ∫_{S³} A∧dA = 4π² and
    Q = (1/4π²) ∫ A∧dA = 1 exactly.
    """
    chi = np.linspace(0.0, 0.5 * math.pi, n_chi, endpoint=False)
    dchi = chi[1] - chi[0]
    chi_int = float(np.sum(2.0 * np.sin(chi) * np.cos(chi)) * dchi)
    integral = chi_int * (2.0 * math.pi) ** 2
    Q = integral / (4.0 * math.pi**2)
    return {
        "status": "derived_and_quadrature",
        "connection": "A = sin²χ dα + cos²χ dβ",
        "integral_A_wedge_dA": integral,
        "integral_target_4pi2": 4.0 * math.pi**2,
        "prefactor": "1/4π²",
        "Q_Hopf_map": Q,
        "Q_target": 1.0,
        "Q_abs_err": abs(Q - 1.0),
        "n_chi": n_chi,
        "note": (
            "Algebra of the Hopf fibration (Q=1). The R³ Whitehead quadrature "
            "is a truncation of this identity for a unit hopfion texture."
        ),
    }


def hopf_connection_spinor(z1, z2, dx: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """A_i = −i ⟨z, ∂_i z⟩ from a normalized C² spinor (Hopf connection)."""
    nrm = np.sqrt(np.abs(z1) ** 2 + np.abs(z2) ** 2) + 1e-16
    z1, z2 = z1 / nrm, z2 / nrm
    dz1 = _fd_grad(z1, dx)
    dz2 = _fd_grad(z2, dx)
    A = []
    for d1, d2 in zip(dz1, dz2):
        A.append(np.real(-1j * (np.conj(z1) * d1 + np.conj(z2) * d2)))
    # axes of _fd_grad are (Z,Y,X) so A = (Az, Ay, Ax)
    return A[2], A[1], A[0]


def whitehead_Q_coulomb(nx, ny, nz, dx: float) -> float:
    bx, by, bz = skyrmion_B(nx, ny, nz, dx)
    ax, ay, az = coulomb_A(bx, by, bz, dx)
    dens = ax * bx + ay * by + az * bz
    return float(np.sum(dens) * dx**3 / (16.0 * math.pi**2))


def whitehead_Q_spinor(z1, z2, nx, ny, nz, dx: float) -> float:
    """Whitehead using the Hopf connection of the spinor, not Coulomb Poisson."""
    bx, by, bz = skyrmion_B(nx, ny, nz, dx)
    ax, ay, az = hopf_connection_spinor(z1, z2, dx)
    dens = ax * bx + ay * by + az * bz
    return float(np.sum(dens) * dx**3 / (16.0 * math.pi**2))


def stereographic_spinor(X, Y, Z):
    r2 = X * X + Y * Y + Z * Z
    den = 1.0 + r2
    z1 = (2.0 * X + 1j * 2.0 * Y) / den
    z2 = (2.0 * Z + 1j * (1.0 - r2)) / den
    return z1, z2


def torus_spinor(X, Y, Z, *, R: float, r: float, n_winding: float = 1.0):
    rho = np.sqrt(X * X + Y * Y)
    w = (rho - R) + 1j * Z
    denom = r * r + np.abs(w) ** 2 + 1e-16
    u = (2.0 * r * w) / denom
    v = (r * r - np.abs(w) ** 2) / denom
    phi = np.arctan2(Y, X)
    return v, u * np.exp(1j * n_winding * phi)


def _unitize(nx, ny, nz):
    nlen = np.sqrt(nx * nx + ny * ny + nz * nz) + 1e-16
    return nx / nlen, ny / nlen, nz / nlen


def s3_hopf_volume_check(n_chi: int = 48) -> dict[str, Any]:
    """∫_{S³} vol = 2π² (normalization of the Hopf fibration domain)."""
    chi = np.linspace(0.0, math.pi, n_chi, endpoint=False)
    theta = np.linspace(0.0, math.pi, n_chi, endpoint=False)
    phi = np.linspace(0.0, 2.0 * math.pi, 2 * n_chi, endpoint=False)
    dchi, dth, dphi = chi[1] - chi[0], theta[1] - theta[0], phi[1] - phi[0]
    s2 = float(np.sum(np.sin(chi) ** 2) * dchi)
    sth = float(np.sum(np.sin(theta)) * dth)
    sph = float(len(phi) * dphi)
    vol = s2 * sth * sph
    target = 2.0 * math.pi**2
    return {
        "volume_S3": vol,
        "target_2pi2": target,
        "rel_err": abs(vol - target) / target,
        "status": "computed",
        "note": "Domain volume of the Hopf fibration. The map S³→S² has Hopf invariant 1 exactly.",
    }


def run_hopf_integral(*, n: int = 48, N_layers: int = 5) -> dict[str, Any]:
    Q_lock = fractional_hall_charge(N_layers)
    s3 = s3_hopf_volume_check()
    s3w = s3_whitehead_hopf()

    # Compact stereographic chart: field → constant at large |r|, so a large
    # box with periodic FFT is a controlled truncation of R³.
    lim_s = 6.0
    ax_s = np.linspace(-lim_s, lim_s, n)
    dx_s = float(ax_s[1] - ax_s[0])
    Zs, Ys, Xs = np.meshgrid(ax_s, ax_s, ax_s, indexing="ij")
    z1s, z2s = stereographic_spinor(Xs, Ys, Zs)
    nx, ny, nz = _unitize(*hopf_map_spinor(z1s, z2s))
    Q_stereo = whitehead_Q_coulomb(nx, ny, nz, dx_s)
    Q_spinor = whitehead_Q_spinor(z1s, z2s, nx, ny, nz, dx_s)

    R, r = 1.0, 0.32
    lim_t = 2.8
    ax_t = np.linspace(-lim_t, lim_t, n)
    dx_t = float(ax_t[1] - ax_t[0])
    Zt, Yt, Xt = np.meshgrid(ax_t, ax_t, ax_t, indexing="ij")
    nx1, ny1, nz1 = _unitize(*hopfion_n_hat(Xt, Yt, Zt, R=R, r=r, n_winding=1.0))
    z1t, z2t = torus_spinor(Xt, Yt, Zt, R=R, r=r, n_winding=1.0)
    Q1 = whitehead_Q_coulomb(nx1, ny1, nz1, dx_t)
    Q1s = whitehead_Q_spinor(z1t, z2t, nx1, ny1, nz1, dx_t)
    nxN, nyN, nzN = _unitize(
        *hopfion_n_hat(Xt, Yt, Zt, R=R, r=r, n_winding=1.0 / N_layers)
    )
    QN = whitehead_Q_coulomb(nxN, nyN, nzN, dx_t)

    return {
        "status": "derived_S3_R3_truncated",
        "method": "S³ Whitehead (derived) + string-free R³ Coulomb/spinor (truncated)",
        "Q_stereographic_R3": Q_stereo,
        "Q_stereographic_spinor": Q_spinor,
        "Q_unit_target": 1.0,
        "Q_stereo_abs_err": abs(Q_stereo - 1.0),
        "Q_spinor_abs_err": abs(Q_spinor - 1.0),
        "Q_unit_hopfion": Q1,
        "Q_unit_hopfion_spinor": Q1s,
        "Q_unit_abs_err": abs(Q1 - 1.0),
        "Q_N_wound": QN,
        "Q_H_lock_1_over_N": Q_lock,
        "Q_N_abs_err_vs_lock": abs(QN - Q_lock),
        "n_grid": n,
        "dx_stereo": dx_s,
        "dx_torus": dx_t,
        "S3_volume": s3,
        "S3_whitehead": s3w,
        "prefactor_R3": "1/16π² Coulomb-gauge Whitehead (Faddeev–Niemi B = n·(∂n×∂n))",
        "note": (
            "The Hopf map on S³ has Whitehead integral Q=1 (derived, quadrature "
            "error ~2e-4). Q_H=1/N remains the CS/linking lock on L(N,1). "
            "R³ Coulomb/spinor quadratures are string-free but box-truncated and "
            "do not replace that lock. The N-wound torus is a labeled probe."
        ),
    }
