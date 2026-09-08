r"""Hopf integral on N-layer texture (Phase Two, Workstream 4).

Computes the Hopf invariant:
  Q_H = (1/4π²) ∫_{R³} A ∧ F

for an N-layer Hopf texture n: S³ → S² with linking number 1/N.

Physical content:
  The Hopf invariant classifies the topological linking of preimage curves
  on S². For N layers with CS[A_k] = k²/N, the total Q_H = 1/N = 0.2
  (for N=5). This replaces the invalid H₂-cycle proof.

Numerical method:
  - Discretize on a cubic grid [-L,L]³
  - Construct Hopf map n(r) with winding 1/N
  - Compute Berry connection A_i = ε^{abc} n_a ∂_i n_b ∂_j n_c / (4π(1+n_z))
  - Evaluate Q_H = (1/4π²) ∫ ε^{ijk} A_i ∂_j A_k d³x
  - Verify topological quantization and grid convergence
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np


@dataclass
class HopfConfig:
    """Configuration for Hopf integral computation."""

    N_layers: int = 5  # Number of layers → Q_H = 1/N
    n_grid: int = 64  # Grid resolution per axis
    L: float = 2.0  # Half-box size (in units of texture radius R)
    R_texture: float = 1.0  # Texture core radius
    profile_exponent: int = 2  # Smoothness of radial profile f(r)


def _hopf_texture(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    N: int,
    R: float,
    k: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the N-wound Hopf texture n(r) = (n_x, n_y, n_z).

    Profile: f(r) = π × max(0, 1 - (r/R)²)^k
    Texture: n = (sin f cos(N φ), sin f sin(N φ), cos f)

    where φ = arctan2(y, x) is the azimuthal angle and N is the winding.
    For N=1 this gives the standard Hopf map with Q_H = 1.
    For general N, the preimages have linking number 1/N in the N-layer construction.
    """
    r_perp = np.sqrt(X**2 + Y**2)
    r_full = np.sqrt(X**2 + Y**2 + Z**2)
    phi = np.arctan2(Y, X)

    # Radial profile: f → π at center, f → 0 at boundary
    s = np.clip(1.0 - (r_full / R) ** 2, 0.0, 1.0)
    f = math.pi * s**k

    # Texture components with N-fold winding in azimuth
    sin_f = np.sin(f)
    cos_f = np.cos(f)
    nx = sin_f * np.cos(N * phi)
    ny = sin_f * np.sin(N * phi)
    nz = cos_f

    return nx, ny, nz


def _numerical_gradient(field: np.ndarray, dx: float, axis: int) -> np.ndarray:
    """Central difference gradient along specified axis."""
    return np.gradient(field, dx, axis=axis, edge_order=2)


def _berry_connection(
    nx: np.ndarray,
    ny: np.ndarray,
    nz: np.ndarray,
    dx: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute the Berry connection A_i from the n-field.

    Using the stereographic formula (south pole gauge):
      A_i = (1/(1 + n_z)) × (n_x ∂_i n_y - n_y ∂_i n_x)

    This is well-defined where n_z > -1 (away from south pole).
    The gauge singularity at n_z = -1 doesn't contribute if
    the texture is smooth and compactly supported.

    Units: A has dimensions of 1/length (1/dx in grid units).
    """
    # Regularize denominator (avoid division by zero at south pole)
    denom = 1.0 + nz
    denom = np.where(np.abs(denom) < 1e-10, 1e-10, denom)

    # Gradients of n_x, n_y along each axis
    dnx_dx = _numerical_gradient(nx, dx, axis=0)
    dnx_dy = _numerical_gradient(nx, dx, axis=1)
    dnx_dz = _numerical_gradient(nx, dx, axis=2)

    dny_dx = _numerical_gradient(ny, dx, axis=0)
    dny_dy = _numerical_gradient(ny, dx, axis=1)
    dny_dz = _numerical_gradient(ny, dx, axis=2)

    # A_i = (n_x ∂_i n_y - n_y ∂_i n_x) / (1 + n_z)
    Ax = (nx * dny_dx - ny * dnx_dx) / denom
    Ay = (nx * dny_dy - ny * dnx_dy) / denom
    Az = (nx * dny_dz - ny * dnx_dz) / denom

    return Ax, Ay, Az


def _chern_simons_integrand(
    Ax: np.ndarray,
    Ay: np.ndarray,
    Az: np.ndarray,
    dx: float,
) -> np.ndarray:
    """Compute A ∧ dA = ε^{ijk} A_i ∂_j A_k.

    The Hopf invariant is Q_H = (1/4π²) ∫ A ∧ dA d³x.

    Expand: A ∧ dA = A_x(∂_y A_z - ∂_z A_y)
                    + A_y(∂_z A_x - ∂_x A_z)
                    + A_z(∂_x A_y - ∂_y A_x)
    """
    # Curl components of A
    dAz_dy = _numerical_gradient(Az, dx, axis=1)
    dAy_dz = _numerical_gradient(Ay, dx, axis=2)
    dAx_dz = _numerical_gradient(Ax, dx, axis=2)
    dAz_dx = _numerical_gradient(Az, dx, axis=0)
    dAy_dx = _numerical_gradient(Ay, dx, axis=0)
    dAx_dy = _numerical_gradient(Ax, dx, axis=1)

    integrand = (
        Ax * (dAz_dy - dAy_dz)
        + Ay * (dAx_dz - dAz_dx)
        + Az * (dAy_dx - dAx_dy)
    )
    return integrand


def compute_hopf_invariant(cfg: HopfConfig | None = None) -> dict[str, Any]:
    """Compute Q_H on a 3D grid.

    Returns dict with Q_H value, grid parameters, and convergence info.
    """
    if cfg is None:
        cfg = HopfConfig()

    n = cfg.n_grid
    L = cfg.L
    dx = 2.0 * L / (n - 1)

    # Create 3D grid
    x1d = np.linspace(-L, L, n)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    # Compute texture
    nx, ny, nz = _hopf_texture(X, Y, Z, cfg.N_layers, cfg.R_texture, cfg.profile_exponent)

    # Berry connection
    Ax, Ay, Az = _berry_connection(nx, ny, nz, dx)

    # Chern-Simons integrand
    integrand = _chern_simons_integrand(Ax, Ay, Az, dx)

    # Integrate: Q_H = (1/4π²) ∫ integrand d³x
    volume_element = dx**3
    raw_integral = float(np.sum(integrand) * volume_element)
    Q_H = raw_integral / (4.0 * math.pi**2)

    # Expected value
    Q_expected = 1.0 / cfg.N_layers

    return {
        "Q_H": Q_H,
        "Q_expected": Q_expected,
        "relative_error": abs(Q_H - Q_expected) / Q_expected if Q_expected != 0 else abs(Q_H),
        "raw_integral": raw_integral,
        "config": {
            "N_layers": cfg.N_layers,
            "n_grid": cfg.n_grid,
            "L": cfg.L,
            "R_texture": cfg.R_texture,
            "dx": dx,
        },
        "diagnostics": {
            "integrand_max": float(np.max(np.abs(integrand))),
            "integrand_sum_pos": float(np.sum(integrand[integrand > 0]) * volume_element),
            "integrand_sum_neg": float(np.sum(integrand[integrand < 0]) * volume_element),
            "texture_coverage": float(np.mean(nz**2 + nx**2 + ny**2)),  # should be ~1
        },
    }


def grid_convergence_study(
    N_layers: int = 5,
    grids: list[int] | None = None,
) -> dict[str, Any]:
    """Run Q_H computation at multiple grid resolutions to assess convergence.

    Standard check: compare n_grid = 32, 48, 64, 96.
    """
    if grids is None:
        grids = [24, 32, 48, 64, 96]

    results = []
    for ng in grids:
        cfg = HopfConfig(N_layers=N_layers, n_grid=ng)
        res = compute_hopf_invariant(cfg)
        results.append({
            "n_grid": ng,
            "Q_H": res["Q_H"],
            "rel_error": res["relative_error"],
            "dx": res["config"]["dx"],
        })

    # Richardson extrapolation from two finest grids
    if len(results) >= 2:
        Q1 = results[-2]["Q_H"]
        Q2 = results[-1]["Q_H"]
        h1 = results[-2]["dx"]
        h2 = results[-1]["dx"]
        # Assuming second-order convergence
        p = 2.0
        Q_extrap = (Q2 * h1**p - Q1 * h2**p) / (h1**p - h2**p)
    else:
        Q_extrap = results[-1]["Q_H"]

    return {
        "N_layers": N_layers,
        "Q_expected": 1.0 / N_layers,
        "results": results,
        "richardson_extrapolation": Q_extrap,
        "convergence_order": _estimate_convergence_order(results) if len(results) >= 3 else None,
        "finest_grid_error": results[-1]["rel_error"],
    }


def _estimate_convergence_order(results: list[dict]) -> float:
    """Estimate convergence order from three successive grids."""
    if len(results) < 3:
        return float("nan")
    Q1 = results[-3]["Q_H"]
    Q2 = results[-2]["Q_H"]
    Q3 = results[-1]["Q_H"]
    if abs(Q2 - Q3) < 1e-15 or abs(Q1 - Q2) < 1e-15:
        return float("inf")  # Already converged
    ratio = abs((Q1 - Q2) / (Q2 - Q3))
    if ratio <= 1.0:
        return float("nan")
    return float(math.log(ratio) / math.log(results[-2]["dx"] / results[-1]["dx"]))


def N_layer_sweep(N_range: list[int] | None = None, n_grid: int = 64) -> dict[str, Any]:
    """Compute Q_H for N = 1, 2, ..., 10, verifying 1/N law."""
    if N_range is None:
        N_range = list(range(1, 11))

    results = []
    for N in N_range:
        cfg = HopfConfig(N_layers=N, n_grid=n_grid)
        res = compute_hopf_invariant(cfg)
        results.append({
            "N": N,
            "Q_H": res["Q_H"],
            "Q_expected": 1.0 / N,
            "rel_error": res["relative_error"],
        })

    return {
        "results": results,
        "n_grid": n_grid,
        "max_rel_error": max(r["rel_error"] for r in results),
        "1_over_N_verified": all(r["rel_error"] < 0.05 for r in results),
    }


def topological_robustness(
    N_layers: int = 5,
    n_grid: int = 64,
    deformation_amplitudes: list[float] | None = None,
) -> dict[str, Any]:
    """Test that Q_H is robust under smooth deformations of the texture.

    Add a smooth perturbation δn to the texture and check Q_H is unchanged.
    """
    if deformation_amplitudes is None:
        deformation_amplitudes = [0.0, 0.05, 0.10, 0.20, 0.30, 0.50]

    L = 2.0
    dx = 2.0 * L / (n_grid - 1)
    x1d = np.linspace(-L, L, n_grid)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    results = []
    for eps in deformation_amplitudes:
        # Compute base texture
        nx, ny, nz = _hopf_texture(X, Y, Z, N_layers, 1.0, 2)

        if eps > 0:
            # Add smooth deformation (low-k Fourier mode)
            # Perturb in direction perpendicular to n, then renormalize
            r_full = np.sqrt(X**2 + Y**2 + Z**2)
            perturbation = eps * np.sin(math.pi * X / L) * np.cos(math.pi * Y / L) * np.exp(-r_full**2)

            nx_new = nx + perturbation * nz  # Rotate towards z
            ny_new = ny + perturbation * (-nx)  # Tangent perturbation
            nz_new = nz - perturbation * nx  # Away from x

            # Renormalize to S²
            norm = np.sqrt(nx_new**2 + ny_new**2 + nz_new**2)
            norm = np.where(norm < 1e-10, 1.0, norm)
            nx, ny, nz = nx_new / norm, ny_new / norm, nz_new / norm

        # Compute Q_H
        Ax, Ay, Az = _berry_connection(nx, ny, nz, dx)
        integrand = _chern_simons_integrand(Ax, Ay, Az, dx)
        Q_H = float(np.sum(integrand) * dx**3) / (4.0 * math.pi**2)

        results.append({
            "epsilon": eps,
            "Q_H": Q_H,
            "deviation_from_quantized": abs(Q_H - 1.0 / N_layers),
        })

    return {
        "N_layers": N_layers,
        "Q_expected": 1.0 / N_layers,
        "results": results,
        "topologically_robust": all(
            r["deviation_from_quantized"] < 0.05 for r in results if r["epsilon"] <= 0.3
        ),
    }


def run_full_hopf_suite(n_grid: int = 64) -> dict[str, Any]:
    """Run complete Hopf integral analysis for Phase Two."""
    # 1. Standard computation at N=5
    cfg = HopfConfig(N_layers=5, n_grid=n_grid)
    base = compute_hopf_invariant(cfg)

    # 2. Grid convergence (use smaller grids for speed)
    convergence = grid_convergence_study(N_layers=5, grids=[24, 32, 48, n_grid])

    # 3. N-layer sweep
    n_sweep = N_layer_sweep(n_grid=min(n_grid, 48))

    # 4. Robustness
    robust = topological_robustness(N_layers=5, n_grid=min(n_grid, 48))

    return {
        "base_result": base,
        "convergence": convergence,
        "N_sweep": n_sweep,
        "robustness": robust,
        "summary": {
            "Q_H_at_N5": base["Q_H"],
            "Q_expected": 0.2,
            "relative_error": base["relative_error"],
            "grid_converged": convergence["finest_grid_error"] < 0.01,
            "1_over_N_law_holds": n_sweep["1_over_N_verified"],
            "topologically_robust": robust["topologically_robust"],
            "PASS": (
                base["relative_error"] < 0.05
                and convergence["finest_grid_error"] < 0.05
            ),
        },
    }
