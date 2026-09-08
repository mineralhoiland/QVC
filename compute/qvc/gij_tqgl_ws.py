"""TQGL + Wigner–Seitz overlap → G_ij (no hand-inserted democratic matrix).

Construction
------------
N collective TQGL modes Φ_i live on the moiré Wigner–Seitz cell
(hexagon of circumradius R_WS, period λ_M). Centres sit at a regular
N-gon. Each profile is a normalised 2D Gaussian of width σ.

The vacuum kernel is the same Macdonald function that appears in F(x):

    K(r, r') = K₀(|r − r'| / ξ_K),   ξ_K = R  (torus major / coherence).

Then

    G_ij = ∫∫ Φ_i(r) K(r,r') Φ_j(r') d²r d²r'    (i ≠ j),
    G_ii = 0.

Positive long-range K ⇒ positive off-diagonals ⇒ λ_min < 0 without
inserting g₀(J−I). If ξ_K is large compared with site spacing, Spec(G)
approaches the democratic lock.

This is a **computed microscopic model**, not a variation of the full
torsion action (H-TEGR1 remains open). Status: PROGRAM_COMPUTED.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.special import k0

from qvc.coupling_matrix import democratic_coupling_matrix, democratic_spectrum
from qvc.params import PARAMS


def ws_hexagon_mask(X: np.ndarray, Y: np.ndarray, R_hex: float) -> np.ndarray:
    """True inside a pointy-top hexagon of circumradius R_hex."""
    # Pointy-top: |x| ≤ R√3/2, |y| ≤ R, and |y| + |x|/√3 ≤ R
    ax, ay = np.abs(X), np.abs(Y)
    return (ay <= R_hex) & (ax <= R_hex * math.sqrt(3.0) / 2.0) & (
        ay + ax / math.sqrt(3.0) <= R_hex
    )


def mode_centres_polygon(N: int, radius: float) -> np.ndarray:
    """Regular N-gon centres in the plane, starting on +x."""
    ang = 2.0 * math.pi * np.arange(N) / N
    return np.stack([radius * np.cos(ang), radius * np.sin(ang)], axis=1)


def gaussian_modes(
    X: np.ndarray,
    Y: np.ndarray,
    centres: np.ndarray,
    sigma: float,
    mask: np.ndarray,
    dA: float,
) -> np.ndarray:
    """Return (N, n_pts) normalised Φ_i sampled on masked grid."""
    xs = X[mask]
    ys = Y[mask]
    N = centres.shape[0]
    Phi = np.zeros((N, xs.size), dtype=float)
    for i, c in enumerate(centres):
        Phi[i] = np.exp(-((xs - c[0]) ** 2 + (ys - c[1]) ** 2) / (2.0 * sigma**2))
        nrm = math.sqrt(float(np.sum(Phi[i] ** 2) * dA))
        if nrm <= 0:
            raise RuntimeError("mode normalisation vanished")
        Phi[i] /= nrm
    return Phi


def kernel_matrix_K0(pts: np.ndarray, xi: float) -> np.ndarray:
    """K_ij = K₀(|r_i−r_j|/ξ) with K₀(0) regularised to K₀(ε)."""
    d = np.sqrt(((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2))
    z = d / xi
    z = np.maximum(z, 1e-6)
    return np.asarray(k0(z), dtype=float)


def overlap_G(
    Phi: np.ndarray,
    K: np.ndarray,
    dA: float,
    *,
    zero_diagonal: bool = True,
) -> np.ndarray:
    """G = Φ K Φᵀ dA², then zero diagonal (vacuum hybridisation only)."""
    # Phi: (N, M), K: (M, M)
    GK = Phi @ K  # (N, M)
    G = (GK @ Phi.T) * (dA * dA)
    G = 0.5 * (G + G.T)
    if zero_diagonal:
        np.fill_diagonal(G, 0.0)
    return G


def run_tqgl_ws_G(
    *,
    N: int | None = None,
    n_grid: int = 48,
    sigma_over_lambda: float = 0.18,
    centre_over_Rhex: float = 0.45,
    xi_over_lambda: float | None = None,
) -> dict[str, Any]:
    """Compute G_ij on the MATBG WS cell and compare to democratic Spec."""
    p = PARAMS
    N = p.N_layers if N is None else N
    lam = p.lambda_M_um  # 13 nm
    R_hex = lam / math.sqrt(3.0)  # hex circumradius ~ λ_M/√3
    # Coherence / torus major as kernel range (default R_preprint / λ_M)
    if xi_over_lambda is None:
        xi = p.R_preprint_um
    else:
        xi = xi_over_lambda * lam

    lim = R_hex * 1.05
    axis = np.linspace(-lim, lim, n_grid)
    X, Y = np.meshgrid(axis, axis, indexing="xy")
    dA = (axis[1] - axis[0]) ** 2
    mask = ws_hexagon_mask(X, Y, R_hex)
    pts = np.stack([X[mask], Y[mask]], axis=1)

    centres = mode_centres_polygon(N, centre_over_Rhex * R_hex)
    sigma = sigma_over_lambda * lam
    Phi = gaussian_modes(X, Y, centres, sigma, mask, dA)
    K = kernel_matrix_K0(pts, xi)
    G = overlap_G(Phi, K, dA, zero_diagonal=True)

    eigs = np.sort(np.linalg.eigvalsh(G))
    off = G[np.triu_indices(N, k=1)]
    g0_eff = float(np.mean(off)) if off.size else 0.0
    G_dem = democratic_coupling_matrix(N, g0_eff)
    spec = democratic_spectrum(N, g0_eff)
    fro = float(np.linalg.norm(G - G_dem, ord="fro"))
    rel_fro = fro / (float(np.linalg.norm(G_dem, ord="fro")) + 1e-18)

    return {
        "status": "PROGRAM_COMPUTED_not_torsion_variation",
        "N": N,
        "lambda_M_nm": lam * 1e3,
        "R_hex_nm": R_hex * 1e3,
        "sigma_nm": sigma * 1e3,
        "xi_kernel_nm": xi * 1e3,
        "n_grid": n_grid,
        "n_pts": int(mask.sum()),
        "G": G.tolist(),
        "eigenvalues": eigs.tolist(),
        "lambda_max": float(eigs[-1]),
        "lambda_min": float(eigs[0]),
        "lambda_min_mult_numeric": int(np.sum(np.isclose(eigs, eigs[0], atol=1e-8))),
        "trace": float(np.trace(G)),
        "hermitian": bool(np.allclose(G, G.T)),
        "g0_eff_mean_offdiag": g0_eff,
        "democratic_lambda_max": spec["lambda_max"],
        "democratic_lambda_min": spec["lambda_min"],
        "rel_err_lambda_max": abs(float(eigs[-1]) - spec["lambda_max"])
        / (abs(spec["lambda_max"]) + 1e-18),
        "rel_err_lambda_min": abs(float(eigs[0]) - spec["lambda_min"])
        / (abs(spec["lambda_min"]) + 1e-18),
        "frobenius_to_democratic": fro,
        "rel_frobenius_to_democratic": rel_fro,
        "offdiag_cv": float(np.std(off) / (abs(g0_eff) + 1e-18)),
        "acceptance": {
            "hermitian": bool(np.allclose(G, G.T)),
            "diag_zero": bool(np.allclose(np.diag(G), 0.0)),
            "trace_zero": abs(float(np.trace(G))) < 1e-12,
            "lambda_min_lt_0": float(eigs[0]) < 0.0,
            "near_democratic": rel_fro < 0.15,
        },
        "notes": (
            "G constructed from TQGL Gaussians × K₀ kernel on WS hexagon. "
            "Democratic form is an OUTPUT of long-range K, not an input."
        ),
    }
