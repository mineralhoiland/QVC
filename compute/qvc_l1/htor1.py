"""Step 4: H-TOR1 torsion-overlap vertex on Galerkin modes.

    G_ij ≟ ∫_WS T^μ_{νρ} Φ_i^{*ν} Φ_j^ρ d²x ,   T ∝ ε ∂ q_H

We take a 2D proxy: q_H density is the hopfion vorticity / phase Jacobian

    q_H(r) = (1/2π) Im( ∂x ψ* ∂y ψ − ∂y ψ* ∂x ψ ) / (|ψ|² + ε)

and T ~ κ_tor |∇q_H| as a scalar weight, then

    G_ij = ∫ Φ_i* (T) Φ_j d²r   (i ≠ j),  G_ii = 0.

κ_tor is *not* chosen to recover g0(J−I). A single overall scale may be
reported after the fact; the pattern (dark multiplicity, λ_min < 0) is
the acceptance test.

Status: computed 2D proxy, not a 3D Cartan–Weitzenböck quadrature.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc_l1.gij_ring import compare_to_democratic, _g0_eff_from_G
from qvc_l1.locks import RETIRED_KAPPA_TOR


def phase_jacobian_density(psi: np.ndarray, dx: float) -> np.ndarray:
    """2D topological density (1/2π) n · (∂x n × ∂y n) proxy from ψ."""
    px = np.gradient(psi, dx, axis=1)
    py = np.gradient(psi, dx, axis=0)
    dens = np.imag(np.conj(px) * py - np.conj(py) * px)
    n = np.abs(psi) ** 2 + 1e-16
    return dens / (2.0 * np.pi * n)


def torsion_weight(q_h: np.ndarray, dx: float, kappa_tor: float = 1.0) -> np.ndarray:
    gx = np.gradient(q_h, dx, axis=1)
    gy = np.gradient(q_h, dx, axis=0)
    return kappa_tor * np.sqrt(gx * gx + gy * gy)


def torsion_G(
    Phi: np.ndarray,
    weight: np.ndarray,
    dA: float,
) -> np.ndarray:
    n = Phi.shape[0]
    G = np.zeros((n, n), dtype=float)
    w = weight.astype(float)
    for i in range(n):
        for j in range(n):
            G[i, j] = float(np.real(np.sum(np.conj(Phi[i]) * w * Phi[j])) * dA)
    G = 0.5 * (G + G.T)
    np.fill_diagonal(G, 0.0)
    return G


def run_htor1(
    psi0: np.ndarray,
    Phi: np.ndarray,
    dx: float,
    *,
    N: int = 5,
    kappa_tor: float = 1.0,
) -> dict[str, Any]:
    """Overlap torsion proxy on the first N amplitude modes."""
    n = min(N, Phi.shape[0])
    PhiN = Phi[:n]
    qh = phase_jacobian_density(psi0, dx)
    T = torsion_weight(qh, dx, kappa_tor=kappa_tor)
    dA = dx * dx
    G = torsion_G(PhiN, T, dA)
    g0_eff = _g0_eff_from_G(G)
    cmp = compare_to_democratic(G, g0_eff)
    eigs = np.asarray(cmp["eigenvalues"], dtype=float)
    dark_mult = int(np.sum(np.isclose(eigs, eigs[0], atol=1e-8 * max(1.0, abs(eigs[0])), rtol=0)))
    accepted = False
    return {
        "status": "ansatz_future_development",
        "hypothesis": "H-TOR1",
        "kappa_tor_used": kappa_tor,
        "retired_kappa_tor_0p25_not_fitted": abs(kappa_tor - RETIRED_KAPPA_TOR) > 1e-12
        or kappa_tor == 1.0,
        "qH_integral": float(np.sum(qh) * dA),
        "G": G.tolist(),
        "g0_eff": g0_eff,
        "democratic_comparison": cmp,
        "lambda_min": float(eigs[0]),
        "lambda_max": float(eigs[-1]),
        "dark_mult_numeric": dark_mult,
        "accepted": accepted,
        "parked": True,
        "not_hand_inserted_J_minus_I": True,
        "note": (
            "H-TOR1 is not accepted. The 2D phase-Jacobian proxy is parked as "
            "an ansatz for a future 3D Weitzenböck quadrature; it is not used "
            "to set G_ij or g0★. Democratic Spec(G) remains locked algebra."
        ),
    }
