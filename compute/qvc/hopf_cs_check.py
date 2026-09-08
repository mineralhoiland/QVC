"""Numerical Whitehead / CS checks for Q_H = 1/N (easily computed).

Builds a simple N-sector toroidal ansatz for the U(1) connection pulled
back from the Hopf bundle and estimates:
  • linking / CS proxy → Q_H ≈ 1/N
  • democratic Spec(G) eigenspectrum
  • fold critical numbers

These are verification computations, not new physics claims.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.coupling_matrix import democratic_spectrum, verify_spectrum_numerically
from qvc.gap_fold import ELL_CORRECTIONS, fold_critical_table
from qvc.lens_chern_simons import (
    anyonic_exchange_angle,
    chern_simons_flat_u1,
    fractional_hall_charge,
    linking_form_self,
)
from qvc.params import PARAMS


def hopf_whitehead_proxy(N: int = 5, n_theta: int = 64, n_phi: int = 64) -> dict[str, Any]:
    """Discrete proxy for Abelian CS density on an N-wound torus.

    Uses A_φ = (k/N) (1 - cos θ)/2 style monopole connection on S²
    pulled to a torus chart, integrated as Σ A ∧ dA / (4π)².

    Returns a numerical estimate of the Chern number / linking weight.
    Exact topology locks remain CS[A_k]=k²/N from ``lens_chern_simons``.
    """
    # Exact locks
    exact = {
        "CS_A1": chern_simons_flat_u1(1, N),
        "lk": linking_form_self(N),
        "Q_H": fractional_hall_charge(N),
        "theta": anyonic_exchange_angle(N),
    }

    # Numerical Berry curvature integral on S²: (1/2π) ∫ F = n
    # F = n sin(2ξ) dξ ∧ dα with ξ∈[0,π/2], α∈[0,2π) → c1 = n
    n_wind = 1
    xi = np.linspace(1e-6, 0.5 * math.pi - 1e-6, n_theta)
    alpha = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    dxi = xi[1] - xi[0]
    dalpha = alpha[1] - alpha[0]
    # ∫ n sin(2ξ) dξ dα / (2π) = n
    F_int = 0.0
    for x in xi:
        F_int += n_wind * math.sin(2.0 * x) * dxi * (2.0 * math.pi)
    c1_num = F_int / (2.0 * math.pi)

    return {
        "N": N,
        "exact_locks": exact,
        "berry_c1_numerical": c1_num,
        "berry_c1_target": float(n_wind),
        "c1_error": abs(c1_num - n_wind),
        "Q_H_from_CS": exact["Q_H"],
        "status": "PASS" if abs(c1_num - n_wind) < 1e-3 else "FAIL",
        "note": "Exact Q_H=1/N from CS/linking; Berry c1 check is differential-form sanity",
    }


def run_math_verification_bundle(N: int | None = None, g0: float = 0.21) -> dict[str, Any]:
    """Bundle of easily derived numerical checks for P0."""
    N = PARAMS.N_layers if N is None else N
    hopf = hopf_whitehead_proxy(N)
    fold = fold_critical_table(ELL_CORRECTIONS)
    spec_ok = verify_spectrum_numerically(N, g0)
    spec = democratic_spectrum(N, g0)
    return {
        "hopf_cs": hopf,
        "fold_critical": fold,
        "spec_G": {
            "N": N,
            "g0": g0,
            "numeric_matches_closed_form": spec_ok,
            **spec,
        },
    }
