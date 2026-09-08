"""Democratic vacuum coupling matrix G and dark (catalyzon) subspace.

CLOSED-FORM LOCK (Tier C)
-------------------------
For G = g₀ (J_N - I_N):
    λ_max = (N-1) g₀   (multiplicity 1)   — bright / in-phase
    λ_min = -g₀        (multiplicity N-1) — dark / catalyzon

Preprint claims λ_min = -4 g₀ (N=5) and supplemental -5 g₀ are FALSE
for the democratic matrix.

Catalyzon gap shift: Δ_cat = Δ_ex - g₀
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np


def democratic_coupling_matrix(N: int, g0: float) -> np.ndarray:
    """G = g0 (J_N - I_N). Off-diagonal vacuum coupling, no self-coupling."""
    if N < 2:
        raise ValueError("N must be >= 2")
    return g0 * (np.ones((N, N)) - np.eye(N))


def democratic_spectrum(N: int, g0: float) -> dict[str, Any]:
    """Exact Spec(G) for the democratic matrix."""
    return {
        "lambda_max": (N - 1) * g0,
        "lambda_max_mult": 1,
        "lambda_min": -g0,
        "lambda_min_mult": N - 1,
        "catalyzon_gap_shift": -g0,
    }


def catalyzon_dark_basis(N: int) -> np.ndarray:
    """Helmert orthonormal basis of the sum-zero (dark) subspace."""
    basis = []
    for k in range(1, N):
        v = np.zeros(N)
        v[:k] = 1.0
        v[k] = -float(k)
        v /= math.sqrt(k * (k + 1))
        basis.append(v)
    return np.asarray(basis)


def catalyzon_gap_shift(Delta_ex: float, g0: float) -> float:
    """Δ_cat = Δ_ex - g0  (from λ_min = -g0 on the dark subspace)."""
    return Delta_ex - g0


def verify_spectrum_numerically(N: int, g0: float, tol: float = 1e-10) -> bool:
    """Return True if numpy eigenspectrum matches closed-form Spec(G)."""
    G = democratic_coupling_matrix(N, g0)
    eigs = np.sort(np.linalg.eigvalsh(G))
    spec = democratic_spectrum(N, g0)
    dark = eigs[:-1]
    bright = eigs[-1]
    return (
        abs(bright - spec["lambda_max"]) <= tol
        and np.allclose(dark, spec["lambda_min"], atol=tol)
        and abs(float(np.trace(G))) <= tol
    )
