"""Live gap estimators that *must* use the GPE field ψ.

The withdrawn solver prescribed Δ(t) from a three-stage schedule
(0.600 → 0.389 → 0.200 → 0.032) and ignored ψ. That function is not
reimplemented. Time is not an argument of the estimator.
"""

from __future__ import annotations

import numpy as np


def ring_weight(X: np.ndarray, Y: np.ndarray, R_um: float, a_H_um: float) -> np.ndarray:
    """Gaussian weight on the torus ring (major R, core width a_H)."""
    rho = np.sqrt(X * X + Y * Y)
    d = np.abs(rho - R_um)
    return np.exp(-d * d / (2.0 * a_H_um * a_H_um))


def mean_density_weighted(psi: np.ndarray, weight: np.ndarray) -> float:
    """Ring-weighted mean of |ψ|². Uses psi."""
    n = np.abs(psi) ** 2
    wsum = float(np.sum(weight))
    if wsum <= 0.0:
        raise ValueError("weight must have positive mass")
    return float(np.sum(n * weight) / wsum)


def delta_live_from_psi(
    psi: np.ndarray,
    *,
    weight: np.ndarray,
    n_ref: float,
    Delta0_meV: float,
) -> float:
    """Primary live-gap estimator: Δ_live = Δ₀ √(n_w / n_ref).

    n_w is the ring-weighted mean of |ψ|². ``n_ref`` is the same
    average on the initial field (so Δ_live(0) = Δ₀ if n_w(0)=n_ref).

    This function uses ``psi``. It does not accept a time argument and
    cannot reproduce the withdrawn 0.389 / 0.032 schedule.
    """
    if n_ref <= 0.0:
        raise ValueError("n_ref must be positive")
    n_w = mean_density_weighted(psi, weight)
    return float(Delta0_meV * np.sqrt(n_w / n_ref))


def mu_nl_from_psi(
    psi: np.ndarray,
    *,
    weight: np.ndarray,
    alpha: float,
    beta: float,
) -> float:
    """Density-weighted cubic-quintic potential ⟨α|ψ|² + β|ψ|⁴⟩_n.

    Secondary chemical-potential proxy. Uses psi. Not a gap schedule.
    """
    n = np.abs(psi) ** 2
    v = alpha * n + beta * n * n
    mass = float(np.sum(n * weight))
    if mass <= 0.0:
        return float("nan")
    return float(np.sum(v * n * weight) / mass)


def coherence_g1(psi: np.ndarray) -> float:
    """Global first-order coherence g^{(1)} = |⟨ψ⟩|² / ⟨|ψ|²⟩."""
    n_mean = float(np.mean(np.abs(psi) ** 2))
    if n_mean <= 0.0:
        return 0.0
    return float(min(np.abs(np.mean(psi)) ** 2 / n_mean, 1.0))
