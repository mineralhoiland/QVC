"""Independent-mode quadrature enhancement of zero-point amplitude.

CLOSED-FORM LOCK
----------------
For n independent Gaussian ZPF modes:
    A_eff = √n · A_single

With running healing a_H(Δ) = ħ c_s / Δ and geometric F:
    A_eff(Δ) = √n_modes / (4π²) · (1/a_H) · F(2π a_H / R)
"""

from __future__ import annotations

import math

from qvc.vacuum_freeze import bessel_mode_sum


def A_eff_mode_enhancement(A_single: float, n_modes: int) -> float:
    """A_eff = √n_modes · A_single."""
    if n_modes < 1:
        raise ValueError("n_modes must be >= 1")
    return math.sqrt(n_modes) * A_single


def A_zpf_eff(Delta_meV: float, hbar_cs: float, R_um: float, n_modes: int) -> float:
    """Effective ZPF amplitude with mode quadrature and running healing length.

    Units: µm^{-1} if ħ c_s in meV·µm and Δ in meV.
    """
    if Delta_meV <= 0:
        return float("nan")
    a_H = hbar_cs / Delta_meV
    x = 2.0 * math.pi * a_H / R_um
    F = bessel_mode_sum(x)
    A_single = (1.0 / (4.0 * math.pi**2)) * (1.0 / a_H) * F
    return A_eff_mode_enhancement(A_single, n_modes)


def mode_enhancement_ratio(n_modes: int) -> float:
    """√n_modes — exact quadrature factor for independent modes."""
    if n_modes < 1:
        raise ValueError("n_modes must be >= 1")
    return math.sqrt(n_modes)
