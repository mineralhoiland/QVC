"""Acoustic Hawking R_WS matching helper.

LABEL: MATCHING CONDITION — NOT A DUALITY PROOF
----------------------------------------------
T_H = ħ c_s / (π k_B R_WS)  is an analogue-gravity identity.
Inverting for R_WS given an observed T* is a *matching* step,
not a derivation of QED/gravity duality.
"""

from __future__ import annotations

import math

from qvc.params import PARAMS


def acoustic_hawking_T_K(
    hbar_cs_meV_um: float, R_WS_um: float, kB: float | None = None
) -> float:
    """T_H [K] from κ = 2 c_s / R_WS and T = ħκ/(2π k_B)."""
    kB = PARAMS.kB_meV_per_K if kB is None else kB
    return hbar_cs_meV_um / (math.pi * kB * R_WS_um)


def acoustic_hawking_T_mK(
    hbar_cs_meV_um: float, R_WS_um: float, kB: float | None = None
) -> float:
    """T_H in mK."""
    return 1e3 * acoustic_hawking_T_K(hbar_cs_meV_um, R_WS_um, kB)


def R_WS_from_T_star(
    hbar_cs_meV_um: float, T_star_K: float, kB: float | None = None
) -> float:
    """Invert T_H = T* for R_WS [µm]. Matching condition, not a derivation."""
    kB = PARAMS.kB_meV_per_K if kB is None else kB
    return hbar_cs_meV_um / (math.pi * kB * T_star_K)
