"""Retardation factor from causal phonon propagator — Tier C algebra lock.

CLOSED-FORM LOCK
----------------
η_ret = 2 Δ₀ ξ / (ħ c_s)

where ξ is the relevant coherence / healing length.
This is an algebraic identity once Δ₀, ξ, ħ c_s are chosen;
it does not by itself fix those parameters.
"""

from __future__ import annotations


def eta_ret(Delta0: float, xi: float, hbar_cs: float) -> float:
    """η_ret = 2 Δ₀ ξ / (ħ c_s).

    Parameters
    ----------
    Delta0 : gap scale (meV)
    xi : length scale (µm)
    hbar_cs : ħ c_s (meV·µm)
    """
    if hbar_cs == 0:
        raise ValueError("hbar_cs must be nonzero")
    return 2.0 * Delta0 * xi / hbar_cs


def eta_ret_from_healing(Delta0: float, hbar_cs: float) -> float:
    """With ξ = a_H = ħ c_s / Δ₀ one obtains η_ret = 2 (identity)."""
    return eta_ret(Delta0, hbar_cs / Delta0, hbar_cs)
