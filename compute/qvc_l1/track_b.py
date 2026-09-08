"""Track B: continuum/Hubbard *matching* of GL coefficients.

A Schrieffer–Wolff projection of the Bistritzer–MacDonald model is not
performed. Track B records the unique dimensional identifications from
the locked GL/GPE scales (U, t, a_H, λ_M) and the ratios they imply.

Healing identification (Track A, Option A*):

    κ_GL = Δ₀ a_H² = ħ² / (2 m*)     ⇒    t_heal ≡ κ_GL / a_H² = Δ₀ .

If the tight-binding lattice constant were the moiré period instead,

    t_λ ≡ κ_GL / λ_M² = Δ₀ (a_H / λ_M)² .

On-site matching U ↔ Δ₀ is a gap-scale identification, not the Coulomb
U of MATBG (tens of meV). Bandwidth proxy W = ħ c_s / λ_M uses the
locked phonon speed as a moiré velocity stand-in. Holographic γ_QM is
not imported.

Status: matched / ansatz. Does not replace Track A.
"""

from __future__ import annotations

import math
from typing import Any

from qvc_l1.gl_coefficients import kappa_GL_mev_um2
from qvc_l1.locks import HOLOGRAPHIC_GAMMA_QM_MEV_NM2, L1Locks, build_l1_locks

LAMBDA_M_UM = 0.013  # MATBG magic-angle moiré period ~ 13 nm


def run_track_b(locks: L1Locks | None = None) -> dict[str, Any]:
    locks = locks or build_l1_locks()
    kappa = kappa_GL_mev_um2(locks)
    a_H = locks.a_H_um
    lam_M = LAMBDA_M_UM
    t_heal = kappa / (a_H**2)  # = Δ0 by construction
    t_lambda = kappa / (lam_M**2)
    U_gap = locks.Delta0_meV
    W = locks.hbar_cs_meV_um / lam_M
    # m* from κ_GL = ħ²/(2m*). ħ in meV·ps, convert length µm.
    # κ [meV·µm²] = ħ²/(2m*) ⇒ m* = ħ² / (2 κ).
    hbar = locks.hbar_meV_ps
    m_star_meV_ps2_um2 = (hbar * hbar) / (2.0 * kappa)
    xi_over_lambda = a_H / lam_M
    n0_heal = 1.0 / (math.pi * a_H * a_H)
    a_cell = 0.5 * math.sqrt(3.0) * lam_M * lam_M
    n0_moire = 1.0 / a_cell
    return {
        "status": "matched_not_hubbard",
        "U_gap_meV": U_gap,
        "t_heal_meV": t_heal,
        "t_over_U": t_heal / U_gap,
        "t_lambda_M_meV": t_lambda,
        "t_lambda_over_U": t_lambda / U_gap,
        "lambda_M_um": lam_M,
        "a_H_um": a_H,
        "a_H_over_lambda_M": xi_over_lambda,
        "W_moire_proxy_meV": W,
        "U_over_W_gap_scale": U_gap / W,
        "m_star_meV_ps2_um2": m_star_meV_ps2_um2,
        "n0_healing_um2": n0_heal,
        "n0_moire_um2": n0_moire,
        "alpha_GL_meV": -locks.Delta0_meV,
        "kappa_GL_meV_um2": kappa,
        "formulae": {
            "kappa_GL": "Δ₀ a_H²",
            "t_heal": "κ_GL / a_H² = Δ₀",
            "t_lambda": "κ_GL / λ_M² = Δ₀ (a_H/λ_M)²",
            "m_star": "ħ² / (2 κ_GL)",
            "W": "ħ c_s / λ_M",
            "U": "Δ₀  (gap-scale matching, not Coulomb U)",
        },
        "holographic_gamma_QM_not_used": HOLOGRAPHIC_GAMMA_QM_MEV_NM2,
        "not_BM_Schrieffer_Wolff": True,
        "note": (
            "t_heal/U = 1 is the Option A* healing freeze, not a Hubbard "
            "derivation. t_λ/U = (a_H/λ_M)² is the alternate lattice matching. "
            "U/W with U=Δ₀ is O(0.1) and is not the MATBG Coulomb U/W→∞. "
            "Track A remains the L1 coefficient lock."
        ),
    }
