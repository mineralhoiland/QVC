"""Locked kinematics: m*, drive frequencies, two velocities, VK bound.

Package ħ c_s and graphene c_LA are distinct. Do not identify them.
Withdrawn E_Cas=0.244 is not used for the VK constraint.
"""

from __future__ import annotations

import math
from typing import Any

from qvc.tqgl.gpe import kinetic_coeff_mev_um2
from qvc.tqgl.locks import HBAR_MEV_PS, FrozenTQGL

from qvc_l1.locks import L1Locks

# Graphene longitudinal acoustic speed (phonon_folding.MoirePhononModel.c_LA).
# This locates q* in the moiré Brillouin zone. It is not ħ c_s.
C_LA_M_S = 21400.0

# ħ²/(2 m_e) in (meV, µm): 3.8099821 eV Å² → 3.8099821e-5 meV µm².
KAPPA_ELECTRON_MEV_UM2 = 3.8099821e-5

# 1 µm/ps = 1e6 m/s
UM_PER_PS_TO_M_PER_S = 1.0e6


def locked_kinematics(tqgl: FrozenTQGL, l1: L1Locks) -> dict[str, Any]:
    """Algebra of Option A* lengths and ħ c_s. Status: derived."""
    kappa = kinetic_coeff_mev_um2(tqgl)
    hbar = HBAR_MEV_PS
    m_over_me = KAPPA_ELECTRON_MEV_UM2 / kappa
    a_H = float(tqgl.r_um)
    hw_aH = float(tqgl.hbar_cs_meV_um) / a_H
    c_s_um_ps = float(tqgl.hbar_cs_meV_um) / hbar
    c_s_m_s = c_s_um_ps * UM_PER_PS_TO_M_PER_S
    lam_M_nm = 13.0
    K_M = 4.0 * math.pi / (3.0 * lam_M_nm)

    def resonance(delta_meV: float) -> dict[str, float]:
        omega_ps = 2.0 * delta_meV / hbar  # rad/ps
        nu_THz = omega_ps / (2.0 * math.pi)
        omega_rad_s = omega_ps * 1.0e12
        q_LA_nm = (omega_rad_s / C_LA_M_S) * 1.0e-9
        lam_LA_nm = 2.0 * math.pi / (omega_rad_s / C_LA_M_S) * 1.0e9
        q_cs_um = 2.0 * delta_meV / float(tqgl.hbar_cs_meV_um)
        lam_cs_nm = 1.0e3 * 2.0 * math.pi / q_cs_um
        return {
            "omega_ps_inv": omega_ps,
            "nu_THz": nu_THz,
            "q_LA_nm_inv": q_LA_nm,
            "lambda_LA_nm": lam_LA_nm,
            "q_over_K_M": q_LA_nm / K_M,
            "q_cs_um_inv": q_cs_um,
            "lambda_cs_nm": lam_cs_nm,
        }

    res0 = resonance(float(tqgl.Delta0_meV))
    res_star = resonance(float(tqgl.Delta_star_meV))
    n_m34 = float(tqgl.N_layers) ** (-0.75)
    C_vk = float(tqgl.E_vac_meV) / n_m34
    return {
        "kappa_GL_meV_um2": kappa,
        "m_star_over_m_e": m_over_me,
        "hbar_omega_ph_aH_meV": hw_aH,
        "c_s_m_s": c_s_m_s,
        "c_LA_m_s": C_LA_M_S,
        "K_M_nm_inv": K_M,
        "drive_Delta0": res0,
        "drive_Deltastar": res_star,
        "N_to_minus_3_4": n_m34,
        "C_VK_min_meV": C_vk,
        "E_vac_meV": float(tqgl.E_vac_meV),
        "not_E_cas_0244": True,
        "velocities_not_identified": abs(c_s_m_s - C_LA_M_S) > 1.0e4,
        "note": (
            "ħ c_s is the Casimir/gap-map/retardation/Goldstone lock. "
            "c_LA is graphene longitudinal acoustic and locates q* in the MBZ. "
            "VK uses E_vac, not withdrawn 0.244 meV."
        ),
    }
