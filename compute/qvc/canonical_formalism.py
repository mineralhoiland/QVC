"""
QVC Canonical Formalism — closed-form locks (Tier C) facade.

Single import surface for locked equations. Implementation lives in
dedicated modules; this file re-exports and documents the lock table.

Physical platform: N-layer TEVC / MATBG-like moiré excitonic condensate.
"""

from __future__ import annotations

from qvc.acoustic_matching import (
    R_WS_from_T_star,
    acoustic_hawking_T_K,
    acoustic_hawking_T_mK,
)
from qvc.coupling_matrix import (
    catalyzon_dark_basis,
    catalyzon_gap_shift,
    democratic_coupling_matrix,
    democratic_spectrum,
)
from qvc.gap_fold import (
    ell_from_delta_c,
    fold_condition,
    fold_square_root_scaling,
    gap_rhs,
    tree_level_ell,
)
from qvc.lens_chern_simons import (
    anyonic_exchange_angle,
    chern_simons_flat_u1,
    fractional_hall_charge,
    h2_vanishes_note,
    linking_form_self,
)
from qvc.mode_enhancement import (
    A_eff_mode_enhancement,
    A_zpf_eff,
    mode_enhancement_ratio,
)
from qvc.params import PARAMS, QVCParams
from qvc.retardation import eta_ret, eta_ret_from_healing
from qvc.vacuum_freeze import (
    bessel_mode_sum,
    casimir_preprint_formula,
    jacobi_theta3_alternating,
    vacuum_self_energy_corrected,
)

# ---------------------------------------------------------------------------
# Lock registry — provenance for each closed-form identity
# ---------------------------------------------------------------------------

CLOSED_FORM_LOCKS: dict[str, dict[str, str]] = {
    "democratic_spectrum": {
        "status": "PROVEN",
        "statement": "λ_max=(N-1)g0, λ_min=-g0 (mult N-1) for G=g0(J-I)",
        "module": "coupling_matrix",
    },
    "helmert_dark_basis": {
        "status": "PROVEN",
        "statement": "Helmert rows span sum-zero subspace; orthonormal",
        "module": "coupling_matrix",
    },
    "catalyzon_gap_shift": {
        "status": "PROVEN",
        "statement": "Δ_cat = Δ_ex - g0",
        "module": "coupling_matrix",
    },
    "fold_bifurcation": {
        "status": "PROVEN",
        "statement": "δ=1+α(lnδ+ℓ); δ_c(1-lnδ_c-ℓ)=1; α_c=δ_c; √ scaling",
        "module": "gap_fold",
    },
    "tree_level_ell": {
        "status": "PROVEN",
        "statement": "ℓ_tree = ln(1/2)+γ_E for log arg=1/2",
        "module": "gap_fold",
    },
    "chern_simons_linking": {
        "status": "PROVEN",
        "statement": "CS=k²/N, lk=1/N, Q_H=1/N; H2=0 invalidates 2-cycle",
        "module": "lens_chern_simons",
    },
    "anyonic_theta": {
        "status": "PROGRAM_PREDICTION",
        "statement": "θ=π/N from linking form (abelian anyon program)",
        "module": "lens_chern_simons",
    },
    "vacuum_form": {
        "status": "PROVEN_FORM_GEOMETRY_OPEN",
        "statement": "E_vac=(1/4π²)(ħc_s/r)F(2πr/R); freeze geometry before meV",
        "module": "vacuum_freeze",
    },
    "mode_quadrature": {
        "status": "PROVEN",
        "statement": "A_eff = √n A_single for independent Gaussian ZPF",
        "module": "mode_enhancement",
    },
    "retardation": {
        "status": "PROVEN_ALGEBRA",
        "statement": "η_ret = 2 Δ0 ξ / (ħ c_s)",
        "module": "retardation",
    },
    "acoustic_R_WS": {
        "status": "MATCHING_NOT_PROOF",
        "statement": "R_WS from T* via analogue T_H; not a duality derivation",
        "module": "acoustic_matching",
    },
}

TORSION_COUPLING_PROGRAM = r"""
Proposed (not yet numerically closed) first-principles kernel:

  S = S_FN[n] + S_sf[ρ_s, θ] + S_int[n, ρ_s, θ]

  J^μ_top = (1/8π²) ε^{μνρσ} F_{νρ} A_σ
  q_H     = ∂_μ J^μ_top = (1/4π²) F ∧ A

  T^μ_{νρ} = (κ_tor / 4π²) ε^μ_{νρσ} ∂^σ q_H

  G_ij = ∫_{WS} T^μ_{νρ} Φ_i^{*ν} Φ_j^ρ d³x

Acceptance: Hermitian G, G_ii=0, Tr G=0, λ_min<0 WITHOUT inserting
democratic form by hand. Compare Spec(G) to democratic exact result.
"""

__all__ = [
    "PARAMS",
    "QVCParams",
    "CLOSED_FORM_LOCKS",
    "TORSION_COUPLING_PROGRAM",
    "democratic_coupling_matrix",
    "democratic_spectrum",
    "catalyzon_dark_basis",
    "catalyzon_gap_shift",
    "gap_rhs",
    "fold_condition",
    "fold_square_root_scaling",
    "tree_level_ell",
    "ell_from_delta_c",
    "chern_simons_flat_u1",
    "linking_form_self",
    "fractional_hall_charge",
    "anyonic_exchange_angle",
    "h2_vanishes_note",
    "bessel_mode_sum",
    "casimir_preprint_formula",
    "vacuum_self_energy_corrected",
    "jacobi_theta3_alternating",
    "A_zpf_eff",
    "A_eff_mode_enhancement",
    "mode_enhancement_ratio",
    "eta_ret",
    "eta_ret_from_healing",
    "acoustic_hawking_T_K",
    "acoustic_hawking_T_mK",
    "R_WS_from_T_star",
]
