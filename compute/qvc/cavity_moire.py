"""Moiré cavity factor F_moiré applied to locked A_geom.

Option A* static lock does *not* include F_moiré:
    α = λ E_vac / Δ₀,    A_geom = E_vac / (ħ c_s).

The cavity section keeps F_moiré as a labeled enhancement of the same
geometric amplitude:

    A_eff = F_moiré A_geom,
    α_cav = λ F_moiré E_vac / Δ₀,
    λ_fold^(F) = α_c Δ₀ / (F_moiré E_vac) = λ_fold / F_moiré.

F_moiré(N) is the self-cavity formula that evaluates to 1.119 at N=5,
not the stale A_ZPF = 0.018 / 0.0201 pair.
"""

from __future__ import annotations

import math
from typing import Any

from qvc.bridge import derive_bridge
from qvc.gap_fold import fold_condition
from qvc.gap_reduction import lambda_fold_critical, static_at_lambda
from qvc.params import PARAMS

ALPHA_FS = 1.0 / 137.0
U_OVER_W_CAVITY = 10.0


def F_moire(N: int, *, U_over_W: float = U_OVER_W_CAVITY, alpha_fs: float = ALPHA_FS) -> float:
    """Self-cavity enhancement. N=5 returns F ≈ 1.119."""
    if N < 2:
        return 1.0
    r_eff2 = (math.pi * alpha_fs) ** 2 * U_over_W
    Q = 10.0 + 2.0 * (N - 2)
    return math.sqrt(1.0 + (N - 2) * r_eff2 * Q)


def apply_F_moire_to_A_geom(*, N: int = 5, lambda_star_frac: float = 0.90) -> dict[str, Any]:
    """Re-check λ_fold after multiplying A_geom by F_moiré(N)."""
    p = PARAMS
    b = derive_bridge(lam=1.0)
    E_vac = float(b["E_vac_meV"])
    A_geom = float(b["A_geom_um_inv"])
    ell = float(b["ell_frozen"])
    Delta0 = p.Delta0_meV
    delta_c, alpha_c = fold_condition(ell)
    F = F_moire(N)
    A_eff = F * A_geom
    lam_fold = lambda_fold_critical(E_vac, alpha_c, Delta0)
    lam_fold_F = lam_fold / F
    lam_star = lambda_star_frac * lam_fold
    lam_star_F = lambda_star_frac * lam_fold_F
    at_star_unenhanced = static_at_lambda(lam_star, E_vac, ell, Delta0)
    # Cavity-enhanced map at the *unenhanced* working λ★
    at_star_with_F = static_at_lambda(lam_star * F, E_vac, ell, Delta0)
    at_star_F = static_at_lambda(lam_star_F, E_vac, ell, Delta0)
    withdrawn_018 = 0.018
    withdrawn_0201 = 0.0201
    return {
        "status": "labeled_cavity_enhancement",
        "not_the_option_A_star_static_lock": True,
        "N": N,
        "F_moire": F,
        "A_geom_um_inv": A_geom,
        "A_eff_um_inv": A_eff,
        "E_vac_meV": E_vac,
        "lambda_fold": lam_fold,
        "lambda_fold_F": lam_fold_F,
        "lambda_star": lam_star,
        "lambda_star_F": lam_star_F,
        "star_unenhanced_folded_out": bool(at_star_unenhanced["folded_out"]),
        "star_times_F_folded_out": bool(at_star_with_F["folded_out"]),
        "retuned_star_F_folded_out": bool(at_star_F["folded_out"]),
        "retuned_Delta_star_F_meV": at_star_F["Delta_star_meV"],
        "physical_branch_exists_for_lambda_le_lambda_fold_F": lam_fold_F > 0.0,
        "withdrawn_A_ZPF_0p018_not_used": abs(A_geom - withdrawn_018) > 0.5,
        "withdrawn_A_cav_0p0201_not_used": abs(A_eff - withdrawn_0201) > 0.5,
        "note": (
            "Option A* keeps α = λ E_vac/Δ₀. F_moiré multiplies A_geom in the "
            "cavity section only. If λ★_A* · F_moiré > λ_fold the enhanced "
            "static map is past the fold; retune to λ★_F = 0.90 λ_fold/F."
        ),
    }
