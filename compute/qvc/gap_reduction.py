"""Recompute gap-reduction stages and derive participation λ.

Channels (do not stack as one percentage)
----------------------------------------
1. Static fold:  Δ*(λ) from δ=1+α(lnδ+ℓ),  α=λ E_vac/Δ₀.
   Physical branch exists only for λ ≤ λ_fold = α_c Δ₀ / E_vac.
   Maximum static reduction is 1−δ_c ≈ 81.8%, not 92%.

2. Catalyzon Spec(G): Δ_cat = Δ_ex − g₀,  g₀ = λ E_vac
   under the identification that the same participation feeds off-diagonal
   hybridisation. Separate from the fold loop.

3. Coherent TQGL 92%: preprint Δ(1 ps)=0.032 meV used V_ZPF=−E_Cas with
   E_Cas=0.244 meV and started from Δ_QVC=0.389. There is no split-step
   solver in-repo. Recompute drive scale E_vac/0.244 and prove that
   δ=0.032/0.6 < δ_c cannot be a static physical branch.

Original folklore (withdrawn baselines)
--------------------------------------
  35%  Δ: 0.600 → 0.389   (static loop at 0.244 meV)
  52%  Δ: 0.600 → 0.29    (catalyzon; mixed Δ_QVC−0.1)
  92%  Δ: 0.600 → 0.032   (coherent TQGL 1 ps; also ≈92% of 0.389)
"""

from __future__ import annotations

import math
from typing import Any

from qvc.bridge import (
    alpha_for_target_delta,
    alpha_from_E_vac,
    derive_bridge,
    ell_frozen_cavity,
)
from qvc.gap_fold import fold_condition, solve_fold_branches
from qvc.params import PARAMS
from qvc.self_consistency import catalyzon_dressed_gap


# Preprint folklore (not locks)
E_CAS_CLAIM_MEV = 0.244
DELTA_QVC_CLAIM_MEV = 0.389
DELTA_1PS_CLAIM_MEV = 0.032
G0_PREPRINT_MEV = 0.1  # |λ_min| estimate from E_Cas=0.244
G0_SUITE_MEV = 0.21  # working Spec(G) example used in freeze suite


def _pct(delta_from: float, delta_to: float) -> float:
    return 100.0 * (1.0 - delta_to / delta_from)


def lambda_fold_critical(E_vac: float, alpha_c: float, Delta0: float) -> float:
    """λ_fold = α_c Δ₀ / E_vac — unit participation is fold-out iff λ_fold < 1."""
    if E_vac <= 0:
        raise ValueError("E_vac must be positive")
    return alpha_c * Delta0 / E_vac


def lambda_from_g0(g0: float, E_vac: float) -> float:
    """λ_cat = g₀ / E_vac if g₀ = λ E_vac (same participation)."""
    if E_vac <= 0:
        raise ValueError("E_vac must be positive")
    return g0 / E_vac


def static_at_lambda(lam: float, E_vac: float, ell: float, Delta0: float) -> dict[str, Any]:
    """Evaluate fold branches at a given participation."""
    alpha = alpha_from_E_vac(E_vac, Delta0=Delta0, lam=lam)
    br = solve_fold_branches(alpha, ell)
    phys = br["phys"]
    return {
        "lambda": lam,
        "alpha": alpha,
        "alpha_c": br["alpha_c"],
        "alpha_over_ac": alpha / br["alpha_c"] if br["alpha_c"] else float("nan"),
        "folded_out": phys is None,
        "delta_star": phys,
        "Delta_star_meV": None if phys is None else phys * Delta0,
        "gap_reduction_pct": None if phys is None else 100.0 * (1.0 - phys),
        "delta_unst": br["unst"],
        "Delta_unst_meV": None if br["unst"] is None else br["unst"] * Delta0,
    }


def recompute_gap_reduction() -> dict[str, Any]:
    """Full recomputation vs withdrawn 35/52/92 and derived λ."""
    p = PARAMS
    Delta0 = p.Delta0_meV
    b = derive_bridge(lam=1.0)
    E_vac = b["E_vac_meV"]
    ell = b["ell_frozen"]
    alpha_c = b["alpha_c"]
    delta_c = b["delta_c"]

    lam_fold = lambda_fold_critical(E_vac, alpha_c, Delta0)
    lam_sub = 0.90 * lam_fold  # subcritical operating point
    lam_fold_in = 0.999 * lam_fold  # just inside, so the branch solver returns phys

    # --- Folklore inventory ---
    folklore = {
        "E_Cas_claimed_meV": E_CAS_CLAIM_MEV,
        "stage1_static": {
            "Delta_from": Delta0,
            "Delta_to": DELTA_QVC_CLAIM_MEV,
            "reduction_pct": _pct(Delta0, DELTA_QVC_CLAIM_MEV),
            "baseline": "Delta0 -> Delta_QVC",
        },
        "stage2_catalyzon": {
            "Delta_from": Delta0,
            "Delta_to": DELTA_QVC_CLAIM_MEV - G0_PREPRINT_MEV,
            "reduction_pct": _pct(Delta0, DELTA_QVC_CLAIM_MEV - G0_PREPRINT_MEV),
            "note": "preprint mixed (0.389-0.1)/0.600 ≈ 52%",
        },
        "stage3_tqgl_from_Delta0": {
            "Delta_from": Delta0,
            "Delta_to": DELTA_1PS_CLAIM_MEV,
            "reduction_pct": _pct(Delta0, DELTA_1PS_CLAIM_MEV),
            "note": "quoted 92% from Delta0 but 0.032/0.6 is 94.7%",
        },
        "stage3_tqgl_from_Delta_QVC": {
            "Delta_from": DELTA_QVC_CLAIM_MEV,
            "Delta_to": DELTA_1PS_CLAIM_MEV,
            "reduction_pct": _pct(DELTA_QVC_CLAIM_MEV, DELTA_1PS_CLAIM_MEV),
            "note": "(0.389-0.032)/0.389 ≈ 91.8% — origin of the 92% label",
        },
    }

    # --- Fold channel ---
    fold_unit = static_at_lambda(1.0, E_vac, ell, Delta0)
    fold_crit = static_at_lambda(lam_fold_in, E_vac, ell, Delta0)
    fold_crit["lambda_exact_critical"] = lam_fold
    fold_crit["Delta_c_analytic_meV"] = delta_c * Delta0
    fold_sub = static_at_lambda(lam_sub, E_vac, ell, Delta0)
    # 35% folklore inverted
    delta_35 = DELTA_QVC_CLAIM_MEV / Delta0
    alpha_35 = alpha_for_target_delta(delta_35, ell)
    lam_35 = alpha_35 * Delta0 / E_vac
    fold_35 = static_at_lambda(lam_35, E_vac, ell, Delta0)

    # 92% as a *static* target is below δ_c — prove it
    delta_92 = DELTA_1PS_CLAIM_MEV / Delta0
    delta_92_label = 0.08  # exact 92% remaining
    below_fold = delta_92 < delta_c
    # Unstable-branch α that would have an unstable root near 0.032
    alpha_near_032 = alpha_for_target_delta(max(delta_92, 1e-6), ell)
    lam_unstable_032 = alpha_near_032 * Delta0 / E_vac
    br_032 = solve_fold_branches(alpha_near_032, ell)

    fold_channel = {
        "E_vac_meV": E_vac,
        "ell": ell,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "Delta_c_meV": delta_c * Delta0,
        "max_static_reduction_pct": 100.0 * (1.0 - delta_c),
        "lambda_fold": lam_fold,
        "lambda_subcritical_0p90": lam_sub,
        "lambda_for_withdrawn_35pct": lam_35,
        "unit_lambda": fold_unit,
        "at_fold": fold_crit,
        "at_0p90_fold": fold_sub,
        "at_35pct_lambda": fold_35,
        "delta_1ps_over_Delta0": delta_92,
        "delta_1ps_below_delta_c": below_fold,
        "alpha_for_delta_032": alpha_near_032,
        "lambda_for_unstable_032": lam_unstable_032,
        "branches_at_alpha_032": {
            "phys": br_032["phys"],
            "unst": br_032["unst"],
        },
        "theorem": (
            "Static physical branch cannot reach 92% reduction: min δ is δ_c≈0.182 "
            "(81.8%). Preprint 0.032 meV is below the fold (unstable or dynamical)."
        ),
    }

    # --- Catalyzon channel ---
    g0_scaled = G0_PREPRINT_MEV * (E_vac / E_CAS_CLAIM_MEV)
    lam_cat_preprint_g0 = lambda_from_g0(G0_PREPRINT_MEV, E_vac)
    lam_cat_scaled = lambda_from_g0(g0_scaled, E_vac)  # = 0.1/0.244
    lam_cat_suite = lambda_from_g0(G0_SUITE_MEV, E_vac)

    def cat_row(g0: float, lam: float) -> dict[str, Any]:
        d = catalyzon_dressed_gap(g0, Delta0)
        return {
            **d,
            "lambda": lam,
            "g0_over_Evac": g0 / E_vac,
            "exceeds_unit_participation": lam > 1.0,
            "exceeds_fold": lam > lam_fold,
        }

    catalyzon_channel = {
        "identification": "g0 = lambda * E_vac  (same participation as fold α)",
        "g0_scaled_from_0p1_and_0p244": g0_scaled,
        "at_preprint_g0_0p1": cat_row(G0_PREPRINT_MEV, lam_cat_preprint_g0),
        "at_scaled_g0": cat_row(g0_scaled, lam_cat_scaled),
        "at_suite_g0_0p21": cat_row(G0_SUITE_MEV, lam_cat_suite),
        "at_lambda_fold": cat_row(lam_fold * E_vac, lam_fold),
        "at_lambda_sub": cat_row(lam_sub * E_vac, lam_sub),
        "note": (
            "Catalyzon reduction = g0/Δ₀. At λ_fold this equals α_c≈0.182 (18%), "
            "not 52%. The 52% used Δ_QVC−0.1 over Δ_ex and is withdrawn as a mix."
        ),
    }

    # --- TQGL 92% drive rescaling (not a TDGL rerun) ---
    drive_ratio = E_vac / E_CAS_CLAIM_MEV
    # Phenomenological: remaining fraction ~ exp(−κ E/Δ₀) fit to 0.032 at 0.244
    rem_old = DELTA_1PS_CLAIM_MEV / Delta0
    kappa = -math.log(rem_old) / (E_CAS_CLAIM_MEV / Delta0)
    rem_new = math.exp(-kappa * (E_vac / Delta0))
    tqgl = {
        "status": "DRIVE_RESCALE_NOT_A_TDGL_RERUN",
        "V_ZPF_scale_new_over_old": drive_ratio,
        "no_split_step_solver_in_repo": True,
        "exponential_kappa": kappa,
        "Delta_1ps_phenomenological_meV": rem_new * Delta0,
        "reduction_pct_phenomenological": 100.0 * (1.0 - rem_new),
        "started_from_withdrawn_Delta_QVC": DELTA_QVC_CLAIM_MEV,
        "cannot_start_from_Delta_QVC": fold_unit["folded_out"],
        "note": (
            "92% used V_ZPF=−0.244 meV and Δ_QVC=0.389. With Tier-C E_vac the "
            "drive is 0.527× weaker; static theory has no 0.389 fixed point at λ=1."
        ),
    }

    # --- Combined λ recommendation ---
    # Shared participation: largest λ that stays on the physical fold branch.
    # Catalyzon then uses g0=λ E_vac at that same λ (not added on top of Δ*).
    lam_star = lam_sub
    rec_fold = static_at_lambda(lam_star, E_vac, ell, Delta0)
    rec_cat = catalyzon_dressed_gap(lam_star * E_vac, Delta0)
    recommendation = {
        "lambda_star": lam_star,
        "rule": "λ* = 0.90 λ_fold = 0.90 α_c Δ₀ / E_vac  (subcritical, same λ for both channels)",
        "fold": rec_fold,
        "catalyzon_parallel_not_stacked": rec_cat,
        "do_not_add_percentages": True,
        "lambda_fold": lam_fold,
        "lambda_cat_scaled_g0": lam_cat_scaled,
        "lambda_cat_preprint_0p1": lam_cat_preprint_g0,
        "why_not_92": (
            "92% is below δ_c; no λ on the physical branch produces it. "
            "It is dynamical/unstable, and was computed with withdrawn 0.244."
        ),
    }

    return {
        "Delta0_meV": Delta0,
        "E_vac_meV": E_vac,
        "E_Cas_withdrawn_meV": E_CAS_CLAIM_MEV,
        "folklore": folklore,
        "fold_channel": fold_channel,
        "catalyzon_channel": catalyzon_channel,
        "tqgl_92": tqgl,
        "recommendation": recommendation,
    }


def evaluate_partition(
    lam_alpha: float,
    lam_G: float = 0.0,
    *,
    E_vac: float | None = None,
    Delta0: float | None = None,
    ell: float | None = None,
) -> dict[str, Any]:
    """Split vacuum participation: λ_α → fold, λ_G → catalyzon, λ_tot = λ_α+λ_G.

    Unitarity: λ_tot ≤ 1 means the condensate cannot couple to more than the
    full geometric vacuum energy. Fold-in requires λ_α ≤ λ_fold.
    """
    b = derive_bridge(lam=1.0)
    E_vac = b["E_vac_meV"] if E_vac is None else E_vac
    Delta0 = b["Delta0_meV"] if Delta0 is None else Delta0
    ell = b["ell_frozen"] if ell is None else ell
    static = static_at_lambda(lam_alpha, E_vac, ell, Delta0)
    g0 = lam_G * E_vac
    cat = catalyzon_dressed_gap(g0, Delta0)
    return {
        "lambda_alpha": lam_alpha,
        "lambda_G": lam_G,
        "lambda_tot": lam_alpha + lam_G,
        "unitarity_ok": (lam_alpha + lam_G) <= 1.0 + 1e-12,
        "static": static,
        "catalyzon": cat,
        "bridge": {"E_vac_meV": E_vac, "Delta0_meV": Delta0, "ell": ell},
    }


def run_gap_reduction_suite() -> dict[str, Any]:
    """Partition scenarios used by tests and the λ derivation."""
    full = recompute_gap_reduction()
    b = derive_bridge(lam=1.0)
    E_vac = b["E_vac_meV"]
    lam_fold = full["fold_channel"]["lambda_fold"]
    leftover = max(0.0, 1.0 - 0.999 * lam_fold)
    scenarios = {
        "unit_alpha_foldout": evaluate_partition(1.0, lam_G=0.0),
        "equal_split": evaluate_partition(0.5, lam_G=0.5),
        "fold_saturated_leftover_G": evaluate_partition(0.999 * lam_fold, lam_G=leftover),
    }
    return {
        **full,
        "bridge": b,
        "lambda_alpha_max": lam_fold,
        "g0_working_exceeds_Evac": G0_SUITE_MEV > E_vac,
        "scenarios": scenarios,
    }
