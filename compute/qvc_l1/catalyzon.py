"""Catalyzon g0 from the bridge: g0 = λ E_vac.

Working point λ★: g0★ = λ★ E_vac, Δ_cat = Δ₀ − g0.
At λ_fold: g0 = Δ_c ≈ α_c Δ₀ (~18.2% = α_c).

Parallel channel only — do not stack with fold or live-gap percentages.
Retire g0~0.21 as an insertion. Democratic Spec(G): λ_min = −g0 (mult. N−1).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.coupling_matrix import (  # noqa: E402
    democratic_coupling_matrix,
    democratic_spectrum,
    verify_spectrum_numerically,
)
from qvc.self_consistency import catalyzon_dressed_gap  # noqa: E402

from qvc_l1.locks import RETIRED_G0_MEV, L1Locks


def g0_from_lambda(locks: L1Locks, lam: float) -> float:
    """g0 = λ E_vac. Status: derived (bridge identification)."""
    return float(lam * locks.E_vac_meV)


def catalyzon_at_lambda(locks: L1Locks, lam: float, *, label: str) -> dict[str, Any]:
    g0 = g0_from_lambda(locks, lam)
    dressed = catalyzon_dressed_gap(g0, locks.Delta0_meV)
    spec = democratic_spectrum(locks.N_layers, g0)
    G = democratic_coupling_matrix(locks.N_layers, g0)
    eigs = np.sort(np.linalg.eigvalsh(G))
    return {
        "label": label,
        "lambda": float(lam),
        "g0_meV": g0,
        "Delta_cat_meV": float(dressed["Delta_cat_meV"]),
        "reduction_pct": float(dressed["reduction_pct"]),
        "formula": "g0 = λ E_vac,  Δ_cat = Δ₀ − g0,  λ_min = −g0",
        "status": "derived",
        "spec": {
            "lambda_max_meV": float(spec["lambda_max"]),
            "lambda_max_mult": spec["lambda_max_mult"],
            "lambda_min_meV": float(spec["lambda_min"]),
            "lambda_min_mult": spec["lambda_min_mult"],
            "eigs_numerical_meV": [float(x) for x in eigs],
            "numerical_matches_closed_form": bool(
                verify_spectrum_numerically(locks.N_layers, g0)
            ),
        },
        "parallel_channel_only": True,
        "do_not_stack_with_fold_or_live_gap": True,
        "retired_g0_0p21_not_used": True,
        "retired_g0_meV": RETIRED_G0_MEV,
    }


def catalyzon_card(locks: L1Locks) -> dict[str, Any]:
    """Working-point and fold-point catalyzon rows plus Spec(G) at g0★."""
    star = catalyzon_at_lambda(locks, locks.lambda_star, label="lambda_star")
    fold = catalyzon_at_lambda(locks, locks.lambda_fold, label="lambda_fold")
    return {
        "identification": "g0 = λ E_vac (same participation as the gap-map α)",
        "identification_status": "derived",
        "at_lambda_star": star,
        "at_lambda_fold": fold,
        "g0_star_meV": star["g0_meV"],
        "Delta_cat_star_meV": star["Delta_cat_meV"],
        "g0_fold_meV": fold["g0_meV"],
        "Delta_cat_fold_meV": fold["Delta_cat_meV"],
        "g0_fold_equals_Delta_c": abs(fold["g0_meV"] - locks.Delta_c_meV) < 1e-9,
        "reduction_fold_equals_alpha_c": abs(fold["reduction_pct"] / 100.0 - locks.alpha_c)
        < 1e-9,
        "lambda_min_equals_minus_g0": abs(
            star["spec"]["lambda_min_meV"] + star["g0_meV"]
        )
        < 1e-12,
        "N_layers": locks.N_layers,
        "note": (
            "At λ_fold, g0 = α_c Δ₀ = Δ_c so the catalyzon reduction is α_c "
            "(~18.2%), not 52%. g0~0.21 is a retired insertion and is not derived."
        ),
    }
