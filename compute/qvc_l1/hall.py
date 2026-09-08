"""Hall convention from CS level k=2N and Q_H=1/N. Algebra only.

Locked dictionary (do not reset k or Q_H to match a textbook evaluation):
    k = 2N,   Q_H = 1/N,   k Q_H = 2   (product of two locks, not a pairing theorem)
    σ_xy^{uncat} = (k/2) e²/h = N e²/h
    δσ_xy = e²/(k h) = Q_H e²/(2h) = e²/(2N h)

This (k/2) e²/h evaluation is the QVC dictionary. The usual reading of
(k/4π)∫ a∧da is σ = k e²/h, which would be 2N e²/h at the same k and is
not used to change the locks.

Mean-field a_θ = Q_H/ρ has holonomy 2π Q_H = 2π/N. Textbook flux attachment
is 2π/k = π/N. Both statements are recorded; only the mean-field holonomy
is what the covariant GPE integrates.

Exchange angle in the CS-level convention is θ = π/k = π/(2N).
The linking-form θ=π/N is a different formula and is not mixed into δσ_xy.
Chern–Simons is not inserted into the default GPE EOM.
"""

from __future__ import annotations

import math
from typing import Any

from qvc_l1.locks import L1Locks

# 2019 SI exact
E_CHARGE_C = 1.602176634e-19
H_PLANCK_JS = 6.62607015e-34
E2_OVER_H_S = (E_CHARGE_C * E_CHARGE_C) / H_PLANCK_JS


def hall_convention(locks: L1Locks) -> dict[str, Any]:
    """Unique CS-level Hall staircase. Status: derived (algebra)."""
    n = int(locks.N_layers)
    k = 2 * n
    Q_H = 1.0 / float(n)
    pairing = k * Q_H
    sigma_parent = 0.5 * k * E2_OVER_H_S  # (k/2) e²/h = N e²/h
    d_sigma = E2_OVER_H_S / k  # e²/(k h)
    theta = math.pi / k
    steps = list(range(0, n + 1))
    sigma_xy = [i * d_sigma for i in steps]
    return {
        "N": n,
        "k_CS": k,
        "Q_H": Q_H,
        "pairing_k_QH": pairing,
        "pairing_identity": abs(pairing - 2.0) < 1e-15,
        "formula_parent": "σ_xy = (k/2) e²/h = N e²/h  (dictionary; not textbook k e²/h)",
        "formula_step": "δσ_xy = e²/(k h) = Q_H e²/(2h) = e²/(2N h)",
        "formula_status": "derived",
        "textbook_sigma_would_be": "k e²/h = 2N e²/h at the same k (not used)",
        "mean_field_holonomy": "2π Q_H = 2π/N",
        "textbook_attachment": "2π/k = π/N (not used to reset k)",
        "pairing_is_lock_product": True,
        "e2_over_h_S": E2_OVER_H_S,
        "sigma_xy_parent_S": sigma_parent,
        "delta_sigma_xy_S": d_sigma,
        "theta_CS_rad": theta,
        "theta_CS_formula": "θ = π/k = π/(2N)",
        "linking_theta_pi_over_N_not_mixed": math.pi / n,
        "staircase": {
            "filling_index": steps,
            "sigma_xy_S": sigma_xy,
        },
        "not_in_GPE_EOM": True,
        "note": (
            "Dictionary convention: k=2N, Q_H=1/N, σ_uncat=N e²/h, δσ=e²/(2N h). "
            "Not the unique Wen/Laughlin evaluation of (k/4π)∫a∧da. "
            "Algebra only; no CS current is discretised on the default GPE grid."
        ),
    }
