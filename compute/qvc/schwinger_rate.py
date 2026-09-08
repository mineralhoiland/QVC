"""Locked Schwinger rate used in the paper box.

    Γ = ω₀ (A_eff / A_c) exp(−π A_c / A_eff),
    A_eff / A_c = α★ = λ★ E_vac / Δ₀,
    T★ = E_pair / (k_B π / α★),    E_pair = 2 g0★.

F_moiré is not inserted into this lock. It belongs to the cavity section.
"""

from __future__ import annotations

import math
from typing import Any

from qvc.tqgl.diagnostics import schwinger_twin_channel
from qvc.tqgl.locks import FrozenTQGL, build_frozen_tqgl


def boxed_schwinger_rate(lock: FrozenTQGL | None = None) -> dict[str, Any]:
    """Evaluate the boxed formula and cross-check T*."""
    lock = build_frozen_tqgl() if lock is None else lock
    sch = schwinger_twin_channel(lock)
    alpha_star = lock.lambda_star * lock.E_vac_meV / lock.Delta0_meV
    ratio = float(sch["A_eff_over_Ac_schwinger"])
    omega0 = float(sch["omega0_ps_inv"])
    Gamma = omega0 * ratio * math.exp(-math.pi / ratio)
    T_star = float(sch["T_star_K"])
    E_pair = float(sch["E_pair_meV"])
    T_from_alpha = E_pair / (lock.kB_meV_per_K * math.pi / alpha_star)
    Gamma_Hz = Gamma * 1e12
    return {
        "formula": "Gamma = omega0 * (A_eff/A_c) * exp(-pi A_c/A_eff)",
        "A_eff_over_Ac": ratio,
        "alpha_star": alpha_star,
        "ratio_is_alpha_star": abs(ratio - alpha_star) < 1e-12,
        "omega0_ps_inv": omega0,
        "schwinger_exponent": math.pi / ratio,
        "Gamma_sch_ps_inv": Gamma,
        "Gamma_sch_Hz": Gamma_Hz,
        "E_pair_meV": E_pair,
        "T_star_K": T_star,
        "T_star_from_alpha_K": T_from_alpha,
        "T_star_matches_alpha": abs(T_star - T_from_alpha) < 1e-9,
        "F_moire_not_in_exponent": True,
        "matches_diagnostics": abs(Gamma - float(sch["Gamma_sch_ps_inv"])) < 1e-18,
    }
