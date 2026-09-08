"""Acoustic Hawking R_WS matching and Berry R_mean (ratio, not a radius).

R_WS = ħ c_s / (π k_B T*) is a matching condition, not a duality proof.
The old 0.91 µm is the same formula at withdrawn T*=285 mK.
Berry R_mean≈0.908 is a GOE/GUE ratio, not a radius.
"""

from __future__ import annotations

from typing import Any

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.acoustic_matching import R_WS_from_T_star  # noqa: E402

from qvc_l1.locks import (
    BERRY_R_MEAN_RATIO,
    R_WS_WITHDRAWN_UM,
    TSTAR_WITHDRAWN_K,
    L1Locks,
)


def acoustic_matching(locks: L1Locks) -> dict[str, Any]:
    """Invert T_H = T* for R_WS at the locked Schwinger T*."""
    R_ws = float(
        R_WS_from_T_star(
            locks.hbar_cs_meV_um, locks.T_star_K, kB=locks.kB_meV_per_K
        )
    )
    R_ws_withdrawn = float(
        R_WS_from_T_star(
            locks.hbar_cs_meV_um, TSTAR_WITHDRAWN_K, kB=locks.kB_meV_per_K
        )
    )
    return {
        "formula": "R_WS = ħ c_s / (π k_B T*)",
        "status": "matched",
        "not_a_duality": True,
        "T_star_K": locks.T_star_K,
        "T_star_mK": locks.T_star_mK,
        "R_WS_um": R_ws,
        "R_WS_withdrawn_285mK_um": R_ws_withdrawn,
        "withdrawn_label_um": R_WS_WITHDRAWN_UM,
        "withdrawn_T_mK": 1e3 * TSTAR_WITHDRAWN_K,
        "Berry_R_mean_is_ratio": True,
        "Berry_R_mean": BERRY_R_MEAN_RATIO,
        "Berry_R_mean_status": "computed",
        "Berry_R_mean_note": (
            "Berry R_mean≈0.908 is ⟨|λ_min|⟩_GUE / ⟨|λ_min|⟩_GOE, not a radius. "
            "The numerical proximity to withdrawn R_WS(285 mK)≈0.91 µm is coincidental."
        ),
        "notes": (
            "Analogue Hawking matching, not a holographic dictionary. "
            "Do not identify R_WS with Berry R_mean."
        ),
    }
