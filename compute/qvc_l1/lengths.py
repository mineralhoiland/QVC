"""Three lengths and the retardation factor η_ret = 2 Δ₀ ξ / (ħ c_s).

Conventions are unique and must not be mixed:

- a_H = r = 8 nm = 0.008 µm — Option A* GL/GPE healing freeze
- ξ_ph = ħ c_s / Δ₀ — phonon healing; η_ret(ξ_ph) = 2 identically
- ξ_181 = 0.181 µm — preprint causal-propagator coherence that produced
  η_ret ≈ 3.1. Origin: ħ c_s / Δ_QVC with withdrawn Δ_QVC = 0.389 meV
  (0.07/0.389 ≈ 0.180 µm, rounded to 0.181). Alternate, not the lock.

η_ret(a_H) ≈ 0.137 is not the identity 2. Never box 3.1 as the unique value.
"""

from __future__ import annotations

from typing import Any

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.retardation import eta_ret  # noqa: E402

from qvc_l1.locks import (
    ETA_RET_ALTERNATE,
    FOLKLORE_DELTA_QVC_MEV,
    XI_ALTERNATE_UM,
    L1Locks,
)

LENGTH_NAMES = ("a_H", "xi_ph", "xi_181")


def phonon_healing_um(locks: L1Locks) -> float:
    """ξ_ph = ħ c_s / Δ₀  [µm]. Status: derived."""
    return float(locks.hbar_cs_meV_um / locks.Delta0_meV)


def xi_181_origin_um(locks: L1Locks) -> float:
    """ħ c_s / Δ_QVC with withdrawn Δ_QVC = 0.389 meV.

    This is the source of the preprint 0.181 µm figure, not a lock.
    """
    return float(locks.hbar_cs_meV_um / FOLKLORE_DELTA_QVC_MEV)


def length_dictionary(locks: L1Locks) -> dict[str, Any]:
    """Three-length card. Status tags are part of the contract."""
    a_H = float(locks.a_H_um)
    xi_ph = phonon_healing_um(locks)
    xi_181 = float(XI_ALTERNATE_UM)
    xi_from_withdrawn = xi_181_origin_um(locks)

    eta_aH = float(eta_ret(locks.Delta0_meV, a_H, locks.hbar_cs_meV_um))
    eta_ph = float(eta_ret(locks.Delta0_meV, xi_ph, locks.hbar_cs_meV_um))
    eta_181 = float(eta_ret(locks.Delta0_meV, xi_181, locks.hbar_cs_meV_um))
    eta_from_withdrawn = float(
        eta_ret(locks.Delta0_meV, xi_from_withdrawn, locks.hbar_cs_meV_um)
    )

    return {
        "formula": "η_ret = 2 Δ₀ ξ / (ħ c_s)",
        "formula_status": "locked",
        "lengths": {
            "a_H": {
                "value_um": a_H,
                "value_nm": 1e3 * a_H,
                "eta_ret": eta_aH,
                "status": "locked",
                "origin": "Option A* GL/GPE healing freeze a_H = r (fat-torus tube)",
                "is_lock": True,
            },
            "xi_ph": {
                "value_um": xi_ph,
                "value_nm": 1e3 * xi_ph,
                "eta_ret": eta_ph,
                "status": "derived",
                "origin": "phonon healing ξ_ph = ħ c_s / Δ₀; η_ret = 2 identically",
                "is_lock": False,
                "identity_eta_ret_is_2": abs(eta_ph - 2.0) < 1e-12,
            },
            "xi_181": {
                "value_um": xi_181,
                "value_nm": 1e3 * xi_181,
                "eta_ret": eta_181,
                "status": "alternate",
                "origin": (
                    "Preprint causal-propagator coherence ξ=0.181 µm that produced "
                    "η_ret≈3.1. Numerically ħ c_s/Δ_QVC with withdrawn "
                    f"Δ_QVC={FOLKLORE_DELTA_QVC_MEV} meV "
                    f"(ħ c_s/Δ_QVC={xi_from_withdrawn:.5f} µm). Not the lock."
                ),
                "is_lock": False,
                "eta_ret_boxed_3p1_is_unique": False,
            },
        },
        "eta_ret_a_H": eta_aH,
        "eta_ret_xi_ph": eta_ph,
        "eta_ret_xi_181": eta_181,
        "eta_ret_from_withdrawn_Delta_QVC": eta_from_withdrawn,
        "alternate_eta_ret_label": ETA_RET_ALTERNATE,
        "unique_lock_is_not_3p1": True,
        "a_H_is_not_xi_181": abs(a_H - xi_181) > 0.05,
        "xi_ph_is_not_xi_181": abs(xi_ph - xi_181) > 0.02,
        "notes": (
            "Three lengths exist. Mixing ξ=181 nm into a_H=8 nm is a convention error. "
            "η_ret(a_H)≈0.137 is not the identity 2. η_ret(ξ_ph)=2. "
            "η_ret≈3.1 is the alternate ξ=0.181 evaluation, not the unique lock."
        ),
    }
