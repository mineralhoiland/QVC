"""Canonical physical parameters for QVC / MATBG-like TEVC platform.

Status
------
- Coupling, fold, CS/linking identities: Tier C closed-form locks.
- Option A* fat torus (r, R) = (8, 50) nm is frozen working in
  ``qvc.tqgl.locks``. Candidate B is audit-only.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class QVCParams:
    """Working parameter set used by verification and demos.

    Units
    -----
    Energy: meV. Length: µm unless a function documents SI/nm.
    ``hbar_cs_meV_um`` is the optical/phonon scale from the preprint
    (ħ c_s ≈ 0.07 meV·µm). Geometry (r, R) for Option A* is the fat torus
    (8, 50) nm, frozen working in ``qvc.tqgl.locks``. Candidate B remains
    recorded for audit only.
    """

    J_meV: float = 0.300
    hbar_cs_meV_um: float = 0.070
    N_layers: int = 5
    lambda_M_um: float = 0.013
    gamma_E: float = 0.5772156649015329
    kB_meV_per_K: float = 0.08617333262

    # Candidate A — preprint fat-torus (healing / coherence scale)
    r_preprint_um: float = 0.008  # 8 nm
    R_preprint_um: float = 0.050  # 50 nm

    # Candidate B — Corrections local-curvature / moiré scale
    r_corr_um: float = 0.000335  # 0.335 nm
    R_corr_um: float = 0.0134  # 13.4 nm

    n_modes_default: int = 14

    @property
    def Delta0_meV(self) -> float:
        """Bare excitonic gap scale Δ₀ = 2J."""
        return 2.0 * self.J_meV


PARAMS = QVCParams()
