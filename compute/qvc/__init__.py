"""
QVC compute-proof package (Tier B/C scaffolding).

Closed-form identities live in dedicated modules and are re-exported from
``canonical_formalism``. Vacuum meV numerics are *provisional* until a
geometry is frozen via ``vacuum_freeze.freeze``.
"""

from qvc.params import PARAMS, QVCParams
from qvc.canonical_formalism import (
    CLOSED_FORM_LOCKS,
    democratic_spectrum,
    fold_condition,
    chern_simons_flat_u1,
    linking_form_self,
    fractional_hall_charge,
    vacuum_self_energy_corrected,
    eta_ret,
    A_eff_mode_enhancement,
)

__version__ = "0.1.0"

__all__ = [
    "PARAMS",
    "QVCParams",
    "CLOSED_FORM_LOCKS",
    "democratic_spectrum",
    "fold_condition",
    "chern_simons_flat_u1",
    "linking_form_self",
    "fractional_hall_charge",
    "vacuum_self_energy_corrected",
    "eta_ret",
    "A_eff_mode_enhancement",
    "__version__",
]
