"""BKT C1–C6 computational verification suite for QVC / TEVC vortex RG.

Imports solvers from QVCCompute (``qvc``). This package is ``qvc_bkt`` so it
does not shadow that tree. Submodules: ``stiffness`` (C1), ``vortex_core``
(C2), ``kt_flow`` (C3), ``xi_kt`` (C4), ``model_compare`` (C5),
``phonon_damping`` (C6), ``suite`` (driver).
"""

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc_bkt.locks import BKTLocks, build_bkt_locks  # noqa: E402

__all__ = [
    "BKTLocks",
    "build_bkt_locks",
    "ensure_qvccompute",
]
