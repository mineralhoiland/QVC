"""Honest TQGL numerical suite (v4) for Quantum Vacuum Catalyzation.

Solves a truncated 2D GPE plus a coupled DCE amplitude ODE. Static
comparison is the saddle-node fold at frozen Option A* (fat torus).
Chern–Simons, Berry, and Skyrme terms are *not* in the EOM.

See ``output/tqgl_v4/METHODS.md`` for the equations actually solved.
"""

from qvc.tqgl.locks import FrozenTQGL, build_frozen_tqgl
from qvc.tqgl.live_gap import coherence_g1, delta_live_from_psi

__all__ = [
    "FrozenTQGL",
    "build_frozen_tqgl",
    "coherence_g1",
    "delta_live_from_psi",
]
