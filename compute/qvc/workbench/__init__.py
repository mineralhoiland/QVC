"""qvc.workbench — TEVC / QVC–TDVT 3D crystal and materials workbench.

Design rule (S0 of the simulation plan): every viewport reads its numbers
from :mod:`qvc.workbench.locks` (Option A* freeze built from the existing
``qvc`` solvers), is stamped with an epistemic layer (L1 / L2 / L3), and
refuses withdrawn folklore numbers.

Physics lives in the existing ``qvc`` modules. This package only wires
them to 3D fields, graphs, and exporters (VESTA / ParaView / Panel).

Submodules
----------
locks       S0 — frozen locks, reject list, layer stamps, unit conventions
catalyzon   S3 — democratic Spec(G), Helmert dark basis, WS overlap comparison
barriers    S4 — Mexican hat, fold map, WKB/Arrhenius sketch, A_ZPF reshaping
dynamics    S5 — qvc.tqgl GPE↔DCE replay, 3D torus extrusion, time graphs
topology    S6 — Hall staircase, Berry curvature, preimage linking
analogue    S7 — L2 acoustic metric, sonic surface, torsion proxy, Goldstones
app         S8 — Panel four-viewport workbench (``python -m qvc.workbench``)
"""

from qvc.workbench.locks import (
    DEPRECATED,
    Layer,
    WorkbenchLocks,
    assert_not_deprecated,
    get_locks,
    stamp,
)

__all__ = [
    "DEPRECATED",
    "Layer",
    "WorkbenchLocks",
    "assert_not_deprecated",
    "get_locks",
    "stamp",
]
