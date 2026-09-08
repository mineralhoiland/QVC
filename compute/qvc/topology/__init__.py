"""
qvc.topology — Topological invariants and transport for TEVC.

Tier 3: computes observables from the effective theory.
- Chern number (lattice gauge / Fukui method)
- Hopf charge (linking number on discretized S³→S²)
- Hall staircase σ_xy = e²/(2Nh)
- Berry phase statistics (GUE/GOE random-matrix comparison)
"""

from qvc.topology.chern_hopf_hall import (
    fukui_chern_number,
    hopfion_field,
    hopf_charge,
    hall_staircase,
)
