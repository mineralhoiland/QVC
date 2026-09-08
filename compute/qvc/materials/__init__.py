"""
qvc.materials — Crystallographic layer for TEVC physics.

Tier 1 of the QVC computational framework:
- Moiré supercell construction (pymatgen)
- Bistritzer-MacDonald band structure
- Folded phonon dispersion
- Berry curvature and Chern number
- LAMMPS relaxation interface
- Parameter bridge to effective theory (Tier 2)

Dependencies: numpy (always), pymatgen (for structure I/O),
ase (alternative), phonopy (phonon calculations).
"""

from qvc.materials.lattice import (
    graphene_lattice_vectors,
    graphene_unit_cell,
    honeycomb_positions,
)
from qvc.materials.moire import (
    moire_period,
    commensurate_indices,
    build_moire_bilayer,
    build_tevc_stack,
)
from qvc.materials.bistritzer_macdonald import (
    BistritzerMacDonald,
    matbg_model,
    tevc_model,
    k_path_high_symmetry,
    fukui_chern,
)
