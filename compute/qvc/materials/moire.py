"""
qvc.materials.moire — Moiré supercell construction for TEVC.

Builds twisted multilayer graphene structures at arbitrary twist angles.
For the TEVC platform: N=5 layers with alternating ±θ twist.

The moiré period is:
    λ_M = a / (2 sin(θ/2))

At θ = 1.1°: λ_M ≈ 12.8 nm ≈ 128 Å, containing ~11,000 atoms per bilayer.

This module provides:
1. Analytic moiré geometry (period, vectors, BZ)
2. Commensurate approximant finder (smallest supercell matching twist)
3. Explicit atomic position generation (for LAMMPS/visualization)
4. pymatgen Structure export (for DFT input generation)

The commensurate condition for twisted bilayer graphene:
    cos(θ) = (3m² + 3mr + r²/2) / (3m² + 3mr + r²)
where (m, r) are integers with r = 1 giving the smallest cells.
For θ near 1.1°, m ≈ 30, giving N_atoms ≈ 4(3m² + 3m + 1) ≈ 11,164.

Units: Angstroms, degrees.
"""

from __future__ import annotations

import numpy as np
from typing import Tuple, Optional

from qvc.materials.lattice import (
    A_CC,
    A_LATTICE,
    D_INTERLAYER,
    graphene_lattice_vectors,
    graphene_unit_cell,
    honeycomb_positions,
    rotate_positions,
    rotation_matrix_2d,
)


# === Moiré Geometry (Analytic) ===


def moire_period(theta_deg: float, a: float = A_LATTICE) -> float:
    """Moiré superlattice period λ_M = a / (2 sin(θ/2)).

    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees.
    a : float
        Graphene lattice constant (2.4595 Å).

    Returns
    -------
    lambda_M : float
        Moiré period in Angstroms.

    Examples
    --------
    >>> moire_period(1.1)  # MATBG magic angle
    128.03...  # ≈ 12.8 nm
    """
    theta = np.radians(theta_deg)
    return a / (2 * np.sin(theta / 2))


def moire_reciprocal_constant(theta_deg: float, a: float = A_LATTICE) -> float:
    """Moiré reciprocal lattice constant K_M = 4π / (3 λ_M).

    This is the magnitude of the moiré reciprocal vectors.
    The moiré BZ is a hexagon with this scale.

    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees.

    Returns
    -------
    K_M : float
        In units of 1/Angstrom.
    """
    lam = moire_period(theta_deg, a)
    return 4 * np.pi / (3 * lam)


def moire_lattice_vectors(theta_deg: float, a: float = A_LATTICE) -> np.ndarray:
    """Moiré superlattice vectors (real space).

    The moiré lattice vectors are:
        L₁ = λ_M (1, 0)
        L₂ = λ_M (1/2, √3/2)

    (Same orientation as graphene lattice, scaled by λ_M/a.)

    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees.

    Returns
    -------
    L : ndarray, shape (2, 2)
        Moiré lattice vectors (row vectors), in Angstroms.
    """
    lam = moire_period(theta_deg, a)
    return lam * np.array([
        [1.0, 0.0],
        [0.5, np.sqrt(3) / 2],
    ])


def moire_unit_cell_area(theta_deg: float, a: float = A_LATTICE) -> float:
    """Area of moiré unit cell in Å².

    A_M = (√3/2) λ_M²
    """
    lam = moire_period(theta_deg, a)
    return (np.sqrt(3) / 2) * lam**2


def atoms_per_moire_cell(theta_deg: float, n_layers: int = 2) -> int:
    """Approximate number of atoms in one moiré unit cell.

    For a bilayer: N ≈ 2 × 2 × (λ_M / a)²  (2 layers × 2 atoms/cell × cells)
    More precisely: N = 4 × (3m² + 3m + 1) for commensurate (m, 1).

    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees.
    n_layers : int
        Number of layers (each contributes ~N_bilayer/2 atoms).
    """
    lam = moire_period(theta_deg)
    cells_per_moire = (lam / A_LATTICE) ** 2  # approximate
    atoms_per_layer = 2 * cells_per_moire  # 2 atoms per graphene unit cell
    return int(np.round(atoms_per_layer * n_layers))


# === Commensurate Approximant ===


def commensurate_indices(theta_deg: float, r: int = 1, tol: float = 0.01) -> Tuple[int, int, float]:
    """Find commensurate (m, r) indices closest to target twist angle.

    The commensurate condition:
        cos(θ) = (3m² + 3mr + r²/2) / (3m² + 3mr + r²)

    We search for the m that gives the closest match to theta_deg.

    Parameters
    ----------
    theta_deg : float
        Target twist angle in degrees.
    r : int
        Sublattice index (1 for smallest cells, 3 for 60° family).
    tol : float
        Acceptable angle mismatch in degrees.

    Returns
    -------
    m : int
        Commensurate index.
    n_atoms : int
        Number of atoms in one bilayer moiré cell = 4(3m² + 3mr + r²).
    actual_theta : float
        Actual commensurate angle in degrees.
    """
    theta_target = np.radians(theta_deg)
    cos_target = np.cos(theta_target)

    best_m = None
    best_diff = np.inf

    for m in range(1, 200):
        num = 3 * m**2 + 3 * m * r + r**2 / 2
        den = 3 * m**2 + 3 * m * r + r**2
        cos_theta = num / den
        diff = abs(np.arccos(cos_theta) - theta_target)

        if diff < best_diff:
            best_diff = diff
            best_m = m

    # Compute actual angle
    num = 3 * best_m**2 + 3 * best_m * r + r**2 / 2
    den = 3 * best_m**2 + 3 * best_m * r + r**2
    actual_theta_deg = np.degrees(np.arccos(num / den))

    n_atoms_bilayer = 4 * (3 * best_m**2 + 3 * best_m * r + r**2)

    return best_m, n_atoms_bilayer, actual_theta_deg


# === Explicit Position Generation ===


def build_moire_bilayer(
    theta_deg: float = 1.1,
    n_cells: int = 1,
    a_cc: float = A_CC,
    d_perp: float = D_INTERLAYER,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build atomic positions for twisted bilayer graphene.

    Layer 1 is rotated by +θ/2, layer 2 by -θ/2, giving relative twist θ.

    Parameters
    ----------
    theta_deg : float
        Total twist angle between layers.
    n_cells : int
        Number of moiré cells along each direction (1 = single moiré cell).
    a_cc : float
        C-C bond length.
    d_perp : float
        Interlayer separation.

    Returns
    -------
    positions_1 : ndarray, shape (N, 3)
        Layer 1 atomic positions [x, y, z] in Angstroms.
    positions_2 : ndarray, shape (N, 3)
        Layer 2 atomic positions [x, y, z] in Angstroms.
    moire_vecs : ndarray, shape (2, 2)
        Moiré lattice vectors for periodic boundary conditions.
    """
    # Determine how many graphene unit cells tile one moiré period
    lam = moire_period(theta_deg, a_cc * np.sqrt(3))
    a = a_cc * np.sqrt(3)
    n_per_side = int(np.ceil(lam / a)) * n_cells + 2  # extra for clipping

    # Generate large honeycomb patch
    pos = honeycomb_positions(n_per_side, n_per_side, a_cc, center=True)

    # Rotate layers
    pos_1_2d = rotate_positions(pos, +theta_deg / 2)
    pos_2_2d = rotate_positions(pos, -theta_deg / 2)

    # Clip to moiré cell (hexagonal region of radius λ_M * n_cells)
    # Use a circular approximation for simplicity; exact hexagonal clip available
    r_max = lam * n_cells * 0.6  # conservative clip radius
    mask_1 = np.linalg.norm(pos_1_2d, axis=1) < r_max
    mask_2 = np.linalg.norm(pos_2_2d, axis=1) < r_max

    pos_1_2d = pos_1_2d[mask_1]
    pos_2_2d = pos_2_2d[mask_2]

    # Add z-coordinates
    z1 = np.full(len(pos_1_2d), 0.0)
    z2 = np.full(len(pos_2_2d), d_perp)

    positions_1 = np.column_stack([pos_1_2d, z1])
    positions_2 = np.column_stack([pos_2_2d, z2])

    moire_vecs = moire_lattice_vectors(theta_deg, a)

    return positions_1, positions_2, moire_vecs


def build_tevc_stack(
    theta_deg: float = 1.1,
    n_layers: int = 5,
    d_perp: float = D_INTERLAYER,
    a_cc: float = A_CC,
    alternating: bool = True,
) -> Tuple[list, np.ndarray]:
    """Build full TEVC multilayer structure.

    Generates N layers with alternating twist angles:
        Layer 1: +θ/2
        Layer 2: -θ/2
        Layer 3: +θ/2
        Layer 4: -θ/2
        Layer 5: +θ/2

    This produces adjacent-layer twist of θ, with the full stack
    having moiré period λ_M = a / (2 sin(θ/2)).

    Parameters
    ----------
    theta_deg : float
        Twist angle between adjacent layers (degrees).
    n_layers : int
        Number of graphene layers.
    d_perp : float
        Interlayer spacing (Angstroms).
    a_cc : float
        C-C bond length.
    alternating : bool
        If True, twist angles alternate ±θ/2.
        If False, each layer is twisted by θ relative to the one below.

    Returns
    -------
    layer_positions : list of ndarray
        Each element is shape (N_i, 3) with atomic positions for layer i.
    moire_vecs : ndarray, shape (2, 2)
        Moiré superlattice vectors.
    """
    lam = moire_period(theta_deg, a_cc * np.sqrt(3))
    a = a_cc * np.sqrt(3)
    n_per_side = int(np.ceil(lam / a)) + 2

    # Generate base honeycomb
    pos_base = honeycomb_positions(n_per_side, n_per_side, a_cc, center=True)

    # Clip radius
    r_max = lam * 0.6

    layer_positions = []

    for i in range(n_layers):
        if alternating:
            # Alternating: +θ/2, -θ/2, +θ/2, ...
            angle = (theta_deg / 2) * ((-1) ** i)
        else:
            # Cumulative: each layer twisted by θ from bottom
            angle = i * theta_deg

        pos_2d = rotate_positions(pos_base, angle)

        # Clip to moiré cell
        mask = np.linalg.norm(pos_2d, axis=1) < r_max
        pos_2d = pos_2d[mask]

        # Z coordinate
        z = np.full(len(pos_2d), i * d_perp)
        pos_3d = np.column_stack([pos_2d, z])

        layer_positions.append(pos_3d)

    moire_vecs = moire_lattice_vectors(theta_deg, a)

    return layer_positions, moire_vecs


# === Stacking Registry Analysis ===


def stacking_registry(
    pos_top: np.ndarray,
    pos_bot: np.ndarray,
    a_cc: float = A_CC,
) -> np.ndarray:
    """Compute local stacking registry between two layers.

    The stacking vector d(r) = (dx, dy) gives the relative lateral
    displacement between layers at each point. This determines:
        - AA stacking: d = 0
        - AB stacking: d = (a_CC, 0) or equivalent
        - SP (saddle point): intermediate

    For the TEVC, the AA regions are where the toroidal Casimir
    cavities form (Section 4 of the paper).

    Parameters
    ----------
    pos_top : ndarray, shape (N, 3) or (N, 2)
        Top layer positions.
    pos_bot : ndarray, shape (M, 3) or (M, 2)
        Bottom layer positions.
    a_cc : float
        C-C bond length (for normalization).

    Returns
    -------
    registry_map : ndarray
        Stacking character: 0 = AA, 1 = AB, intermediate = SP domain wall.
        (Requires interpolation onto regular grid for visualization.)
    """
    # This is a simplified version; full implementation uses
    # Voronoi tessellation or grid interpolation
    raise NotImplementedError(
        "Full stacking registry analysis requires grid interpolation. "
        "Use LAMMPS relaxation output for accurate domain maps."
    )


# === pymatgen Interface ===


def to_pymatgen_structure(
    layer_positions: list,
    moire_vecs: np.ndarray,
    d_perp: float = D_INTERLAYER,
    vacuum: float = 20.0,
):
    """Convert TEVC layer positions to pymatgen Structure.

    Requires pymatgen to be installed (optional dependency).

    Parameters
    ----------
    layer_positions : list of ndarray
        Output from build_tevc_stack.
    moire_vecs : ndarray, shape (2, 2)
        Moiré lattice vectors (2D).
    d_perp : float
        Interlayer spacing.
    vacuum : float
        Vacuum spacing above/below slab (Angstroms).

    Returns
    -------
    structure : pymatgen.core.Structure
        Full 3D periodic structure suitable for DFT or LAMMPS input.
    """
    try:
        from pymatgen.core import Structure, Lattice
    except ImportError:
        raise ImportError(
            "pymatgen required for Structure export. "
            "Install with: pip install pymatgen"
        )

    n_layers = len(layer_positions)
    total_height = (n_layers - 1) * d_perp + 2 * vacuum

    # Build 3D lattice: moiré vectors in-plane, vacuum out-of-plane
    a1_3d = np.array([moire_vecs[0, 0], moire_vecs[0, 1], 0.0])
    a2_3d = np.array([moire_vecs[1, 0], moire_vecs[1, 1], 0.0])
    a3_3d = np.array([0.0, 0.0, total_height])

    lattice = Lattice([a1_3d, a2_3d, a3_3d])

    # Collect all positions and species
    all_coords = []
    all_species = []

    for positions in layer_positions:
        for pos in positions:
            # Convert Cartesian to fractional
            cart = np.array([pos[0], pos[1], pos[2] + vacuum])
            all_coords.append(cart)
            all_species.append("C")

    # Create structure with Cartesian coordinates
    structure = Structure(
        lattice,
        all_species,
        all_coords,
        coords_are_cartesian=True,
    )

    return structure


def to_lammps_data(
    layer_positions: list,
    moire_vecs: np.ndarray,
    filename: str = "tevc.data",
    d_perp: float = D_INTERLAYER,
):
    """Write TEVC structure as LAMMPS data file.

    Uses 'atomic' style with atom types = layer index.

    Parameters
    ----------
    layer_positions : list of ndarray
        Output from build_tevc_stack.
    moire_vecs : ndarray, shape (2, 2)
        Moiré lattice vectors.
    filename : str
        Output filename.
    d_perp : float
        Interlayer spacing.
    """
    n_layers = len(layer_positions)
    n_atoms = sum(len(lp) for lp in layer_positions)
    total_height = (n_layers - 1) * d_perp + 40.0  # 20 Å vacuum each side

    # Triclinic box from moiré vectors
    # LAMMPS triclinic: xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    ax, ay = moire_vecs[0]
    bx, by = moire_vecs[1]

    xlo, xhi = 0.0, ax
    ylo, yhi = 0.0, by
    zlo, zhi = -20.0, total_height - 20.0
    xy = bx  # tilt factor

    with open(filename, 'w') as f:
        f.write(f"TEVC {n_layers}-layer moire structure\n\n")
        f.write(f"{n_atoms} atoms\n")
        f.write(f"{n_layers} atom types\n\n")
        f.write(f"{xlo:.6f} {xhi:.6f} xlo xhi\n")
        f.write(f"{ylo:.6f} {yhi:.6f} ylo yhi\n")
        f.write(f"{zlo:.6f} {zhi:.6f} zlo zhi\n")
        f.write(f"{xy:.6f} 0.000000 0.000000 xy xz yz\n\n")
        f.write("Masses\n\n")
        for i in range(1, n_layers + 1):
            f.write(f"{i} 12.011  # C (layer {i})\n")
        f.write("\nAtoms\n\n")

        atom_id = 0
        for layer_idx, positions in enumerate(layer_positions):
            for pos in positions:
                atom_id += 1
                f.write(
                    f"{atom_id} {layer_idx + 1} "
                    f"{pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n"
                )
