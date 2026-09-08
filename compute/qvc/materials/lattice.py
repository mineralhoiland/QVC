"""
qvc.materials.lattice — Graphene lattice primitives.

Provides the fundamental lattice operations for graphene:
- Unit cell construction (2-atom basis on hexagonal lattice)
- Lattice vector generation
- Real-space honeycomb coordinate generation for visualization
- Rotation operations for twisted layers

Units: Angstroms for lengths, degrees for angles.

The graphene lattice constant a = √3 × a_CC where a_CC = 1.42 Å
is the carbon-carbon bond length. The unit cell contains two atoms
(A and B sublattices) at positions:
    τ_A = (0, 0)
    τ_B = (a_CC, 0) = (1/3, 2/3) in fractional coords

This module has NO external dependencies beyond numpy.
"""

from __future__ import annotations

import numpy as np
from typing import Tuple


# === Constants ===

A_CC: float = 1.42  # C-C bond length in Angstroms
A_LATTICE: float = A_CC * np.sqrt(3)  # lattice constant = 2.4595 Å
D_INTERLAYER: float = 3.35  # interlayer spacing in Angstroms (graphite)


def graphene_lattice_vectors(a: float = A_LATTICE) -> np.ndarray:
    """Primitive lattice vectors for graphene (2D, in Angstroms).

    Convention: a₁ along x, a₂ at 60° (hexagonal).

        a₁ = a (1, 0)
        a₂ = a (1/2, √3/2)

    Parameters
    ----------
    a : float
        Lattice constant (default: 2.4595 Å for graphene).

    Returns
    -------
    vecs : ndarray, shape (2, 2)
        Row vectors: vecs[0] = a₁, vecs[1] = a₂.
    """
    return a * np.array([
        [1.0, 0.0],
        [0.5, np.sqrt(3) / 2],
    ])


def reciprocal_vectors(a: float = A_LATTICE) -> np.ndarray:
    """Reciprocal lattice vectors for graphene.

        b₁ = (2π/a)(1, -1/√3)
        b₂ = (2π/a)(0, 2/√3)

    Parameters
    ----------
    a : float
        Real-space lattice constant.

    Returns
    -------
    bvecs : ndarray, shape (2, 2)
        Row vectors: bvecs[0] = b₁, bvecs[1] = b₂.
    """
    return (2 * np.pi / a) * np.array([
        [1.0, -1.0 / np.sqrt(3)],
        [0.0, 2.0 / np.sqrt(3)],
    ])


def graphene_unit_cell(a_cc: float = A_CC) -> Tuple[np.ndarray, np.ndarray]:
    """Two-atom graphene unit cell.

    Parameters
    ----------
    a_cc : float
        C-C bond length.

    Returns
    -------
    lattice_vecs : ndarray, shape (2, 2)
        Primitive lattice vectors (row vectors).
    basis : ndarray, shape (2, 2)
        Atomic positions in Cartesian coords: basis[0] = A site, basis[1] = B site.
    """
    a = a_cc * np.sqrt(3)
    vecs = graphene_lattice_vectors(a)

    # Basis atoms in Cartesian coordinates
    # A sublattice at origin, B sublattice displaced by (a_cc, 0) rotated
    # In fractional coords: A = (1/3, 2/3), B = (2/3, 1/3)
    # Converting to Cartesian:
    tau_A = (1 / 3) * vecs[0] + (2 / 3) * vecs[1]
    tau_B = (2 / 3) * vecs[0] + (1 / 3) * vecs[1]

    basis = np.array([tau_A, tau_B])
    return vecs, basis


def honeycomb_positions(
    n1: int = 10,
    n2: int = 10,
    a_cc: float = A_CC,
    center: bool = True,
) -> np.ndarray:
    """Generate honeycomb lattice positions over n1 × n2 unit cells.

    Parameters
    ----------
    n1, n2 : int
        Number of unit cells along a₁ and a₂ directions.
    a_cc : float
        C-C bond length.
    center : bool
        If True, center the positions around the origin.

    Returns
    -------
    positions : ndarray, shape (2 * n1 * n2, 2)
        All atomic positions in Cartesian coordinates (Angstroms).
    """
    vecs, basis = graphene_unit_cell(a_cc)

    positions = []
    for i in range(n1):
        for j in range(n2):
            R = i * vecs[0] + j * vecs[1]
            positions.append(R + basis[0])
            positions.append(R + basis[1])

    positions = np.array(positions)

    if center:
        positions -= positions.mean(axis=0)

    return positions


def rotate_positions(positions: np.ndarray, theta_deg: float) -> np.ndarray:
    """Rotate 2D positions by angle θ (degrees).

    Parameters
    ----------
    positions : ndarray, shape (N, 2)
        Atomic positions.
    theta_deg : float
        Rotation angle in degrees (counterclockwise positive).

    Returns
    -------
    rotated : ndarray, shape (N, 2)
        Rotated positions.
    """
    theta = np.radians(theta_deg)
    R = np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)],
    ])
    return positions @ R.T


def rotation_matrix_2d(theta_deg: float) -> np.ndarray:
    """2D rotation matrix for angle θ (degrees)."""
    theta = np.radians(theta_deg)
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)],
    ])


def high_symmetry_points_moire(K_M: float) -> dict:
    """High-symmetry points in the moiré Brillouin zone.

    Parameters
    ----------
    K_M : float
        Moiré reciprocal lattice constant |K_M| = 4π/(3λ_M).

    Returns
    -------
    points : dict
        Mapping from label to (kx, ky) coordinates.
    """
    # For hexagonal mBZ:
    # Γ = (0, 0)
    # K_M = K_M * (1, 0)  [one of the K points]
    # M = midpoint of Γ-K edge
    return {
        "Γ": np.array([0.0, 0.0]),
        "K": K_M * np.array([1.0, 0.0]),
        "K'": K_M * np.array([0.5, np.sqrt(3) / 2]),
        "M": K_M * np.array([0.75, np.sqrt(3) / 4]),
    }
