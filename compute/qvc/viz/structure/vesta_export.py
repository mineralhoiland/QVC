"""S1 — Geometry construction and VESTA / OVITO export.

What this module guarantees
---------------------------
* The **magic-angle TEVC** (θ ≈ 1.1°, λ_M ≈ 12.8 nm, 5 layers) is built from
  :func:`qvc.materials.moire.build_tevc_stack`. The 13–22° commensurate
  CIFs in ``QVCCursor/TEVC_structures`` are *proxies* and are never used as
  the TEVC crystal.
* CIF, POSCAR and LAMMPS files are written from the same position arrays.
* The hopfion BCC crystal is written in its proper Im-3m (No. 229) setting
  with one Wyckoff-2a site; dummy He/Ne preimage tori live in a *separate*
  P1 overlay CIF and are labelled as schematic, not chemistry.
* Gaussian CUBE writer converts nm → bohr so VESTA reads the volume at
  the right scale.

Units: Å for atomic structures (CIF/POSCAR/LAMMPS); nm in for CUBE (converted).
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import numpy as np

from qvc.materials.lattice import A_CC, D_INTERLAYER
from qvc.materials.moire import (
    build_tevc_stack,
    moire_lattice_vectors,
    moire_period,
    to_lammps_data,
)

BOHR_PER_NM = 18.897261246
BOHR_PER_ANGSTROM = 1.8897261246

# --------------------------------------------------------------------------
# Low-level writers (no pymatgen dependency)
# --------------------------------------------------------------------------


def _lattice_params(lat: np.ndarray) -> tuple[float, float, float, float, float, float]:
    a, b, c = (float(np.linalg.norm(v)) for v in lat)

    def ang(u: np.ndarray, v: np.ndarray) -> float:
        return math.degrees(
            math.acos(float(np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1, 1)))
        )

    return a, b, c, ang(lat[1], lat[2]), ang(lat[0], lat[2]), ang(lat[0], lat[1])


def write_cif(
    path: str | Path,
    lattice: np.ndarray,
    species: Sequence[str],
    frac_coords: np.ndarray,
    *,
    data_name: str = "structure",
    space_group: str = "P 1",
    sg_number: int = 1,
    comment: str = "",
    labels: Sequence[str] | None = None,
    u_iso: float = 0.005,
) -> Path:
    """Write a minimal CIF (VESTA-readable). ``lattice`` rows are Å vectors."""
    path = Path(path)
    lattice = np.asarray(lattice, dtype=float)
    frac = np.asarray(frac_coords, dtype=float) % 1.0
    a, b, c, al, be, ga = _lattice_params(lattice)
    if labels is None:
        counts: dict[str, int] = {}
        labels = []
        for s in species:
            counts[s] = counts.get(s, 0) + 1
            labels.append(f"{s}{counts[s]}")
    with path.open("w") as f:
        if comment:
            for line in comment.splitlines():
                f.write(f"# {line}\n")
        f.write(f"data_{data_name}\n\n")
        f.write(f"_cell_length_a    {a:.6f}\n_cell_length_b    {b:.6f}\n_cell_length_c    {c:.6f}\n")
        f.write(f"_cell_angle_alpha {al:.4f}\n_cell_angle_beta  {be:.4f}\n_cell_angle_gamma {ga:.4f}\n\n")
        f.write(f"_symmetry_space_group_name_H-M   '{space_group}'\n")
        f.write(f"_symmetry_Int_Tables_number       {sg_number}\n\n")
        f.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
        f.write("_atom_site_U_iso_or_equiv\n_atom_site_occupancy\n")
        for lab, sym, (x, y, z) in zip(labels, species, frac):
            f.write(f"{lab:8s} {sym:4s} {x:10.6f} {y:10.6f} {z:10.6f} {u_iso:.4f} 1.0000\n")
    return path


def write_poscar(
    path: str | Path,
    lattice: np.ndarray,
    species: Sequence[str],
    frac_coords: np.ndarray,
    *,
    comment: str = "TEVC",
) -> Path:
    """Write VASP POSCAR (Direct coordinates); species grouped in order."""
    path = Path(path)
    lattice = np.asarray(lattice, dtype=float)
    frac = np.asarray(frac_coords, dtype=float) % 1.0
    order = sorted(set(species), key=lambda s: species.index(s))
    idx = [i for s in order for i, sp in enumerate(species) if sp == s]
    with path.open("w") as f:
        f.write(f"{comment}\n1.0\n")
        for v in lattice:
            f.write(f"  {v[0]:.10f} {v[1]:.10f} {v[2]:.10f}\n")
        f.write("  " + " ".join(order) + "\n")
        f.write("  " + " ".join(str(species.count(s)) for s in order) + "\n")
        f.write("Direct\n")
        for i in idx:
            f.write(f"  {frac[i,0]:.10f} {frac[i,1]:.10f} {frac[i,2]:.10f}\n")
    return path


def write_cube(
    path: str | Path,
    data: np.ndarray,
    *,
    origin_nm: Sequence[float],
    spacing_nm: float,
    title: str,
    subtitle: str,
    atoms: Iterable[tuple[int, float, float, float]] | None = None,
) -> Path:
    """Write a Gaussian CUBE file (lengths converted nm → bohr for VESTA).

    ``atoms`` are (Z, x_nm, y_nm, z_nm). If none, one dummy atom (Z=0) is
    written at the grid centre so VESTA has a structure to attach the volume to.
    """
    path = Path(path)
    data = np.asarray(data, dtype=float)
    nx, ny, nz = data.shape
    o = np.asarray(origin_nm, dtype=float) * BOHR_PER_NM
    h = spacing_nm * BOHR_PER_NM
    if atoms is None:
        c = o + 0.5 * h * np.array([nx - 1, ny - 1, nz - 1])
        atoms = [(0, c[0] / BOHR_PER_NM, c[1] / BOHR_PER_NM, c[2] / BOHR_PER_NM)]
    atoms = list(atoms)
    with path.open("w") as f:
        f.write(f"{title}\n{subtitle}\n")
        f.write(f"{len(atoms):5d} {o[0]:12.6f} {o[1]:12.6f} {o[2]:12.6f}\n")
        f.write(f"{nx:5d} {h:12.6f} {0.0:12.6f} {0.0:12.6f}\n")
        f.write(f"{ny:5d} {0.0:12.6f} {h:12.6f} {0.0:12.6f}\n")
        f.write(f"{nz:5d} {0.0:12.6f} {0.0:12.6f} {h:12.6f}\n")
        for Z, x, y, z in atoms:
            f.write(
                f"{int(Z):5d} {float(Z):12.6f} {x*BOHR_PER_NM:12.6f} "
                f"{y*BOHR_PER_NM:12.6f} {z*BOHR_PER_NM:12.6f}\n"
            )
        flat = data.reshape(nx * ny, nz)
        for row in flat:
            for k in range(0, nz, 6):
                f.write("".join(f" {v:13.5E}" for v in row[k : k + 6]) + "\n")
    return path


# --------------------------------------------------------------------------
# TEVC magic-angle stack
# --------------------------------------------------------------------------


@dataclass
class TEVCGeometry:
    theta_deg: float
    n_layers: int
    lambda_M_A: float
    moire_vecs: np.ndarray  # (2,2) Å
    lattice: np.ndarray  # (3,3) Å, with vacuum
    layer_positions: list[np.ndarray]  # Å, cartesian (z from 0)
    cart: np.ndarray  # (N,3) Å, shifted by vacuum
    frac: np.ndarray  # (N,3)
    species: list[str]
    layer_index: np.ndarray
    vacuum_A: float
    fragment: dict[str, Any] = field(default_factory=dict)
    provenance: str = (
        "magic-angle TEVC from qvc.materials.moire.build_tevc_stack (unrelaxed, "
        "analytic ±θ/2 alternating stack). TEVC_structures 13–22° CIFs are proxies."
    )

    @property
    def n_atoms(self) -> int:
        return int(self.cart.shape[0])


def build_tevc_geometry(
    theta_deg: float = 1.1,
    n_layers: int = 5,
    *,
    d_perp: float = D_INTERLAYER,
    a_cc: float = A_CC,
    vacuum_A: float = 20.0,
    fragment_radius_nm: float | None = 3.0,
) -> TEVCGeometry:
    """Build the 5-layer TEVC cell and a downsampled interactive fragment.

    The fragment keeps atoms within ``fragment_radius_nm`` of the AA-site
    (moiré origin) — enough to show the AA/AB registry, light enough for
    an in-browser lattice viewport (~26k atoms full cell is too heavy).
    """
    layers, moire_vecs = build_tevc_stack(theta_deg, n_layers, d_perp, a_cc, alternating=True)
    lam = moire_period(theta_deg, a_cc * math.sqrt(3))
    height = (n_layers - 1) * d_perp + 2.0 * vacuum_A
    lattice = np.array(
        [
            [moire_vecs[0, 0], moire_vecs[0, 1], 0.0],
            [moire_vecs[1, 0], moire_vecs[1, 1], 0.0],
            [0.0, 0.0, height],
        ]
    )
    cart = np.vstack(layers).astype(float)
    cart[:, 2] += vacuum_A
    layer_index = np.concatenate([np.full(len(l), i) for i, l in enumerate(layers)])
    frac = cart @ np.linalg.inv(lattice)
    species = ["C"] * cart.shape[0]

    frag: dict[str, Any] = {}
    if fragment_radius_nm is not None:
        rad = fragment_radius_nm * 10.0
        # build_tevc_stack centres the honeycomb patch at the origin = AA site
        allpos = np.vstack(layers)
        mask = np.linalg.norm(allpos[:, :2], axis=1) < rad
        frag = {
            "radius_nm": fragment_radius_nm,
            "cart_A": allpos[mask],
            "layer_index": layer_index[mask],
            "n_atoms": int(mask.sum()),
            "label": "downsampled AA-site fragment for interactive atoms (not the periodic cell)",
        }
    return TEVCGeometry(
        theta_deg=theta_deg,
        n_layers=n_layers,
        lambda_M_A=lam,
        moire_vecs=moire_vecs,
        lattice=lattice,
        layer_positions=layers,
        cart=cart,
        frac=frac,
        species=species,
        layer_index=layer_index,
        vacuum_A=vacuum_A,
        fragment=frag,
    )


def export_tevc(
    outdir: str | Path,
    geom: TEVCGeometry | None = None,
    *,
    theta_deg: float = 1.1,
    n_layers: int = 5,
    write_lammps: bool = True,
    write_full_cif: bool = True,
) -> dict[str, Any]:
    """Write TEVC.cif, POSCAR, tevc.data and a fragment CIF from one array set."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    g = geom or build_tevc_geometry(theta_deg, n_layers)
    comment = (
        f"TEVC {g.n_layers}-layer twisted graphene, theta={g.theta_deg:.3f} deg, "
        f"lambda_M={g.lambda_M_A/10:.2f} nm, {g.n_atoms} C atoms, unrelaxed.\n{g.provenance}"
    )
    out: dict[str, Any] = {
        "theta_deg": g.theta_deg,
        "n_layers": g.n_layers,
        "lambda_M_nm": g.lambda_M_A / 10.0,
        "n_atoms": g.n_atoms,
        "layer": "L1",
        "provenance": g.provenance,
        "files": {},
    }
    if write_full_cif:
        out["files"]["cif"] = str(
            write_cif(
                outdir / "TEVC.cif", g.lattice, g.species, g.frac,
                data_name=f"TEVC_{g.n_layers}L_{g.theta_deg:.2f}deg", comment=comment,
            )
        )
        out["files"]["poscar"] = str(
            write_poscar(outdir / "POSCAR", g.lattice, g.species, g.frac, comment=comment.splitlines()[0])
        )
    if write_lammps:
        p = outdir / "tevc.data"
        to_lammps_data(g.layer_positions, g.moire_vecs, str(p))
        out["files"]["lammps"] = str(p)
    if g.fragment:
        fr = g.fragment
        rad = fr["radius_nm"] * 10.0
        box = np.diag([2.2 * rad, 2.2 * rad, (g.n_layers - 1) * D_INTERLAYER + 2 * g.vacuum_A])
        cart = fr["cart_A"].copy()
        cart[:, :2] += 1.1 * rad
        cart[:, 2] += g.vacuum_A
        frac = cart @ np.linalg.inv(box)
        out["files"]["fragment_cif"] = str(
            write_cif(
                outdir / "TEVC_fragment.cif", box, ["C"] * len(cart), frac,
                data_name="TEVC_AA_fragment",
                comment=f"{fr['label']}; radius {fr['radius_nm']} nm; {fr['n_atoms']} atoms; non-periodic box",
            )
        )
        out["fragment_n_atoms"] = fr["n_atoms"]
    (outdir / "TEVC_export.json").write_text(json.dumps(out, indent=2))
    return out


def tevc_pymatgen_structure(geom: TEVCGeometry):
    """Optional pymatgen Structure (for Crystal Toolkit / spglib)."""
    from pymatgen.core import Lattice, Structure

    return Structure(Lattice(geom.lattice), geom.species, geom.cart, coords_are_cartesian=True)


# --------------------------------------------------------------------------
# Hopfion BCC crystal (L2 schematic) — proper Im-3m + separate overlay
# --------------------------------------------------------------------------


def hopfion_bcc_cif(
    path: str | Path,
    *,
    a_nm: float = 50.0,
    display_scale_A_per_nm: float = 1.0,
) -> Path:
    """BCC hopfion crystal in Im-3m (No. 229): one Wyckoff-2a site.

    VESTA generates the body-centre from the space group. The site symbol
    'Ho' is a *display proxy* for a Q=1 hopfion centre, not holmium.
    Lattice constant a = ξ ≈ 50 nm is written as 50 display-Å (1 Å ↔ 1 nm).
    """
    a = a_nm * display_scale_A_per_nm
    lat = np.diag([a, a, a])
    return write_cif(
        path, lat, ["Ho"], np.array([[0.0, 0.0, 0.0]]),
        data_name="hopfion_crystal_BCC_Im3m",
        space_group="I m -3 m", sg_number=229,
        labels=["Hopf1"],
        comment=(
            "L2 SCHEMATIC: BCC hopfion crystal, target group Im-3m x U(1)_fiber (Paper C). "
            f"a = xi = {a_nm:.0f} nm written as {a:.0f} display-Angstrom. "
            "Site 'Ho' marks a Q=1 hopfion centre; it is not a chemical species. "
            "Wyckoff 2a; body centre generated by symmetry (replaces the dummy P1 file)."
        ),
        u_iso=0.05,
    )


def preimage_overlay_cif(
    path: str | Path,
    curves: Sequence[np.ndarray],
    *,
    a_nm: float = 50.0,
    display_scale_A_per_nm: float = 1.0,
    symbols: Sequence[str] = ("He", "Ne", "Ar", "Kr"),
    centre_frac: Sequence[float] = (0.5, 0.5, 0.5),
    curve_extent_nm: float | None = None,
) -> Path:
    """P1 overlay CIF tracing *computed* preimage curves with dummy sites.

    ``curves`` are (M,3) arrays in nm centred on the hopfion. They are scaled
    so the largest curve spans ~0.6 a and placed about ``centre_frac``.
    Keep this file separate from :func:`hopfion_bcc_cif`; it is a VESTA overlay
    of n⁻¹(p) fibres, not chemistry.
    """
    a = a_nm * display_scale_A_per_nm
    lat = np.diag([a, a, a])
    if curve_extent_nm is None:
        curve_extent_nm = max(float(np.max(np.abs(c))) for c in curves) or 1.0
    scale = 0.30 / curve_extent_nm  # largest radius → 0.3 a
    species: list[str] = []
    frac: list[np.ndarray] = []
    labels: list[str] = []
    for ci, c in enumerate(curves):
        sym = symbols[ci % len(symbols)]
        for k, p in enumerate(np.asarray(c, dtype=float)):
            frac.append(np.asarray(centre_frac) + scale * p)
            species.append(sym)
            labels.append(f"{sym}{ci+1}_{k+1}")
    return write_cif(
        path, lat, species, np.array(frac),
        data_name="hopfion_preimage_overlay",
        space_group="P 1", sg_number=1, labels=labels,
        comment=(
            "L2 SCHEMATIC OVERLAY: dummy sites trace computed preimage curves n^-1(p) "
            "of the Whitehead hopfion (two fibres link once). Not chemistry. "
            "Open together with hopfion_crystal_BCC_Im3m.cif in VESTA."
        ),
        u_iso=0.02,
    )
