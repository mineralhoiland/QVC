r"""S2 — Static 3D fields on a grid containing the fat torus.

Fields
------
n(x)            Whitehead hopfion texture R³ → S² (stereographic Hopf map),
                smoothly windowed to the vacuum n = (0,0,−1) beyond R_cut so
                the FFT Coulomb-gauge inversion is periodic (fix_v2 recipe).
ρ_H = A·B/16π²  Hopf charge density with B_i = ½ε_ijk F_jk,
                F_ij = n·(∂_i n × ∂_j n), A from ∇×A = B in Coulomb gauge.
ε_FN            Faddeev–Niemi energy density c₂|∇n|² + c₄ F_ij F_ij /4.
ρ_Cas           toroidal Casimir/ZPF envelope (Gaussian ring, same shape as
                the GPE ``gauss_ring``) scaled to E_vac.
ρ_s, θ          superfluid density and phase θ = Q_H φ on the ring (2D slice
                extruded along the torus tube) — 2D CS-flux diagnostic → 1/N.

Two distinct topological numbers are always reported side by side:

* **Whitehead Q** on S³ (this module)               → 1 for the Q=1 hopfion
* **2D CS-flux / lens-space lock Q_H = 1/N** (phase winding of θ) → 0.2 at N=5

They are not the same integral and the workbench never conflates them.

Normalisation note: the TEVC_structures cube header Q≈3.81 used 1/4π² on
A·B. With F_ij = n·(∂_i n×∂_j n) (area form, ∫_{S²} = 4π) the Whitehead
integral is Q = (1/16π²)∫A·B d³x — the same number divided by 4 → 0.95 at
80³ — plus a discretisation error that this module tracks by grid refinement.

Units: nm for coordinates (grid spans ±L·R_core), energies in meV via E_vac.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

# --------------------------------------------------------------------------
# Grid
# --------------------------------------------------------------------------


@dataclass
class Grid3D:
    x: np.ndarray  # 1D axis (nm), endpoint=False (periodic for FFT)
    X: np.ndarray
    Y: np.ndarray
    Z: np.ndarray
    dx: float
    L: float  # half box (nm)

    @property
    def n(self) -> int:
        return int(self.x.size)

    @property
    def origin(self) -> np.ndarray:
        return np.array([self.x[0], self.x[0], self.x[0]])


def make_grid(n_grid: int = 64, L_nm: float = 40.0) -> Grid3D:
    x = np.linspace(-L_nm, L_nm, n_grid, endpoint=False)
    dx = float(x[1] - x[0])
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    return Grid3D(x=x, X=X, Y=Y, Z=Z, dx=dx, L=L_nm)


# --------------------------------------------------------------------------
# Whitehead hopfion texture (windowed)
# --------------------------------------------------------------------------


def hopfion_texture(
    grid: Grid3D,
    *,
    core_nm: float = 10.0,
    R_cut_over_core: float = 2.5,
    window_width_over_core: float = 0.4,
    winding: int = 1,
) -> np.ndarray:
    """Stereographic Hopf map with a smooth window to the vacuum.

    ``core_nm`` is the length scale where |r| = 1 in the Hopf map (the
    unit-circle preimage of n = (0,0,−1)?? — no: the circle x²+y²=1, z=0 maps
    to n = (0,0,+1); the z-axis and infinity map to (0,0,−1)). ``winding``
    multiplies the fibre phase (Q = winding² for the stereographic map).

    Returns n with shape (N,N,N,3), |n| = 1 to machine precision.
    """
    X, Y, Z = grid.X / core_nm, grid.Y / core_nm, grid.Z / core_nm
    r2 = X * X + Y * Y + Z * Z
    r = np.sqrt(r2)
    den = 1.0 + r2
    w0, w1, w2, w3 = 2 * X / den, 2 * Y / den, 2 * Z / den, (1 - r2) / den
    z1 = w0 + 1j * w1
    z2 = w2 + 1j * w3
    if winding != 1:
        z1 = np.abs(z1) * np.exp(1j * winding * np.angle(z1))
    nx = 2.0 * np.real(z1 * np.conj(z2))
    ny = 2.0 * np.imag(z1 * np.conj(z2))
    nz = np.abs(z1) ** 2 - np.abs(z2) ** 2

    R_cut = R_cut_over_core
    delta = window_width_over_core
    t = 0.5 * (1.0 + np.tanh((r - R_cut) / delta))
    nx_b = (1 - t) * nx
    ny_b = (1 - t) * ny
    nz_b = (1 - t) * nz - t
    norm = np.sqrt(nx_b**2 + ny_b**2 + nz_b**2)
    norm = np.maximum(norm, 1e-15)
    n = np.stack([nx_b / norm, ny_b / norm, nz_b / norm], axis=-1)
    return n


def _pgrad(f: np.ndarray, dx: float, axis: int) -> np.ndarray:
    """Periodic central difference."""
    return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2.0 * dx)


def berry_curvature_B(n: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray]:
    """Emergent field B_i = ½ ε_ijk F_jk, F_ij = n·(∂_i n × ∂_j n).

    Returns (B (N,N,N,3), dn (3,N,N,N,3) gradients ∂_i n_a).
    """
    dn = np.stack([_pgrad(n, dx, ax) for ax in range(3)], axis=0)  # (i, ..., a)

    def F(i: int, j: int) -> np.ndarray:
        cross = np.cross(dn[i], dn[j])
        return np.einsum("...a,...a->...", n, cross)

    Bx = F(1, 2)
    By = F(2, 0)
    Bz = F(0, 1)
    return np.stack([Bx, By, Bz], axis=-1), dn


def coulomb_gauge_A(B: np.ndarray, dx: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Solve ∇×A = B, ∇·A = 0 by FFT after projecting B divergence-free.

    Returns (A, B_clean, max|∇·B_raw|).
    """
    N = B.shape[0]
    k = 2.0 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY, KZ = np.meshgrid(k, k, k, indexing="ij")
    K2 = KX**2 + KY**2 + KZ**2
    K2s = np.where(K2 > 0, K2, 1.0)
    Bh = [np.fft.fftn(B[..., i]) for i in range(3)]
    kdotB = KX * Bh[0] + KY * Bh[1] + KZ * Bh[2]
    Bc = [Bh[0] - KX * kdotB / K2s, Bh[1] - KY * kdotB / K2s, Bh[2] - KZ * kdotB / K2s]
    for arr in Bc:
        arr[0, 0, 0] = 0.0
    cx = KY * Bc[2] - KZ * Bc[1]
    cy = KZ * Bc[0] - KX * Bc[2]
    cz = KX * Bc[1] - KY * Bc[0]
    Ah = [-1j * cx / K2s, -1j * cy / K2s, -1j * cz / K2s]
    for arr in Ah:
        arr[0, 0, 0] = 0.0
    A = np.stack([np.real(np.fft.ifftn(a)) for a in Ah], axis=-1)
    B_clean = np.stack([np.real(np.fft.ifftn(b)) for b in Bc], axis=-1)
    divB = sum(_pgrad(B[..., i], dx, i) for i in range(3))
    return A, B_clean, float(np.max(np.abs(divB)))


def whitehead_charge(n: np.ndarray, dx: float) -> dict[str, Any]:
    """Whitehead Q = (1/16π²) ∫ A·B d³x in Coulomb gauge, with diagnostics."""
    B, dn = berry_curvature_B(n, dx)
    A, Bc, divB_raw = coulomb_gauge_A(B, dx)
    AdotB = np.einsum("...i,...i->...", A, Bc)
    dV = dx**3
    raw = float(np.sum(AdotB) * dV)
    Q = raw / (16.0 * math.pi**2)
    norm_dev = float(np.max(np.abs(np.linalg.norm(n, axis=-1) - 1.0)))
    return {
        "Q_whitehead": Q,
        "Q_if_4pi2_normalisation": raw / (4.0 * math.pi**2),
        "raw_integral_AdotB": raw,
        "rho_H": AdotB / (16.0 * math.pi**2),
        "A": A,
        "B": Bc,
        "dn": dn,
        "max_divB_raw": divB_raw,
        "norm_deviation": norm_dev,
        "nz_at_boundary": float(n[0, n.shape[1] // 2, n.shape[2] // 2, 2]),
        "definition": "Whitehead integral on S^3 (compactified R^3); distinct from Q_H=1/N lens-space lock",
    }


WARD_VK_CONSTANT = (3.0 / 16.0) ** (3.0 / 8.0)  # ≈ 0.534
WARD_Q1_MINIMUM = 1.22  # Battye–Sutcliffe / Ward relaxed Q=1 hopfion energy


def faddeev_niemi_energy(
    n: np.ndarray,
    dx: float,
    *,
    dn: np.ndarray | None = None,
    length_unit_nm: float = 1.0,
) -> dict[str, Any]:
    """Faddeev–Niemi energy in Ward's normalisation.

        E_W = (1/32π²) ∫ ( ∂_i n·∂_i n + ½ F_ij F_ij ) d³x   (x in units of the core)

    so that the Vakulenko–Kapitansky bound reads E_W ≥ (3/16)^{3/8} |Q|^{3/4}
    ≈ 0.534 |Q|^{3/4} and the relaxed Q=1 hopfion sits at E_W ≈ 1.22.
    ``length_unit_nm`` converts the grid (nm) to core units.
    """
    if dn is None:
        dn = np.stack([_pgrad(n, dx, ax) for ax in range(3)], axis=0)
    grad2 = np.einsum("i...a,i...a->...", dn, dn)  # 1/nm²
    B, _ = berry_curvature_B(n, dx)
    F2 = 2.0 * np.einsum("...i,...i->...", B, B)  # F_ij F_ij = 2 B·B, 1/nm⁴
    s = length_unit_nm
    eps_density = (grad2 * s**2 + 0.5 * F2 * s**4) / (32.0 * math.pi**2)  # per core³
    dV_core = (dx / s) ** 3
    E2 = float(np.sum(grad2 * s**2) * dV_core) / (32.0 * math.pi**2)
    E4 = float(np.sum(0.5 * F2 * s**4) * dV_core) / (32.0 * math.pi**2)
    return {
        "energy_density": eps_density,
        "E_total": E2 + E4,
        "E_quadratic": E2,
        "E_quartic": E4,
        "VK_constant": WARD_VK_CONSTANT,
        "relaxed_Q1_reference": WARD_Q1_MINIMUM,
        "normalisation": "Ward: E=(1/32pi^2) int(|dn|^2 + F^2/2)",
        "quadratic_grad2_field": grad2,
    }


def vk_bound_check(E_total: float, Q: float, VK_C: float = WARD_VK_CONSTANT) -> dict[str, Any]:
    bound = VK_C * abs(Q) ** 0.75
    return {"E": E_total, "VK_bound": bound, "ratio": E_total / bound if bound else float("inf"), "satisfied": E_total >= bound}


# --------------------------------------------------------------------------
# Torus fields: Casimir envelope, superfluid density / phase
# --------------------------------------------------------------------------


def torus_coordinates(grid: Grid3D, R_nm: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(ρ_tube, φ_toroidal, χ_poloidal) about a ring of major radius R in the xy-plane."""
    rho_xy = np.sqrt(grid.X**2 + grid.Y**2)
    phi = np.arctan2(grid.Y, grid.X)
    d = rho_xy - R_nm
    rho_tube = np.sqrt(d * d + grid.Z**2)
    chi = np.arctan2(grid.Z, d)
    return rho_tube, phi, chi


def casimir_ring_envelope(grid: Grid3D, *, R_nm: float, r_nm: float, E_vac_meV: float) -> np.ndarray:
    """ZPF/Casimir envelope −E_vac·exp(−ρ_tube²/2r²) on the fat torus (meV)."""
    rho_tube, _, _ = torus_coordinates(grid, R_nm)
    return -E_vac_meV * np.exp(-(rho_tube**2) / (2.0 * r_nm**2))


def superfluid_torus_fields(
    grid: Grid3D, *, R_nm: float, r_nm: float, Q_H: float, core_depletion: float = 0.10
) -> dict[str, np.ndarray]:
    """ρ_s (core-depleted ring) and phase θ = Q_H φ — 3D extrusion of the GPE slice."""
    rho_tube, phi, _ = torus_coordinates(grid, R_nm)
    ring = np.exp(-(rho_tube**2) / (2.0 * r_nm**2))
    rho_s = (1.0 - core_depletion * ring) ** 2
    theta = Q_H * phi
    return {"rho_s": rho_s, "theta": theta, "ring": ring}


def cs_flux_2d(theta_slice: np.ndarray, x: np.ndarray, R_nm: float) -> dict[str, float]:
    """2D Chern–Simons-flux diagnostic: (1/2π)∮ ∇θ·dl on the ring → Q_H.

    A fractional vortex θ = Q_H φ is not single valued: going once around the
    ring accumulates 2π Q_H and the remainder is returned at a single branch
    cut (the lens-space identification). We report the smooth accumulated
    winding with the cut removed, and the cut jump itself; they must agree.
    This is the lens-space Q_H lock and is *not* the Whitehead integral.
    """
    n_pts = 720
    ang = np.linspace(0, 2 * np.pi, n_pts, endpoint=False)
    px = R_nm * np.cos(ang)
    py = R_nm * np.sin(ang)
    ix = np.interp(px, x, np.arange(x.size))
    iy = np.interp(py, x, np.arange(x.size))
    i0 = np.clip(np.round(ix).astype(int), 0, x.size - 1)
    j0 = np.clip(np.round(iy).astype(int), 0, x.size - 1)
    th = theta_slice[i0, j0]
    dth = np.diff(np.concatenate([th, th[:1]]))
    cut = int(np.argmax(np.abs(dth)))
    smooth = np.delete(dth, cut)
    return {
        "Q_H_flux": float(np.sum(smooth) / (2 * np.pi)),
        "Q_H_from_cut": float(-dth[cut] / (2 * np.pi)),
        "n_samples": float(n_pts),
    }


# --------------------------------------------------------------------------
# Full static-field bundle + acceptance
# --------------------------------------------------------------------------


def compute_static_fields(
    *,
    n_grid: int = 64,
    R_nm: float = 50.0,
    r_nm: float = 8.0,
    core_nm: float | None = None,
    box_over_core: float = 4.0,
    E_vac_meV: float = 0.1287,
    Q_H: float = 0.2,
    q_tol: float = 0.08,
) -> dict[str, Any]:
    """All S2 fields plus acceptance gates.

    Two grids are used deliberately:

    * ``grid`` (hopfion): the Hopf-map core scale is set to the torus major
      radius R (so the n=(0,0,1) fibre is the torus ring) and the box is
      ±``box_over_core``·core so the vacuum window fits. Q depends on dx/core.
    * ``grid_torus`` (Casimir / superfluid): box ±(R + 4r) so the 8 nm tube
      is resolved. The same n-field is *not* resampled there.
    """
    core = R_nm if core_nm is None else core_nm
    grid = make_grid(n_grid, box_over_core * core)
    n = hopfion_texture(grid, core_nm=core)
    wh = whitehead_charge(n, grid.dx)
    fn = faddeev_niemi_energy(n, grid.dx, dn=wh["dn"], length_unit_nm=core)
    vk = vk_bound_check(fn["E_total"], wh["Q_whitehead"], fn["VK_constant"])

    grid_t = make_grid(n_grid, R_nm + 4.0 * r_nm)
    cas = casimir_ring_envelope(grid_t, R_nm=R_nm, r_nm=r_nm, E_vac_meV=E_vac_meV)
    sf = superfluid_torus_fields(grid_t, R_nm=R_nm, r_nm=r_nm, Q_H=Q_H)
    mid = n_grid // 2
    flux = cs_flux_2d(sf["theta"][:, :, mid], grid_t.x, R_nm)

    acceptance = {
        "unit_norm": wh["norm_deviation"] < 1e-8,
        "whitehead_Q_is_1": abs(wh["Q_whitehead"] - 1.0) < q_tol,
        "old_4pi2_header_rejected": abs(wh["Q_if_4pi2_normalisation"] - 1.0) > q_tol,
        "cs_flux_is_1_over_N": abs(flux["Q_H_flux"] - Q_H) < 0.02 and abs(flux["Q_H_from_cut"] - Q_H) < 0.02,
        "vk_bound": vk["satisfied"],
        "boundary_is_vacuum": abs(wh["nz_at_boundary"] + 1.0) < 1e-3,
    }
    acceptance["all_pass"] = all(acceptance.values())
    return {
        "layer": "L2 geometry (n-field, FN energy) + L1 locks (E_vac envelope, Q_H)",
        "grid": grid,
        "grid_torus": grid_t,
        "n": n,
        "rho_H": wh["rho_H"],
        "A": wh["A"],
        "B": wh["B"],
        "energy_density": fn["energy_density"],
        "casimir_envelope_meV": cas,
        "rho_s": sf["rho_s"],
        "theta": sf["theta"],
        "Q_whitehead": wh["Q_whitehead"],
        "Q_if_4pi2_normalisation": wh["Q_if_4pi2_normalisation"],
        "Q_H_flux_2d": flux["Q_H_flux"],
        "Q_H_from_cut": flux["Q_H_from_cut"],
        "Q_H_lock": Q_H,
        "E_FN_ward": fn["E_total"],
        "E_FN_quadratic": fn["E_quadratic"],
        "E_FN_quartic": fn["E_quartic"],
        "E_FN_relaxed_reference": fn["relaxed_Q1_reference"],
        "VK": vk,
        "max_divB_raw": wh["max_divB_raw"],
        "R_nm": R_nm,
        "r_nm": r_nm,
        "core_nm": core,
        "acceptance": acceptance,
        "labels": {
            "Q_whitehead": "Whitehead integral (1/16π²)∫A·B on S³ — expect 1 (unrelaxed stereographic ansatz)",
            "Q_H_flux_2d": "2D CS-flux / lens-space lock (1/2π)∮∇θ with cut removed — expect 1/N",
            "E_FN_ward": "Ward-normalised FN energy; relaxed Q=1 minimum ≈ 1.22; VK bound 0.534|Q|^{3/4}",
        },
    }


def whitehead_convergence(grids=(32, 48, 64), core_nm: float = 50.0, box_over_core: float = 4.0) -> dict[str, Any]:
    rows = []
    for ng in grids:
        g = make_grid(ng, box_over_core * core_nm)
        n = hopfion_texture(g, core_nm=core_nm)
        wh = whitehead_charge(n, g.dx)
        rows.append({"n_grid": ng, "dx_nm": g.dx, "Q": wh["Q_whitehead"]})
    if len(rows) >= 2:
        h1, h2 = rows[-2]["dx_nm"], rows[-1]["dx_nm"]
        Q1, Q2 = rows[-2]["Q"], rows[-1]["Q"]
        rich = (Q2 * h1**2 - Q1 * h2**2) / (h1**2 - h2**2)
    else:
        rich = rows[-1]["Q"]
    return {"rows": rows, "richardson_Q": rich}


# --------------------------------------------------------------------------
# Preimage curves (S6 payoff; computed here because they need n)
# --------------------------------------------------------------------------


def preimage_curve(n: np.ndarray, grid: Grid3D, p: np.ndarray, *, tol_rad: float = 0.12, max_pts: int = 600) -> np.ndarray:
    """Points x with angle(n(x), p) < tol, ordered around the fibre.

    Hopf fibres are circles, so ordering by PCA-plane angle is exact.
    """
    p = np.asarray(p, float)
    p = p / np.linalg.norm(p)
    cosang = np.clip(np.einsum("...a,a->...", n, p), -1, 1)
    mask = np.arccos(cosang) < tol_rad
    pts = np.stack([grid.X[mask], grid.Y[mask], grid.Z[mask]], axis=1)
    if pts.shape[0] < 8:
        return pts
    c = pts.mean(axis=0)
    u, s, vt = np.linalg.svd(pts - c, full_matrices=False)
    e1, e2 = vt[0], vt[1]
    ang = np.arctan2((pts - c) @ e2, (pts - c) @ e1)
    order = np.argsort(ang)
    pts = pts[order]
    if pts.shape[0] > max_pts:
        idx = np.linspace(0, pts.shape[0] - 1, max_pts).astype(int)
        pts = pts[idx]
    return pts


def gauss_linking_number(c1: np.ndarray, c2: np.ndarray) -> float:
    """Gauss double integral for two closed polylines (should be ±1 for Hopf fibres)."""
    if len(c1) < 3 or len(c2) < 3:
        return float("nan")
    d1 = np.roll(c1, -1, axis=0) - c1
    d2 = np.roll(c2, -1, axis=0) - c2
    m1 = 0.5 * (c1 + np.roll(c1, -1, axis=0))
    m2 = 0.5 * (c2 + np.roll(c2, -1, axis=0))
    r = m1[:, None, :] - m2[None, :, :]
    rn = np.linalg.norm(r, axis=-1) ** 3
    rn = np.where(rn < 1e-12, np.inf, rn)
    cross = np.cross(d1[:, None, :], d2[None, :, :])
    val = np.einsum("ijk,ijk->ij", cross, r) / rn
    return float(np.sum(val) / (4.0 * math.pi))


def preimage_pair(n: np.ndarray, grid: Grid3D, *, tol_rad: float = 0.12) -> dict[str, Any]:
    """Two fibres n⁻¹(p₁), n⁻¹(p₂) and their Gauss linking number."""
    p1 = np.array([0.0, 0.0, 1.0])
    p2 = np.array([1.0, 0.0, 0.0])
    c1 = preimage_curve(n, grid, p1, tol_rad=tol_rad)
    c2 = preimage_curve(n, grid, p2, tol_rad=tol_rad)
    lk = gauss_linking_number(c1, c2)
    return {"curve_1": c1, "curve_2": c2, "p1": p1, "p2": p2, "linking_number": lk, "expected": 1.0}


# --------------------------------------------------------------------------
# Export: CUBE (VESTA) + VTK (ParaView) from the same arrays
# --------------------------------------------------------------------------


def export_static_fields(fields: dict[str, Any], outdir: str | Path, *, stem: str = "hopfion") -> dict[str, str]:
    """Write CUBE + VTK image data (+ PVD manifest) from ``compute_static_fields``."""
    from qvc.viz.structure.vesta_export import write_cube

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    g: Grid3D = fields["grid"]
    gt: Grid3D = fields["grid_torus"]
    files: dict[str, str] = {}
    Q = fields["Q_whitehead"]
    header = f"Whitehead Q={Q:.4f} (1/16pi^2 A.B); 2D CS-flux Q_H={fields['Q_H_flux_2d']:.3f}; grid {g.n}^3, dx={g.dx:.3f} nm"
    files["cube_rho_H"] = str(
        write_cube(outdir / f"{stem}_hopf_density.cube", fields["rho_H"], origin_nm=g.origin, spacing_nm=g.dx,
                   title="Hopfion topological charge density rho_H = A.B/16pi^2 [L2 geometry]", subtitle=header)
    )
    files["cube_energy"] = str(
        write_cube(outdir / f"{stem}_FN_energy_density.cube", fields["energy_density"], origin_nm=g.origin, spacing_nm=g.dx,
                   title="Faddeev-Niemi energy density (Ward norm) [L2 geometry]",
                   subtitle=f"E_W={fields['E_FN_ward']:.3f} VK_bound={fields['VK']['VK_bound']:.3f} relaxed_ref=1.22")
    )
    files["cube_nz"] = str(
        write_cube(outdir / f"{stem}_nz.cube", fields["n"][..., 2], origin_nm=g.origin, spacing_nm=g.dx,
                   title="Hopfion n_z component [L2 geometry]", subtitle=header)
    )
    files["cube_casimir"] = str(
        write_cube(outdir / f"{stem}_casimir_envelope.cube", fields["casimir_envelope_meV"], origin_nm=gt.origin, spacing_nm=gt.dx,
                   title="Toroidal Casimir/ZPF envelope -E_vac exp(-rho^2/2r^2) [meV, L1 lock E_vac]",
                   subtitle=f"R={fields['R_nm']:.1f} nm r={fields['r_nm']:.1f} nm E_vac=0.1287 meV (torus grid)")
    )
    files["vti_hopfion"] = str(write_vtk_image(outdir / f"{stem}_fields.vti", g, {
        "n": fields["n"],
        "rho_H": fields["rho_H"],
        "energy_density": fields["energy_density"],
        "A": fields["A"],
        "B": fields["B"],
    }))
    files["vti_torus"] = str(write_vtk_image(outdir / f"{stem}_torus.vti", gt, {
        "casimir_meV": fields["casimir_envelope_meV"],
        "rho_s": fields["rho_s"],
        "theta": fields["theta"],
    }))
    return files


def write_vtk_image(path: str | Path, grid: Grid3D, arrays: dict[str, np.ndarray]) -> Path:
    """Write a .vti ImageData (pyvista if available, else ASCII XML fallback)."""
    path = Path(path)
    try:
        import pyvista as pv

        img = pv.ImageData(dimensions=(grid.n, grid.n, grid.n), spacing=(grid.dx,) * 3, origin=tuple(grid.origin))
        for name, arr in arrays.items():
            if arr.ndim == 4:
                img.point_data[name] = arr.reshape(-1, arr.shape[-1], order="F")
            else:
                img.point_data[name] = arr.ravel(order="F")
        img.save(str(path))
        return path
    except Exception:
        pass
    n = grid.n
    with path.open("w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <ImageData WholeExtent="0 {n-1} 0 {n-1} 0 {n-1}" Origin="{grid.origin[0]} {grid.origin[1]} {grid.origin[2]}" Spacing="{grid.dx} {grid.dx} {grid.dx}">\n')
        f.write(f'    <Piece Extent="0 {n-1} 0 {n-1} 0 {n-1}">\n      <PointData>\n')
        for name, arr in arrays.items():
            ncomp = arr.shape[-1] if arr.ndim == 4 else 1
            flat = arr.reshape(-1, ncomp, order="F") if arr.ndim == 4 else arr.ravel(order="F")[:, None]
            f.write(f'        <DataArray type="Float64" Name="{name}" NumberOfComponents="{ncomp}" format="ascii">\n')
            np.savetxt(f, flat, fmt="%.8e")
            f.write("        </DataArray>\n")
        f.write("      </PointData>\n    </Piece>\n  </ImageData>\n</VTKFile>\n")
    return path


def write_pvd(path: str | Path, files: list[tuple[float, str]]) -> Path:
    """ParaView collection manifest for time series."""
    path = Path(path)
    with path.open("w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n  <Collection>\n')
        for t, fn in files:
            f.write(f'    <DataSet timestep="{t}" group="" part="0" file="{Path(fn).name}"/>\n')
        f.write("  </Collection>\n</VTKFile>\n")
    return path
