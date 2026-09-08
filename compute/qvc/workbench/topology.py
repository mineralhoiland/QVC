"""S6 — Topology and transport (L1 locks + L2 3D payoff).

* Hall staircase σ_xy(B) vs N and hopfion count: step δσ_xy = e²/(2Nh)
  (= e²/10h at N = 5) from :func:`qvc.topology.chern_hopf_hall.hall_staircase`.
* Berry curvature Ω(k) on the moiré BZ from the Bistritzer–MacDonald model
  (Fukui lattice gauge), with the Chern number of the tracked flat band.
* Preimage linking: two fibres n⁻¹(p₁), n⁻¹(p₂) of the Whitehead hopfion
  and their Gauss linking number (→ 1). This is the 3D object the VESTA
  dummy He/Ne sites only hint at.

The workbench names which integer is on screen:
  Whitehead Q (S³ linking, → 1)   vs   Q_H = 1/N lens-space CS/linking lock (→ 0.2).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc.materials.bistritzer_macdonald import BistritzerMacDonald
from qvc.topology.chern_hopf_hall import hall_staircase
from qvc.viz.fields3d import Grid3D, gauss_linking_number, preimage_curve
from qvc.workbench.locks import Layer, get_locks, stamp


def hall_panel(*, N_values: tuple[int, ...] = (3, 4, 5, 6), B_max: float = 6.0, n_B: int = 1200) -> dict[str, Any]:
    L = get_locks()
    B = np.linspace(0.0, B_max, n_B)
    curves = {}
    for N in N_values:
        hs = hall_staircase(B, N_layers=N, Delta_ex_meV=L.Delta0_meV)
        curves[N] = {
            "sigma_xy_e2h": hs["sigma_xy"],
            "sigma_xy_sharp_e2h": hs["sigma_xy_sharp"],
            "step_e2h": hs["step_height"],
            "step_S": hs["step_height_SI"],
            "base_e2h": hs["base_conductance"],
            "nucleation_fields": hs["nucleation_fields"],
        }
    lock_step = curves[L.N]["step_e2h"] if L.N in curves else 1.0 / (2 * L.N)
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, f"Hall staircase: base N e²/h, step e²/(2Nh) per nucleated hopfion → e²/{2*L.N}h at N={L.N}"),
        "B_over_B0": B,
        "by_N": curves,
        "lock_step_e2h": L.hall_step_e2_over_h,
        "lock_step_matches": abs(lock_step - L.hall_step_e2_over_h) < 1e-12,
        "note": "B in units of B₀ = Φ₀/(πR²); thermal smearing k_BT/Δ = 0.02 (model in chern_hopf_hall)",
    }


def berry_panel(
    *,
    theta_deg: float = 1.1,
    n_layers: int = 2,
    k_mesh: int = 16,
    n_shells: int = 3,
    band_index: int = 0,
) -> dict[str, Any]:
    """Berry curvature heat map + Chern number of the tracked flat band (one valley)."""
    bm = BistritzerMacDonald(theta_deg, n_layers=n_layers)
    kx, ky, om = bm.berry_curvature(k_mesh_size=k_mesh, band_index=band_index, n_shells=n_shells)
    C = float(np.sum(om) / (2.0 * np.pi))
    info = bm.get_moire_bz_info()
    return {
        "layer": Layer.L1.value,
        "caption": stamp(
            Layer.L1,
            f"BM Berry curvature Ω(k), θ={theta_deg}°, {n_layers} layers, {k_mesh}² mesh, {n_shells} shells — C = {C:+.2f} (one valley)",
        ),
        "kx": kx,
        "ky": ky,
        "omega": om,
        "chern": C,
        "chern_rounded": int(np.round(C)),
        "theta_deg": theta_deg,
        "n_layers": n_layers,
        "bz_info": info,
        "note": "Fukui lattice gauge on a coarse mesh; |C|=1 per valley per spin for the MATBG flat band",
    }


def preimage_panel(n: np.ndarray, grid: Grid3D, *, n_fibres: int = 4, tol_rad: float = 0.12) -> dict[str, Any]:
    """Several fibres n⁻¹(p_k) and their pairwise Gauss linking numbers."""
    L = get_locks()
    pts = [np.array([0.0, 0.0, 1.0]), np.array([1.0, 0.0, 0.0]), np.array([0.0, 1.0, 0.0]), np.array([-1.0, 0.0, 0.0]),
           np.array([0.0, -1.0, 0.0]), np.array([0.0, 0.0, -1.0])][:n_fibres]
    curves = [preimage_curve(n, grid, p, tol_rad=tol_rad) for p in pts]
    lk = np.full((len(curves), len(curves)), np.nan)
    for i in range(len(curves)):
        for j in range(i + 1, len(curves)):
            lk[i, j] = lk[j, i] = gauss_linking_number(curves[i], curves[j])
    off = lk[np.triu_indices(len(curves), k=1)]
    off = off[np.isfinite(off)]
    return {
        "layer": Layer.L2.value,
        "caption": stamp(
            Layer.L2,
            f"Preimage fibres n⁻¹(p) of the Whitehead hopfion; pairwise Gauss linking = {np.nanmean(np.abs(off)) if off.size else float('nan'):.2f} (expect 1). "
            f"Distinct from the L1 lens-space lock Q_H = 1/N = {L.Q_H:.2f}",
        ),
        "points_S2": pts,
        "curves_nm": curves,
        "linking_matrix": lk,
        "mean_abs_linking": float(np.nanmean(np.abs(off))) if off.size else float("nan"),
        "expected_linking": 1.0,
        "Q_H_lock": L.Q_H,
        "distinction": "Whitehead/Gauss linking on S³ (integer) ≠ Q_H = 1/N (lens-space CS/linking form)",
    }


def topology_report(fields: dict[str, Any] | None = None, *, quick: bool = True) -> dict[str, Any]:
    hall = hall_panel()
    berry = berry_panel(k_mesh=10 if quick else 16, n_shells=2 if quick else 3)
    out: dict[str, Any] = {
        "hall": {
            "lock_step_e2h": hall["lock_step_e2h"],
            "lock_step_matches": hall["lock_step_matches"],
            "steps_by_N": {int(N): v["step_e2h"] for N, v in hall["by_N"].items()},
        },
        "berry": {"chern": berry["chern"], "chern_rounded": berry["chern_rounded"], "n_layers": berry["n_layers"]},
    }
    if fields is not None:
        pp = preimage_panel(fields["n"], fields["grid"])
        out["preimage"] = {
            "mean_abs_linking": pp["mean_abs_linking"],
            "expected": 1.0,
            "linking_matrix": np.nan_to_num(pp["linking_matrix"]).tolist(),
            "Q_H_lock": pp["Q_H_lock"],
            "distinction": pp["distinction"],
        }
    return out
