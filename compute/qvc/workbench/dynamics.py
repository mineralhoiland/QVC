"""S5 — Condensate mechanical phenomena (L1).

Single GPE truth: :mod:`qvc.tqgl` (``FrozenTQGL`` → ``run_coupled_gpe_dce``
and ``run_dce_fold_slaved``). The coarser ``qvc/dynamics/tqgl_solver.py`` demo
is *not* used here.

What is shown
-------------
* 2D truncated GPE replay |ψ|²(x,y,t) on the ring (r,R) = (8,50) nm with the
  locked V_ZPF = −E_vac (A/A₀) and phonon drive — **2D evolution**.
* A **3D torus extrusion** of each frame: the in-plane radial profile around
  the ring is wrapped onto the tube surface — **3D display only**, honestly
  labelled; the IC is a 2D fractional vortex, not a 3D Hopf map.
* Time graphs: live Δ(ψ) (no time argument, from ψ), g⁽¹⁾, A(t)/A_c^fold,
  Lorentzian, and the slow 80 ps fold-slaved DCE (resonant / detuning).

Caption discipline (gap_reduction.py forbids stacking): the drive channel can
**raise** the live gap (USAGE: 0.600 → 0.650 meV), while the static fold
**lowers** Δ to Δ* = 0.2324 meV. They are parallel channels, never compounded.

True 3D GPE / Faddeev–Niemi dynamics is deferred (HPC; Architecture §12).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from qvc.tqgl.coupled import run_coupled_gpe_dce
from qvc.tqgl.dce import DCEConfig, run_dce_fold_slaved
from qvc.workbench.locks import Layer, get_locks, stamp

CAPTION_CHANNELS = (
    "Live Δ(ψ) under V_ZPF + phonon drive is a GPE *drive* channel and may rise above Δ₀ "
    "(demo: 0.600 → 0.650 meV). The static fold lowers Δ to Δ* = 0.2324 meV. "
    "Parallel channels — never one compounded reduction (gap_reduction.py)."
)
CAPTION_3D = "2D evolution, 3D display: the ring profile is wrapped onto the torus tube; the IC is a 2D fractional vortex, not a 3D Hopf map."


def run_replay(
    *,
    n_grid: int = 128,
    t_end_ps: float = 1.0,
    dt_ps: float = 0.002,
    n_frames: int = 20,
    pump_mode: str = "resonant",
    seed: int = 42,
) -> dict[str, Any]:
    """qvc.tqgl coupled GPE↔DCE with frame capture."""
    L = get_locks()
    n_steps = int(round(t_end_ps / dt_ps))
    every = max(1, n_steps // max(n_frames, 1))
    res = run_coupled_gpe_dce(
        L.tqgl, n_grid=n_grid, t_end_ps=t_end_ps, dt_ps=dt_ps, pump_mode=pump_mode, seed=seed, snapshot_every=every
    )
    t = res["t_ps"]
    res["layer"] = Layer.L1.value
    res["caption"] = stamp(Layer.L1, f"qvc.tqgl 2D truncated GPE on {n_grid}², ring (r,R)=({L.r_nm:.0f},{L.R_nm:.0f}) nm, {t_end_ps} ps")
    res["caption_channels"] = CAPTION_CHANNELS
    res["caption_3d"] = CAPTION_3D
    res["A_c_fold_um_inv"] = L.tqgl.A_c_fold_um_inv
    res["A_over_Ac"] = res["A_um_inv"] / L.tqgl.A_c_fold_um_inv
    res["Delta0_meV"] = L.Delta0_meV
    res["Delta_star_meV"] = L.Delta_star_meV
    res["Delta_c_meV"] = L.tqgl.Delta_c_meV
    res["live_gap_rose_above_Delta0"] = bool(np.max(res["Delta_live_meV"]) > L.Delta0_meV)
    res["live_gap_max_meV"] = float(np.max(res["Delta_live_meV"]))
    return res


def run_slow_dce() -> dict[str, Any]:
    """80 ps fold-slaved DCE, resonant and detuned (slow timescale; Δ from fold branch)."""
    L = get_locks()
    res = run_dce_fold_slaved(L.tqgl, DCEConfig(pump_mode="resonant"))
    det = run_dce_fold_slaved(L.tqgl, DCEConfig(pump_mode="detuning"))
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, "80 ps DCE slaved to the static fold Δ(A); A*/A_c^fold is reported, not imposed (0.52 is not a target)"),
        "resonant": res,
        "detuning": det,
        "A_c_fold_um_inv": L.tqgl.A_c_fold_um_inv,
    }


# --------------------------------------------------------------------------
# 3D torus extrusion (display only)
# --------------------------------------------------------------------------


def extrude_to_torus(
    density: np.ndarray,
    x_um: np.ndarray,
    y_um: np.ndarray,
    *,
    R_um: float,
    r_um: float,
    n_phi: int = 120,
    n_chi: int = 40,
    tube_scale: float = 1.0,
) -> dict[str, np.ndarray]:
    """Wrap the 2D ring neighbourhood onto a torus surface.

    For toroidal angle φ and poloidal angle χ the surface point is
    (R + ρ cosχ)(cosφ, sinφ), ρ sinχ with ρ = r·tube_scale, coloured by the
    in-plane density sampled at radius R + ρ cosχ along direction φ.
    This is a display map of 2D data; nothing is evolved in 3D.
    """
    phi = np.linspace(0, 2 * np.pi, n_phi)
    chi = np.linspace(0, 2 * np.pi, n_chi)
    PHI, CHI = np.meshgrid(phi, chi, indexing="ij")
    rho = r_um * tube_scale
    rad = R_um + rho * np.cos(CHI)
    X = rad * np.cos(PHI)
    Y = rad * np.sin(PHI)
    Z = rho * np.sin(CHI)
    # sample density at in-plane point (rad cosφ, rad sinφ)
    ix = np.interp(X.ravel(), x_um, np.arange(x_um.size))
    iy = np.interp(Y.ravel(), y_um, np.arange(y_um.size))
    i0 = np.clip(np.round(ix).astype(int), 0, x_um.size - 1)
    j0 = np.clip(np.round(iy).astype(int), 0, y_um.size - 1)
    # density is indexed [row=y, col=x] (meshgrid indexing="xy")
    C = density[j0, i0].reshape(X.shape)
    return {"X_um": X, "Y_um": Y, "Z_um": Z, "color": C}


def replay_frames_3d(res: dict[str, Any], **kw) -> list[dict[str, Any]]:
    L = get_locks()
    frames = []
    for snap in res["snapshots"]:
        tor = extrude_to_torus(snap["density"], res["x_um"], res["y_um"], R_um=L.tqgl.R_um, r_um=L.tqgl.r_um, **kw)
        frames.append({"t_ps": snap["t_ps"], **tor, "density_2d": snap["density"], "phase_2d": snap["phase"]})
    return frames


# --------------------------------------------------------------------------
# ParaView export of the time series
# --------------------------------------------------------------------------


def export_replay_paraview(res: dict[str, Any], outdir: str | Path, *, stem: str = "gpe_replay") -> dict[str, Any]:
    """Write one .vts per frame (2D grid + torus surface) and a .pvd manifest."""
    from qvc.viz.fields3d import write_pvd

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    x, y = res["x_um"], res["y_um"]
    files_2d: list[tuple[float, str]] = []
    files_3d: list[tuple[float, str]] = []
    frames = replay_frames_3d(res)
    try:
        import pyvista as pv
    except Exception:  # pragma: no cover
        pv = None
    for k, fr in enumerate(frames):
        t = fr["t_ps"]
        f2 = outdir / f"{stem}_2d_{k:03d}.vts"
        f3 = outdir / f"{stem}_torus_{k:03d}.vts"
        if pv is not None:
            X, Y = np.meshgrid(x, y, indexing="xy")
            g = pv.StructuredGrid(X[..., None], Y[..., None], np.zeros_like(X)[..., None])
            g.point_data["density"] = fr["density_2d"].ravel(order="F")
            g.point_data["phase"] = fr["phase_2d"].ravel(order="F")
            g.save(str(f2))
            s = pv.StructuredGrid(fr["X_um"][..., None], fr["Y_um"][..., None], fr["Z_um"][..., None])
            s.point_data["density"] = fr["color"].ravel(order="F")
            s.save(str(f3))
        else:
            _write_vts_ascii(f2, np.meshgrid(x, y, indexing="xy") + (None,), {"density": fr["density_2d"], "phase": fr["phase_2d"]})
            _write_vts_ascii(f3, (fr["X_um"], fr["Y_um"], fr["Z_um"]), {"density": fr["color"]})
        files_2d.append((t, str(f2)))
        files_3d.append((t, str(f3)))
    pvd2 = write_pvd(outdir / f"{stem}_2d.pvd", files_2d)
    pvd3 = write_pvd(outdir / f"{stem}_torus.pvd", files_3d)
    # time series CSV for external plotting
    csv = outdir / f"{stem}_timeseries.csv"
    np.savetxt(
        csv,
        np.column_stack([res["t_ps"], res["Delta_live_meV"], res["Delta_fold_of_A_meV"], res["g1"], res["A_um_inv"], res["A_over_Ac"]]),
        delimiter=",",
        header="t_ps,Delta_live_meV,Delta_fold_of_A_meV,g1,A_um_inv,A_over_Ac_fold",
        comments="",
    )
    return {"pvd_2d": str(pvd2), "pvd_torus": str(pvd3), "timeseries_csv": str(csv), "n_frames": len(frames), "caption_3d": CAPTION_3D}


def _write_vts_ascii(path: Path, XYZ, arrays: dict[str, np.ndarray]) -> None:
    X, Y, Z = XYZ
    if Z is None:
        Z = np.zeros_like(X)
    nx, ny = X.shape
    with path.open("w") as f:
        f.write('<?xml version="1.0"?>\n<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'  <StructuredGrid WholeExtent="0 {nx-1} 0 {ny-1} 0 0">\n    <Piece Extent="0 {nx-1} 0 {ny-1} 0 0">\n      <PointData>\n')
        for name, arr in arrays.items():
            f.write(f'        <DataArray type="Float64" Name="{name}" format="ascii">\n')
            np.savetxt(f, arr.ravel(order="F")[None, :], fmt="%.8e")
            f.write("        </DataArray>\n")
        f.write('      </PointData>\n      <Points>\n        <DataArray type="Float64" NumberOfComponents="3" format="ascii">\n')
        pts = np.column_stack([X.ravel(order="F"), Y.ravel(order="F"), Z.ravel(order="F")])
        np.savetxt(f, pts, fmt="%.8e")
        f.write("        </DataArray>\n      </Points>\n    </Piece>\n  </StructuredGrid>\n</VTKFile>\n")


def dynamics_report(res: dict[str, Any], slow: dict[str, Any] | None = None) -> dict[str, Any]:
    L = get_locks()
    out = {
        "n_grid": res["n_grid"],
        "t_end_ps": res["t_end_ps"],
        "Delta_live_initial_meV": float(res["Delta_live_meV"][0]),
        "Delta_live_final_meV": res["Delta_live_final_meV"],
        "Delta_live_max_meV": res["live_gap_max_meV"],
        "Delta_live_min_meV": res["Delta_live_min_meV"],
        "live_gap_rose_above_Delta0": res["live_gap_rose_above_Delta0"],
        "Delta_star_static_fold_meV": L.Delta_star_meV,
        "g1_final": res["g1_final"],
        "A_final_over_Ac_fold": res["A_final_over_Ac_fold"],
        "norm_rel_drift": res["norm_rel_drift"],
        "truncation": res["truncation"],
        "caption_channels": CAPTION_CHANNELS,
        "caption_3d": CAPTION_3D,
        "n_frames": len(res["snapshots"]),
    }
    if slow is not None:
        out["dce_slow"] = {
            "resonant_A_star_over_Ac": slow["resonant"]["A_star_over_Ac"],
            "detuning_A_star_over_Ac": slow["detuning"]["A_star_over_Ac"],
            "resonant_Delta_star_dce_meV": slow["resonant"]["Delta_star_dce_meV"],
            "runaway_clipped": slow["resonant"]["runaway_clipped"],
            "note": slow["resonant"]["note"],
        }
    return out
