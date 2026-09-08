#!/usr/bin/env python
"""Batch: build → solve → export the TEVC / QVC–TDVT crystal suite.

    python scripts/run_crystal_suite.py                    # full run → output/crystal_suite/
    python scripts/run_crystal_suite.py --quick            # small grids, for CI / smoke
    python scripts/run_crystal_suite.py --out /tmp/suite   # custom folder
    python scripts/run_crystal_suite.py --html             # also save a static workbench HTML

Everything the Panel workbench shows is produced here headlessly, so VESTA /
ParaView users never need the server:

  VESTA      TEVC.cif, POSCAR, TEVC_fragment.cif, hopfion_crystal_BCC_Im3m.cif,
             hopfion_preimage_overlay.cif, hopfion_*.cube
  OVITO      tevc.data (LAMMPS)
  ParaView   hopfion_*.vti (static), gpe_replay_*.pvd (+ .vts / .vti frames)
  graphs     *.json reports + PNG (if kaleido is installed) for fold, Spec(G),
             Mexican hat, rates, Hall, Berry, Goldstone, acoustic slice
  locks      workbench_locks.json (frozen numbers + reject list), acceptance.json

Numbers are read from qvc.workbench.locks (S0). If any acceptance gate fails the
exit code is 1 — a rendering that disagrees with the locks is a bug.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _json_default(o):
    if isinstance(o, np.ndarray):
        return o.tolist() if o.size <= 64 else f"<ndarray {o.shape}>"
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    return str(o)


def _dump(path: Path, obj) -> None:
    path.write_text(json.dumps(obj, indent=2, default=_json_default))


def _try_png(fig, path: Path) -> bool:
    try:
        fig.write_image(str(path), scale=1.5)
        return True
    except Exception as exc:  # kaleido missing / no chrome
        print(f"    [png skipped] {path.name}: {type(exc).__name__}")
        return False


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default=str(ROOT / "output" / "crystal_suite"))
    ap.add_argument("--quick", action="store_true", help="small grids (fields 32³, GPE 48², Berry 8×8)")
    ap.add_argument("--fields-grid", type=int, default=None)
    ap.add_argument("--gpe-grid", type=int, default=None)
    ap.add_argument("--gpe-t-ps", type=float, default=None)
    ap.add_argument("--gpe-frames", type=int, default=12)
    ap.add_argument("--berry-mesh", type=int, default=None)
    ap.add_argument("--theta", type=float, default=1.10, help="twist angle (deg)")
    ap.add_argument("--no-png", action="store_true")
    ap.add_argument("--html", action="store_true", help="also save a static Panel HTML of the workbench")
    ap.add_argument("--l2", action="store_true", default=True, help="include the L2 analogue-geometry products (default on)")
    args = ap.parse_args(argv)

    fields_grid = args.fields_grid or (32 if args.quick else 48)
    gpe_grid = args.gpe_grid or (48 if args.quick else 64)
    gpe_t = args.gpe_t_ps or (0.3 if args.quick else 0.5)
    berry_mesh = args.berry_mesh or (8 if args.quick else 12)

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    T0 = time.time()
    report: dict = {"config": vars(args) | {"fields_grid": fields_grid, "gpe_grid": gpe_grid, "gpe_t_ps": gpe_t, "berry_mesh": berry_mesh}}

    # ------------------------------------------------------------------ S0
    from qvc.workbench.locks import DEPRECATED, check_against_phase_two, get_locks

    print("[S0] locks")
    L = get_locks()
    s0 = check_against_phase_two()
    _dump(out / "workbench_locks.json", {"locks": L.as_dict(), "deprecated_reject_list": DEPRECATED, "phase_two_check": s0})
    print(f"    Option A*: E_vac={L.E_vac_meV:.4f} meV, Δ*={L.Delta_star_meV:.4f} meV, g₀*={L.g0_star_meV:.4f}, Q_H={L.Q_H}, phase-two check: {s0['all_pass']}")
    gates = {"S0_locks_match_phase_two": s0["all_pass"]}

    # ------------------------------------------------------------------ S1
    from qvc.viz.structure import build_tevc_geometry, export_tevc, hopfion_bcc_cif, preimage_overlay_cif

    print("[S1] geometry → CIF / POSCAR / LAMMPS")
    geom = build_tevc_geometry(args.theta, L.N, fragment_radius_nm=2.0)
    s1 = export_tevc(out, geom)
    hopfion_bcc_cif(out / "hopfion_crystal_BCC_Im3m.cif")
    moire_nm = geom.lambda_M_A / 10.0
    print(f"    TEVC: θ={args.theta:.2f}°, {L.N} layers, {s1['n_atoms']} atoms, λ_M={moire_nm:.2f} nm → {len(s1['files'])} files")
    report["S1"] = {"n_atoms": s1["n_atoms"], "moire_nm": moire_nm, "files": s1["files"]}

    # ------------------------------------------------------------------ S2
    from qvc.viz.fields3d import compute_static_fields, export_static_fields, whitehead_convergence

    print(f"[S2] hopfion / Casimir fields at {fields_grid}³ → CUBE + VTI")
    f = compute_static_fields(n_grid=fields_grid, R_nm=L.R_nm, r_nm=L.r_nm, E_vac_meV=L.E_vac_meV, Q_H=L.Q_H)
    s2files = export_static_fields(f, out)
    conv = whitehead_convergence(grids=(32, 48, 64), core_nm=L.R_nm)
    Q_rich = conv["richardson_Q"]
    print(
        f"    Whitehead Q={f['Q_whitehead']:.3f} at {fields_grid}³ (Richardson over 32/48/64 → {Q_rich:.3f}; "
        f"old 1/4π² header would say {f['Q_if_4pi2_normalisation']:.2f} — rejected); "
        f"CS-flux Q_H={f['Q_H_flux_2d']:.3f}; E_FN(Ward)={f['E_FN_ward']:.3f}, VK ratio={f['VK']['ratio']:.2f}"
    )
    acc = dict(f["acceptance"])
    acc["whitehead_Q_is_1_at_this_grid"] = acc.pop("whitehead_Q_is_1")
    acc["whitehead_Q_is_1_richardson"] = abs(Q_rich - 1.0) < 0.05
    grid_independent = {k: v for k, v in acc.items() if k not in ("whitehead_Q_is_1_at_this_grid", "all_pass")}
    acc["all_pass"] = all(grid_independent.values())
    report["S2"] = {k: f[k] for k in ("Q_whitehead", "Q_if_4pi2_normalisation", "Q_H_flux_2d", "Q_H_from_cut", "E_FN_ward", "VK", "labels")} | {
        "acceptance": acc, "whitehead_convergence": conv, "files": s2files}
    gates["S2_fields"] = acc["all_pass"]

    # ------------------------------------------------------------------ S6 (needs n-field)
    from qvc.workbench import analogue, barriers, catalyzon, dynamics, topology

    print("[S6] topology: Hall staircase, Berry curvature, preimage linking")
    hall = topology.hall_panel()
    berry = topology.berry_panel(theta_deg=args.theta, k_mesh=berry_mesh, n_shells=2)
    pre = topology.preimage_panel(f["n"], f["grid"], n_fibres=4)
    preimage_overlay_cif(out / "hopfion_preimage_overlay.cif", pre["curves_nm"][:2])
    print(f"    Hall step at N={L.N}: {hall['by_N'][L.N]['step_e2h']:.3f} e²/h (lock {L.hall_step_e2_over_h:.3f}, match={hall['lock_step_matches']}); "
          f"preimage pairwise linking = {pre['mean_abs_linking']:.2f} (expect 1); Chern(BM, {berry_mesh}²) = {berry['chern']:+.2f}")
    report["S6"] = {
        "hall": {"lock_step_e2h": hall["lock_step_e2h"], "lock_step_matches": hall["lock_step_matches"],
                 "steps_by_N": {N: c["step_e2h"] for N, c in hall["by_N"].items()}, "note": hall["note"]},
        "berry": {k: v for k, v in berry.items() if not isinstance(v, np.ndarray)},
        "preimage": {k: v for k, v in pre.items() if k not in ("curves_nm",)},
    }
    gates["S6_hall_step_matches_lock"] = bool(hall["lock_step_matches"])
    # Gauss linking of grid-traced fibres converges ~ dx/core; 0.15 at 48³ (dx = 8.3 nm), looser on coarse grids
    lk_tol = 0.15 * max(1.0, 48.0 / fields_grid) ** 2
    gates["S6_preimage_linking_is_1"] = abs(pre["mean_abs_linking"] - 1.0) < lk_tol
    report["S6"]["preimage"]["tolerance_used"] = lk_tol

    # ------------------------------------------------------------------ S3
    print("[S3] catalyzon spectrum")
    s3 = catalyzon.catalyzon_report()
    report["S3"] = s3
    dem = s3["democratic"]
    gates["S3_lambda_min_is_minus_g0"] = abs(dem["lambda_min_meV"] + dem["g0_meV"]) < 1e-9 and dem["lambda_min_multiplicity"] == L.N - 1
    print(f"    Spec(G): λ_max={L.lambda_max_meV:.4f}, λ_min={L.lambda_min_meV:.4f} (×{L.lambda_min_multiplicity}); Δ_cat={L.Delta_cat_meV:.4f} meV ({L.catalyzon_reduction_pct:.1f}%, parallel to fold)")

    # ------------------------------------------------------------------ S4
    print("[S4] barriers: fold map, Mexican hat, rates, U_eff(A)")
    s4 = barriers.barriers_report()
    report["S4"] = s4
    fold = barriers.fold_map()
    gates["S4_fold_alpha_c"] = abs(fold["alpha_c"] - L.tqgl.alpha_c) < 1e-6
    print(f"    fold at α_c=δ_c={fold['alpha_c']:.4f}; λ*={L.lambda_star:.4f} → Δ*={L.Delta_star_meV:.4f} meV")

    # ------------------------------------------------------------------ S5
    print(f"[S5] GPE↔DCE replay {gpe_grid}² × {gpe_t} ps → .pvd")
    res = dynamics.run_replay(n_grid=gpe_grid, t_end_ps=gpe_t, n_frames=args.gpe_frames)
    slow = dynamics.run_slow_dce()
    ex = dynamics.export_replay_paraview(res, out)
    s5 = dynamics.dynamics_report(res, slow)
    report["S5"] = s5 | {"export": ex}
    print(f"    {ex['n_frames']} frames → {Path(ex['pvd_torus']).name}, {Path(ex['pvd_2d']).name}; "
          f"live Δ(ψ) end = {s5['Delta_live_final_meV']:.4f} meV (drive channel); static fold Δ* = {L.Delta_star_meV:.4f} meV; norm drift {s5['norm_rel_drift']:.1e}")
    gates["S5_norm_conserved"] = abs(s5["norm_rel_drift"]) < 1e-3

    # ------------------------------------------------------------------ S7 (L2)
    if args.l2:
        print("[S7] L2 analogue geometry (diagnostic only)")
        s7 = analogue.analogue_report(f["rho_H"])
        report["S7_L2"] = s7
        ac, acd = s7["acoustic_vortex_only"], s7["acoustic_with_drain_1p2cs"]
        print(f"    acoustic slice: vortex-only sonic radius {ac['rho_h_nm']:.3g} nm (inside core={ac['horizon_inside_core']} → no horizon); "
              f"with drain 1.2c_s: ρ_h={acd['rho_h_nm']:.3g} nm, T_H={acd['T_H_mK']:.3g} mK ≠ T*={ac['T_star_mK']} mK; "
              f"H-EIN1 TT overlap = {s7['h_ein1']['overlap_plus']:.1e}; H-TOR1 {s7['h_tor1']}")
        gates["S7_TA_is_not_graviton"] = not s7["h_ein1"]["is_graviton_like"]

    # ------------------------------------------------------------------ graphs
    if not args.no_png:
        print("[graphs] PNG via plotly/kaleido")
        from qvc.workbench import app

        figs = {
            "fold_map": app.fig_fold(L.lambda_star),
            "specG": app.fig_specG(L.N, L.g0_star_meV),
            "mexican_hat": app.fig_mexican(0.9),
            "rates": app.fig_rates(),
            "hall_staircase": app.fig_hall(L.N),
            "goldstone_L2": app.fig_goldstone(),
            "acoustic_L2": app.fig_acoustic(1.2),
            "dynamics_timeseries": app.fig_dynamics_timeseries(res, slow),
        }
        n_ok = sum(_try_png(fig, out / f"{name}.png") for name, fig in figs.items())
        print(f"    {n_ok}/{len(figs)} PNGs written")

    if args.html:
        print("[html] static workbench snapshot")
        from qvc.workbench import app

        tpl = app.build_app(fields_grid=fields_grid, gpe_grid=gpe_grid, gpe_t_ps=gpe_t, gpe_frames=args.gpe_frames)
        tpl.save(str(out / "workbench.html"), embed=False)
        print(f"    {out / 'workbench.html'}")

    # ------------------------------------------------------------------ acceptance
    gates["all_pass"] = all(gates.values())
    report["acceptance"] = gates
    report["elapsed_s"] = round(time.time() - T0, 1)
    _dump(out / "acceptance.json", gates)
    _dump(out / "suite_report.json", report)
    print(f"\nacceptance: {gates}\nwrote → {out}  ({report['elapsed_s']} s)")
    return 0 if gates["all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
