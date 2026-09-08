"""Driver: run B1 (G_ij) and B2 (Lindblad), write summary.json + figures."""

from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np

from qvc_l1.figures import write_all_figures
from qvc_l1.gij_ring import run_b1
from qvc_l1.lindblad import run_b2
from qvc_l1.locks import L1Locks, build_l1_locks, folklore_values

DEFAULT_OUT = Path("/Users/mineralhoiland/Code/QVCCursor/notebooks/gij_lindblad/output")
MIRROR_OUT = Path("/Users/mineralhoiland/Code/QVCCursor/output/gij_lindblad")


def _jsonify(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, np.ndarray):
        if obj.size > 2500 and obj.ndim >= 2:
            return {
                "shape": list(obj.shape),
                "omitted": True,
                "note": "array stored as figure / kept in-memory only",
            }
        return _jsonify(obj.tolist())
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        x = float(obj)
        return x if math.isfinite(x) else None
    if obj is None or isinstance(obj, str):
        return obj
    return str(obj)


def _slim_b1(b1: dict[str, Any]) -> dict[str, Any]:
    out = dict(b1)
    for key in ("ws_gaussian", "ring_fourier"):
        blk = dict(out.get(key) or {})
        for drop in (
            "wannier_abs_0",
            "fourier_abs_stack",
            "envelope",
            "x_um",
            "y_um",
            "profiles",
        ):
            blk.pop(drop, None)
        if "gpe_quadratic" in blk:
            gq = dict(blk["gpe_quadratic"])
            # keep G (5×5) and spectrum; drop nothing essential
            blk["gpe_quadratic"] = gq
        out[key] = blk
    return out


def _slim_b2(b2: dict[str, Any]) -> dict[str, Any]:
    def slim_run(r: dict[str, Any]) -> dict[str, Any]:
        keep = {
            k: r[k]
            for k in r
            if k
            not in (
                # keep t, Delta, g1 for plots already written; downsample for JSON
            )
        }
        t = np.asarray(r["t_ps"])
        if t.size > 80:
            idx = np.linspace(0, t.size - 1, 80).astype(int)
            keep["t_ps"] = t[idx]
            keep["Delta_live_meV"] = np.asarray(r["Delta_live_meV"])[idx]
            keep["g1"] = np.asarray(r["g1"])[idx]
            keep["norm"] = np.asarray(r["norm"])[idx]
            keep["A_um_inv"] = np.asarray(r["A_um_inv"])[idx]
        return keep

    out = {k: v for k, v in b2.items() if k not in ("runs_1ps", "runs_20ps")}
    out["runs_1ps"] = {name: slim_run(r) for name, r in b2["runs_1ps"].items()}
    out["runs_20ps"] = {name: slim_run(r) for name, r in b2["runs_20ps"].items()}
    return out


def key_numbers(locks: L1Locks, b1: dict[str, Any], b2: dict[str, Any]) -> dict[str, Any]:
    ws = b1["ws_gaussian"]
    ring = b1["ring_fourier"]
    return {
        "r_nm": locks.r_nm,
        "R_nm": locks.R_nm,
        "hbar_cs_meV_um": locks.hbar_cs_meV_um,
        "Delta0_meV": locks.Delta0_meV,
        "E_vac_meV": locks.E_vac_meV,
        "lambda_star": locks.lambda_star,
        "g0_star_meV": locks.g0_star_meV,
        "N_layers": locks.N_layers,
        "ws": {
            "n_grid": ws["n_grid"],
            "eigenvalues": ws["eigenvalues"],
            "lambda_max": ws["lambda_max"],
            "lambda_min": ws["lambda_min"],
            "g0_eff_kernel": ws["g0_eff_kernel"],
            "rel_frobenius_to_democratic": ws["rel_frobenius_to_democratic"],
            "near_democratic": ws["acceptance"]["near_democratic"],
        },
        "ring": {
            "n_grid": ring["n_grid"],
            "eigenvalues": ring["eigenvalues"],
            "lambda_max": ring["lambda_max"],
            "lambda_min": ring["lambda_min"],
            "g0_eff_kernel": ring["g0_eff_kernel"],
            "rel_frobenius_to_democratic": ring["rel_frobenius_to_democratic"],
            "near_democratic": ring["acceptance"]["near_democratic"],
        },
        "lindblad_table_1ps": b2["table_1ps"],
        "lindblad_table_20ps": b2["table_20ps"],
        "lindblad_channel_split_1ps": b2["channel_split_1ps"],
        "coherent_1ps": b2["this_run_coherent_1ps"],
        "published_coherent_baseline": b2["published_coherent_baseline"],
    }


def render_markdown(bundle: dict[str, Any]) -> str:
    kn = bundle["key_numbers"]
    ws, ring = kn["ws"], kn["ring"]
    lines = [
        "# L1 dictionary — B1 $G_{ij}$ and B2 Lindblad GPE",
        "",
        "Computed kernel model + truncated GPE with phenomenological dissipation.",
        "Not torsion, not a linearized TQGL theorem, not Hubbard / CS-in-EOM / H-TOR1 / holography.",
        "Hopfion ≠ 2D phase vortex ≠ CS fluxon. CS/Berry/Skyrme truncated.",
        "",
        "## Locks (Option A*)",
        "",
        f"- $(r,R)=({kn['r_nm']:.0f},{kn['R_nm']:.0f})$ nm, $\\hbar c_s={kn['hbar_cs_meV_um']}$ meV·µm, $\\Delta_0={kn['Delta0_meV']:.3f}$ meV",
        f"- $E_{{\\mathrm{{vac}}}}={kn['E_vac_meV']:.4f}$ meV, $\\lambda^\\star={kn['lambda_star']:.4f}$",
        f"- **BRIDGE algebra:** $g_0^\\star=\\lambda^\\star E_{{\\mathrm{{vac}}}}={kn['g0_star_meV']:.4f}$ meV",
        "- Democratic Spec (algebra): $\\lambda_{\\max}=4g_0$, $\\lambda_{\\min}=-g_0$ (mult 4). **Not inserted into $G_{ij}$.**",
        "",
        "## B1  $G_{ij}$  (computed kernel, K0 units ≠ meV)",
        "",
        "### WS Gaussian (finer grid)",
        "",
        f"- grid {ws['n_grid']}, $g_{{0,\\mathrm{{eff}}}}={ws['g0_eff_kernel']:.6g}$ (kernel)",
        f"- Spec $= {np.array2string(np.asarray(ws['eigenvalues']), precision=5)}",
        f"- $\\lambda_{{\\min}}={ws['lambda_min']:.5g}$, $\\lambda_{{\\max}}={ws['lambda_max']:.5g}$",
        f"- rel Frobenius to democratic at $g_{{0,\\mathrm{{eff}}}}$: **{ws['rel_frobenius_to_democratic']:.4f}**",
        f"- near-democratic (rel Fro < 0.15): {ws['near_democratic']}",
        "",
        "### Ring Fourier / Wannier",
        "",
        f"- grid {ring['n_grid']}, $g_{{0,\\mathrm{{eff}}}}={ring['g0_eff_kernel']:.6g}$ (kernel)",
        f"- Spec $= {np.array2string(np.asarray(ring['eigenvalues']), precision=5)}",
        f"- $\\lambda_{{\\min}}={ring['lambda_min']:.5g}$, $\\lambda_{{\\max}}={ring['lambda_max']:.5g}$",
        f"- rel Frobenius to democratic at $g_{{0,\\mathrm{{eff}}}}$: **{ring['rel_frobenius_to_democratic']:.4f}**",
        f"- near-democratic: {ring['near_democratic']}",
        "",
        "$g_{0,\\mathrm{eff}}$ is **not** forced equal to $g_0^\\star$. Different units, different objects.",
        "",
        "## B2  Lindblad scan",
        "",
        "Coherent 1 ps target (published truncated GPE+DCE): "
        f"$\\Delta_{{\\mathrm{{live}}}}\\approx {kn['published_coherent_baseline']['Delta_live_1ps_meV']:.4f}$ meV, "
        f"$g^{{(1)}}\\approx {kn['published_coherent_baseline']['g1_1ps']:.4f}$.",
        "",
        "### Table, 1 ps (coupled DCE)",
        "",
        "| label | channel | $\\tau_{\\mathrm{dec}}$ (ps) | $\\Delta_{\\mathrm{live}}$ (meV) | $g^{(1)}$ |",
        "|---|---|---:|---:|---:|",
    ]
    for row in kn["lindblad_table_1ps"] + kn["lindblad_channel_split_1ps"]:
        tau = "—" if row["tau_dec_ps"] is None else f"{row['tau_dec_ps']:.0f}"
        dlt = "—" if row["Delta_live_meV"] is None else f"{row['Delta_live_meV']:.4f}"
        g1 = "—" if row["g1"] is None else f"{row['g1']:.4f}"
        lines.append(f"| {row['label']} | {row['channel']} | {tau} | {dlt} | {g1} |")
    lines.extend(
        [
            "",
            "### Table, 20 ps window (V_ZPF frozen, no DCE)",
            "",
            "| label | channel | $\\tau_{\\mathrm{dec}}$ (ps) | $\\Delta_{\\mathrm{live}}$ (meV) | $g^{(1)}$ |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in kn["lindblad_table_20ps"]:
        tau = "—" if row["tau_dec_ps"] is None else f"{row['tau_dec_ps']:.0f}"
        dlt = "—" if row["Delta_live_meV"] is None else f"{row['Delta_live_meV']:.4f}"
        g1 = "—" if row["g1"] is None else f"{row['g1']:.4f}"
        lines.append(f"| {row['label']} | {row['channel']} | {tau} | {dlt} | {g1} |")
    lines.extend(
        [
            "",
            "## Computed vs ansatz",
            "",
            "- **Computed:** K0 overlaps; Spec(G); Δ_live(ψ); g^{(1)}(ψ); Lindblad trajectories.",
            "- **Locked algebra:** democratic (4g0, −g0×4); g0★=λ★ E_vac.",
            "- **Phenomenological:** γ=1/τ_dec; local phase kicks.",
            "- **Not:** J−I insertion, torsion theorem, 92%, prescribed Δ(t).",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def run_all(
    *,
    out_dir: Path | str | None = None,
    make_figures: bool = True,
    ws_n_grid: int = 80,
    ring_n_grid: int = 96,
    n_grid_1ps: int = 128,
    n_grid_20ps: int = 96,
    dt_ps: float = 0.002,
    t_short_ps: float = 1.0,
    t_long_ps: float = 20.0,
    include_channel_split: bool = True,
) -> dict[str, Any]:
    """Execute B1 then B2, save artifacts."""
    out = Path(out_dir) if out_dir is not None else DEFAULT_OUT
    out.mkdir(parents=True, exist_ok=True)
    locks = build_l1_locks()
    b1 = run_b1(locks, ws_n_grid=ws_n_grid, ring_n_grid=ring_n_grid)
    b2 = run_b2(
        locks,
        n_grid_1ps=n_grid_1ps,
        n_grid_20ps=n_grid_20ps,
        dt_ps=dt_ps,
        t_short_ps=t_short_ps,
        t_long_ps=t_long_ps,
        include_channel_split=include_channel_split,
    )
    kn = key_numbers(locks, b1, b2)
    bundle = {
        "suite": "L1 dictionary  B1 G_ij + B2 Lindblad GPE",
        "locks": locks.as_dict(),
        "B1": b1,
        "B2": b2,
        "key_numbers": kn,
        "computed_vs_ansatz": {
            "B1": b1["computed_vs_ansatz"],
            "B2": b2["computed_vs_ansatz"],
        },
        "folklore_not_results": {
            "values": list(folklore_values()),
            "do_not_quote_as_physics": True,
        },
    }
    figs: list[str] = []
    if make_figures:
        figs = write_all_figures(bundle, out)
        bundle["figures"] = figs

    slim = {
        "suite": bundle["suite"],
        "locks": bundle["locks"],
        "B1": _slim_b1(b1),
        "B2": _slim_b2(b2),
        "key_numbers": kn,
        "computed_vs_ansatz": bundle["computed_vs_ansatz"],
        "folklore_not_results": bundle["folklore_not_results"],
        "figures": figs,
    }
    (out / "summary.json").write_text(json.dumps(_jsonify(slim), indent=2))
    (out / "summary.md").write_text(render_markdown(bundle))
    MIRROR_OUT.mkdir(parents=True, exist_ok=True)
    for item in out.iterdir():
        dest = MIRROR_OUT / item.name
        if item.is_dir():
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)
    return bundle


if __name__ == "__main__":
    import argparse
    import os

    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-qvc-l1")

    p = argparse.ArgumentParser(description="L1 G_ij + Lindblad suite")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--quick", action="store_true", help="coarse grids / short windows")
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args()
    kw: dict[str, Any]
    if args.quick:
        kw = dict(
            ws_n_grid=32,
            ring_n_grid=48,
            n_grid_1ps=32,
            n_grid_20ps=32,
            dt_ps=0.004,
            t_short_ps=0.08,
            t_long_ps=0.20,
            include_channel_split=False,
        )
    else:
        kw = {}
    bundle = run_all(out_dir=args.out, make_figures=not args.no_figures, **kw)
    kn = bundle["key_numbers"]
    print("WS eigs", kn["ws"]["eigenvalues"], "relFro", kn["ws"]["rel_frobenius_to_democratic"])
    print("ring eigs", kn["ring"]["eigenvalues"], "relFro", kn["ring"]["rel_frobenius_to_democratic"])
    print("g0_star_meV", kn["g0_star_meV"], "g0_eff_ws", kn["ws"]["g0_eff_kernel"], "g0_eff_ring", kn["ring"]["g0_eff_kernel"])
    print("1ps table")
    for row in kn["lindblad_table_1ps"]:
        print(" ", row["label"], row["Delta_live_meV"], row["g1"])
    print("wrote", args.out)
