#!/usr/bin/env python3
"""Derive bridge, freeze fat-torus working standard, compute TQGL G_ij, plot.

Usage:
  .venv/bin/python scripts/run_bridge_and_gij.py
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qvc.bridge import derive_bridge, frozen_fat_torus_spec
from qvc.gap_fold import fold_curve, solve_fold_branches
from qvc.gij_tqgl_ws import run_tqgl_ws_G
from qvc.params import PARAMS
from qvc.vacuum_freeze import freeze

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 11,
        "legend.fontsize": 9,
        "figure.dpi": 120,
        "savefig.dpi": 200,
        "axes.grid": True,
        "grid.alpha": 0.3,
    }
)


def _fmt(x, d=5):
    if x is None:
        return "—"
    if isinstance(x, float) and x != x:
        return "nan"
    if isinstance(x, bool):
        return str(x)
    return f"{x:.{d}g}"


def make_figures(bridge: dict, gij: dict, curve: dict, out: Path) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    paths = []

    # 1. Fold branches + operating point
    fig, ax = plt.subplots(figsize=(6.2, 4.2), constrained_layout=True)
    a = np.array(curve["alpha_over_ac"])
    phys = np.array([np.nan if v is None else v for v in curve["delta_phys"]])
    unst = np.array([np.nan if v is None else v for v in curve["delta_unst"]])
    ax.plot(a, phys, color="#1f4e79", lw=2, label="physical branch")
    ax.plot(a, unst, color="#9c2a2a", lw=1.5, ls="--", label="unstable branch")
    ax.axvline(1.0, color="#555", ls=":", lw=1, label="fold α=α_c")
    ax.axvline(
        bridge["alpha_over_ac"],
        color="#c47b17",
        ls="-",
        lw=1.4,
        label=f"λ=1 operating α/α_c={bridge['alpha_over_ac']:.3f}",
    )
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$\delta=\Delta/\Delta_0$")
    ax.set_title("Frozen-cavity fold map (fat torus, derived ℓ)")
    ax.set_xlim(0, 1.25)
    ax.set_ylim(0, 1.05)
    ax.legend(loc="lower left")
    p = out / "fold_operating_point.png"
    fig.savefig(p)
    plt.close(fig)
    paths.append(str(p))

    # 2. F(x) fat torus
    from qvc.vacuum_freeze import bessel_mode_sum

    xs = np.linspace(0.2, 2.4, 80)
    Fs = [bessel_mode_sum(float(x)) for x in xs]
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    ax.plot(xs, Fs, color="#1f4e79", lw=2, label=r"$F(x)=\sum K_0(nx)$")
    ax.axvline(bridge["x"], color="#c47b17", ls="--", label=f"fat torus x={bridge['x']:.3f}")
    ax.scatter([bridge["x"]], [bridge["F"]], color="#c47b17", zorder=5)
    ax.axhline(0.596, color="#888", ls=":", label="withdrawn F=0.596")
    ax.set_xlabel(r"$x=2\pi r/R$")
    ax.set_ylabel(r"$F(x)$")
    ax.set_title("Bessel mode sum vs aspect ratio")
    ax.legend()
    p = out / "F_fat_torus.png"
    fig.savefig(p)
    plt.close(fig)
    paths.append(str(p))

    # 3. G heatmap
    G = np.array(gij["G"])
    fig, ax = plt.subplots(figsize=(5.2, 4.4), constrained_layout=True)
    im = ax.imshow(G, cmap="coolwarm", origin="upper")
    ax.set_xticks(range(G.shape[0]))
    ax.set_yticks(range(G.shape[0]))
    ax.set_xlabel("mode j")
    ax.set_ylabel("mode i")
    ax.set_title(r"TQGL–WS $G_{ij}$ (diag forced 0)")
    fig.colorbar(im, ax=ax, label="overlap (arb.)")
    p = out / "Gij_heatmap.png"
    fig.savefig(p)
    plt.close(fig)
    paths.append(str(p))

    # 4. Spectrum vs democratic
    fig, ax = plt.subplots(figsize=(6.2, 4.0), constrained_layout=True)
    eigs = np.array(gij["eigenvalues"])
    dem = np.array(
        [gij["democratic_lambda_min"]] * (len(eigs) - 1) + [gij["democratic_lambda_max"]]
    )
    idx = np.arange(1, len(eigs) + 1)
    ax.scatter(idx, eigs, s=60, color="#1f4e79", label="TQGL–WS numeric", zorder=3)
    ax.plot(idx, dem, "o--", color="#c47b17", label="democratic lock at ⟨G_ij⟩")
    ax.axhline(0, color="#444", lw=0.8)
    ax.set_xlabel("eigenvalue index (ascending)")
    ax.set_ylabel(r"$\lambda$")
    ax.set_title("Spec(G): computed vs democratic")
    ax.legend()
    p = out / "Gij_spectrum.png"
    fig.savefig(p)
    plt.close(fig)
    paths.append(str(p))

    return paths


def write_markdown(report: dict, path: Path, fig_paths: list[str]) -> None:
    b = report["bridge"]
    g = report["gij"]
    lines = [
        "# Bridge freeze + TQGL $G_{ij}$ report",
        "",
        "## Epistemic status (do not flatten)",
        "",
        "| Object | Status |",
        "|--------|--------|",
        "| Gap-equation algebra → α, ℓ(r) | **PROVEN** |",
        "| $A=E_{\\mathrm{vac}}/(\\hbar c_s)$, $\\alpha=\\lambda E_{\\mathrm{vac}}/\\Delta_0$ | **PROVEN** at $\\lambda=1$ identification |",
        "| Fat-torus $(r,R)$, $\\hbar c_s=0.07$ | **FROZEN_WORKING** |",
        "| $E_{\\mathrm{vac}}$ meV | **COMPUTED** (replaces 0.244) |",
        "| $\\Delta^*$ at $\\lambda=1$ | **FOLD THEOREM**: folded out |",
        "| $G_{ij}$ from TQGL–WS $K_0$ | **PROGRAM_COMPUTED** (not full torsion variation) |",
        "",
        "## Derived bridge",
        "",
        f"- $\\ell(r)={_fmt(b['ell_frozen'],6)}$ vs Corrections $\\ell={_fmt(b['ell_corrections'],6)}$ "
        f"(rel. err. {_fmt(100*b['ell_rel_err_vs_corrections'],3)}%)",
        f"- $F={_fmt(b['F'],5)}$ (withdrawn claim $F=0.596$)",
        f"- $E_{{\\mathrm{{vac}}}}={_fmt(b['E_vac_meV'],5)}$ meV (withdrawn $0.244$ meV)",
        f"- $A_{{\\mathrm{{geom}}}}={_fmt(b['A_geom_um_inv'],4)}\\,\\mu\\mathrm{{m}}^{{-1}}$, "
        f"$A R={_fmt(b['A_times_R'],4)}$ (withdrawn $A\\xi=0.296$)",
        f"- $\\alpha={_fmt(b['alpha'],5)}$, $\\alpha_c={_fmt(b['alpha_c'],5)}$, "
        f"$\\alpha/\\alpha_c={_fmt(b['alpha_over_ac'],4)}$ → folded_out={b['folded_out']}",
        f"- ℏ$c_s$ at fold: {_fmt(b['hbar_cs_at_fold_meV_um'],5)} meV·µm",
        f"- $E_{{\\mathrm{{vac}}}}$ that would reproduce withdrawn 35%: {_fmt(b['E_vac_for_withdrawn_35pct_meV'],4)} meV",
        "",
        "## $G_{ij}$ acceptance",
        "",
        f"- Hermitian {g['acceptance']['hermitian']}; diag 0 {g['acceptance']['diag_zero']}; "
        f"$\\lambda_{{\\min}}<0$ {g['acceptance']['lambda_min_lt_0']}; near-democratic {g['acceptance']['near_democratic']}",
        f"- $\\lambda_{{\\min}}={_fmt(g['lambda_min'],5)}$, $\\lambda_{{\\max}}={_fmt(g['lambda_max'],5)}$",
        f"- $g_{{0,\\mathrm{{eff}}}}={_fmt(g['g0_eff_mean_offdiag'],5)}$; rel. Frobenius {_fmt(g['rel_frobenius_to_democratic'],4)}",
        "",
        "Figures:",
        "",
    ]
    for fp in fig_paths:
        lines.append(f"- `{fp}`")
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    print("=== Derived bridge + TQGL G_ij ===\n")
    os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

    bridge = derive_bridge()
    gij = run_tqgl_ws_G(n_grid=24)
    curve = fold_curve(ell=bridge["ell_frozen"], n=60)

    out_dir = Path(os.environ.get("QVC_OUTPUT_DIR", str(ROOT / "output" / "bridge_freeze")))
    fig_dir = out_dir / "figures"
    out_dir.mkdir(parents=True, exist_ok=True)
    locks_dir = Path(os.environ.get("QVC_LOCKS_DIR", str(out_dir / "locks")))

    spec = frozen_fat_torus_spec()
    try:
        lock = freeze(spec, lock_name="frozen_working_fat_torus.json", locks_dir=locks_dir)
        lock_disp = str(lock)
    except OSError as exc:
        lock_disp = f"(lock write skipped: {exc})"

    figs = make_figures(bridge, gij, curve, fig_dir)
    report = {
        "bridge": bridge,
        "gij": {k: v for k, v in gij.items() if k != "G"} | {"G": gij["G"]},
        "lockfile": lock_disp,
        "figures": figs,
    }
    json_path = out_dir / "bridge_gij_report.json"
    md_path = out_dir / "bridge_gij_report.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")
    write_markdown(report, md_path, figs)

    print("Bridge:")
    print(f"  ell_frozen={bridge['ell_frozen']:.6f}  F={bridge['F']:.5f}  E_vac={bridge['E_vac_meV']:.5f} meV")
    print(f"  alpha/ac={bridge['alpha_over_ac']:.4f}  folded_out={bridge['folded_out']}")
    print("G_ij acceptance:", gij["acceptance"])
    print("Wrote", json_path)
    print("Wrote", md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
