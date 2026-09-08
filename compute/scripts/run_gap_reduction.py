#!/usr/bin/env python3
"""Recompute 35/52/92 gap reductions and derive λ from fold + catalyzon.

Usage:
  MPLBACKEND=Agg .venv/bin/python scripts/run_gap_reduction.py
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

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qvc.gap_reduction import recompute_gap_reduction

plt.rcParams.update(
    {
        "font.size": 10,
        "axes.titlesize": 11,
        "figure.dpi": 120,
        "savefig.dpi": 200,
        "axes.grid": True,
        "grid.alpha": 0.3,
    }
)


def _fmt(x, d=5):
    if x is None:
        return "—"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, float) and x != x:
        return "nan"
    return f"{x:.{d}g}"


def make_figure(rep: dict, path: Path) -> None:
    rec = rep["recommendation"]
    fold = rep["fold_channel"]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8), constrained_layout=True)

    ax = axes[0]
    labels = ["35% folklore", "52% folklore", "92% (of Δ_QVC)", "92% (of Δ₀)", "max static (fold)"]
    vals = [
        rep["folklore"]["stage1_static"]["reduction_pct"],
        rep["folklore"]["stage2_catalyzon"]["reduction_pct"],
        rep["folklore"]["stage3_tqgl_from_Delta_QVC"]["reduction_pct"],
        rep["folklore"]["stage3_tqgl_from_Delta0"]["reduction_pct"],
        fold["max_static_reduction_pct"],
    ]
    colors = ["#8aa0b4", "#8aa0b4", "#8aa0b4", "#8aa0b4", "#1f4e79"]
    ax.barh(labels, vals, color=colors)
    ax.set_xlabel("gap reduction (%)")
    ax.set_title("Folklore vs max static fold")
    ax.set_xlim(0, 100)

    ax = axes[1]
    lams = [
        ("λ_fold", fold["lambda_fold"]),
        ("λ* (0.90 fold)", rec["lambda_star"]),
        ("λ_35", fold["lambda_for_withdrawn_35pct"]),
        ("λ_cat (g0=0.1)", rec["lambda_cat_preprint_0p1"]),
        ("λ=1 (unit)", 1.0),
    ]
    ax.bar([n[0] for n in lams], [n[1] for n in lams], color="#1f4e79")
    ax.axhline(fold["lambda_fold"], color="#c47b17", ls="--", lw=1, label="fold-out boundary")
    ax.set_ylabel("participation λ")
    ax.set_title("Derived participation")
    ax.tick_params(axis="x", rotation=25)
    ax.legend()

    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path)
    plt.close(fig)


def write_md(rep: dict, path: Path, fig: Path) -> None:
    f = rep["fold_channel"]
    c = rep["catalyzon_channel"]
    t = rep["tqgl_92"]
    r = rep["recommendation"]
    folk = rep["folklore"]
    lines = [
        "# Gap reduction recomputation and λ",
        "",
        "Channels are **parallel**, not stacked. 0.244 meV is withdrawn.",
        "",
        "## Folklore (withdrawn)",
        "",
        f"- 35% static: 0.600 → 0.389 ({_fmt(folk['stage1_static']['reduction_pct'],3)}%)",
        f"- 52% catalyzon mix: 0.600 → 0.289 ({_fmt(folk['stage2_catalyzon']['reduction_pct'],3)}%)",
        f"- 92% of Δ_QVC: 0.389 → 0.032 ({_fmt(folk['stage3_tqgl_from_Delta_QVC']['reduction_pct'],3)}%)",
        f"- vs Δ₀: 0.600 → 0.032 ({_fmt(folk['stage3_tqgl_from_Delta0']['reduction_pct'],3)}%, not 92%)",
        "",
        "## Fold channel (proven)",
        "",
        f"- E_vac = {_fmt(f['E_vac_meV'])} meV, ℓ = {_fmt(f['ell'],6)}, α_c = {_fmt(f['alpha_c'],5)}",
        f"- λ_fold = α_c Δ₀ / E_vac = **{_fmt(f['lambda_fold'],4)}**  (λ=1 is fold-out)",
        f"- Max static reduction = **{_fmt(f['max_static_reduction_pct'],3)}%** at Δ_c = {_fmt(f['Delta_c_meV'])} meV",
        f"- λ for withdrawn 35%: {_fmt(f['lambda_for_withdrawn_35pct'],4)}",
        f"- δ(1 ps)=0.032/0.6 below δ_c: **{f['delta_1ps_below_delta_c']}**",
        f"- {f['theorem']}",
        "",
        "## Catalyzon channel",
        "",
        f"- Identification g₀ = λ E_vac. At λ_fold, g₀ = {_fmt(c['at_lambda_fold']['g0_meV'])} meV, "
        f"Δ_cat = {_fmt(c['at_lambda_fold']['Delta_cat_meV'])} meV "
        f"({_fmt(c['at_lambda_fold']['reduction_pct'],3)}%).",
        f"- Scaled g₀ from 0.1×(E_vac/0.244) = {_fmt(c['g0_scaled_from_0p1_and_0p244'])} meV "
        f"→ λ_cat = {_fmt(c['at_scaled_g0']['lambda'],4)}.",
        f"- Preprint g₀=0.1 / E_vac → λ = {_fmt(c['at_preprint_g0_0p1']['lambda'],4)} "
        f"(exceeds unit: {c['at_preprint_g0_0p1']['exceeds_unit_participation']}).",
        "",
        "## TQGL 92%",
        "",
        f"- Status: `{t['status']}`",
        f"- Drive scale E_vac/0.244 = {_fmt(t['V_ZPF_scale_new_over_old'],4)}",
        f"- Phenomenological remaining Δ(1 ps) ≈ {_fmt(t['Delta_1ps_phenomenological_meV'])} meV "
        f"({_fmt(t['reduction_pct_phenomenological'],3)}%) — not a split-step rerun.",
        "",
        "## Recommended λ*",
        "",
        f"- **λ* = {_fmt(r['lambda_star'],4)}** = 0.90 λ_fold",
        f"- Fold Δ* = {_fmt(r['fold']['Delta_star_meV'])} meV "
        f"({_fmt(r['fold']['gap_reduction_pct'],3)}%)",
        f"- Parallel catalyzon Δ_cat = {_fmt(r['catalyzon_parallel_not_stacked']['Delta_cat_meV'])} meV "
        f"({_fmt(r['catalyzon_parallel_not_stacked']['reduction_pct'],3)}%)",
        f"- {r['why_not_92']}",
        "",
        f"Figure: `{fig}`",
        "",
    ]
    path.write_text("\n".join(lines))


def main() -> int:
    print("=== Gap reduction + λ participation ===\n")
    rep = recompute_gap_reduction()
    out = Path(os.environ.get("QVC_OUTPUT_DIR", str(ROOT / "output" / "gap_reduction")))
    out.mkdir(parents=True, exist_ok=True)
    fig = out / "gap_reduction_lambda.png"
    make_figure(rep, fig)
    json_path = out / "gap_reduction_report.json"
    md_path = out / "gap_reduction_report.md"
    json_path.write_text(json.dumps(rep, indent=2) + "\n")
    write_md(rep, md_path, fig)
    rec = rep["recommendation"]
    print(f"λ_fold={rep['fold_channel']['lambda_fold']:.4f}  λ*={rec['lambda_star']:.4f}")
    print(f"max static reduction={rep['fold_channel']['max_static_reduction_pct']:.2f}%")
    print(f"Δ* at λ*={rec['fold']['Delta_star_meV']}")
    print(f"92% below fold: {rep['fold_channel']['delta_1ps_below_delta_c']}")
    print("Wrote", json_path)
    print("Wrote", md_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
