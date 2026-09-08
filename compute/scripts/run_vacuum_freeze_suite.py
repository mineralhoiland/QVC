#!/usr/bin/env python3
"""Vacuum freeze computational suite + rewritten self-consistency table.

Runs:
  1. Fat-torus / geometry vacuum diagnostics (Tier-C vs old vs PFA)
  2. Inverse freeze diagnostics (what ℏc_s or r would hit 0.244 meV)
  3. Self-consistency table under locked fold (ℓ≈−2.791)
  4. Fold critical appendix numbers
  5. Spec(G) + Hopf/CS easy numerical checks
  6. Refresh provisional lockfile

Does NOT promote status to frozen — prints open decisions.

Usage:
  .venv/bin/python scripts/run_vacuum_freeze_suite.py
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qvc.gap_fold import ELL_CORRECTIONS, fold_critical_table, fold_curve
from qvc.hopf_cs_check import run_math_verification_bundle
from qvc.self_consistency import catalyzon_dressed_gap, self_consistency_table
from qvc.vacuum_freeze import freeze, provisional_corrections_self_energy_spec
from qvc.vacuum_suite import run_vacuum_freeze_suite

def _fmt(x: float | None, digits: int = 6) -> str:
    if x is None:
        return "—"
    if isinstance(x, float) and (x != x):  # NaN
        return "nan"
    return f"{x:.{digits}g}"


def write_markdown(report: dict, path: Path) -> None:
    sc = report["self_consistency"]
    fold = report["fold_critical"]
    vac = report["vacuum_suite"]
    mathb = report["math_verification"]

    lines = [
        "# Vacuum Freeze Suite Report",
        "",
        "**Status:** all vacuum meV values remain **PROVISIONAL** until a lock",
        "is promoted from `provisional` → `frozen` (see `docs/VACUUM_FREEZE_SPEC.md`).",
        "",
        "## 1. Fold critical table (locked)",
        "",
        "| Quantity | Value | Status |",
        "|----------|-------|--------|",
        f"| ℓ (Corrections) | {_fmt(fold['ell'], 4)} | locked |",
        f"| ℓ_tree = ln(1/2)+γ_E | {_fmt(fold['ell_tree'], 6)} | locked |",
        f"| δ_c | {_fmt(fold['delta_c'], 6)} | locked |",
        f"| α_c = δ_c | {_fmt(fold['alpha_c'], 6)} | locked |",
        f"| Δ₀ = 2J (meV) | {_fmt(fold['Delta0_meV'], 4)} | input |",
        f"| Δ_c = δ_c Δ₀ (meV) | {_fmt(fold['Delta_c_meV'], 5)} | locked |",
        f"| Near-fold scaling | `{fold['scaling']}` | locked |",
        f"| BKT | `{fold['bkt_status']}` | contingent |",
        "",
        "## 2. Self-consistency table (rewritten)",
        "",
        "Primary bridge (fixed geometry): α = √n · E_vac / Δ₀ → fold branch "
        f"δ(α; ℓ={ELL_CORRECTIONS}).",
        "Running-healing A_ZPF(Δ) is a separate diagnostic (often exponentially small when a_H ≫ R).",
        "",
        "Mechanisms kept separate: freeze loop ≠ Spec(G) catalyzon shift ≠ TQGL dynamics.",
        "",
        "| Geometry | x | F | E_vac Tier-C (meV) | Δ* (meV) | δ* | gap ↓ % | α/α_c | conv |",
        "|----------|---|---|-------------------|----------|-----|---------|-------|------|",
    ]
    for row in sc:
        lines.append(
            "| {geometry_label} | {x} | {F} | {E} | {D} | {d} | {g} | {a} | {c} |".format(
                geometry_label=row["geometry_label"],
                x=_fmt(row["x"], 4),
                F=_fmt(row["F"], 5),
                E=_fmt(row["E_vac_tierC_meV"], 5),
                D=_fmt(row["Delta_star_meV"], 5),
                d=_fmt(row["delta_star"], 5),
                g=_fmt(row["gap_reduction_pct"], 3),
                a=_fmt(row["alpha_over_ac"], 4),
                c="yes" if row["converged"] else "no/folded",
            )
        )

    if report.get("self_consistency_running_healing"):
        lines += [
            "",
            "### Running-healing comparison (diagnostic)",
            "",
            "| Geometry | Δ* | α/α_c | conv |",
            "|----------|-----|-------|------|",
        ]
        for row in report["self_consistency_running_healing"]:
            lines.append(
                f"| {row['geometry_label']} | {_fmt(row['Delta_star_meV'],5)} | "
                f"{_fmt(row['alpha_over_ac'],4)} | "
                f"{'yes' if row['converged'] else 'no/folded'} |"
            )

    cat = report["catalyzon"]
    lines += [
        "",
        "### Catalyzon channel (separate from freeze loop)",
        "",
        f"| Δ_ex | g₀ | Δ_cat=Δ_ex−g₀ | reduction |",
        f"|------|----|---------------|-----------|",
        f"| {_fmt(cat['Delta_ex_meV'],4)} | {_fmt(cat['g0_meV'],4)} | {_fmt(cat['Delta_cat_meV'],4)} | {_fmt(cat['reduction_pct'],3)}% |",
        "",
        f"_Mechanism:_ {cat['mechanism']}",
        "",
        "## 3. Fat-torus vacuum diagnostics",
        "",
    ]
    fat = vac["fat_torus"]
    lines += [
        f"- Regime: **{fat['regime']}** with x={_fmt(fat['x'],4)}, F={_fmt(fat['F'],5)}",
        f"- E_TierC = {_fmt(fat['E_tierC_meV'],5)} meV (provisional)",
        f"- E_old = {_fmt(fat['E_old_meV'],5)} meV — audit: `{fat['claimed_0p244_audit']['audit']}`",
        f"- E_old/E_PFA ≈ {_fmt(fat['E_old_over_E_pfa'],3)}; E_TierC/E_PFA ≈ {_fmt(fat['E_tierC_over_E_pfa'],3)}",
        f"- Jacobi θ₃ (moire α) = {_fmt(fat['jacobi_theta3_at_moire_alpha'],4)} ≠ F (PASS withdrawal)",
        "",
        "### Fold-in bounds (bridge rigor)",
        "",
    ]
    bounds = vac.get("fold_in_bounds", {})
    if bounds:
        lines += [
            f"- α_c = {_fmt(bounds.get('alpha_c'), 6)}; Δ₀ = {_fmt(bounds.get('Delta0_meV'), 4)} meV; "
            f"n = {bounds.get('n_modes')}",
            f"- E_vac ≤ **{_fmt(bounds.get('E_vac_max_bare_meV'), 5)}** meV (bare) for fold-in",
            f"- E_vac ≤ **{_fmt(bounds.get('E_vac_max_enhanced_meV'), 5)}** meV (√n enhanced) for fold-in",
            "",
        ]
    if vac.get("bridge_scan"):
        lines += [
            "| Geometry | E_TierC | bare α/α_c | bare Δ* | √n α/α_c | √n Δ* |",
            "|----------|---------|------------|---------|----------|-------|",
        ]
        for b in vac["bridge_scan"]:
            lines.append(
                "| {g} | {E} | {ba} | {bd} | {ea} | {ed} |".format(
                    g=b["geometry_label"],
                    E=_fmt(b["E_vac_tierC_meV"], 5),
                    ba=_fmt(b["bare"]["alpha_over_ac"], 4),
                    bd=_fmt(b["bare"]["Delta_star_meV"], 5),
                    ea=_fmt(b["enhanced_sqrt_n"]["alpha_over_ac"], 4),
                    ed=_fmt(b["enhanced_sqrt_n"]["Delta_star_meV"], 5),
                )
            )
        lines.append("")
    fd = vac.get("fat_torus_F_diagnostics")
    if fd:
        lines += [
            "### F(x) fat-torus diagnostics",
            "",
            f"- F_exact = {_fmt(fd['F_exact'], 5)}; F_n1 = {_fmt(fd['F_n1_only'], 5)} "
            f"(rel err {_fmt(fd['rel_err_n1_vs_exact'], 3)})",
            f"- Large-x asym ≈ {_fmt(fd['F_large_x_asym'], 5)}",
            "",
        ]
    lines += [
        "### Geometry comparison",
        "",
        "| Geometry | r (nm) | R (nm) | x | F | E_TierC | E_old | old→0.244? |",
        "|----------|--------|--------|---|---|---------|-------|------------|",
    ]
    for g in vac["geometry_table"]:
        lines.append(
            f"| {g['geometry']} | {_fmt(g['r_nm'],4)} | {_fmt(g['R_nm'],4)} | "
            f"{_fmt(g['x'],4)} | {_fmt(g['F'],5)} | {_fmt(g['E_tierC_meV'],5)} | "
            f"{_fmt(g['E_old_meV'],5)} | {g['old_matches_0p244']} |"
        )

    inv = vac["inverse_hbar_cs_for_0p244_tierC_fat"]
    invr = vac["inverse_r_for_0p244_tierC_fat"]
    lines += [
        "",
        "## 4. Inverse freeze diagnostics (not a freeze)",
        "",
        "What parameters would make Tier-C hit the provisional 0.244 meV *target* on the fat torus?",
        "",
        f"- Required ℏc_s = **{_fmt(inv['hbar_cs_required_meV_um'],5)}** meV·µm "
        f"(package uses {_fmt(inv['package_hbar_cs'],4)}; ratio {_fmt(inv['ratio_vs_package'],3)})",
        f"- Or, at package ℏc_s and R=50 nm: r* ≈ **{_fmt(invr['r_star_nm'],4)}** nm "
        f"(x*={_fmt(invr['x_star'],4)})",
        "",
        "## 5. Easy math verification bundle",
        "",
        f"- Spec(G) numeric match: **{mathb['spec_G']['numeric_matches_closed_form']}** "
        f"(λ_max={_fmt(mathb['spec_G']['lambda_max'],4)}, λ_min={_fmt(mathb['spec_G']['lambda_min'],4)} ×{mathb['spec_G']['lambda_min_mult']})",
        f"- Berry c₁ numerical: {_fmt(mathb['hopf_cs']['berry_c1_numerical'],6)} "
        f"(target {mathb['hopf_cs']['berry_c1_target']}) → {mathb['hopf_cs']['status']}",
        f"- Q_H = {_fmt(mathb['hopf_cs']['exact_locks']['Q_H'],4)} (CS/linking lock)",
        "",
        "## 6. Open freeze decisions",
        "",
    ]
    for d in vac["freeze_open_decisions"]:
        lines.append(f"- [ ] {d}")

    lines += [
        "",
        f"Provisional lockfile: `{report['lockfile']}`",
        "",
        "## Provenance",
        "",
        "- `qvc/vacuum_suite.py`, `qvc/self_consistency.py`, `qvc/gap_fold.py`",
        "- `docs/VACUUM_FREEZE_SPEC.md`, `docs/CLOSED_FORM_LOCKS.md`",
        "",
    ]
    path.write_text("\n".join(lines) + "\n")


def main() -> int:
    print("=== QVC Vacuum Freeze Computational Suite ===\n")

    vacuum = run_vacuum_freeze_suite()
    sc = self_consistency_table(mode="fixed_geometry")
    sc_run = self_consistency_table(mode="running_healing")
    fold = fold_critical_table()
    curve = fold_curve(n=40)
    mathb = run_math_verification_bundle()
    cat = catalyzon_dressed_gap(g0=0.21)

    out_dir = Path(os.environ.get("QVC_OUTPUT_DIR", str(ROOT / "output")))
    out_dir.mkdir(parents=True, exist_ok=True)
    locks_dir = Path(os.environ.get("QVC_LOCKS_DIR", str(ROOT / "qvc" / "locks")))
    try:
        lock = freeze(
            provisional_corrections_self_energy_spec(),
            lock_name="provisional_corrections_self_energy.json",
            locks_dir=locks_dir,
        )
        lock_disp = str(lock)
    except OSError as exc:
        lock_disp = f"(lock write skipped: {exc})"
        lock = None

    report = {
        "vacuum_suite": vacuum,
        "self_consistency": sc,
        "self_consistency_running_healing": sc_run,
        "fold_critical": fold,
        "fold_curve_sample": {
            "alpha_over_ac": curve["alpha_over_ac"],
            "delta_phys": curve["delta_phys"],
            "delta_unst": curve["delta_unst"],
        },
        "math_verification": mathb,
        "catalyzon": cat,
        "lockfile": lock_disp,
        "epistemic": {
            "vacuum_meV": "provisional",
            "fold": "locked",
            "Spec_G": "locked_phenomenological",
            "Q_H": "locked",
            "BKT": "contingent_vortex_RG",
        },
    }

    json_path = out_dir / "vacuum_freeze_suite_report.json"
    md_path = out_dir / "vacuum_freeze_suite_report.md"
    json_path.write_text(json.dumps(report, indent=2, default=str) + "\n")
    write_markdown(report, md_path)

    print("Fold critical:")
    print(json.dumps(fold, indent=2))
    print("\nSelf-consistency table:")
    for row in sc:
        print(
            f"  {row['geometry_label']}: E={row['E_vac_tierC_meV']:.6g} meV, "
            f"Δ*={row['Delta_star_meV']}, α/αc={row['alpha_over_ac']}, "
            f"conv={row['converged']}"
        )
    print(f"\nWrote {json_path}")
    print(f"Wrote {md_path}")
    print(f"Lockfile {lock_disp}")
    print("\nWARNING: meV values are PROVISIONAL — not a freeze promotion.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
