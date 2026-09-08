#!/usr/bin/env python3
"""Run the honest TQGL v4 GPE+DCE suite and write output/tqgl_v4/.

Usage:
  MPLBACKEND=Agg .venv/bin/python scripts/run_tqgl_v4.py
  MPLBACKEND=Agg .venv/bin/python scripts/run_tqgl_v4.py --quick
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl")

import matplotlib

matplotlib.use("Agg")
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qvc.gap_fold import fold_curve
from qvc.tqgl.coupled import run_two_timescale_suite
from qvc.tqgl.figures import write_figures
from qvc.tqgl.locks import (
    DELTA_1PS_FOLKLORE_MEV,
    DELTA_QVC_FOLKLORE_MEV,
    E_CAS_WITHDRAWN_MEV,
    build_frozen_tqgl,
)

MIRROR = Path("/Users/mineralhoiland/Code/QVCCursor/output/tqgl_v4")


def _jsonable(x: Any) -> Any:
    if isinstance(x, dict):
        return {str(k): _jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        if x.size > 4000 and x.ndim >= 2:
            return {"shape": list(x.shape), "omitted": True, "note": "array stored as CSV/figure"}
        return [_jsonable(v) for v in x.tolist()]
    if isinstance(x, (np.floating, float)):
        v = float(x)
        if math.isnan(v) or math.isinf(v):
            return None
        return v
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    return x


def _csv(path: Path, header: str, cols: list[np.ndarray]) -> None:
    arr = np.column_stack(cols)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(path, arr, delimiter=",", header=header, comments="")


def write_csvs(bundle: dict[str, Any], outdir: Path) -> list[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    c = bundle["coupled_gpe_dce"]
    res = bundle["dce_fold_resonant"]
    det = bundle["dce_fold_detuning"]
    sch = bundle["schwinger"]
    lock = bundle["lock"]
    paths: list[str] = []

    p = outdir / "coupled_timeseries.csv"
    _csv(
        p,
        "t_ps,Delta_live_meV,Delta_fold_of_A_meV,g1,A_um_inv,lorentzian,mu_nl_meV",
        [
            c["t_ps"],
            c["Delta_live_meV"],
            np.nan_to_num(c["Delta_fold_of_A_meV"], nan=np.nan),
            c["g1"],
            c["A_um_inv"],
            c["lorentzian"],
            c["mu_nl_meV"],
        ],
    )
    paths.append(str(p))

    p = outdir / "dce_resonant.csv"
    _csv(
        p,
        "t_ps,A_um_inv,Delta_fold_meV,lorentzian,folded_out",
        [res["t_ps"], res["A_um_inv"], res["Delta_fold_meV"], res["lorentzian"], res["folded_out"].astype(float)],
    )
    paths.append(str(p))

    p = outdir / "dce_detuning.csv"
    _csv(
        p,
        "t_ps,A_um_inv,Delta_fold_meV,lorentzian,folded_out",
        [det["t_ps"], det["A_um_inv"], det["Delta_fold_meV"], det["lorentzian"], det["folded_out"].astype(float)],
    )
    paths.append(str(p))

    fc = fold_curve(lock["ell"], n=150)
    p = outdir / "fold_curve.csv"
    phys = np.array([np.nan if x is None else x for x in fc["delta_phys"]])
    unst = np.array([np.nan if x is None else x for x in fc["delta_unst"]])
    _csv(p, "alpha_over_ac,delta_phys,delta_unst", [np.array(fc["alpha_over_ac"]), phys, unst])
    paths.append(str(p))

    p = outdir / "schwinger.csv"
    _csv(
        p,
        "A_over_Ac,Gamma_sch_ps_inv,T_K,Gamma_thermal,Gamma_total",
        [sch["ratio_sweep"], sch["Gamma_sch_sweep"], sch["T_K"], sch["Gamma_thermal"], sch["Gamma_total"]],
    )
    paths.append(str(p))

    berry = bundle["berry"]
    p = outdir / "berry_lambda_min.csv"
    _csv(
        p,
        "lambda_min_GOE,lambda_min_GUE",
        [berry["lambda_min_GOE"], berry["lambda_min_GUE"]],
    )
    paths.append(str(p))
    return paths


def write_methods(lock: dict[str, Any], bundle: dict[str, Any], path: Path) -> None:
    c = bundle["coupled_gpe_dce"]
    res = bundle["dce_fold_resonant"]
    text = f"""# TQGL v4 methods (equations actually solved)

Status: truncated 2D GPE + DCE ODE. Option A* fat torus (r,R)=({lock['r_nm']:.0f},{lock['R_nm']:.0f}) nm.

## Truncation (not in the EOM)

Chern–Simons, Berry connection, and Skyrme/Hopf terms are **not** discretised.
Democratic Spec(G) is a static lock (λ_min=−g₀, multiplicity N−1), not a GPE term.
The 1 ps window does **not** include Lindblad dissipation.

## GPE (split-step Fourier)

$$
i\\hbar\\partial_t\\psi = \\Bigl[-\\tfrac{{\\hbar^2}}{{2m^*}}\\nabla^2 + \\alpha|\\psi|^2 + \\beta|\\psi|^4 + V_{{\\rm ZPF}}(\\mathbf r,t) + V_{{\\rm ph}}(\\mathbf r,t)\\Bigr]\\psi
$$

- Box: {c['L_um']*1e3:.0f} nm, grid {c['n_grid']}×{c['n_grid']}, dt={c['dt_ps']*1e3:.1f} fs, T={c['t_end_ps']} ps.
- ℏ²/(2m*)=Δ₀ a_H² with a_H=r={lock['r_nm']:.0f} nm; α=−Δ₀, β=Δ₀/n₀, n₀=1.
- Hopfion IC: 10% core depletion, phase Q_H θ with Q_H=1/N, noise 3%.
- $V_{{\\rm ZPF}}=-E_{{\\rm vac}}\\,(A(t)/A(0))\\,\\mathrm{{Gaussian}}_{{\\rm ring}}$. Drive is locked $E_{{\\rm vac}}={lock['E_vac_meV']:.4f}$ meV, **not** 0.244. Scale vs withdrawn Casimir: {lock['drive_scale_vs_old_0244']:.3f}.
- $V_{{\\rm ph}}=0.1\\Delta_0\\sin(\\omega_0 t)\\,\\mathrm{{Gaussian}}_{{\\rm ring}}$, $\\omega_0=2\\Delta_0/\\hbar$.

## Live gap (uses ψ; no time schedule)

$$
\\Delta_{{\\rm live}}=\\Delta_0\\sqrt{{n_w/n_{{\\rm ref}}}},\\qquad n_w=\\langle|\\psi|^2\\rangle_{{\\rm ring}}.
$$

Compared to the static physical fold branch Δ(A) at the running DCE amplitude, and to Δ★(λ★)={lock['Delta_star_meV']:.4f} meV. If Δ_live<Δ_c={lock['Delta_c_meV']:.4f} meV the run is flagged unphysical vs the static fold. There is **no** 0.600→0.389→0.032 prescription. Do not quote 92%.

## DCE ODE (RK4)

$$
\\dot A=\\bigl[\\Lambda(\\Delta)\\,\\mathcal L(\\delta\\omega/\\Gamma)-\\Lambda_{{\\rm loss}}\\bigr]A,\\quad
\\Lambda(\\Delta)=(m\\Delta/2\\hbar)\\,Q/(2\\pi),\\quad \\mathcal L(x)=1/(1+x^2).
$$

- Coupled 1 ps resonant pump tracks the IC gap (ω/2 = 2Δ₀/ℏ). Fold-slaved 80 ps resonant pump tracks Δ★; detuning-feedback uses +2Γ.
- Λ_loss = Γ_cav(Δ★) (cavity linewidth), **not** fitted to A*/A_c=0.52.
- A_c = α_c Δ₀/(ℏ c_s) = {lock['A_c_fold_um_inv']:.4f} µm⁻¹ (fold, not BKT 1.2233).
- Coupled 1 ps: Δ = Δ_live(ψ). Long window (~80 ps): Δ slaved to the physical fold branch Δ(A). Two-timescale justification: Γ_cav⁻¹ ~ tens of ps ≫ 1 ps GPE window.

Fold-slaved resonant A*/A_c (computed, not imposed): {res['A_star_over_Ac']:.4f}.

## Schwinger / Berry (diagnostics only)

E_pair = 2|λ_min| with |λ_min|=g₀=λ★ E_vac (democratic). Withdrawn E_pair=0.172 is not used unless that number is re-derived. Berry GOE/GUE Monte Carlo is a random-matrix probe, **separate** from democratic G.
"""
    path.write_text(text, encoding="utf-8")


def write_paper_numbers(lock: dict[str, Any], bundle: dict[str, Any], path: Path) -> None:
    c = bundle["coupled_gpe_dce"]
    res = bundle["dce_fold_resonant"]
    det = bundle["dce_fold_detuning"]
    sch = bundle["schwinger"]
    berry = bundle["berry"]
    text = f"""# Numbers for the paper (TQGL v4, Option A*)

Do **not** quote 92%, Δ(1 ps)=0.032 meV, Δ_QVC=0.389 meV, or E_Cas=0.244 meV as results.
Withdrawn 0.244 appears only as the drive-ratio denominator E_vac/0.244={lock['drive_scale_vs_old_0244']:.3f}.

## Frozen locks (reuse bridge / fold / λ)

| Quantity | Value |
|----------|-------|
| (r, R) | ({lock['r_nm']:.0f}, {lock['R_nm']:.0f}) nm |
| E_vac | {lock['E_vac_meV']:.4f} meV |
| F(x) | {lock['F']:.5f} |
| α_c = δ_c | {lock['alpha_c']:.4f} |
| Δ_c | {lock['Delta_c_meV']:.4f} meV |
| Δ₀ | {lock['Delta0_meV']:.3f} meV |
| λ_fold | {lock['lambda_fold']:.4f} |
| λ★ = 0.90 λ_fold | {lock['lambda_star']:.4f} |
| Δ★(λ★) | {lock['Delta_star_meV']:.4f} meV ({lock['gap_reduction_star_pct']:.1f}% static) |
| max static reduction | {lock['max_static_reduction_pct']:.1f}% (not 92%) |
| A_c^{{fold}} | {lock['A_c_fold_um_inv']:.4f} µm⁻¹ |
| A(λ★)/A_c^{{fold}} | {lock['A_star_um_inv']/lock['A_c_fold_um_inv']:.4f} |

## Coupled GPE (1 ps, live gap from ψ)

| Quantity | Value |
|----------|-------|
| Δ_live(0) | {c['Delta_live_meV'][0]:.4f} meV (=Δ₀ by n_ref) |
| Δ_live(1 ps) | {c['Delta_live_final_meV']:.4f} meV |
| min Δ_live | {c['Delta_live_min_meV']:.4f} meV |
| live < Δ_c ? | {c['any_live_below_delta_c']} |
| g¹(1 ps) | {c['g1_final']:.4f} |
| A(1 ps)/A_c^{{fold}} | {c['A_final_over_Ac_fold']:.4f} |
| ‖ψ‖ drift | {c['norm_rel_drift']:.2e} |

## DCE (fold-slaved, 80 ps) — A*/A_c is computed, not 0.52

| Case | A*/A_c^{{fold}} | Δ_ss (meV) | folded out | L(0) |
|------|-----------------|------------|------------|------|
| resonant | {res['A_star_over_Ac']:.4f} | {res['Delta_star_dce_meV']:.4f} | {res['any_folded_out']} | {res['L0']:.4f} |
| detuning | {det['A_star_over_Ac']:.4f} | {det['Delta_star_dce_meV']:.4f} | {det['any_folded_out']} | {det['L0']:.4f} |

## Schwinger / Berry (not in EOM)

| Quantity | Value |
|----------|-------|
| g₀ = λ★ E_vac | {sch['g0_democratic_meV']:.4f} meV |
| |λ_min| democratic | {abs(sch['lambda_min_democratic_meV']):.4f} meV |
| E_pair = 2|λ_min| | {sch['E_pair_meV']:.4f} meV (not 0.172) |
| A_eff/A_{{c,Sch}} | {sch['A_eff_over_Ac_schwinger']:.4f} |
| T* | {sch['T_star_K']*1e3:.1f} mK |
| Berry R_mean | {berry['R_mean']:.3f} (not forced to 1.43) |
"""
    path.write_text(text, encoding="utf-8")


def results_payload(lock: dict[str, Any], bundle: dict[str, Any], fig_paths: list[str], csv_paths: list[str]) -> dict[str, Any]:
    c = bundle["coupled_gpe_dce"]
    res = bundle["dce_fold_resonant"]
    det = bundle["dce_fold_detuning"]
    sch = bundle["schwinger"]
    berry = bundle["berry"]
    return {
        "status": "PROGRAM_COMPUTED_truncated_GPE_DCE",
        "geometry": "option_A_star_fat_torus",
        "locks": {
            "r_nm": lock["r_nm"],
            "R_nm": lock["R_nm"],
            "E_vac_meV": lock["E_vac_meV"],
            "F": lock["F"],
            "alpha_c": lock["alpha_c"],
            "delta_c": lock["delta_c"],
            "Delta_c_meV": lock["Delta_c_meV"],
            "Delta0_meV": lock["Delta0_meV"],
            "lambda_fold": lock["lambda_fold"],
            "lambda_star": lock["lambda_star"],
            "Delta_star_meV": lock["Delta_star_meV"],
            "delta_star": lock["delta_star"],
            "gap_reduction_star_pct": lock["gap_reduction_star_pct"],
            "max_static_reduction_pct": lock["max_static_reduction_pct"],
            "A_c_fold_um_inv": lock["A_c_fold_um_inv"],
            "A_star_um_inv": lock["A_star_um_inv"],
            "drive_scale_vs_old_0244": lock["drive_scale_vs_old_0244"],
        },
        "withdrawn_folklore_not_results": {
            "E_Cas_meV": E_CAS_WITHDRAWN_MEV,
            "Delta_QVC_meV": DELTA_QVC_FOLKLORE_MEV,
            "Delta_1ps_meV": DELTA_1PS_FOLKLORE_MEV,
            "do_not_quote_92_pct": True,
            "used_only_as_drive_denominator": True,
        },
        "coupled_gpe_1ps": {
            "n_grid": c["n_grid"],
            "dt_ps": c["dt_ps"],
            "Delta_live_0_meV": float(c["Delta_live_meV"][0]),
            "Delta_live_1ps_meV": c["Delta_live_final_meV"],
            "Delta_live_min_meV": c["Delta_live_min_meV"],
            "live_below_delta_c": c["any_live_below_delta_c"],
            "g1_final": c["g1_final"],
            "A_final_over_Ac_fold": c["A_final_over_Ac_fold"],
            "norm_rel_drift": c["norm_rel_drift"],
            "live_gap_uses_psi": True,
            "no_time_schedule": True,
        },
        "dce": {
            "resonant": {
                "A_star_over_Ac_fold": res["A_star_over_Ac"],
                "A_star_um_inv": res["A_star_um_inv"],
                "Delta_star_dce_meV": res["Delta_star_dce_meV"],
                "any_folded_out": res["any_folded_out"],
                "runaway_clipped": res["runaway_clipped"],
                "L0": res["L0"],
                "Lorentzian_identically_zero": res["Lorentzian_identically_zero"],
            },
            "detuning": {
                "A_star_over_Ac_fold": det["A_star_over_Ac"],
                "A_star_um_inv": det["A_star_um_inv"],
                "Delta_star_dce_meV": det["Delta_star_dce_meV"],
                "any_folded_out": det["any_folded_out"],
                "runaway_clipped": det["runaway_clipped"],
                "L0": det["L0"],
            },
            "A_c_is_fold_not_BKT": True,
            "A_c_not_backcomputed_from_0p52": True,
            "paper_ratio_0p52_not_imposed": True,
        },
        "schwinger": {
            "g0_democratic_meV": sch["g0_democratic_meV"],
            "lambda_min_democratic_meV": sch["lambda_min_democratic_meV"],
            "E_pair_meV": sch["E_pair_meV"],
            "A_eff_over_Ac_schwinger": sch["A_eff_over_Ac_schwinger"],
            "T_star_K": sch["T_star_K"],
            "Gamma_sch_ps_inv": sch["Gamma_sch_ps_inv"],
            "E_pair_not_0p172": sch["E_pair_not_0p172_unless_rederived"],
        },
        "berry": {
            "R_mean": berry["R_mean"],
            "R_median": berry["R_median"],
            "mean_abs_GOE": berry["mean_abs_GOE"],
            "mean_abs_GUE": berry["mean_abs_GUE"],
            "separate_from_democratic_G": True,
        },
        "truncation": c["truncation"],
        "two_timescale": bundle["two_timescale"],
        "figures": fig_paths,
        "csv": csv_paths,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="32² / short DCE for smoke tests")
    args = ap.parse_args()

    outdir = ROOT / "output" / "tqgl_v4"
    outdir.mkdir(parents=True, exist_ok=True)
    figdir = outdir / "figures"
    lock_obj = build_frozen_tqgl(lambda_star_frac=0.90)
    lock = lock_obj.as_dict()

    if args.quick:
        from qvc.tqgl.coupled import run_coupled_gpe_dce
        from qvc.tqgl.dce import DCEConfig, run_dce_fold_slaved
        from qvc.tqgl.diagnostics import berry_goe_gue, schwinger_twin_channel

        coupled = run_coupled_gpe_dce(
            lock_obj, n_grid=32, t_end_ps=0.05, dt_ps=0.002, pump_mode="resonant"
        )
        bundle = {
            "lock": lock,
            "coupled_gpe_dce": coupled,
            "dce_fold_resonant": run_dce_fold_slaved(
                lock_obj, DCEConfig(pump_mode="resonant", t_end_ps=8.0, dt_ps=0.05)
            ),
            "dce_fold_detuning": run_dce_fold_slaved(
                lock_obj, DCEConfig(pump_mode="detuning", t_end_ps=8.0, dt_ps=0.05)
            ),
            "schwinger": schwinger_twin_channel(lock_obj),
            "berry": berry_goe_gue(n_samples=128, g_bar=abs(lock_obj.lambda_star * lock_obj.E_vac_meV)),
            "two_timescale": {"gpe_window_ps": 0.05, "dce_window_ps": 8.0, "justification": "quick smoke"},
        }
    else:
        bundle = run_two_timescale_suite(lock_obj, n_grid=128, t_gpe_ps=1.0, dt_gpe_ps=0.002, n_berry=2000)

    csv_paths = write_csvs(bundle, outdir / "csv")
    fig_paths = write_figures(bundle, figdir)
    write_methods(lock, bundle, outdir / "METHODS.md")
    write_paper_numbers(lock, bundle, outdir / "PAPER_NUMBERS.md")
    payload = results_payload(lock, bundle, fig_paths, csv_paths)
    (outdir / "results.json").write_text(json.dumps(_jsonable(payload), indent=2), encoding="utf-8")

    MIRROR.mkdir(parents=True, exist_ok=True)
    for item in outdir.iterdir():
        dest = MIRROR / item.name
        if item.is_dir():
            if dest.exists():
                shutil.rmtree(dest)
            shutil.copytree(item, dest)
        else:
            shutil.copy2(item, dest)

    print(f"wrote {outdir}")
    print(f"mirrored {MIRROR}")
    print(f"E_vac={lock['E_vac_meV']:.4f} meV  Delta_star={lock['Delta_star_meV']:.4f} meV")
    print(f"Delta_live(1ps)={payload['coupled_gpe_1ps']['Delta_live_1ps_meV']:.4f} meV")
    print(f"A*/Ac resonant={payload['dce']['resonant']['A_star_over_Ac_fold']:.4f}")
    print(f"A*/Ac detuning={payload['dce']['detuning']['A_star_over_Ac_fold']:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
