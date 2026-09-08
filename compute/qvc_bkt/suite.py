"""Driver: run C1–C6 in order, write figures + JSON/markdown summary."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from qvc_bkt.locks import BKTLocks, build_bkt_locks
from qvc_bkt.stiffness import run_c1
from qvc_bkt.vortex_core import run_c2
from qvc_bkt.kt_flow import competing_scale_scan, run_c3
from qvc_bkt.xi_kt import build_vortex_gap_curve, run_c4
from qvc_bkt.model_compare import compare_models, run_c5
from qvc_bkt.phonon_damping import run_c6
from qvc_bkt.figures import write_all_figures

DEFAULT_OUT = Path("/Users/mineralhoiland/Code/QVCCursor/notebooks/bkt_suite/output")


def _jsonify(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonify(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, np.ndarray):
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


def _drop_profiles(bundle: dict[str, Any]) -> dict[str, Any]:
    """JSON summary without large radial profiles (kept on disk as figures)."""
    out = dict(bundle)
    c2 = dict(out.get("C2") or {})
    c2.pop("profile", None)
    out["C2"] = c2
    # flows keep ell/K/y — downsample
    c3 = dict(out.get("C3") or {})
    flows = []
    for st in c3.get("flows_Tstar") or []:
        fl = dict(st.get("flow") or {})
        for key in ("ell", "K", "y"):
            arr = np.asarray(fl.get(key, []))
            if arr.size > 40:
                idx = np.linspace(0, arr.size - 1, 40).astype(int)
                fl[key] = arr[idx]
        flows.append({**st, "flow": fl})
    c3["flows_Tstar"] = flows
    out["C3"] = c3
    return out


def run_all(
    *,
    out_dir: Path | str | None = None,
    n_alpha: int = 80,
    make_figures: bool = True,
) -> dict[str, Any]:
    """Execute C1→C6, save artifacts, return the in-memory bundle."""
    out = Path(out_dir) if out_dir is not None else DEFAULT_OUT
    out.mkdir(parents=True, exist_ok=True)
    locks = build_bkt_locks()

    c1 = run_c1(locks, n_alpha=n_alpha)
    # Use GP conventional UV κ for the C2 logarithm (geometric floor is not in GP).
    kappa_gp_uv = c1["variants"]["conv_only"]["kappa_uv_meV"]
    c2 = run_c2(locks, kappa_uv_meV=kappa_gp_uv)
    c3 = run_c3(c1, c2, locks)
    c4 = run_c4(c3, c2, locks)

    # Counterfactual vortex Δ at T=1 K, conventional-only, if T* does not unbind.
    if not c4["primary_Tstar"]["unbinding_occurs"]:
        sc1 = competing_scale_scan(
            c1["variants"]["conv_only"]["rows"],
            locks=locks,
            T_K=1.0,
            E_core_over_kappa=float(c3["E_core_over_kappa_used"]),
            ell_ir=float(c3["scan_Tstar_primary"]["ell_ir"]),
        )
        curve = build_vortex_gap_curve(
            sc1, locks=locks, a_c_um=float(c2["a_c_um"]), ell_ir=sc1["ell_ir"]
        )
        c4["counterfactual_T1K_conv_only"] = {
            "scan": {
                k: sc1[k]
                for k in (
                    "T_K",
                    "kt_before_fold",
                    "fold_wins",
                    "alpha_kt_bare_over_ac",
                    "K_c",
                    "cs",
                )
            },
            "curve": curve,
            "label": "counterfactual: T=1 K, conventional κ only, no geometric floor",
        }
        if curve["unbinding_occurs"]:
            rows = curve["rows"]
            tA = np.array([r["t_A"] for r in rows], dtype=float)
            dv = np.array([r["Delta_vortex_meV"] for r in rows], dtype=float) / locks.Delta0_meV
            finite = np.array([r["finite_unbinding"] for r in rows], dtype=bool) & np.isfinite(dv)
            if np.count_nonzero(finite) >= 8:
                c4["counterfactual_T1K_conv_only"]["c5_preview"] = compare_models(
                    tA[finite], dv[finite], 0.0, label="vortex_Delta_T1K_conv_only"
                )

    c5 = run_c5(c4, locks)
    c6 = run_c6(c1, c3, locks)

    bundle = {
        "suite": "QVC / TEVC vortex RG  C1–C6",
        "locks": locks.as_dict(),
        "C1": c1,
        "C2": c2,
        "C3": c3,
        "C4": c4,
        "C5": c5,
        "C6": c6,
        "computed_vs_assumed": computed_vs_assumed(c1, c2, c3, c4, c5, c6, locks),
        "h_vrg1": h_vrg1_statement(c3, c1),
        "key_numbers": key_numbers(locks, c1, c2, c3, c4, c5, c6),
    }

    slim = _drop_profiles(bundle)
    (out / "summary.json").write_text(json.dumps(_jsonify(slim), indent=2))
    (out / "summary.md").write_text(render_markdown(bundle))
    figs: list[str] = []
    if make_figures:
        figs = write_all_figures(bundle, out)
        bundle["figures"] = figs
        slim["figures"] = figs
        (out / "summary.json").write_text(json.dumps(_jsonify(slim), indent=2))
    return bundle


def computed_vs_assumed(
    c1: dict, c2: dict, c3: dict, c4: dict, c5: dict, c6: dict, locks: BKTLocks
) -> dict[str, Any]:
    return {
        "C1": {
            "computed": [
                "physical-branch δ_+(α) from locked gap map",
                "⟨ρ²⟩=δ² via live-gap identification",
                "uniform phase-twist consistency of 2 n0 κ_GL ⟨ρ²⟩",
                "hopfion ring-weighted ⟨|ψ|²⟩ (amplitude sector only)",
            ],
            "derived": [
                "two-sector identity κ=2 κ_GL ⟨ρ²⟩+κ_geom (algebra)",
                "κ_GL=Δ₀ a_H² from GPE kinetic match",
                "K=κ/T_eff in Kosterlitz convention",
            ],
            "ansatz": [
                "n0=1/(π a_H²) healing-disk density (primary)",
                "Peotta–Törmä κ_geom=(C_eff/2π)Δ₀ frozen geometric floor",
                "alternate moiré-cell n0 and C_eff=Q_H variants",
            ],
            "not_a_theorem": True,
        },
        "C2": {
            "computed": [
                "radial GP winding-1 profile",
                "E_v on a disk",
                "half-density core radius",
                "E_core = E_v − πκ ln(R/a_c) given the cutoff",
            ],
            "ansatz": ["a_c=ξ healing-length cutoff used in the IR logarithm unless noted"],
            "derived": ["GP energy functional matching the truncated EOM"],
        },
        "C3": {
            "computed": [
                "Kosterlitz ODE at frozen α",
                "physical-branch competing-scale scan",
                "T_eff grid winner table",
            ],
            "derived": ["H-VRG1 criterion (KT surface before fold)", "y=exp(−E_core/T_eff)"],
            "truncated": ["Chern–Simons (default)"],
            "labeled_option": ["anyon K_c=(2/π)(1−1/N)²"],
            "not_used": ["microscopic β_α", "T_BKT^bare input"],
        },
        "C4": {
            "computed": ["ξ_KT from flow ℓ* when unbinding occurs"],
            "conjecture": ["Δ∼ħ c_s/ξ_KT as a gap identification"],
            "locked": ["finite Δ_c of the amplitude fold is not replaced"],
        },
        "C5": {
            "computed": [
                "OLS slopes and Gaussian log-likelihood on ≳1.5 decades of t_A",
                "locked fold branch vs √t_A and vs essential singularity",
            ],
            "disclaimer": c5["disclaimer"],
        },
        "C6": {
            "computed": ["AHNS ratio (ħ/τ_v)/(πκ) along the branch given the estimates"],
            "labeled_estimate": [
                "ħ/τ_v ∼ g_ph=0.1Δ₀ from existing GPE phonon potential",
                "ħ/τ_v ∼ 2π λ_ep kT with λ_ep=1 Gao–Khalaf-type OOM",
            ],
            "locked_algebra": ["η_ret = 2 Δ₀ ξ/(ħ c_s)"],
        },
        "locks_status": "Option A* fat torus from qvc.tqgl.locks / schwinger_twin_channel",
        "Q_H": locks.Q_H,
        "lambda_min_equals_minus_g0": abs(locks.lambda_min_meV + locks.g0_meV) < 1e-12,
    }


def h_vrg1_statement(c3: dict, c1: dict) -> dict[str, Any]:
    h = c3["h_vrg1_Tstar_primary"]
    teff = c3["T_eff_grid"]
    return {
        "at_Tstar_primary_closure": h,
        "primary_closure": c1["primary_closure"],
        "Tstar_conv_only_kt_before_fold": c3["scan_Tstar_conv_only"]["kt_before_fold"],
        "Tstar_anyon_kt_before_fold": c3["scan_Tstar_anyon"]["kt_before_fold"],
        "T_eff_grid": teff,
        "reading": (
            "H-VRG1 is a competing-scale test on a computed-or-ansatz κ(α), "
            "not a BKT theorem. Option A* keeps α frozen. CS is truncated "
            "in the default scan."
        ),
    }


def key_numbers(locks: BKTLocks, c1, c2, c3, c4, c5, c6) -> dict[str, Any]:
    p = c1["primary"]
    return {
        "r_nm": 8.0,
        "R_nm": 50.0,
        "hbar_cs_meV_um": locks.hbar_cs_meV_um,
        "F": locks.F,
        "E_vac_meV": locks.E_vac_meV,
        "alpha_c": locks.alpha_c,
        "delta_c": locks.delta_c,
        "Delta_c_meV": locks.Delta_c_meV,
        "Delta0_meV": locks.Delta0_meV,
        "lambda_star": locks.lambda_star,
        "Delta_star_meV": locks.Delta_star_meV,
        "T_star_K": locks.T_star_K,
        "T_star_mK": 1e3 * locks.T_star_K,
        "E_pair_meV": locks.E_pair_meV,
        "Q_H": locks.Q_H,
        "lambda_min_meV": locks.lambda_min_meV,
        "g0_meV": locks.g0_meV,
        "kT_star_meV": locks.kT_star_meV,
        "kappa_GL_meV_um2": p["kappa_GL_meV_um2"],
        "n0_um2": p["n0_um2"],
        "kappa_uv_meV": p["kappa_uv_meV"],
        "kappa_fold_meV": p["kappa_fold_meV"],
        "K_uv_Tstar": p["K_uv_Tstar"],
        "K_fold_Tstar": p["K_fold_Tstar"],
        "K_c_bosonic": locks.k_c_bosonic,
        "K_c_anyon": locks.k_c_anyon,
        "E_core_over_kappa": c2["E_core_over_kappa"],
        "a_c_um": c2["a_c_um"],
        "a_c_label": c2["a_c_label"],
        "xi_um": c2["radial"]["xi_um"],
        "kt_before_fold_Tstar_primary": c3["scan_Tstar_primary"]["kt_before_fold"],
        "kt_before_fold_Tstar_conv_only": c3["scan_Tstar_conv_only"]["kt_before_fold"],
        "C5_fold_slope_free": c5["fold_comparison"]["sqrt"]["slope_free"],
        "C5_delta_ll_sqrt_minus_essential": c5["fold_comparison"]["delta_ll_sqrt_minus_essential"],
        "C6_jump_vs_crossover": c6["jump_vs_crossover"],
        "twist_rel_err": p["twist_uniform"]["rel_err"],
    }


def render_markdown(bundle: dict[str, Any]) -> str:
    kn = bundle["key_numbers"]
    h = bundle["h_vrg1"]
    cva = bundle["computed_vs_assumed"]
    c3 = bundle["C3"]
    lines = [
        "# QVC / TEVC vortex RG — C1–C6 summary",
        "",
        "BKT is treated as a physically interesting vortex sector. C1 is partly ansatz;",
        "this suite is **not** a BKT theorem. Hopfion ≠ 2D phase vortex ≠ CS fluxon.",
        "Chern–Simons is truncated unless the anyon $K_c$ option is on. Option A* freezes $\\alpha$.",
        "",
        "## Locked numbers (Option A*, fat torus)",
        "",
        f"- $(r,R)=(8,50)$ nm, $\\hbar c_s={kn['hbar_cs_meV_um']}$ meV·µm, $F={kn['F']:.5f}$, $E_{{\\mathrm{{vac}}}}={kn['E_vac_meV']:.4f}$ meV",
        f"- $\\alpha_c=\\delta_c={kn['alpha_c']:.4f}$, $\\Delta_c={kn['Delta_c_meV']:.4f}$ meV, $\\Delta_0={kn['Delta0_meV']:.3f}$ meV",
        f"- $\\lambda^\\star={kn['lambda_star']:.4f}$, $\\Delta^\\star={kn['Delta_star_meV']:.4f}$ meV, $T^*={kn['T_star_mK']:.1f}$ mK, $E_{{\\mathrm{{pair}}}}={kn['E_pair_meV']:.4f}$ meV",
        f"- Democratic $\\lambda_{{\\min}}=-g_0={kn['lambda_min_meV']:.4f}$ meV, $Q_H=1/N={kn['Q_H']}$",
        "",
        "## Competing-scale test (H-VRG1)",
        "",
        f"- Primary closure (healing-disk $n_0$ + frozen $\\kappa_{{\\mathrm{{geom}}}}$) at $T^*$: **KT before fold = {kn['kt_before_fold_Tstar_primary']}**",
        f"- Conventional-only (no geometric floor) at $T^*$: **KT before fold = {kn['kt_before_fold_Tstar_conv_only']}**"
        + (
            f" (IR plasma at $\\alpha/\\alpha_c\\approx {c3['scan_Tstar_conv_only']['alpha_kt_ir_over_ac']:.3f}$, finite-$y$ separatrix; bare $K$ stays above $2/\\pi$)"
            if c3["scan_Tstar_conv_only"]["alpha_kt_ir_over_ac"] is not None
            and c3["scan_Tstar_conv_only"]["alpha_kt_bare_over_ac"] is None
            else ""
        ),
        f"- $K(\\alpha=0)={kn['K_uv_Tstar']:.3f}$, $K(\\alpha_c)={kn['K_fold_Tstar']:.3f}$, $K_c=2/\\pi={kn['K_c_bosonic']:.4f}$",
        f"- $\\kappa_{{\\mathrm{{UV}}}}={kn['kappa_uv_meV']:.4f}$ meV, $\\kappa_{{\\mathrm{{fold}}}}={kn['kappa_fold_meV']:.4f}$ meV",
        "",
        h["reading"],
        "",
        "## C2 core",
        "",
        f"- $E_{{\\mathrm{{core}}}}/\\kappa={kn['E_core_over_kappa']:.4f}$ (cutoff: {kn['a_c_label']})",
        f"- $a_c={kn['a_c_um']:.6f}$ µm, $\\xi={kn['xi_um']:.6f}$ µm",
        "",
        "## C5 fold-branch comparison (not a BKT disproof)",
        "",
        f"- free slope of $\\ln(\\delta-\\delta_c)$ vs $\\ln t_A$: {kn['C5_fold_slope_free']:.4f} (fold theorem: $1/2$)",
        f"- $\\Delta\\ell\\ell$ (√ minus essential): {kn['C5_delta_ll_sqrt_minus_essential']:.3f}",
        "",
        "## C6 H-VRG2",
        "",
        f"- {kn['C6_jump_vs_crossover']}",
        "",
        "## Computed vs assumed",
        "",
    ]
    for key in ("C1", "C2", "C3", "C4", "C5", "C6"):
        block = cva[key]
        lines.append(f"### {key}")
        for lab in ("computed", "derived", "ansatz", "conjecture", "labeled_estimate", "truncated"):
            if lab in block:
                lines.append(f"- **{lab}:** " + "; ".join(block[lab]))
        lines.append("")
    lines.extend(
        [
            "Withdrawn Casimir / gap-schedule folklore is not used as physics.",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser(description="QVC BKT C1–C6 suite")
    p.add_argument("--out", type=Path, default=DEFAULT_OUT)
    p.add_argument("--n-alpha", type=int, default=80)
    p.add_argument("--no-figures", action="store_true")
    args = p.parse_args()
    bundle = run_all(out_dir=args.out, n_alpha=args.n_alpha, make_figures=not args.no_figures)
    kn = bundle["key_numbers"]
    print(json.dumps({k: kn[k] for k in kn if k != "locks"}, indent=2, default=str))
    print("wrote", args.out)
