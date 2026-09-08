"""Track A L1 parameter dictionary driver: algebra card, figures, summaries.

This module is independent of the B1/B2 ``qvc_l1.suite`` (G_ij / Lindblad).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from qvc_l1.acoustic import acoustic_matching
from qvc_l1.catalyzon import catalyzon_card
from qvc_l1.dictionary_figures import write_all_figures
from qvc_l1.gl_coefficients import gl_coefficients
from qvc_l1.hall import hall_convention
from qvc_l1.lengths import length_dictionary
from qvc_l1.locks import (
    FOLKLORE_DELTA_1PS_MEV,
    FOLKLORE_DELTA_QVC_MEV,
    FOLKLORE_E_CAS_MEV,
    G1_1PS,
    LIVE_GAP_1PS_MEV,
    RETIRED_G0_MEV,
    RETIRED_KAPPA_TOR,
    L1Locks,
    build_l1_locks,
)
from qvc_l1.mexican_hat import mexican_hat

DEFAULT_OUT = Path("/Users/mineralhoiland/Code/QVCCursor/notebooks/l1_dictionary/output")
MIRROR_OUT = Path("/Users/mineralhoiland/Code/QVCCursor/output/l1_dictionary")


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


def _drop_arrays(bundle: dict[str, Any]) -> dict[str, Any]:
    """JSON card without dense plot arrays."""
    out = dict(bundle)
    hat = dict(out.get("mexican_hat") or {})
    for key in ("rho", "V_bare", "V_gpe_sign", "V_catalyzon_dressed"):
        hat.pop(key, None)
    out["mexican_hat"] = hat
    return out


def status_ledger(bundle: dict[str, Any]) -> dict[str, list[str]]:
    """What is locked / derived / matched / ansatz / computed."""
    return {
        "locked": [
            "fat torus (r,R)=(8,50) nm, ħ c_s=0.07 meV·µm, Δ₀=0.600 meV",
            "a_H = r (Option A* GL/GPE healing freeze)",
            "η_ret = 2 Δ₀ ξ / (ħ c_s) as algebra (ξ remains a choice)",
            "λ★ = 0.90 λ_fold, N=5, Q_H=1/N, k=2N",
            "F, E_vac, α_c=δ_c, Δ_c from the frozen-cavity bridge",
        ],
        "derived": [
            "ξ_ph = ħ c_s / Δ₀ and η_ret(ξ_ph)=2",
            "η_ret(a_H)=2 Δ₀ a_H/(ħ c_s)≈0.137",
            "κ_GL = Δ₀ a_H²",
            "α = λ E_vac / Δ₀ (gap-map coupling)",
            "g0 = λ E_vac, Δ_cat = Δ₀ − g0, Spec(G): λ_min=−g0",
            "λ_fold = α_c Δ₀ / E_vac",
            "Hall: k=2N, δσ_xy = e²/(k h) = Q_H e²/(2h)",
            "pairing k Q_H = 2",
        ],
        "matched": [
            "α_GL = −Δ₀ (Landau mass |α_GL|=Δ₀)",
            "β = Δ₀ / n_ref so the Mexican-hat vacuum sits on hopfion n_ref",
            "R_WS = ħ c_s/(π k_B T*) at T* (analogue Hawking matching, not duality)",
        ],
        "ansatz": [
            "κ_S ≈ κ_GL a_H² = Δ₀ a_H⁴ (Derrick balance on the locked tube)",
        ],
        "computed": [
            "n_ref from hopfion_initial (ring-weighted |ψ|²)",
            "Mexican-hat well depths and catalyzon barrier reduction",
            "T*, E_pair from Schwinger twin-channel at g0★",
            "live gap 0.6305 meV / g(1)=0.8832 coherent GPE baseline (not rerun here)",
            "Berry R_mean≈0.908 (ratio, not a radius)",
        ],
        "alternate_not_lock": [
            "ξ=0.181 µm → η_ret≈3.1 (preprint coherence; ħ c_s/Δ_QVC withdrawn)",
            "R_WS(285 mK)≈0.91 µm (same matching formula, withdrawn T*)",
        ],
        "retired": [
            "g0~0.21 as a derived coupling",
            "κ_tor≈0.25 as an insertion",
            "cartoon α0=−0.5, α_QVC=−0.35, 42% barrier drop",
            "holographic γ_QM≈3.73 meV·nm² as κ_GL",
            "folklore 0.244 / 0.389 / 0.032 / 92%",
        ],
        "out_of_scope": [
            "Hubbard microscopic derivation of GL coefficients",
            "Chern–Simons in the GPE EOM",
            "H-TOR1 torsion vertex",
            "holographic dictionary",
        ],
    }


def parameter_card(locks: L1Locks, bundle: dict[str, Any]) -> dict[str, Any]:
    gl = bundle["gl"]
    cat = bundle["catalyzon"]
    hat = bundle["mexican_hat"]
    lengths = bundle["lengths"]
    hall = bundle["hall"]
    ac = bundle["acoustic"]
    return {
        "r_nm": 1e3 * locks.r_um,
        "R_nm": 1e3 * locks.R_um,
        "hbar_cs_meV_um": locks.hbar_cs_meV_um,
        "Delta0_meV": locks.Delta0_meV,
        "F": locks.F,
        "E_vac_meV": locks.E_vac_meV,
        "alpha_c": locks.alpha_c,
        "delta_c": locks.delta_c,
        "Delta_c_meV": locks.Delta_c_meV,
        "lambda_fold": locks.lambda_fold,
        "lambda_star": locks.lambda_star,
        "Delta_star_meV": locks.Delta_star_meV,
        "T_star_mK": locks.T_star_mK,
        "E_pair_meV": locks.E_pair_meV,
        "N": locks.N_layers,
        "Q_H": locks.Q_H,
        "k_CS": hall["k_CS"],
        "a_H_um": lengths["lengths"]["a_H"]["value_um"],
        "xi_ph_um": lengths["lengths"]["xi_ph"]["value_um"],
        "xi_181_um": lengths["lengths"]["xi_181"]["value_um"],
        "eta_ret_a_H": lengths["eta_ret_a_H"],
        "eta_ret_xi_ph": lengths["eta_ret_xi_ph"],
        "eta_ret_xi_181": lengths["eta_ret_xi_181"],
        "kappa_GL_meV_um2": gl["kappa_GL"]["value_meV_um2"],
        "alpha_GL_meV": gl["alpha_GL"]["value_meV"],
        "alpha_gap_map_star": gl["alpha_gap_map_star"]["alpha"],
        "n_ref": gl["n_ref"]["n_ref"],
        "beta_meV": gl["beta"]["beta_meV"],
        "kappa_S_meV_um4": gl["kappa_S"]["kappa_S_meV_um4"],
        "g0_star_meV": cat["g0_star_meV"],
        "Delta_cat_star_meV": cat["Delta_cat_star_meV"],
        "g0_star_reduction_pct": cat["at_lambda_star"]["reduction_pct"],
        "g0_fold_meV": cat["g0_fold_meV"],
        "Delta_cat_fold_meV": cat["Delta_cat_fold_meV"],
        "mexican_hat_barrier_reduction_pct": hat["catalyzon_dressed"]["barrier_reduction_pct"],
        "mexican_hat_gpe_depth_change_frac": hat["gpe_sign_shift"]["depth_change_frac"],
        "delta_sigma_xy_S": hall["delta_sigma_xy_S"],
        "sigma_xy_parent_S": hall["sigma_xy_parent_S"],
        "R_WS_um": ac["R_WS_um"],
        "Berry_R_mean": ac["Berry_R_mean"],
        "Delta_live_1ps_meV": LIVE_GAP_1PS_MEV,
        "g1_1ps": G1_1PS,
        "retired_g0_meV": RETIRED_G0_MEV,
        "retired_kappa_tor": RETIRED_KAPPA_TOR,
        "folklore_not_results": {
            "E_Cas_meV": FOLKLORE_E_CAS_MEV,
            "Delta_QVC_meV": FOLKLORE_DELTA_QVC_MEV,
            "Delta_1ps_meV": FOLKLORE_DELTA_1PS_MEV,
        },
    }


def render_markdown(bundle: dict[str, Any]) -> str:
    locks: L1Locks = bundle["locks_obj"]
    p = bundle["parameter_card"]
    gl = bundle["gl"]
    hat = bundle["mexican_hat"]
    cat = bundle["catalyzon"]
    ledger = bundle["status_ledger"]
    lines = [
        "# Track A L1 parameter dictionary",
        "",
        "Fat-torus Option A* working point. Three lengths are distinct. "
        "Catalyzon is a parallel channel (not stacked with fold or live-gap). "
        "Hall algebra uses CS level $k=2N$ with $Q_H=1/N$.",
        "",
        "## Parameter card",
        "",
        f"- $(r,R)=({p['r_nm']:.0f},{p['R_nm']:.0f})$ nm, "
        f"$\\hbar c_s={p['hbar_cs_meV_um']}$ meV·µm, $\\Delta_0={p['Delta0_meV']:.3f}$ meV  "
        f"**[locked]**",
        f"- $F={p['F']:.5f}$, $E_{{\\mathrm{{vac}}}}={p['E_vac_meV']:.4f}$ meV, "
        f"$\\alpha_c=\\delta_c={p['alpha_c']:.4f}$, $\\Delta_c={p['Delta_c_meV']:.4f}$ meV  "
        f"**[locked / derived]**",
        f"- $\\lambda_{{\\mathrm{{fold}}}}={p['lambda_fold']:.4f}$, "
        f"$\\lambda^\\star={p['lambda_star']:.4f}$, $\\Delta^\\star={p['Delta_star_meV']:.4f}$ meV  "
        f"**[derived / locked]**",
        f"- $T^*={p['T_star_mK']:.1f}$ mK, $E_{{\\mathrm{{pair}}}}={p['E_pair_meV']:.4f}$ meV, "
        f"$N={p['N']}$, $Q_H={p['Q_H']}$, $k={p['k_CS']}$  **[computed / locked]**",
        "",
        "### Lengths",
        "",
        f"- $a_H=r={p['a_H_um']:.3f}$ µm $=8$ nm, $\\eta_{{\\mathrm{{ret}}}}(a_H)={p['eta_ret_a_H']:.3f}$  **[locked]**",
        f"- $\\xi_{{\\mathrm{{ph}}}}=\\hbar c_s/\\Delta_0={p['xi_ph_um']:.4f}$ µm, "
        f"$\\eta_{{\\mathrm{{ret}}}}(\\xi_{{\\mathrm{{ph}}}})={p['eta_ret_xi_ph']:.0f}$ identically  **[derived]**",
        f"- $\\xi_{{181}}={p['xi_181_um']:.3f}$ µm (alternate; preprint coherence "
        f"$\\hbar c_s/\\Delta_{{\\mathrm{{QVC}}}}$), $\\eta_{{\\mathrm{{ret}}}}={p['eta_ret_xi_181']:.2f}$  "
        f"**[alternate, not the lock]**",
        "",
        "### GL coefficients (matched, not Hubbard)",
        "",
        f"- $\\alpha=\\lambda^\\star E_{{\\mathrm{{vac}}}}/\\Delta_0={p['alpha_gap_map_star']:.4f}$  **[derived]**",
        f"- $\\alpha_{{\\mathrm{{GL}}}}=-\\Delta_0={p['alpha_GL_meV']:.3f}$ meV  **[matched]**",
        f"- $\\kappa_{{\\mathrm{{GL}}}}=\\Delta_0 a_H^2={p['kappa_GL_meV_um2']:.4e}$ meV·µm$^2$  **[derived]**",
        f"- $\\lambda=\\lambda^\\star={p['lambda_star']:.4f}$  **[locked]**",
        f"- $n_{{\\mathrm{{ref}}}}={p['n_ref']:.5f}$ on hopfion_initial  **[computed]**",
        f"- $\\beta=\\Delta_0/n_{{\\mathrm{{ref}}}}={p['beta_meV']:.4f}$ meV  **[matched]**",
        f"- $\\kappa_S\\approx\\Delta_0 a_H^4={p['kappa_S_meV_um4']:.4e}$ meV·µm$^4$  **[ansatz]**",
        "",
        "### Catalyzon (parallel channel)",
        "",
        f"- $g_0^\\star=\\lambda^\\star E_{{\\mathrm{{vac}}}}={p['g0_star_meV']:.4f}$ meV, "
        f"$\\Delta_{{\\mathrm{{cat}}}}={p['Delta_cat_star_meV']:.4f}$ meV "
        f"({p['g0_star_reduction_pct']:.1f}%)  **[derived]**",
        f"- at $\\lambda_{{\\mathrm{{fold}}}}$: $g_0=\\Delta_c={p['g0_fold_meV']:.4f}$ meV "
        f"({100.0 * locks.alpha_c:.1f}%=$\\alpha_c$)  **[derived]**",
        f"- Spec$(G)$: $\\lambda_{{\\min}}=-g_0^\\star$, multiplicity $N-1={locks.N_layers - 1}$  **[locked algebra]**",
        "",
        "### Mexican hat",
        "",
        f"- GPE-sign well-depth change: {100.0 * hat['gpe_sign_shift']['depth_change_frac']:.1f}% "
        f"(attractive $V_{{\\mathrm{{ZPF}}}}$, well deepens)  **[computed]**",
        f"- Catalyzon-dressed barrier reduction: "
        f"{p['mexican_hat_barrier_reduction_pct']:.1f}% "
        f"$=1-(\\Delta_{{\\mathrm{{cat}}}}/\\Delta_0)^2$ (replaces cartoon 42%)  **[computed]**",
        "",
        "### Hall and matching",
        "",
        f"- $\\delta\\sigma_{{xy}}=e^2/(k h)={p['delta_sigma_xy_S']:.4e}$ S, "
        f"parent $\\sigma_{{xy}}=(k/2)e^2/h={p['sigma_xy_parent_S']:.4e}$ S  **[derived]**",
        f"- $R_{{\\mathrm{{WS}}}}=\\hbar c_s/(\\pi k_B T^*)={p['R_WS_um']:.2f}$ µm  **[matched]**",
        f"- Berry $R_{{\\mathrm{{mean}}}}={p['Berry_R_mean']:.3f}$ is a ratio, not a radius  **[computed]**",
        "",
        "### Coherent GPE baseline (Lindblad is another agent)",
        "",
        f"- $\\Delta_{{\\mathrm{{live}}}}(1\\,\\mathrm{{ps}})={p['Delta_live_1ps_meV']:.4f}$ meV, "
        f"$g^{{(1)}}={p['g1_1ps']:.4f}$  **[computed]**",
        "",
        "## Status ledger",
        "",
    ]
    for tag in (
        "locked",
        "derived",
        "matched",
        "ansatz",
        "computed",
        "alternate_not_lock",
        "retired",
        "out_of_scope",
    ):
        lines.append(f"### {tag}")
        for item in ledger[tag]:
            lines.append(f"- {item}")
        lines.append("")
    lines.extend(
        [
            "## Notes",
            "",
            f"- Democratic $\\lambda_{{\\min}}=-g_0$ check: "
            f"{cat['lambda_min_equals_minus_g0']}.",
            f"- $n_{{\\mathrm{{ref}}}}$ formula: {gl['n_ref']['formula']}.",
            f"- $\\kappa_S$ ansatz: {gl['kappa_S']['ansatz_label']}.",
            "- Do not stack catalyzon, fold, and live-gap percentages.",
            "- Do not box $\\eta_{\\mathrm{ret}}=3.1$ as the unique lock.",
            f"- Retired insertions: $g_0\\sim{p['retired_g0_meV']}$ meV, "
            f"$\\kappa_{{\\mathrm{{tor}}}}\\approx{p['retired_kappa_tor']}$.",
            "",
        ]
    )
    return "\n".join(lines) + "\n"


def build_dictionary(*, lambda_star_frac: float = 0.90) -> dict[str, Any]:
    locks = build_l1_locks(lambda_star_frac=lambda_star_frac)
    gl = gl_coefficients(locks)
    cat = catalyzon_card(locks)
    hat = mexican_hat(locks, beta_meV=gl["beta"]["beta_meV"], n_ref=gl["n_ref"]["n_ref"])
    lengths = length_dictionary(locks)
    hall = hall_convention(locks)
    ac = acoustic_matching(locks)
    bundle: dict[str, Any] = {
        "locks": locks.as_dict(),
        "locks_obj": locks,
        "gl": gl,
        "catalyzon": cat,
        "mexican_hat": hat,
        "lengths": lengths,
        "hall": hall,
        "acoustic": ac,
        "live_gap_baseline": {
            "Delta_live_1ps_meV": LIVE_GAP_1PS_MEV,
            "g1_1ps": G1_1PS,
            "status": "computed",
            "note": "Coherent truncated GPE/DCE baseline. Lindblad is another agent.",
        },
    }
    bundle["parameter_card"] = parameter_card(locks, bundle)
    bundle["status_ledger"] = status_ledger(bundle)
    return bundle


def _write_outputs(slim: dict[str, Any], md: str, out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    (out / "summary.json").write_text(json.dumps(_jsonify(slim), indent=2))
    (out / "summary.md").write_text(md)


def run_all(
    *,
    out_dir: Path | str | None = None,
    make_figures: bool = True,
    mirror: bool = True,
) -> dict[str, Any]:
    """Build the dictionary, write summary.json / summary.md and figures."""
    out = Path(out_dir) if out_dir is not None else DEFAULT_OUT
    out.mkdir(parents=True, exist_ok=True)
    bundle = build_dictionary()
    slim = _drop_arrays(bundle)
    slim.pop("locks_obj", None)
    md = render_markdown(bundle)
    _write_outputs(slim, md, out)
    figs: list[str] = []
    if make_figures:
        figs = write_all_figures(bundle, out)
        bundle["figures"] = figs
        slim["figures"] = figs
        _write_outputs(slim, md, out)
    if mirror and out.resolve() != MIRROR_OUT.resolve():
        MIRROR_OUT.mkdir(parents=True, exist_ok=True)
        _write_outputs(slim, md, MIRROR_OUT)
        if make_figures:
            write_all_figures(bundle, MIRROR_OUT)
    return bundle
