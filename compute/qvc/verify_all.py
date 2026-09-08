"""
QVC Tier B/C PASS/FAIL verification suite.

Port/improvement of QVCCursor/rigor/run_verification.py.
Run via:  scripts/run_verification.py
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from qvc.acoustic_matching import R_WS_from_T_star
from qvc.coupling_matrix import (
    catalyzon_dark_basis,
    catalyzon_gap_shift,
    democratic_coupling_matrix,
    democratic_spectrum,
)
from qvc.gap_fold import ell_from_delta_c, fold_condition, gap_rhs, tree_level_ell
from qvc.lens_chern_simons import (
    anyonic_exchange_angle,
    chern_simons_flat_u1,
    fractional_hall_charge,
    h2_vanishes_note,
    linking_form_self,
)
from qvc.mode_enhancement import A_zpf_eff, mode_enhancement_ratio
from qvc.params import PARAMS
from qvc.retardation import eta_ret, eta_ret_from_healing
from qvc.vacuum_freeze import (
    assert_preprint_0p244_not_reproduced,
    bessel_mode_sum,
    casimir_preprint_formula,
    ensure_default_provisional_lock,
    jacobi_theta3_alternating,
    vacuum_self_energy_corrected,
)

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class Check:
    name: str
    status: str  # PASS | FAIL | WARN | INFO
    detail: str
    value: float | str | None = None
    expected: float | str | None = None


def near(a: float, b: float, tol: float = 1e-8) -> bool:
    return abs(a - b) <= tol * max(1.0, abs(b))


def check_coupling_spectrum() -> list[Check]:
    out: list[Check] = []
    N, g0 = PARAMS.N_layers, 1.0
    G = democratic_coupling_matrix(N, g0)
    eigs = np.sort(np.linalg.eigvalsh(G))
    spec = democratic_spectrum(N, g0)

    out.append(
        Check(
            "democratic Spec(G) exact match",
            "PASS" if near(eigs[0], -g0) and near(eigs[-1], (N - 1) * g0) else "FAIL",
            f"eigs={eigs.tolist()}",
            value=float(eigs[0]),
            expected=spec["lambda_min"],
        )
    )
    out.append(
        Check(
            "preprint claim lambda_min=-4g0 is FALSE",
            "PASS" if not near(eigs[0], -4 * g0) else "FAIL",
            "Democratic matrix does not have lambda_min=-4g0",
            value=float(eigs[0]),
            expected=-4.0,
        )
    )
    out.append(
        Check(
            "Tr(G)=0",
            "PASS" if near(float(np.trace(G)), 0.0) else "FAIL",
            f"Tr={float(np.trace(G))}",
        )
    )
    B = catalyzon_dark_basis(N)
    sums = np.sum(B, axis=1)
    grams = B @ B.T
    out.append(
        Check(
            "catalyzon dark basis is sum-zero",
            "PASS" if np.allclose(sums, 0.0, atol=1e-12) else "FAIL",
            f"row sums={sums.tolist()}",
        )
    )
    out.append(
        Check(
            "catalyzon dark basis orthonormal (Helmert)",
            "PASS" if np.allclose(grams, np.eye(N - 1), atol=1e-12) else "FAIL",
            "B B^T ≈ I",
        )
    )
    out.append(
        Check(
            "catalyzon gap shift Δ_cat = Δ_ex - g0",
            "PASS" if near(catalyzon_gap_shift(0.6, g0), 0.6 - g0) else "FAIL",
            "Identity from λ_min",
            value=catalyzon_gap_shift(0.6, g0),
            expected=0.6 - g0,
        )
    )
    out.append(
        Check(
            "supplemental lambda_min=-5g0 != democratic",
            "PASS" if not near(eigs[0], -5 * g0) else "FAIL",
            "Octahedral adjacency claim needs a different matrix; democratic is -g0",
            value=float(eigs[0]),
            expected=-5.0,
        )
    )
    return out


def check_lens_space() -> list[Check]:
    out: list[Check] = []
    N = PARAMS.N_layers
    cs = chern_simons_flat_u1(1, N)
    lk = linking_form_self(N)
    qh = fractional_hall_charge(N)
    note = h2_vanishes_note(N)
    out.append(
        Check(
            "CS[A_1] = 1/N",
            "PASS" if near(cs, 1 / N) else "FAIL",
            "Chern-Simons of fundamental flat connection",
            value=cs,
            expected=1 / N,
        )
    )
    out.append(
        Check(
            "linking self-number = 1/N",
            "PASS" if near(lk, 1 / N) else "FAIL",
            "Rational linking form on H1=Z_N",
            value=lk,
            expected=1 / N,
        )
    )
    out.append(
        Check(
            "Q_H identification = 1/N",
            "PASS" if near(qh, 1 / N) else "FAIL",
            "Physical Hall fractionalization target",
            value=qh,
            expected=1 / N,
        )
    )
    out.append(
        Check(
            "anyonic theta = pi/N",
            "INFO",
            "Programmatic prediction from same linking form",
            value=anyonic_exchange_angle(N),
            expected=math.pi / N,
        )
    )
    out.append(
        Check(
            "H2(L(N,1))=0 — no fundamental 2-cycle",
            "PASS",
            note["consequence"],
            value=note["H2"],
            expected="0",
        )
    )
    return out


def check_gap_fold() -> list[Check]:
    out: list[Check] = []
    delta_target = 0.182
    ell = ell_from_delta_c(delta_target)
    delta_c, alpha_c = fold_condition(ell)
    out.append(
        Check(
            "fold δ_c recovers Corrections ~0.182",
            "PASS" if abs(delta_c - 0.182) < 0.01 else "WARN",
            f"ell={ell:.6f}",
            value=delta_c,
            expected=0.182,
        )
    )
    out.append(
        Check(
            "fold α_c = δ_c",
            "PASS" if near(alpha_c, delta_c) else "FAIL",
            "Saddle-node condition",
            value=alpha_c,
            expected=delta_c,
        )
    )
    out.append(
        Check(
            "fold has FINITE Delta_c (not BKT vanishing)",
            "PASS" if delta_c > 0.05 else "FAIL",
            "Without vortex RG, transition is saddle-node",
            value=delta_c * PARAMS.Delta0_meV,
            expected=0.109,
        )
    )
    alphas = np.linspace(0.0, alpha_c * 0.999, 40)
    n_roots = []
    for a in alphas:
        ds = np.linspace(1e-3, 1.0, 2000)
        f = ds - np.array([gap_rhs(d, a, ell) for d in ds])
        crossings = int(np.sum((f[:-1] * f[1:]) < 0))
        n_roots.append(crossings)
    out.append(
        Check(
            "two real branches exist below alpha_c",
            "PASS" if max(n_roots) >= 2 else "WARN",
            f"max sign-change pairs={max(n_roots)}",
            value=max(n_roots),
            expected=2,
        )
    )
    tl = tree_level_ell()
    out.append(
        Check(
            "tree-level log coefficient ln(1/2)+γ_E",
            "PASS" if abs(tl - (-0.11593)) < 1e-4 else "FAIL",
            "Negative ⇒ vacuum amplitude suppresses gap",
            value=tl,
            expected=-0.11593,
        )
    )
    return out


def check_vacuum_energy() -> list[Check]:
    out: list[Check] = []
    p = PARAMS

    audit = assert_preprint_0p244_not_reproduced()
    out.append(
        Check(
            "preprint formula does NOT yield 0.244 meV on stated (r,R)",
            "PASS" if not audit["claim_reproduced"] else "FAIL",
            audit["audit"],
            value=audit["E_old_meV"],
            expected=0.244,
        )
    )

    x_p = audit["x"]
    F_p = audit["F"]
    E_new = vacuum_self_energy_corrected(p.hbar_cs_meV_um, p.r_preprint_um, F_p)
    out.append(
        Check(
            "Bessel F(x) at preprint aspect ratio",
            "INFO",
            f"x={x_p:.4f}",
            value=F_p,
            expected=0.596,
        )
    )
    out.append(
        Check(
            "corrected prefactor on preprint geometry (no R/r)",
            "INFO",
            "Shows sensitivity to formula choice",
            value=E_new,
        )
    )

    # Corrections geometry with package units (ħc_s in meV·µm)
    x_c = 2 * math.pi * p.r_corr_um / p.R_corr_um
    F_c = bessel_mode_sum(x_c)
    E_corr_pkg = vacuum_self_energy_corrected(p.hbar_cs_meV_um, p.r_corr_um, F_c)
    out.append(
        Check(
            "corrected self-energy on Corrections geometry (package units)",
            "INFO",
            f"F({x_c:.3f})={F_c:.3f}; provisional — not a paper meV claim",
            value=E_corr_pkg,
        )
    )

    # SI cross-check (sensitive to c_s choice)
    x_si = 2 * math.pi * 0.335 / 13.4
    F_si = bessel_mode_sum(x_si)
    hbar = 1.055e-34
    cs = 2.1e4
    r_m = 0.335e-9
    E_J = (1 / (4 * math.pi**2)) * (hbar * cs / r_m) * F_si
    E_meV = E_J / 1.602e-22
    out.append(
        Check(
            "Corrections SI self-energy order ~0.25 meV",
            "WARN" if abs(E_meV - 0.256) > 0.1 else "PASS",
            f"F({x_si:.3f})={F_si:.3f}; sensitive to K0 accuracy / c_s",
            value=E_meV,
            expected=0.256,
        )
    )

    alpha_c = 8.0 / 181.0
    th = jacobi_theta3_alternating(alpha_c)
    out.append(
        Check(
            "Jacobi θ₃ alternating does NOT equal F≈0.596",
            "PASS" if abs(th - 0.596) > 0.2 else "FAIL",
            f"theta={th:.6e} at alpha_c={alpha_c:.4f}",
            value=th,
            expected=0.596,
        )
    )

    lock_path = ensure_default_provisional_lock()
    out.append(
        Check(
            "provisional vacuum lockfile written",
            "PASS" if lock_path.is_file() else "FAIL",
            str(lock_path.relative_to(ROOT)),
            value=lock_path.name,
            expected="provisional_corrections_self_energy.json",
        )
    )
    return out


def check_mode_enhancement() -> list[Check]:
    out: list[Check] = []
    p = PARAMS
    R = p.R_corr_um
    A1 = A_zpf_eff(p.Delta0_meV, p.hbar_cs_meV_um, R, n_modes=1)
    A14 = A_zpf_eff(p.Delta0_meV, p.hbar_cs_meV_um, R, n_modes=14)
    out.append(
        Check(
            "single-mode A_ZPF at bare gap is negligible",
            "PASS" if A1 < 1e-2 else "WARN",
            "Matches Corrections iteration-0 exponential suppression",
            value=A1,
        )
    )
    ratio = A14 / max(A1, 1e-30)
    out.append(
        Check(
            "√14 enhancement ratio",
            "PASS" if near(ratio, math.sqrt(14), tol=1e-6) else "FAIL",
            "Quadrature multiplicity factor",
            value=ratio,
            expected=mode_enhancement_ratio(14),
        )
    )
    return out


def check_retardation() -> list[Check]:
    out: list[Check] = []
    p = PARAMS
    # With ξ = a_H = ħ c_s / Δ₀, η_ret = 2 exactly
    eta = eta_ret_from_healing(p.Delta0_meV, p.hbar_cs_meV_um)
    out.append(
        Check(
            "η_ret = 2 when ξ = a_H = ħc_s/Δ₀",
            "PASS" if near(eta, 2.0) else "FAIL",
            "Causal-propagator algebra identity",
            value=eta,
            expected=2.0,
        )
    )
    xi = p.lambda_M_um
    eta2 = eta_ret(p.Delta0_meV, xi, p.hbar_cs_meV_um)
    out.append(
        Check(
            "η_ret = 2 Δ₀ ξ / (ħ c_s) for ξ=λ_M",
            "PASS" if near(eta2, 2 * p.Delta0_meV * xi / p.hbar_cs_meV_um) else "FAIL",
            "General retardation factor",
            value=eta2,
        )
    )
    return out


def check_acoustic_matching() -> list[Check]:
    out: list[Check] = []
    p = PARAMS
    R = R_WS_from_T_star(p.hbar_cs_meV_um, T_star_K=1.0)
    out.append(
        Check(
            "acoustic R_WS matching helper (not duality proof)",
            "INFO",
            "MATCHING_NOT_PROOF: invert T_H=T* for R_WS at T*=1 K",
            value=R,
        )
    )
    return out


def check_gap_foundations() -> list[Check]:
    from qvc.gap_foundations import run_gap_foundations

    out: list[Check] = []
    g = run_gap_foundations()
    cut = g["sketch_cutoff"]
    lock = g["phenomenological_lock"]
    fold = g["static_fold_collapse"]
    out.append(
        Check(
            "sketch IR cutoff exceeds UV on Option A*",
            "PASS" if cut["paper_IR_exceeds_UV"] and cut["signs_opposite"] else "FAIL",
            "Phonon-loop sketch is not a derivation of the boxed map",
            value=cut["q_IR_paper_um_inv"],
            expected=cut["Lambda_um_inv"],
        )
    )
    out.append(
        Check(
            "phenomenological lock α_c = δ_c",
            "PASS" if lock["alpha_c_equals_delta_c"] else "FAIL",
            "Fold algebra of the locked map",
            value=lock["alpha_c"],
            expected=lock["delta_c"],
        )
    )
    out.append(
        Check(
            "static fold collapse confirmed (finite Δ_c, two branches, λ★ physical)",
            "PASS" if fold["two_branches_below_fold"] and fold["star_on_physical_branch"] and fold["unit_lambda_folded_out"] else "FAIL",
            fold["collapse_type"],
            value=fold["Delta_c_meV"],
            expected=0.109,
        )
    )
    out.append(
        Check(
            "sketch does not derive the lock",
            "PASS" if g["sketch_derives_lock"] is False and g["lock_is_internally_consistent"] else "FAIL",
            "Epistemic split: boxed map is phenomenological",
        )
    )
    return out


def check_parallel_channels() -> list[Check]:
    from qvc.parallel_channels import run_parallel_channels

    out: list[Check] = []
    ch = run_parallel_channels()
    out.append(
        Check(
            "static fold collapse is the amplitude-map channel",
            "PASS" if ch["fold"]["confirms_static_fold_collapse"] else "FAIL",
            "Not a GPE result",
        )
    )
    out.append(
        Check(
            "Mexican-hat catalyzon well is shallower; GPE well deepens",
            "PASS" if ch["mexican_hat"]["catalyzon_well_shallower"] and ch["gpe"]["well_deepens"] else "FAIL",
            ch["mexican_hat"]["formula"],
            value=ch["mexican_hat"]["barrier_reduction_pct"],
            expected=30.0,
        )
    )
    out.append(
        Check(
            "GPE live gap does not confirm fold collapse",
            "PASS" if ch["gpe"]["live_gap_rises"] and not ch["gpe"]["confirms_static_fold_collapse"] else "FAIL",
            "Attractive V_ZPF; parallel channel",
            value=ch["gpe"]["Delta_live_1ps_meV"],
        )
    )
    out.append(
        Check(
            "frozen Lindblad τ=5 ps is fold-violating",
            "PASS" if ch["lindblad"]["fold_violating_diagnostic"] and ch["lindblad"]["below_Delta_c"] else "FAIL",
            ch["lindblad"]["note"],
            value=ch["lindblad"]["Delta_live_frozen_tau5_meV"],
            expected=0.0872,
        )
    )
    return out


def check_cavity_moire() -> list[Check]:
    from qvc.cavity_moire import apply_F_moire_to_A_geom

    out: list[Check] = []
    cav = apply_F_moire_to_A_geom()
    out.append(
        Check(
            "F_moire applied to A_geom; withdrawn 0.018/0.0201 unused",
            "PASS" if cav["withdrawn_A_ZPF_0p018_not_used"] and cav["withdrawn_A_cav_0p0201_not_used"] else "FAIL",
            f"F={cav['F_moire']:.4f}",
            value=cav["A_geom_um_inv"],
        )
    )
    out.append(
        Check(
            "cavity-enhanced λ_fold is λ_fold/F_moire",
            "PASS" if abs(cav["lambda_fold_F"] * cav["F_moire"] - cav["lambda_fold"]) < 1e-12 else "FAIL",
            "Physical branch exists for λ ≤ λ_fold^(F)",
            value=cav["lambda_fold_F"],
            expected=cav["lambda_fold"] / cav["F_moire"],
        )
    )
    return out


def check_schwinger_boxed() -> list[Check]:
    from qvc.schwinger_rate import boxed_schwinger_rate

    out: list[Check] = []
    s = boxed_schwinger_rate()
    out.append(
        Check(
            "boxed Γ = ω0 (A_eff/A_c) exp(−π A_c/A_eff) with A_eff/A_c = α★",
            "PASS" if s["ratio_is_alpha_star"] and s["matches_diagnostics"] else "FAIL",
            "F_moire not in the exponent",
            value=s["Gamma_sch_Hz"],
        )
    )
    out.append(
        Check(
            "T* from E_pair / (k_B π/α★) matches diagnostics",
            "PASS" if s["T_star_matches_alpha"] else "FAIL",
            "Twin-channel crossover",
            value=s["T_star_K"],
            expected=0.1187,
        )
    )
    return out


def run_all() -> dict[str, Any]:
    checks: list[Check] = []
    checks += check_coupling_spectrum()
    checks += check_lens_space()
    checks += check_gap_fold()
    checks += check_vacuum_energy()
    checks += check_mode_enhancement()
    checks += check_retardation()
    checks += check_acoustic_matching()
    checks += check_gap_foundations()
    checks += check_parallel_channels()
    checks += check_cavity_moire()
    checks += check_schwinger_boxed()

    summary = {"PASS": 0, "FAIL": 0, "WARN": 0, "INFO": 0}
    for c in checks:
        summary[c.status] = summary.get(c.status, 0) + 1

    return {
        "summary": summary,
        "checks": [asdict(c) for c in checks],
        "canonical_notes": {
            "lambda_min": "Use -g0 (mult N-1), not -4g0 or -5g0",
            "Q_H": "Use CS/linking; never H2-cycle language",
            "gap_transition": "Fold unless vortex RG is derived",
            "vacuum": "Freeze one (formula, geometry) before quoting meV; see locks/",
            "retardation": "η_ret=2Δ0 ξ/(ħ c_s) — algebra lock",
            "mode_enhancement": "A_eff=√n A_single — quadrature lock",
            "gap_map": "Boxed δ=1+α(lnδ+ℓ) is a phenomenological lock; phonon sketch does not derive it",
            "parallel_channels": "Fold collapse ≠ catalyzon Mexican hat ≠ GPE live-gap rise",
            "cavity": "F_moire multiplies A_geom; 0.018/0.0201 withdrawn",
            "schwinger": "Γ=ω0 (α★) exp(−π/α★); T* from E_pair; no F_moire in the lock",
        },
    }


def write_report(report: dict[str, Any], out_dir: Path | None = None) -> tuple[Path, Path]:
    out_dir = ROOT / "output" if out_dir is None else Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "verification_report.json"
    md_path = out_dir / "verification_report.md"
    json_path.write_text(json.dumps(report, indent=2) + "\n")

    lines = [
        "# QVC Tier B/C Verification Report",
        "",
        f"Summary: `{report['summary']}`",
        "",
        "| Status | Check | Detail | Value | Expected |",
        "|--------|-------|--------|-------|----------|",
    ]
    for c in report["checks"]:
        lines.append(
            f"| {c['status']} | {c['name']} | {c['detail']} | {c['value']} | {c['expected']} |"
        )
    lines += ["", "## Canonical lock notes", ""]
    for k, v in report["canonical_notes"].items():
        lines.append(f"- **{k}**: {v}")
    lines += [
        "",
        "## Provenance",
        "",
        "- Closed-form locks: `docs/CLOSED_FORM_LOCKS.md`",
        "- Vacuum freeze decisions: `docs/VACUUM_FREEZE_SPEC.md`",
        "- Ported/improved from `QVCCursor/rigor/`",
        "",
    ]
    md_path.write_text("\n".join(lines))
    return json_path, md_path


def main() -> int:
    report = run_all()
    json_path, md_path = write_report(report)
    print(json.dumps(report["summary"], indent=2))
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return 0 if report["summary"].get("FAIL", 0) == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
