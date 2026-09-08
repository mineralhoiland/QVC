"""Driver for catalyzon steps 1–5 and L1 closures.

Writes JSON + figures. Arrays from the linearization are stripped from the
card (Phi, grids stay in memory only).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from qvc_l1.c1_metric import run_c1_metric
from qvc_l1.cs_eom import run_cs_eom
from qvc_l1.hcat1 import run_hcat1
from qvc_l1.hopf_integral import run_hopf_integral
from qvc_l1.htor1 import run_htor1
from qvc_l1.linearize import run_linearization
from qvc_l1.locks import build_l1_locks
from qvc_l1.spectral_cat import run_spectral_cat
from qvc_l1.track_b import run_track_b
from qvc_l1.villain_c7 import run_villain_c7


def _jsonify(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): _jsonify(v) for k, v in obj.items() if k not in _DROP}
    if isinstance(obj, (list, tuple)):
        return [_jsonify(v) for v in obj]
    if isinstance(obj, np.ndarray):
        if obj.size > 64:
            return {"ndarray_shape": list(obj.shape), "dtype": str(obj.dtype)}
        return obj.tolist()
    if isinstance(obj, (np.floating, np.integer)):
        return float(obj) if isinstance(obj, np.floating) else int(obj)
    if isinstance(obj, (float, int, str, bool)) or obj is None:
        return obj
    return str(obj)


_DROP = {"Phi", "psi0", "X", "Y", "theta"}


def _review(card: dict[str, Any]) -> dict[str, Any]:
    """Scientific-critical-thinking style claim grades (not a preprint rewrite)."""
    lin = card["linearize"]
    hfn1 = lin["H_FN1"]
    hcat = card["hcat1"]
    htor = card["htor1"]
    cs = card["cs_eom"]
    c1 = card["c1_metric"]
    hopf = card["hopf"]
    vill = card["villain_c7"]
    grades = {
        "H-FN1_soft_modes": {
            "grade": "computed_galerkin",
            "accepted": hfn1["accepted"],
            "evidence": hfn1["eigenvalues_meV"][:8],
            "limit": "Galerkin 2N basis, not continuum Goldstone count",
        },
        "G_ij_from_hessian": {
            "grade": "computed",
            "g0_eff_meV": lin["g0_eff_hessian_meV"],
            "rel_frobenius": lin["democratic_comparison"]["rel_frobenius_to_democratic"],
            "limit": "Amplitude-block Hessian off-diagonals; not torsion",
        },
        "H-CAT1": {
            "grade": "computed_ode",
            "accepted": hcat["accepted"],
            "ratio": hcat["ratio_tau_dark_over_bright"],
            "limit": "γ∝ω² ansatz on Galerkin eigenvalues",
        },
        "H-TOR1": {
            "grade": "ansatz_future",
            "accepted": False,
            "parked": True,
            "rel_frobenius": htor["democratic_comparison"]["rel_frobenius_to_democratic"],
            "limit": "2D proxy parked; future 3D Weitzenböck, not used for G_ij",
        },
        "CS_EOM": {
            "grade": "covariant_stepper",
            "a_in_stepper": cs.get("a_in_stepper", False),
            "flux_rel_err": cs["flux_rel_err"],
            "j_cancellation_ratio": cs.get("j_cancellation_ratio"),
            "limit": "background a in kinetic EOM; Hall remains algebra",
            "analytic_attachment_rms": cs.get("analytic_attachment_rms"),
        },
        "C1_quantum_metric": {
            "grade": "toy_band_physical_branch",
            "chern": c1["quantum_metric_toy"]["chern_numerical"],
            "kappa_geom_meV": c1["quantum_metric_toy"]["kappa_geom_meV"],
            "KT_hit_Tstar": c1["branch_frozen_toy"]["KT_hit_before_fold_at_Tstar"],
            "K_fold_Tstar": c1["branch_frozen_toy"]["K_fold_Tstar"],
            "limit": "QWZ C=1, not MATBG continuum; n0 ansatz",
        },
        "Track_B": {
            "grade": "matched_not_hubbard",
            "t_over_U": card["track_b"]["t_over_U"],
            "t_lambda_over_U": card["track_b"]["t_lambda_over_U"],
            "U_over_W": card["track_b"]["U_over_W_gap_scale"],
        },
        "Villain_C7": {
            "grade": "checkerboard_xy_scan",
            "eta_uv_fit": vill["mc_uv"]["eta_log_fit"],
            "eta_KT_fit": vill["mc_near_KT"]["eta_log_fit"],
            "limit": "Finite-L XY, not thermodynamic Villain",
        },
        "Hopf_integral": {
            "grade": "S3_whitehead_derived",
            "Q_S3": hopf["S3_whitehead"]["Q_Hopf_map"],
            "Q_stereo": hopf["Q_stereographic_R3"],
            "Q_spinor": hopf.get("Q_stereographic_spinor"),
            "Q1": hopf["Q_unit_hopfion"],
            "QN": hopf["Q_N_wound"],
            "lock": hopf["Q_H_lock_1_over_N"],
            "limit": "string-free FFT; lock remains CS/linking 1/N",
        },
    }
    return grades


def run_closure(*, out_dir: Path | None = None) -> dict[str, Any]:
    locks = build_l1_locks()
    lin = run_linearization(locks, n_grid=48)
    eigs = np.asarray(lin["hessian_eigs_meV"], dtype=float)
    hcat = run_hcat1(eigs, N=locks.N_layers)
    htor = run_htor1(lin["psi0"], lin["Phi"], lin["grid_dx"], N=locks.N_layers)
    spec = run_spectral_cat(
        eigs, Delta0_meV=locks.Delta0_meV, g0_star_meV=locks.g0_star_meV
    )
    cs = run_cs_eom(locks, n_grid=48, t_end_ps=0.15, dt_ps=0.003)
    c1 = run_c1_metric(locks)
    tb = run_track_b(locks)
    T_meV = locks.kT_star_meV
    kappa_uv = float(c1["kappa_UV_with_toy_geom_meV"])
    vill = run_villain_c7(kappa_meV=kappa_uv, T_eff_meV=T_meV, L=20)
    hopf = run_hopf_integral(n=36, N_layers=locks.N_layers)
    card = {
        "locks": {
            "g0_star_meV": locks.g0_star_meV,
            "Delta0_meV": locks.Delta0_meV,
            "E_vac_meV": locks.E_vac_meV,
            "lambda_star": locks.lambda_star,
            "Q_H": locks.Q_H,
            "N": locks.N_layers,
        },
        "linearize": lin,
        "hcat1": hcat,
        "htor1": htor,
        "spectral_cat": spec,
        "cs_eom": cs,
        "c1_metric": c1,
        "track_b": tb,
        "villain_c7": vill,
        "hopf": hopf,
    }
    card["review"] = _review(card)
    payload = _jsonify(card)

    out_dir = out_dir or Path(__file__).resolve().parents[1] / "notebooks" / "l1_closure" / "output"
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "summary.json").write_text(json.dumps(payload, indent=2, default=str))
    _figures(card, out_dir)
    _markdown(payload, out_dir)
    mirror = Path(__file__).resolve().parents[1] / "output" / "l1_closure"
    mirror.mkdir(parents=True, exist_ok=True)
    (mirror / "summary.json").write_text((out_dir / "summary.json").read_text())
    tex_figs = Path(__file__).resolve().parents[1] / "tex" / "QVC_arXiv_v4" / "closure_figs"
    if tex_figs.parent.is_dir():
        tex_figs.mkdir(parents=True, exist_ok=True)
        for png in out_dir.glob("fig_*.png"):
            (tex_figs / png.name).write_bytes(png.read_bytes())
    return payload


def _figures(card: dict[str, Any], out_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    eigs = np.asarray(card["linearize"]["hessian_eigs_meV"], dtype=float)
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.axhline(0.0, color="#888", lw=0.8)
    ax.plot(np.arange(1, eigs.size + 1), eigs, "o", color="#1f4e79")
    ax.set_xlabel("mode index")
    ax.set_ylabel(r"Hessian eigenvalue (meV)")
    ax.set_title("Galerkin TQGL Hessian about the hopfion")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_hessian_eigs.png", dpi=140)
    plt.close(fig)

    spec = card["spectral_cat"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(spec["omega_meV"], spec["A_cat"], color="#b85c38")
    ax.set_xlabel(r"$\omega$ (meV)")
    ax.set_ylabel(r"$A_{\mathrm{cat}}(\omega)$")
    ax.set_title("Catalyzon mode spectral function")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_A_cat.png", dpi=140)
    plt.close(fig)

    hcat = card["hcat1"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(hcat["t_ps"], hcat["weight_dark"], label="dark", color="#1f4e79")
    ax.plot(hcat["t_ps"], hcat["weight_bright"], label="bright", color="#b85c38")
    ax.set_xlabel("t (ps)")
    ax.set_ylabel("mode weight")
    ax.legend()
    ax.set_title("H-CAT1 dark vs bright envelopes")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_hcat1.png", dpi=140)
    plt.close(fig)

    vill = card["villain_c7"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(vill["mc_uv"]["r"], vill["mc_uv"]["g_r"], "o-", label="UV K", color="#1f4e79")
    ax.plot(
        vill["mc_near_KT"]["r"],
        vill["mc_near_KT"]["g_r"],
        "s-",
        label=r"K=2/π",
        color="#b85c38",
    )
    ax.set_xlabel("r (sites)")
    ax.set_ylabel(r"$g(r)$")
    ax.set_yscale("log")
    ax.legend()
    ax.set_title(r"C7 XY $g(r)$ (checkerboard Metropolis)")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_gr.png", dpi=140)
    plt.close(fig)

    scan = vill["eta_scan"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    Ks = [row["K"] for row in scan["rows"]]
    eta_fit = [row["eta_log_fit"] for row in scan["rows"]]
    eta_sw = [row["eta_spinwave"] for row in scan["rows"]]
    ax.plot(Ks, eta_fit, "o-", color="#1f4e79", label=r"MC $\eta$ (log fit)")
    ax.plot(Ks, eta_sw, "s--", color="#b85c38", label=r"spin-wave $1/(2\pi K)$")
    ax.axhline(0.25, color="#888", lw=0.8, label=r"$\eta=1/4$")
    ax.axvline(2.0 / np.pi, color="#888", ls=":", lw=0.8)
    ax.set_xlabel(r"$K$")
    ax.set_ylabel(r"$\eta$")
    ax.legend()
    ax.set_title(r"C7 $\eta(K)$ versus spin-wave")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_eta_K.png", dpi=140)
    plt.close(fig)

    br = card["c1_metric"]["branch_frozen_toy"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(br["alpha_over_ac"], br["kappa_meV"], color="#1f4e79")
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$\kappa(\alpha)$ (meV)")
    ax.set_title(r"C1 physical-branch stiffness (QWZ toy $\kappa_{\mathrm{geom}}$)")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_kappa_branch.png", dpi=140)
    plt.close(fig)

    G = np.asarray(card["htor1"]["G"], dtype=float)
    fig, ax = plt.subplots(figsize=(4.4, 3.8))
    im = ax.imshow(G, cmap="RdBu_r")
    fig.colorbar(im, ax=ax, fraction=0.046)
    ax.set_title("H-TOR1 proxy $G_{ij}$ (parked)")
    fig.tight_layout()
    fig.savefig(out_dir / "fig_htor1_G.png", dpi=140)
    plt.close(fig)


def _markdown(payload: dict[str, Any], out_dir: Path) -> None:
    rev = payload["review"]
    lines = ["# L1 closure card", "", "## Review grades", ""]
    for k, v in rev.items():
        lines.append(f"- **{k}**: {v}")
    lines += ["", "## Locks", json.dumps(payload["locks"], indent=2)]
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    run_closure()
