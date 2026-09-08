"""Matplotlib figures for the BKT C1–C6 suite."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-qvc-bkt")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update(
    {
        "figure.dpi": 120,
        "savefig.dpi": 160,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "axes.grid": True,
        "grid.alpha": 0.35,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#000000"]


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_c1_kappa(c1: dict[str, Any], out: Path) -> Path:
    p = c1["primary"]
    rows = p["rows"]
    x = np.array([r["alpha_over_ac"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(x, [r["kappa_conv_meV"] for r in rows], color=C[0], label=r"$\kappa_{\mathrm{conv}}$ (derived $|$ $n_0$ ansatz)")
    ax.plot(x, [r["kappa_geom_meV"] for r in rows], color=C[2], ls="--", label=r"$\kappa_{\mathrm{geom}}$ (Peotta–Törmä ansatz)")
    ax.plot(x, [r["kappa_meV"] for r in rows], color=C[6], lw=2.0, label=r"$\kappa=\kappa_{\mathrm{conv}}+\kappa_{\mathrm{geom}}$")
    conv = c1["variants"]["conv_only"]["rows"]
    ax.plot(
        [r["alpha_over_ac"] for r in conv],
        [r["kappa_meV"] for r in conv],
        color=C[1],
        ls=":",
        label="conventional only",
    )
    ax.axvline(0.90, color=C[3], ls="--", alpha=0.7, label=r"$\lambda^\star$")
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$\kappa$ (meV)")
    ax.set_title("C1  phase stiffness on the physical branch")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c1_kappa_alpha.png")


def fig_c1_K(c1: dict[str, Any], c3: dict[str, Any], out: Path) -> Path:
    p = c1["primary"]
    rows = p["rows"]
    x = np.array([r["alpha_over_ac"] for r in rows])
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(x, [r["K_Tstar"] for r in rows], color=C[0], lw=2, label=r"$K(\alpha)$ at $T^*$ (primary)")
    conv = c1["variants"]["conv_only"]["rows"]
    ax.plot(
        [r["alpha_over_ac"] for r in conv],
        [r["K_Tstar"] for r in conv],
        color=C[1],
        ls="--",
        label=r"$K$ conventional only",
    )
    ax.axhline(p["K_c_bosonic"], color=C[6], ls=":", label=r"$K_c=2/\pi$ (truncated)")
    ax.axhline(c3["scan_Tstar_anyon"]["K_c"], color=C[4], ls=":", label=r"$K_c^{\mathrm{anyon}}$ (CS restored)")
    ax.axvline(0.90, color=C[3], ls="--", alpha=0.7)
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$K=\kappa/T_{\mathrm{eff}}$")
    ax.set_title("C1/C3  Kosterlitz $K$ vs frozen $\\alpha$")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c1_K_alpha.png")


def fig_c2_profile(c2: dict[str, Any], out: Path) -> Path:
    r = np.asarray(c2["profile"]["r_um"])
    n = np.asarray(c2["profile"]["n"])
    xi = c2["radial"]["xi_um"]
    fig, ax = plt.subplots(figsize=(6.2, 3.6))
    ax.plot(r / xi, n, color=C[0], lw=2, label=r"$|\psi|^2$ (winding-1 GP vortex)")
    ax.axvline(c2["core"]["a_c_half_um"] / xi, color=C[1], ls="--", label=r"computed $a_{1/2}/\xi$")
    ax.axvline(1.0, color=C[2], ls=":", label=r"healing $a_c=\xi$ (ansatz cutoff)")
    ax.set_xlabel(r"$r/\xi$")
    ax.set_ylabel(r"$n(r)$")
    ax.set_xlim(0, 12)
    ax.set_title("C2  2D phase vortex (not a hopfion)")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c2_vortex_profile.png")


def fig_c3_flows(c3: dict[str, Any], out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    colors = {"uv": C[0], "star": C[3], "nearfold": C[1]}
    for st in c3["flows_Tstar"]:
        fl = st["flow"]
        ax.plot(fl["ell"], fl["y"], color=colors.get(st["name"], C[6]), lw=2, label=fr"{st['name']} $\alpha/\alpha_c={st['alpha_over_ac']:.2f}$")
    ax.axhline(0.4, color=C[6], ls=":", alpha=0.6, label=r"$y_{\mathrm{cut}}$")
    ax.set_xlabel(r"KT scale $\ell=\ln(L/a_c)$")
    ax.set_ylabel(r"fugacity $y$")
    ax.set_yscale("log")
    ax.set_title("C3  frozen-$\\alpha$ Kosterlitz flow at $T^*$")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c3_kt_flows.png")


def fig_c3_competing(c3: dict[str, Any], out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    sc = c3["scan_Tstar_primary"]
    x = [r["alpha_over_ac"] for r in sc["rows"]]
    ax.plot(x, [r["K_bare"] for r in sc["rows"]], color=C[0], lw=2, label=r"primary $K_{\mathrm{bare}}(T^*)$")
    scc = c3["scan_Tstar_conv_only"]
    ax.plot(
        [r["alpha_over_ac"] for r in scc["rows"]],
        [r["K_bare"] for r in scc["rows"]],
        color=C[1],
        ls="--",
        label="conventional only",
    )
    ax.axhline(sc["K_c"], color=C[6], ls=":", label=r"$K_c=2/\pi$")
    ax.axhline(c3["scan_Tstar_anyon"]["K_c"], color=C[4], ls=":", label=r"$K_c^{\mathrm{anyon}}$")
    ax.axvline(0.90, color=C[3], ls="--", alpha=0.7, label=r"$\lambda^\star$")
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$K$")
    ax.set_title("C3  competing-scale scan (frozen $\\alpha$, Option A*)")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c3_competing_scales.png")


def fig_c4_xi(c4: dict[str, Any], out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    rows = c4["primary_Tstar"]["rows"]
    x = [r["alpha_over_ac"] for r in rows]
    ax.plot(x, [r["Delta_fold_meV"] for r in rows], color=C[0], lw=2, label=r"locked fold $\Delta_+(\alpha)$")
    ax.axhline(c4["vs_locked_Delta_c"]["Delta_c_locked_meV"], color=C[6], ls=":", label=r"locked $\Delta_c$")
    dv = np.array([r["Delta_vortex_meV"] for r in rows], dtype=float)
    if np.any(np.isfinite(dv)):
        ax.plot(x, dv, color=C[1], ls="--", label=r"vortex $\hbar c_s/\xi_{\mathrm{KT}}$ (if unbound)")
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$\Delta$ (meV)")
    ax.set_title("C4  vortex-implied gap vs locked fold")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c4_xi_delta.png")


def fig_c5_compare(c5: dict[str, Any], out: Path) -> Path:
    cur = c5["curves"]
    tA = np.asarray(cur["t_A"])
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.plot(tA, cur["delta_phys"], color=C[0], lw=2, label=r"locked $\delta_+(t_A)$")
    ax.plot(tA, cur["delta_sqrt_local"], color=C[2], ls="--", label=r"local $\delta_c+B\sqrt{t_A}$")
    ax.plot(tA, cur["delta_bkt_toy"], color=C[1], ls=":", label=r"toy $e^{-C/\sqrt{t_A}}$ (not a QVC result)")
    ax.set_xscale("log")
    ax.set_xlabel(r"$t_A=(\alpha_c-\alpha)/\alpha_c$")
    ax.set_ylabel(r"$\delta$")
    ax.set_title("C5  fold map vs essential-singularity shape")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c5_model_compare.png")


def fig_c6_ahns(c6: dict[str, Any], out: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    g = c6["scan_GPE_rows"]
    e = c6["scan_Eliashberg_rows"]
    ax.plot(
        [r["alpha_over_ac"] for r in g],
        [r["ratio"] for r in g],
        color=C[0],
        lw=2,
        label=r"GPE $g_{\mathrm{ph}}$ estimate",
    )
    ax.plot(
        [r["alpha_over_ac"] for r in e],
        [r["ratio"] for r in e],
        color=C[1],
        ls="--",
        label=r"Eliashberg $\lambda_{\mathrm{ep}}=1$ at $T^*$",
    )
    ax.axhline(1.0, color=C[6], ls=":", label=r"AHNS: $\hbar/\tau_v=\pi\kappa$")
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$(\hbar/\tau_v)/(\pi\kappa)$")
    ax.set_yscale("log")
    ax.set_title("C6  overdamped criterion along the physical branch")
    ax.legend(loc="best")
    return _save(fig, out / "fig_c6_ahns.png")


def write_all_figures(bundle: dict[str, Any], out_dir: Path) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    c1, c2, c3, c4, c5, c6 = (bundle[k] for k in ("C1", "C2", "C3", "C4", "C5", "C6"))
    paths = [
        fig_c1_kappa(c1, out_dir),
        fig_c1_K(c1, c3, out_dir),
        fig_c2_profile(c2, out_dir),
        fig_c3_flows(c3, out_dir),
        fig_c3_competing(c3, out_dir),
        fig_c4_xi(c4, out_dir),
        fig_c5_compare(c5, out_dir),
        fig_c6_ahns(c6, out_dir),
    ]
    return [str(p) for p in paths]
