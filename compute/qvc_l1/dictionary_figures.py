"""Matplotlib figures for the Track A L1 dictionary."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-qvc-l1-dict")
os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

mpl.rcParams.update(
    {
        "figure.dpi": 120,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "legend.fontsize": 8,
        "axes.grid": True,
        "grid.alpha": 0.35,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "serif",
    }
)

C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#000000"]


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def fig_lengths(lengths: dict[str, Any], out: Path) -> Path:
    """Three lengths and the three η_ret evaluations."""
    names = ["a_H", r"$\xi_{\mathrm{ph}}$", r"$\xi_{181}$"]
    keys = ["a_H", "xi_ph", "xi_181"]
    um = [lengths["lengths"][k]["value_um"] for k in keys]
    eta = [lengths["lengths"][k]["eta_ret"] for k in keys]
    colors = [C[0], C[2], C[1]]

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.5))
    ax = axes[0]
    bars = ax.bar(names, um, color=colors, width=0.62)
    ax.set_ylabel(r"length (µm)")
    ax.set_title("Three lengths (do not mix)")
    for b, v in zip(bars, um):
        ax.text(
            b.get_x() + b.get_width() / 2.0,
            v + 0.004,
            f"{v:.4f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_ylim(0, max(um) * 1.22)

    ax = axes[1]
    bars = ax.bar(names, eta, color=colors, width=0.62)
    ax.axhline(2.0, color=C[6], ls=":", lw=1.2, label=r"identity $\eta_{\mathrm{ret}}=2$")
    ax.axhline(3.1, color=C[1], ls="--", lw=1.0, alpha=0.7, label=r"alternate $3.1$ (not the lock)")
    ax.set_ylabel(r"$\eta_{\mathrm{ret}}=2\Delta_0\xi/(\hbar c_s)$")
    ax.set_title(r"$\eta_{\mathrm{ret}}$ is not unique")
    ax.legend(loc="upper left", frameon=False)
    for b, v in zip(bars, eta):
        ax.text(
            b.get_x() + b.get_width() / 2.0,
            v + 0.08,
            f"{v:.3f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    ax.set_ylim(0, max(eta) * 1.35)
    fig.suptitle(
        r"Length dictionary: $a_H=r$, $\xi_{\mathrm{ph}}=\hbar c_s/\Delta_0$, $\xi_{181}$ alternate",
        fontsize=11,
    )
    return _save(fig, out / "fig_lengths.png")


def fig_mexican_hat(hat: dict[str, Any], out: Path) -> Path:
    rho = np.asarray(hat["rho"])
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(rho, hat["V_bare"], color=C[0], lw=2.0, label=r"bare $\alpha_{\mathrm{GL}}=-\Delta_0$")
    ax.plot(
        rho,
        hat["V_gpe_sign"],
        color=C[3],
        lw=1.6,
        ls="--",
        label=r"GPE sign $-\lambda^\star E_{\mathrm{vac}}\rho^2$ (deepens)",
    )
    ax.plot(
        rho,
        hat["V_catalyzon_dressed"],
        color=C[1],
        lw=2.0,
        label=r"catalyzon $|\alpha|\to\Delta_{\mathrm{cat}}$ (parallel)",
    )
    ax.axhline(0.0, color=C[6], lw=0.7)
    red = hat["catalyzon_dressed"]["barrier_reduction_pct"]
    ax.fill_between(
        rho,
        hat["V_catalyzon_dressed"],
        hat["V_bare"],
        where=(np.asarray(hat["V_bare"]) < np.asarray(hat["V_catalyzon_dressed"])),
        color=C[1],
        alpha=0.18,
        label=rf"barrier reduction {red:.1f}% (replaces 42%)",
    )
    ax.set_xlabel(r"$\rho=|\psi|$  (GPE axis)")
    ax.set_ylabel(r"$V(\rho)$ (meV)")
    ax.set_title("Mexican hat on the hopfion $n_{\\mathrm{ref}}$ axis")
    ax.legend(loc="best", frameon=False)
    return _save(fig, out / "fig_mexican_hat.png")


def fig_spec_G(cat: dict[str, Any], out: Path) -> Path:
    star = cat["at_lambda_star"]
    eigs = np.asarray(star["spec"]["eigs_numerical_meV"])
    g0 = star["g0_meV"]
    n = cat["N_layers"]
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    ax.hlines(eigs, xmin=0.2, xmax=0.8, colors=C[0], lw=2.4)
    ax.scatter(np.full_like(eigs, 0.5), eigs, c=C[0], s=36, zorder=3)
    ax.axhline(-g0, color=C[1], ls="--", lw=1.2, label=rf"$\lambda_{{\min}}=-g_0^\star={-g0:.4f}$ meV")
    ax.axhline(
        (n - 1) * g0,
        color=C[2],
        ls=":",
        lw=1.2,
        label=rf"$\lambda_{{\max}}=(N-1)g_0^\star={(n-1)*g0:.4f}$ meV",
    )
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.set_ylabel(r"$\mathrm{Spec}(G)$ (meV)")
    ax.set_title(rf"Democratic Spec$(G)$ at $g_0^\star={g0:.4f}$ meV  ($N={n}$)")
    ax.legend(loc="best", frameon=False)
    ax.text(
        0.02,
        0.04,
        r"dark $\times(N-1)$;  $g_0\sim 0.21$ is not derived",
        transform=ax.transAxes,
        fontsize=8,
        color=C[6],
    )
    return _save(fig, out / "fig_spec_G.png")


def write_all_figures(bundle: dict[str, Any], out: Path) -> list[str]:
    out.mkdir(parents=True, exist_ok=True)
    paths = [
        fig_lengths(bundle["lengths"], out),
        fig_mexican_hat(bundle["mexican_hat"], out),
        fig_spec_G(bundle["catalyzon"], out),
    ]
    return [str(p) for p in paths]
