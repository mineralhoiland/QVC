"""Matplotlib figures for B1 G_ij and B2 Lindblad scans."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl-qvc-l1")
os.environ.setdefault("MPLBACKEND", "Agg")

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

# Wong colorblind-safe palette
C = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#000000"]


def _save(fig: plt.Figure, path: Path) -> Path:
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return path


def _stick_spectrum(ax: plt.Axes, eigs: np.ndarray, g0: float, color: str, label: str) -> None:
    eigs = np.asarray(eigs, dtype=float)
    ax.vlines(eigs, 0.0, 1.0, colors=color, linewidth=2.0, label=label)
    ax.plot(eigs, np.ones_like(eigs), "o", color=color, ms=6)
    if abs(g0) > 0:
        ax.axvline(-g0, color=C[6], ls=":", lw=1.2, alpha=0.85)
        ax.axvline(4.0 * g0, color=C[6], ls="--", lw=1.2, alpha=0.85)


def fig_b1_heatmaps(b1: dict[str, Any], out: Path) -> Path:
    ws = b1["ws_gaussian"]
    ring = b1["ring_fourier"]
    G_ws = np.asarray(ws["G"], dtype=float)
    G_r = np.asarray(ring["G"], dtype=float)
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.8), constrained_layout=True)
    vmax = max(np.max(np.abs(G_ws)), np.max(np.abs(G_r)), 1e-18)
    for ax, G, title in (
        (axes[0], G_ws, r"WS Gaussian $G_{ij}$ (K$_0$)"),
        (axes[1], G_r, r"ring Wannier $G_{ij}$ (K$_0$)"),
    ):
        im = ax.imshow(G, cmap="coolwarm", vmin=-vmax, vmax=vmax, origin="upper")
        n = G.shape[0]
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xlabel(r"$j$")
        ax.set_ylabel(r"$i$")
        ax.set_title(title)
        ax.grid(False)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(r"Computed kernel $G$  (not $g_0(J-I)$; $G_{ii}=0$)", fontsize=11)
    return _save(fig, out / "fig_b1_G_heatmap.png")


def fig_b1_spectrum(b1: dict[str, Any], out: Path) -> Path:
    ws = b1["ws_gaussian"]
    ring = b1["ring_fourier"]
    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), constrained_layout=True, sharey=True)
    for ax, blk, title in (
        (axes[0], ws, "WS Gaussian"),
        (axes[1], ring, "ring Wannier"),
    ):
        eigs = np.asarray(blk["eigenvalues"], dtype=float)
        g0 = float(blk["g0_eff_kernel"])
        _stick_spectrum(ax, eigs, g0, C[0], r"computed Spec$(G)$")
        ax.set_xlabel(r"$\lambda$ (kernel units)")
        ax.set_title(title)
        ax.set_ylim(0.0, 1.35)
        ax.set_yticks([])
        ax.legend(loc="upper left", fontsize=7)
        ax.text(
            0.02,
            0.08,
            rf"$g_{{0,\mathrm{{eff}}}}={g0:.4g}$"
            "\n"
            rf"rel Fro={blk['rel_frobenius_to_democratic']:.3f}",
            transform=ax.transAxes,
            fontsize=8,
            va="bottom",
        )
    axes[0].set_ylabel("weight")
    fig.suptitle(
        r"Stick spectrum vs democratic $(-\,g_{0,\mathrm{eff}},\,4 g_{0,\mathrm{eff}})$  [algebra lock, not inserted]",
        fontsize=11,
    )
    return _save(fig, out / "fig_b1_spectrum.png")


def fig_b1_ring_modes(b1: dict[str, Any], out: Path) -> Path:
    ring = b1["ring_fourier"]
    th = np.asarray(ring["profiles"]["theta"])
    freal = np.asarray(ring["profiles"].get("fourier_real", ring["profiles"]["fourier_amp"]))
    wamp = np.asarray(ring["profiles"]["wannier_amp"])
    w0 = np.asarray(ring["wannier_abs_0"])
    x = np.asarray(ring["x_um"]) * 1e3
    y = np.asarray(ring["y_um"]) * 1e3

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.5), constrained_layout=True)
    ax = axes[0]
    for n in range(freal.shape[0]):
        ax.plot(th, freal[n], color=C[n % len(C)], lw=1.6, label=fr"$n={n}$")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$\mathrm{Re}\,\Phi_n$ on $\rho=R$")
    ax.set_title("Fourier harmonics")
    ax.legend(ncol=2, fontsize=7)

    ax = axes[1]
    for j in range(wamp.shape[0]):
        ax.plot(th, wamp[j], color=C[j % len(C)], lw=1.6, label=fr"$w_{j}$")
    ax.set_xlabel(r"$\theta$")
    ax.set_ylabel(r"$|w_j|$ on $\rho=R$")
    ax.set_title("Wannier packets")
    ax.legend(ncol=2, fontsize=7)

    ax = axes[2]
    im = ax.imshow(
        w0,
        origin="lower",
        extent=(x[0], x[-1], y[0], y[-1]),
        cmap="viridis",
        aspect="equal",
    )
    ax.set_xlabel(r"$x$ (nm)")
    ax.set_ylabel(r"$y$ (nm)")
    ax.set_title(r"$|w_0|$ (torus packet)")
    ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle("Ring modes from hopfion envelope × angular harmonics $n=0..4$", fontsize=11)
    return _save(fig, out / "fig_b1_ring_modes.png")


def fig_b2_delta(b2: dict[str, Any], out: Path, *, window: str = "1ps") -> Path:
    key = "runs_1ps" if window == "1ps" else "runs_20ps"
    runs = b2[key]
    fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
    order = ["coherent"] + [k for k in runs if k.startswith("tau_")]
    for i, name in enumerate(order):
        r = runs[name]
        ls = "-" if name == "coherent" else "--"
        ax.plot(
            r["t_ps"],
            r["Delta_live_meV"],
            color=C[i % len(C)],
            ls=ls,
            lw=1.8,
            label=name.replace("tau_", r"$\tau_{\mathrm{dec}}=$") + (" ps" if name.startswith("tau_") else ""),
        )
    ax.axhline(0.600, color=C[6], ls=":", alpha=0.6, label=r"$\Delta_0$")
    ax.set_xlabel(r"$t$ (ps)")
    ax.set_ylabel(r"$\Delta_{\mathrm{live}}$ (meV)")
    ax.set_title(rf"Live gap vs $\tau_{{\mathrm{{dec}}}}$  ({window} window)")
    ax.legend(loc="best", fontsize=7)
    return _save(fig, out / f"fig_b2_delta_live_{window}.png")


def fig_b2_g1(b2: dict[str, Any], out: Path, *, window: str = "1ps") -> Path:
    key = "runs_1ps" if window == "1ps" else "runs_20ps"
    runs = b2[key]
    fig, ax = plt.subplots(figsize=(6.4, 3.8), constrained_layout=True)
    order = ["coherent"] + [k for k in runs if k.startswith("tau_")]
    for i, name in enumerate(order):
        r = runs[name]
        ls = "-" if name == "coherent" else "--"
        ax.plot(
            r["t_ps"],
            r["g1"],
            color=C[i % len(C)],
            ls=ls,
            lw=1.8,
            label=name.replace("tau_", r"$\tau_{\mathrm{dec}}=$") + (" ps" if name.startswith("tau_") else ""),
        )
    ax.set_xlabel(r"$t$ (ps)")
    ax.set_ylabel(r"$g^{(1)}$")
    ax.set_ylim(0.0, 1.02)
    ax.set_title(rf"Coherence $g^{{(1)}}$ vs $\tau_{{\mathrm{{dec}}}}$  ({window} window)")
    ax.legend(loc="best", fontsize=7)
    return _save(fig, out / f"fig_b2_g1_{window}.png")


def write_all_figures(bundle: dict[str, Any], out_dir: Path) -> list[str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    if "B1" in bundle:
        paths.extend(
            [
                fig_b1_heatmaps(bundle["B1"], out_dir),
                fig_b1_spectrum(bundle["B1"], out_dir),
                fig_b1_ring_modes(bundle["B1"], out_dir),
            ]
        )
    if "B2" in bundle:
        paths.extend(
            [
                fig_b2_delta(bundle["B2"], out_dir, window="1ps"),
                fig_b2_g1(bundle["B2"], out_dir, window="1ps"),
                fig_b2_delta(bundle["B2"], out_dir, window="20ps"),
                fig_b2_g1(bundle["B2"], out_dir, window="20ps"),
            ]
        )
    return [str(p) for p in paths]
