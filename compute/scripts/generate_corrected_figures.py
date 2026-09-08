#!/usr/bin/env python3
"""
Regenerate QVC / Supplemental figures with corrected mathematics.

Locked content:
  - Democratic Spec(G): λ_max=(N-1)g0, λ_min=-g0 (mult. N-1)
  - Gap map: saddle-node fold (square-root); BKT only contingent
  - Q_H = 1/N via CS/linking (H2=0)

Outputs PDFs (+ PNG previews) into tex/figures/.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from qvc.coupling_matrix import (  # noqa: E402
    catalyzon_dark_basis,
    democratic_coupling_matrix,
    democratic_spectrum,
)
from qvc.gap_fold import ell_from_delta_c, fold_condition, gap_rhs  # noqa: E402
from qvc.lens_chern_simons import fractional_hall_charge  # noqa: E402
from qvc.params import PARAMS  # noqa: E402

OUTDIR = ROOT / "tex" / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "legend.fontsize": 8,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.linewidth": 0.8,
    }
)

# Colorblind-safe (avoid purple-on-white AI cliché; use teal / coral / slate)
C_TEAL = "#0D7377"
C_CORAL = "#C44536"
C_SLATE = "#2C3E50"
C_GOLD = "#B8860B"
C_BLUE = "#1B4F72"
C_GRAY = "#7F8C8D"
C_DARK = "#1A1A1A"


def save(fig: plt.Figure, stem: str) -> None:
    pdf = OUTDIR / f"{stem}.pdf"
    png = OUTDIR / f"{stem}.png"
    fig.savefig(pdf)
    fig.savefig(png)
    plt.close(fig)
    print(f"  wrote {pdf.name}")


# ---------------------------------------------------------------------------
# 1. Gap fold bifurcation (replaces misleading BKT-shaped external plot usage)
# ---------------------------------------------------------------------------
def fig_gap_fold_collapse() -> None:
    """Physical branch of δ=1+α(lnδ+ℓ) ending in a saddle-node fold."""
    delta_c_target = 0.182
    ell = ell_from_delta_c(delta_c_target)
    delta_c, alpha_c = fold_condition(ell)

    alphas = np.linspace(1e-4, alpha_c * 0.999, 120)
    delta_phys = []
    delta_unst = []
    for a in alphas:
        # find roots of f(δ)=δ-(1+a(lnδ+ell))
        ds = np.linspace(1e-4, 1.0, 4000)
        f = ds - (1.0 + a * (np.log(ds) + ell))
        roots = []
        for i in range(len(f) - 1):
            if f[i] == 0:
                roots.append(ds[i])
            elif f[i] * f[i + 1] < 0:
                # linear interpolate
                t = -f[i] / (f[i + 1] - f[i])
                roots.append(ds[i] + t * (ds[i + 1] - ds[i]))
        roots = sorted(set(round(r, 8) for r in roots), reverse=True)
        delta_phys.append(roots[0] if roots else np.nan)
        delta_unst.append(roots[1] if len(roots) > 1 else np.nan)

    a_norm = alphas / alpha_c
    fig, ax = plt.subplots(figsize=(5.2, 3.6))
    ax.plot(a_norm, delta_phys, color=C_TEAL, lw=2.0, label=r"physical branch $\delta_+$")
    ax.plot(
        a_norm,
        delta_unst,
        color=C_CORAL,
        lw=1.5,
        ls="--",
        label=r"unstable branch $\delta_-$",
    )
    ax.plot(1.0, delta_c, "o", color=C_CORAL, ms=7, zorder=5, label=r"fold $(\alpha_c,\delta_c)$")
    # operating point ~ A/A_c = 0.45 from preprint narrative
    # interpolate physical δ at 0.45
    op = 0.45
    dop = np.interp(op, a_norm, delta_phys)
    ax.plot(op, dop, "s", color=C_GOLD, ms=7, label="MATBG operating point")

    # contingent BKT sketch near criticality (not a solution of the fold map)
    a_bkt = np.linspace(0.85, 0.999, 80)
    bkt = np.exp(-1.2 / np.sqrt(np.maximum(1e-6, 1 - a_bkt)))
    ax.plot(a_bkt, bkt, color=C_GRAY, ls=":", lw=1.2, label="contingent BKT form (vortex RG)")

    ax.axvline(1.0, color=C_CORAL, ls=":", lw=0.9, alpha=0.7)
    ax.set_xlabel(r"$A_{\mathrm{ZPF}} / A_c$")
    ax.set_ylabel(r"$\delta = \Delta / \Delta_{\mathrm{ex}}$")
    ax.set_xlim(0, 1.08)
    ax.set_ylim(0, 1.05)
    ax.set_title("Saddle-node fold of the self-consistent gap map")
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.text(
        0.02,
        -0.02,
        rf"$\delta_c\approx{delta_c:.3f}$, $\ell\approx{ell:.3f}$; near fold: "
        r"$\delta-\delta_c\sim\pm\sqrt{\alpha_c-\alpha}$",
        fontsize=7,
        color=C_GRAY,
    )
    save(fig, "gap_fold_collapse")


# ---------------------------------------------------------------------------
# 2. Five-mode coupling matrix + correct spectrum
# ---------------------------------------------------------------------------
def fig_five_mode_coupling_matrix() -> None:
    N, g0 = 5, 1.0
    G = democratic_coupling_matrix(N, g0)
    eigs = np.sort(np.linalg.eigvalsh(G))
    spec = democratic_spectrum(N, g0)
    labels = ["Bog", "Jos", "Plas", "Hopf", "Ac"]

    fig = plt.figure(figsize=(8.5, 3.6))
    gs = GridSpec(1, 2, width_ratios=[1.1, 1.0], wspace=0.35)

    ax0 = fig.add_subplot(gs[0])
    im = ax0.imshow(G, cmap="RdBu_r", vmin=-g0, vmax=(N - 1) * g0)
    ax0.set_xticks(range(N))
    ax0.set_yticks(range(N))
    ax0.set_xticklabels(labels)
    ax0.set_yticklabels(labels)
    for i in range(N):
        for j in range(N):
            ax0.text(j, i, f"{G[i, j]:.0f}", ha="center", va="center", fontsize=8, color=C_DARK)
    ax0.set_title(r"Democratic $G=g_0(J-I)$")
    cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
    cbar.set_label(r"$G_{ij}/g_0$")

    ax1 = fig.add_subplot(gs[1])
    colors = [C_CORAL] * (N - 1) + [C_TEAL]
    # plot sorted: dark then bright
    x = np.arange(1, N + 1)
    ax1.axhline(0, color=C_GRAY, lw=0.8)
    ax1.scatter(x[:-1], eigs[:-1], c=C_CORAL, s=60, zorder=3, label=rf"$\lambda_{{\min}}=-g_0$ (×{N-1})")
    ax1.scatter([x[-1]], [eigs[-1]], c=C_TEAL, s=70, zorder=3, label=rf"$\lambda_{{\max}}=(N-1)g_0$")
    for xi, yi in zip(x, eigs):
        ax1.plot([xi, xi], [0, yi], color=C_GRAY, lw=1.0, zorder=1)
    ax1.set_xticks(x)
    ax1.set_xticklabels([f"{i}" for i in x])
    ax1.set_xlabel("mode index (sorted)")
    ax1.set_ylabel(r"$\lambda / g_0$")
    ax1.set_title("Exact spectrum (catalyzon = dark subspace)")
    ax1.set_ylim(-1.5, 4.5)
    ax1.legend(frameon=False, loc="upper left")
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    fig.suptitle(
        rf"Locked: $\lambda_{{\min}}={spec['lambda_min']:.0f}\,g_0$ (mult.\ {N-1}), "
        rf"not $-4g_0$ or $-5g_0$",
        fontsize=9,
        y=1.02,
    )
    save(fig, "five_mode_coupling_matrix")


# ---------------------------------------------------------------------------
# 3. Catalyzon formation: dispersions + eigenvalue splitting (corrected)
# ---------------------------------------------------------------------------
def fig_catalyzon_formation_3d() -> None:
    fig = plt.figure(figsize=(9.5, 4.0))
    gs = GridSpec(1, 2, wspace=0.28)

    # Left: schematic dispersions meeting at resonance
    ax0 = fig.add_subplot(gs[0], projection="3d")
    k = np.linspace(0, 0.12, 80)
    k_star = 0.056
    omega0 = 2.0  # arbitrary units; mark resonance plane
    # five schematic branches approaching ω0 at k*
    branches = [
        (omega0 + 0.15 * (k - k_star) ** 2 * 80, C_BLUE, "Bogoliubov amp."),
        (omega0 + 0.08 * np.cos(40 * (k - k_star)), C_TEAL, "Josephson"),
        (omega0 + 0.12 * (k**0.5 - k_star**0.5), C_GOLD, "Plasmon"),
        (omega0 + 0.1 * np.sin(25 * (k - k_star)), "#6C3483", "Hopfion"),
        (PARAMS.hbar_cs_meV_um * k / 0.07 * 1.2, C_CORAL, "Acoustic"),  # linear-ish
    ]
    # normalize acoustic to hit omega0 at k*
    branches[-1] = (omega0 * (k / k_star), C_CORAL, "Acoustic")
    for y, c, lab in branches:
        ax0.plot(k, y, np.zeros_like(k), color=c, lw=1.5, label=lab)
    # resonance plane
    Kk, Ww = np.meshgrid(np.linspace(0, 0.12, 2), np.linspace(omega0, omega0, 2))
    ax0.plot_surface(Kk, Ww, np.zeros_like(Kk), alpha=0.15, color=C_GRAY)
    ax0.scatter([k_star], [omega0], [0], c=C_CORAL, s=40, zorder=5)
    ax0.set_xlabel(r"$k$ (nm$^{-1}$)")
    ax0.set_ylabel(r"$\omega$ (arb.)")
    ax0.set_zlabel("")
    ax0.set_title(r"Resonance at $k^*\approx 0.056\,\mathrm{nm}^{-1}$")
    ax0.view_init(elev=22, azim=-60)
    ax0.legend(loc="upper left", fontsize=6)

    # Right: eigenvalues vs coupling — dark multiplet at -g0, bright at 4g0
    ax1 = fig.add_subplot(gs[1])
    g = np.linspace(0, 1.0, 100)
    ax1.fill_between(g, -g, -g, color=C_CORAL)  # placeholder
    for _ in range(4):
        ax1.plot(g, -g, color=C_CORAL, lw=1.8, alpha=0.85)
    ax1.plot(g, 4 * g, color=C_TEAL, lw=2.0, label=r"bright $\lambda_{\max}=4g_0$")
    ax1.plot(g, -g, color=C_CORAL, lw=2.0, label=r"dark $\lambda_{\min}=-g_0$ (×4)")
    g_op = 0.21  # meV scale from preprint estimate / g0~0.21
    ax1.axvline(g_op, color=C_GOLD, ls="--", lw=1.0, label=r"$g_0\sim 0.21\,\mathrm{meV}$")
    ax1.axhline(0, color=C_GRAY, lw=0.7)
    ax1.set_xlabel(r"$g_0$ (meV)")
    ax1.set_ylabel(r"$\lambda$ (meV)")
    ax1.set_title("Eigenvalue splitting (corrected spectrum)")
    ax1.legend(frameon=False, fontsize=7)
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-1.2, 4.5)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    fig.suptitle("Catalyzon formation: resonance + dark-subspace hybridization", fontsize=11)
    save(fig, "catalyzon_formation_3d")


# ---------------------------------------------------------------------------
# 4. Dark-mode / Helmert basis visualization
# ---------------------------------------------------------------------------
def fig_catalyzon_dark_modes() -> None:
    N = 5
    B = catalyzon_dark_basis(N)
    fig, axes = plt.subplots(1, N - 1, figsize=(9, 2.4), sharey=True)
    layers = np.arange(1, N + 1)
    for ax, vec, k in zip(axes, B, range(1, N)):
        ax.bar(layers, vec, color=C_TEAL, edgecolor=C_SLATE, width=0.7)
        ax.axhline(0, color=C_GRAY, lw=0.8)
        ax.set_title(rf"$e^{{({k})}}$, $\sum_i v_i={vec.sum():.1e}$", fontsize=8)
        ax.set_xlabel("layer")
        ax.set_xticks(layers)
    axes[0].set_ylabel("amplitude")
    fig.suptitle("Catalyzon dark-mode (Helmert) basis — sum-zero subspace", fontsize=10)
    fig.tight_layout()
    save(fig, "catalyzon_dark_modes")


# ---------------------------------------------------------------------------
# 5. Toroidal geometry (2D schematic with corrected labels)
# ---------------------------------------------------------------------------
def fig_toroidal_geometry() -> None:
    fig, ax = plt.subplots(figsize=(5, 4))
    R, r = 2.2, 0.7
    theta = np.linspace(0, 2 * np.pi, 400)
    # outer/inner envelopes of torus tube in 2D projection
    ax.plot((R + r) * np.cos(theta), (R + r) * np.sin(theta), color=C_BLUE, lw=1.5)
    ax.plot((R - r) * np.cos(theta), (R - r) * np.sin(theta), color=C_BLUE, lw=1.5)
    ax.plot(R * np.cos(theta), R * np.sin(theta), color=C_CORAL, lw=1.8, label="core centerline")
    # minor circle
    phi = np.linspace(0, 2 * np.pi, 200)
    ax.plot(R + r * np.cos(phi), r * np.sin(phi), color=C_TEAL, lw=1.5)
    ax.annotate(
        "",
        xy=(R + r, 0),
        xytext=(R, 0),
        arrowprops=dict(arrowstyle="<->", color=C_TEAL, lw=1.2),
    )
    ax.text(R + r / 2, 0.15, r"$r\sim a_H$", color=C_TEAL, ha="center", fontsize=9)
    ax.annotate(
        "",
        xy=(R, 0),
        xytext=(0, 0),
        arrowprops=dict(arrowstyle="<->", color=C_CORAL, lw=1.2),
    )
    ax.text(R / 2, -0.35, r"$R\sim\xi$", color=C_CORAL, ha="center", fontsize=9)
    ax.text(0, R + r + 0.35, r"$Q_H=1/N=1/5$", ha="center", color=C_SLATE, fontsize=10)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title("Toroidal hopfion Casimir cavity (schematic)")
    ax.legend(frameon=False, loc="lower right")
    save(fig, "toroidal_geometry_3d")


# ---------------------------------------------------------------------------
# 6. Lens space / fractional charge schematic
# ---------------------------------------------------------------------------
def fig_lens_space_qh() -> None:
    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    N = 5
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False)
    for i, a in enumerate(angles):
        wedge = plt.matplotlib.patches.Wedge(
            (0, 0),
            1.0,
            np.degrees(a),
            np.degrees(a + 2 * np.pi / N),
            width=0.45,
            facecolor=plt.cm.Blues(0.3 + 0.1 * i),
            edgecolor=C_SLATE,
            lw=1.0,
        )
        ax.add_patch(wedge)
        mid = a + np.pi / N
        ax.text(
            0.75 * np.cos(mid),
            0.75 * np.sin(mid),
            rf"$k={i}$",
            ha="center",
            va="center",
            fontsize=8,
        )
    ax.text(0, 0, r"$L(5,1)$" + "\n" + r"$Q_H=1/5$", ha="center", va="center", fontsize=11)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.3, 1.3)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(r"Lens space sectors: $\mathrm{CS}[A_k]=k^2/N$, $Q_H=1/N$ (no $H_2$ cycle)")
    save(fig, "lens_space_qh")


# ---------------------------------------------------------------------------
# 7. Vakulenko–Kapitansky bound
# ---------------------------------------------------------------------------
def fig_vk_bound() -> None:
    N = np.linspace(1, 12, 200)
    bound = N ** (-0.75)
    fig, ax = plt.subplots(figsize=(4.5, 3.2))
    ax.plot(N, bound, color=C_TEAL, lw=2.0, label=r"$N^{-3/4}$ bound")
    ax.plot(5, 5 ** (-0.75), "o", color=C_CORAL, ms=7, label="MATBG ($N=5$)")
    ax.axhline(0.244, color="#6C3483", ls="--", lw=1.2, label=r"$E_{\mathrm{Cas}}=0.244\,\mathrm{meV}/C$ scale")
    ax.set_xlabel(r"$N$ (layers)")
    ax.set_ylabel(r"$E_{\min}/C$")
    ax.set_title("Vakulenko–Kapitansky energy bound")
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "vk_bound")


# ---------------------------------------------------------------------------
# 8. Phase diagram: fold vs contingent BKT
# ---------------------------------------------------------------------------
def fig_gap_phase_diagram() -> None:
    fig, ax = plt.subplots(figsize=(5.0, 3.4))
    ax.axvspan(0, 1.0, color=C_TEAL, alpha=0.12, label="gapped (fold map)")
    ax.axvspan(1.0, 1.25, color=C_CORAL, alpha=0.15, label="no real gapped solution")
    ax.axvline(1.0, color=C_CORAL, lw=2)
    ax.annotate(
        "saddle-node fold",
        xy=(1.0, 0.5),
        xytext=(0.55, 0.75),
        arrowprops=dict(arrowstyle="->", color=C_SLATE),
        fontsize=9,
    )
    ax.annotate(
        "BKT only if\nvortex RG",
        xy=(0.95, 0.2),
        xytext=(0.35, 0.25),
        arrowprops=dict(arrowstyle="->", color=C_GRAY),
        fontsize=8,
        color=C_GRAY,
    )
    ax.plot(0.45, 0.65, "s", color=C_GOLD, ms=8)
    ax.text(0.47, 0.68, "MATBG op.", fontsize=8, color=C_GOLD)
    ax.set_xlim(0, 1.25)
    ax.set_ylim(0, 1)
    ax.set_xlabel(r"$A_{\mathrm{ZPF}}/A_c$")
    ax.set_ylabel(r"$\Delta/\Delta_{\mathrm{ex}}$ (schematic)")
    ax.set_title("Gap collapse: fold (locked) vs contingent BKT")
    ax.legend(frameon=False, loc="upper right", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "gap_phase_diagram")


# ---------------------------------------------------------------------------
# 9. Hall staircase schematic with Q_H=1/N
# ---------------------------------------------------------------------------
def fig_hall_staircase() -> None:
    N = 5
    qh = fractional_hall_charge(N)
    steps = np.arange(0, N + 1)
    # schematic σ_xy steps of height ∝ 1/(2N) in preprint units narrative
    height = 1.0 / (2 * N)
    sigma = steps * height
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    ax.step(steps, sigma, where="post", color=C_TEAL, lw=2.0)
    ax.scatter(steps, sigma, c=C_CORAL, s=30, zorder=3)
    ax.set_xlabel("sector index / filling proxy")
    ax.set_ylabel(r"$\sigma_{xy}$ (arb. units)")
    ax.set_title(rf"Hall staircase schematic — $Q_H=1/N={qh:.2f}$, step $\propto 1/(2N)$")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    save(fig, "hall_staircase")



# ---------------------------------------------------------------------------
# 10. Cavity QVC hypotheses panel (corrected spectrum / fold language)
# ---------------------------------------------------------------------------
def fig_cavity_qvc_hypotheses() -> None:
    N = np.arange(1, 9)
    F_moire = 1.0 + 0.35 * (N - 1)
    A_over_Ac = np.clip(0.12 * F_moire, 0, 0.95)

    fig = plt.figure(figsize=(9.5, 7.2))
    gs = GridSpec(2, 2, hspace=0.35, wspace=0.32)

    ax = fig.add_subplot(gs[0, 0])
    ax.plot(N, F_moire, "o-", color=C_TEAL, lw=1.8, label=r"$\mathcal{F}_{\mathrm{moir\'e}}$")
    ax.plot(N, A_over_Ac, "s--", color=C_CORAL, lw=1.5, label=r"$A_{\mathrm{ZPF}}/A_c$")
    ax.axhline(1.0, color=C_GRAY, ls=":", lw=1.0)
    ax.set_xlabel(r"$N$ layers")
    ax.set_ylabel("enhancement / amplitude")
    ax.set_title("(a) Moiré self-cavity enhancement")
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = fig.add_subplot(gs[0, 1])
    mu = np.linspace(0, 1, 500)
    sigma = np.zeros_like(mu)
    for k in range(0, 11):
        sigma[mu >= k / 10] = k / (2 * 5)
    ax.step(mu, sigma, where="post", color=C_TEAL, lw=1.8)
    ax.set_xlabel(r"$\mu_{\mathrm{eff}}$ (arb.)")
    ax.set_ylabel(r"$\sigma_{xy}$ (arb.)")
    ax.set_title(r"(b) Hall staircase ($Q_H=1/N$; not BKT)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = fig.add_subplot(gs[1, 0])
    wc = np.linspace(0.5, 3.5, 200)
    g0 = 0.21
    for i in range(4):
        ax.plot(wc, -g0 + 0.02 * i + 0.01 * np.sin(3 * wc), color=C_CORAL, lw=1.2, alpha=0.8)
    omega0 = 2.0
    split = 0.35
    upper = 0.5 * (wc + omega0) + 0.5 * np.sqrt((wc - omega0) ** 2 + 4 * split**2)
    lower = 0.5 * (wc + omega0) - 0.5 * np.sqrt((wc - omega0) ** 2 + 4 * split**2)
    ax.plot(wc, upper, color=C_TEAL, lw=2.0, label=r"bright–cavity $+$")
    ax.plot(wc, lower, color=C_BLUE, lw=2.0, label=r"bright–cavity $-$")
    ax.plot(wc, -g0 * np.ones_like(wc), color=C_CORAL, lw=2.0, label=r"dark $\lambda=-g_0$ (×4)")
    ax.set_xlabel(r"$\omega_c$ (arb.)")
    ax.set_ylabel(r"$\lambda$ (meV arb.)")
    ax.set_title(r"(c) Spectrum: $\lambda_{\min}=-g_0$, not $-4g_0$")
    ax.legend(frameon=False, fontsize=6, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax = fig.add_subplot(gs[1, 1])
    a = np.linspace(0.05, 0.99, 200)
    barrier = np.sqrt(np.maximum(0, 1 - a))
    rate = np.exp(-8.0 * barrier)
    ax.plot(a, barrier, color=C_TEAL, lw=2.0, label="barrier (fold map)")
    ax.plot(a, rate / rate.max(), color=C_CORAL, lw=1.8, label="rate (norm.)")
    ax.axvline(0.45, color=C_GOLD, ls="--", lw=1.0, label="MATBG op.")
    ax.set_xlabel(r"$A_{\mathrm{ZPF}}/A_c$")
    ax.set_ylabel("normalized")
    ax.set_title("(d) Vacuum-catalyzed transition (fold)")
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("Cavity-enhanced QVC predictions (corrected spectrum / fold)", fontsize=11, y=0.98)
    save(fig, "cavity_qvc_hypotheses")


# ---------------------------------------------------------------------------
# 11. Devil staircase + quantum atmosphere (corrected labels)
# ---------------------------------------------------------------------------
def fig_devil_staircase_atmosphere() -> None:
    fig = plt.figure(figsize=(9.5, 4.0))
    gs = GridSpec(1, 2, wspace=0.28)

    ax0 = fig.add_subplot(gs[0], projection="3d")
    mu = np.linspace(0, 1, 80)
    az = np.linspace(0.1, 0.9, 40)
    Mu, Az = np.meshgrid(mu, az)
    Sigma = np.zeros_like(Mu)
    for k in range(0, 11):
        Sigma = np.where(Mu >= k / 10, k / 10.0, Sigma)
    Sigma = Sigma * (1.0 - 0.25 * Az)
    ax0.plot_surface(Mu, Az, Sigma, cmap="viridis", linewidth=0, antialiased=True, alpha=0.9)
    ax0.set_xlabel(r"$\mu_{\mathrm{eff}}$")
    ax0.set_ylabel(r"$A_{\mathrm{ZPF}}/A_c$")
    ax0.set_zlabel(r"$\sigma_{xy}$")
    ax0.set_title("(a) Hall staircase vs vacuum coupling")
    ax0.view_init(elev=28, azim=-55)

    ax1 = fig.add_subplot(gs[1])
    r = np.linspace(0.5, 40, 200)
    lam_M = 13.0
    for C, c in zip([1, 2, 3, 5], [C_TEAL, C_BLUE, C_GOLD, C_CORAL]):
        dE = C * np.exp(-r / lam_M) / r
        ax1.plot(r, dE / dE.max(), color=c, lw=1.8, label=rf"$C={C}$")
    ax1.axvline(lam_M, color=C_GRAY, ls=":", lw=1.0)
    ax1.text(lam_M + 0.5, 0.85, r"$\lambda_M\approx 13\,\mathrm{nm}$", fontsize=8, color=C_GRAY)
    ax1.set_xlabel(r"$r$ (nm)")
    ax1.set_ylabel(r"$\Delta E_{\mathrm{QF}}$ (norm.)")
    ax1.set_title("(b) Quantum atmosphere (moiré screening)")
    ax1.legend(frameon=False, fontsize=7)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    fig.suptitle("Devil's staircase atmosphere (schematic; $Q_H=1/N$)", fontsize=11)
    save(fig, "devil_staircase_atmosphere")


# ---------------------------------------------------------------------------
# 12. Mexican-hat / catalyzon barrier (3D surfaces)
# ---------------------------------------------------------------------------
def fig_mexican_hat_3d_catalyzon() -> None:
    fig = plt.figure(figsize=(10.5, 3.6))
    gs = GridSpec(1, 3, wspace=0.15)

    def pot(X, Y, alpha, beta=1.0):
        rho2 = X**2 + Y**2
        return alpha * rho2 + 0.5 * beta * rho2**2

    xs = np.linspace(-1.4, 1.4, 80)
    X, Y = np.meshgrid(xs, xs)

    for idx, (alpha, title) in enumerate(
        [
            (-0.5, r"(a) bare $\alpha_0=-0.5$"),
            (-0.35, r"(b) QVC $\alpha_{\mathrm{QVC}}=-0.35$"),
        ]
    ):
        ax = fig.add_subplot(gs[idx], projection="3d")
        Z = pot(X, Y, alpha)
        ax.plot_surface(X, Y, Z, cmap="coolwarm", linewidth=0, antialiased=True, alpha=0.92)
        ax.set_xlabel(r"$\mathrm{Re}\,\Psi$")
        ax.set_ylabel(r"$\mathrm{Im}\,\Psi$")
        ax.set_zlabel(r"$V$")
        ax.set_title(title, fontsize=9)
        ax.view_init(elev=28, azim=-60)

    ax = fig.add_subplot(gs[2])
    rho = np.linspace(0, 1.4, 200)
    v0 = -0.5 * rho**2 + 0.5 * rho**4
    vq = -0.35 * rho**2 + 0.5 * rho**4
    b0 = 0.5**2 / 2
    bq = 0.35**2 / 2
    ax.plot(rho, v0, color=C_BLUE, lw=2.0, label="bare")
    ax.plot(rho, vq, color=C_CORAL, lw=2.0, label="QVC")
    ax.fill_between(rho, vq, v0, where=(v0 > vq), color=C_GOLD, alpha=0.25, label="barrier Δ")
    ax.axhline(0, color=C_GRAY, lw=0.7)
    ax.set_xlabel(r"$|\Psi|$")
    ax.set_ylabel(r"$V(|\Psi|)$")
    ax.set_title(rf"(c) radial cut (~{(1-bq/b0)*100:.0f}% barrier drop)", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.suptitle("Mexican-hat reshaping under QVC (schematic)", fontsize=11)
    save(fig, "mexican_hat_3d_catalyzon")


def main() -> None:
    print(f"Writing corrected figures to {OUTDIR}")
    fig_gap_fold_collapse()
    fig_five_mode_coupling_matrix()
    fig_catalyzon_formation_3d()
    fig_catalyzon_dark_modes()
    fig_toroidal_geometry()
    fig_lens_space_qh()
    fig_vk_bound()
    fig_gap_phase_diagram()
    fig_hall_staircase()
    fig_cavity_qvc_hypotheses()
    fig_devil_staircase_atmosphere()
    fig_mexican_hat_3d_catalyzon()
    print("Done.")


if __name__ == "__main__":
    main()
