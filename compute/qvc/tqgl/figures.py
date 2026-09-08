"""TQGL v4 figures: PDF + PNG at 300 dpi."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib.pyplot as plt
import numpy as np

C_TEAL = "#0D7377"
C_CORAL = "#C44536"
C_SLATE = "#2C3E50"
C_GOLD = "#B8860B"
C_BLUE = "#1B4F72"
C_GRAY = "#7F8C8D"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.size": 9,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "legend.fontsize": 8,
        "figure.dpi": 120,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "axes.linewidth": 0.8,
    }
)


def _save(fig: plt.Figure, outdir: Path, stem: str) -> list[str]:
    outdir.mkdir(parents=True, exist_ok=True)
    pdf = outdir / f"{stem}.pdf"
    png = outdir / f"{stem}.png"
    fig.savefig(pdf, dpi=300)
    fig.savefig(png, dpi=300)
    plt.close(fig)
    return [str(pdf), str(png)]


def write_figures(bundle: dict[str, Any], outdir: Path) -> list[str]:
    """Write all suite figures. Returns saved paths."""
    from qvc.gap_fold import fold_curve
    from qvc.tqgl.locks import FrozenTQGL

    lock_d = bundle["lock"]
    lock = FrozenTQGL(**lock_d)
    c = bundle["coupled_gpe_dce"]
    res = bundle["dce_fold_resonant"]
    det = bundle["dce_fold_detuning"]
    sch = bundle["schwinger"]
    berry = bundle["berry"]
    saved: list[str] = []

    # 1. gap vs time (live vs fold)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(c["t_ps"], c["Delta_live_meV"], color=C_TEAL, lw=1.8, label=r"$\Delta_{\mathrm{live}}(\psi)$")
    ax.plot(c["t_ps"], c["Delta_fold_of_A_meV"], color=C_GOLD, lw=1.4, ls="--", label=r"fold $\Delta(A(t))$")
    ax.axhline(lock.Delta_star_meV, color=C_BLUE, ls=":", lw=1.0, label=rf"$\Delta^\star={lock.Delta_star_meV:.3f}$ meV")
    ax.axhline(lock.Delta_c_meV, color=C_CORAL, ls=":", lw=1.0, label=rf"$\Delta_c={lock.Delta_c_meV:.3f}$ meV")
    ax.set_xlabel(r"$t$ (ps)")
    ax.set_ylabel(r"$\Delta$ (meV)")
    ax.set_title("Live gap from $\\psi$ vs static fold")
    ax.legend(frameon=False)
    ax.set_xlim(0, c["t_ps"][-1])
    saved += _save(fig, outdir, "gap_vs_time")

    # 2. g1(t)
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot(c["t_ps"], c["g1"], color=C_SLATE, lw=1.8)
    ax.set_xlabel(r"$t$ (ps)")
    ax.set_ylabel(r"$g^{(1)}$")
    ax.set_title("Global first-order coherence")
    ax.set_ylim(0, 1.05)
    saved += _save(fig, outdir, "g1_vs_time")

    # 3. |ψ|² heatmap
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), constrained_layout=True)
    extent = [c["x_um"][0] * 1e3, c["x_um"][-1] * 1e3, c["y_um"][0] * 1e3, c["y_um"][-1] * 1e3]
    for ax, dens, title in (
        (axes[0], c["density_initial"], r"$|\psi|^2$ at $t=0$"),
        (axes[1], c["density_final"], r"$|\psi|^2$ at $t=1$ ps"),
    ):
        im = ax.imshow(dens, origin="lower", extent=extent, cmap="magma", aspect="equal")
        ax.set_xlabel(r"$x$ (nm)")
        ax.set_ylabel(r"$y$ (nm)")
        ax.set_title(title)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    saved += _save(fig, outdir, "density_heatmap")

    # 4. A(t) long DCE
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(res["t_ps"], res["A_um_inv"], color=C_TEAL, lw=1.7, label="resonant (fold-slaved)")
    ax.plot(det["t_ps"], det["A_um_inv"], color=C_CORAL, lw=1.5, label="detuning-feedback")
    ax.axhline(lock.A_c_fold_um_inv, color=C_GOLD, ls="--", lw=1.0, label=r"$A_c^{\mathrm{fold}}$")
    ax.axhline(lock.A_star_um_inv, color=C_BLUE, ls=":", lw=1.0, label=r"$A(\lambda^\star)$")
    ax.set_xlabel(r"$t$ (ps)")
    ax.set_ylabel(r"$A$ ($\mu\mathrm{m}^{-1}$)")
    ax.set_title(r"DCE amplitude (not forced to $A^\ast/A_c=0.52$)")
    ax.legend(frameon=False)
    saved += _save(fig, outdir, "A_vs_time")

    # 5. Δ(A)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(res["A_um_inv"], res["Delta_fold_meV"], color=C_TEAL, lw=1.6, label="resonant")
    ax.plot(det["A_um_inv"], det["Delta_fold_meV"], color=C_CORAL, lw=1.4, ls="--", label="detuning")
    ax.axhline(lock.Delta_c_meV, color=C_GOLD, ls=":", lw=1.0, label=r"$\Delta_c$")
    ax.axvline(lock.A_c_fold_um_inv, color=C_GRAY, ls=":", lw=0.9)
    ax.set_xlabel(r"$A$ ($\mu\mathrm{m}^{-1}$)")
    ax.set_ylabel(r"$\Delta(A)$ (meV)")
    ax.set_title("Fold-slaved gap vs amplitude")
    ax.legend(frameon=False)
    saved += _save(fig, outdir, "Delta_of_A")

    # 6. Lorentzian
    fig, ax = plt.subplots(figsize=(5.2, 3.2))
    ax.plot(res["t_ps"], res["lorentzian"], color=C_TEAL, lw=1.6, label="resonant")
    ax.plot(det["t_ps"], det["lorentzian"], color=C_CORAL, lw=1.4, label="detuning")
    ax.set_xlabel(r"$t$ (ps)")
    ax.set_ylabel(r"$\mathcal{L}(\delta\omega/\Gamma)$")
    ax.set_title("Cavity Lorentzian (not parked at 0)")
    ax.set_ylim(-0.05, 1.05)
    ax.legend(frameon=False)
    saved += _save(fig, outdir, "lorentzian")

    # 7. (Δ, A) portrait
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.plot(c["A_um_inv"], c["Delta_live_meV"], color=C_TEAL, lw=1.8, label=r"GPE $\Delta_{\mathrm{live}}$")
    ax.plot(res["A_um_inv"], res["Delta_fold_meV"], color=C_GOLD, lw=1.3, ls="--", label="DCE fold resonant")
    ax.plot(lock.A_star_um_inv, lock.Delta_star_meV, "s", color=C_BLUE, ms=7, label=r"$(\lambda^\star,\Delta^\star)$")
    ax.axhline(lock.Delta_c_meV, color=C_CORAL, ls=":", lw=1.0)
    ax.set_xlabel(r"$A$ ($\mu\mathrm{m}^{-1}$)")
    ax.set_ylabel(r"$\Delta$ (meV)")
    ax.set_title(r"$(\Delta,A)$ portrait")
    ax.legend(frameon=False)
    saved += _save(fig, outdir, "portrait_Delta_A")

    # 8. fold curve
    fc = fold_curve(lock.ell, n=120)
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    a = np.array(fc["alpha_over_ac"], dtype=float)
    phys = np.array([np.nan if p is None else p for p in fc["delta_phys"]], dtype=float)
    unst = np.array([np.nan if p is None else p for p in fc["delta_unst"]], dtype=float)
    ax.plot(a, phys, color=C_TEAL, lw=2.0, label="physical")
    ax.plot(a, unst, color=C_CORAL, lw=1.4, ls="--", label="unstable")
    ax.plot(1.0, lock.delta_c, "o", color=C_CORAL, ms=7, label=r"fold $(\alpha_c,\delta_c)$")
    ax.plot(lock.alpha_star_over_ac, lock.delta_star, "s", color=C_GOLD, ms=7, label=r"$\lambda^\star$")
    ax.set_xlabel(r"$\alpha/\alpha_c$")
    ax.set_ylabel(r"$\delta=\Delta/\Delta_0$")
    ax.set_title("Saddle-node fold (not BKT)")
    ax.legend(frameon=False)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 1.05)
    saved += _save(fig, outdir, "fold_curve")

    # 9. Schwinger
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.semilogy(sch["ratio_sweep"], sch["Gamma_sch_sweep"], color=C_SLATE, lw=1.6)
    ax.axvline(sch["A_eff_over_Ac_schwinger"], color=C_TEAL, ls="--", lw=1.0, label=r"$A^\star/A_{c,\mathrm{Sch}}$")
    ax.set_xlabel(r"$A_{\mathrm{eff}}/A_{c,\mathrm{Sch}}$")
    ax.set_ylabel(r"$\Gamma_{\mathrm{Sch}}$ (ps$^{-1}$)")
    ax.set_title("Schwinger rate (new $E_{\\mathrm{vac}}$, democratic $|\\lambda_{\\min}|$)")
    ax.legend(frameon=False)
    saved += _save(fig, outdir, "schwinger_rate")

    # 10. Berry histograms
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    ax.hist(np.abs(berry["lambda_min_GOE"]), bins=40, density=True, alpha=0.55, color=C_GRAY, label="GOE")
    ax.hist(np.abs(berry["lambda_min_GUE"]), bins=40, density=True, alpha=0.55, color=C_TEAL, label="GUE")
    ax.set_xlabel(r"$|\lambda_{\min}|$ (meV)")
    ax.set_ylabel("density")
    ax.set_title(rf"Berry MC (not democratic $G$): $R={berry['R_mean']:.3f}$")
    ax.legend(frameon=False)
    saved += _save(fig, outdir, "berry_goe_gue")

    return saved
