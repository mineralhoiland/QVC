"""Publication-quality matplotlib figures for TEVC material structure.

Generates figures for the TEVC paper showing:
  1. Moire lattice with twist and stacking domains
  2. Band structure E(k) along high-symmetry path
  3. Phase diagram (gap fold map)
  4. Hall staircase sigma_xy(B)
  5. Mexican hat potential V(rho)

Style: PRL single-column format (3.375 x 2.5 inches), 8pt body text,
colorblind-safe palette, vector + raster output.

Dependencies: matplotlib (wrapped in try/except for importability without it).
Data sources: qvc.materials, qvc.gap_fold, qvc.topology.chern_hopf_hall.
"""

from __future__ import annotations

import numpy as np
from typing import Optional

# Deferred matplotlib import: module is importable without matplotlib,
# but functions raise RuntimeError if called without it.
try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend for file output
    import matplotlib.pyplot as plt
    from matplotlib.patches import RegularPolygon, FancyArrowPatch
    from matplotlib.collections import LineCollection
    from mpl_toolkits.mplot3d import Axes3D
    HAS_MPL = True
except ImportError:
    HAS_MPL = False


def _require_matplotlib():
    """Raise informative error if matplotlib is unavailable."""
    if not HAS_MPL:
        raise RuntimeError(
            "matplotlib is required for figure generation. "
            "Install with: pip install matplotlib"
        )


# === Style Configuration ===

# Colorblind-safe palette (Wong, Nature Methods 2011 + custom)
PALETTE = {
    "blue": "#0C5DA5",
    "red": "#FF2C00",
    "green": "#00B945",
    "orange": "#FF9500",
    "purple": "#845B97",
    "gray": "#7F7F7F",
    "black": "#2B2B2B",
}

# PRL single-column dimensions
FIG_WIDTH = 3.375  # inches
FIG_HEIGHT = 2.5   # inches
FONT_BODY = 8      # pt
FONT_LABEL = 9     # pt
LINE_WIDTH = 1.0   # pt


def _apply_style():
    """Set matplotlib rcParams for publication style."""
    plt.rcParams.update({
        "font.size": FONT_BODY,
        "axes.labelsize": FONT_LABEL,
        "axes.titlesize": FONT_LABEL,
        "xtick.labelsize": FONT_BODY,
        "ytick.labelsize": FONT_BODY,
        "legend.fontsize": FONT_BODY - 1,
        "lines.linewidth": LINE_WIDTH,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4,
        "ytick.minor.width": 0.4,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "mathtext.fontset": "dejavusans",
    })


def _save_figure(fig, output: str):
    """Save figure as both PDF (vector) and PNG (raster)."""
    # Save PDF
    if output.endswith(".png"):
        pdf_path = output.replace(".png", ".pdf")
    elif output.endswith(".pdf"):
        pdf_path = output
    else:
        pdf_path = output + ".pdf"

    png_path = pdf_path.replace(".pdf", ".png")

    fig.savefig(pdf_path, format="pdf")
    fig.savefig(png_path, format="png", dpi=300)
    plt.close(fig)


# =============================================================================
# 1. Moire Lattice
# =============================================================================


def plot_moire_lattice(
    theta_deg: float = 1.1,
    n_layers: int = 5,
    output: str = "moire_lattice.png",
):
    """Plot twisted bilayer graphene moire pattern with stacking domains.

    Shows two overlaid honeycomb lattices with relative twist theta,
    color-coding the AA sites (Casimir cavity locations), AB domains,
    and domain walls (where the excitonic gap is maximal).

    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees. Default: 1.1 (magic angle).
    n_layers : int
        Number of layers in the TEVC stack. Annotated on figure.
    output : str
        Output filename (saves both .png and .pdf).
    """
    _require_matplotlib()
    _apply_style()

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_WIDTH))  # square for lattice

    # Moire period
    a_cc = 1.42  # Angstroms
    a = a_cc * np.sqrt(3)  # lattice constant
    theta_rad = np.radians(theta_deg)
    lambda_M = a / (2 * np.sin(theta_rad / 2))  # Angstroms
    lambda_M_nm = lambda_M / 10.0

    # Generate honeycomb lattice points for visualization
    # Show a region ~ 2 moire cells across
    n_cells = int(np.ceil(2.5 * lambda_M / a))
    coords = []
    for i in range(-n_cells, n_cells + 1):
        for j in range(-n_cells, n_cells + 1):
            # A sublattice
            x = i * a + j * a * 0.5
            y = j * a * np.sqrt(3) / 2
            coords.append([x, y])
            # B sublattice
            coords.append([x + a_cc, y])

    coords = np.array(coords)

    # Clip to circular region
    r_clip = 1.3 * lambda_M
    mask = np.linalg.norm(coords, axis=1) < r_clip
    coords = coords[mask]

    # Rotate: layer 1 by +theta/2, layer 2 by -theta/2
    def rotate(pts, angle_rad):
        c, s = np.cos(angle_rad), np.sin(angle_rad)
        R = np.array([[c, -s], [s, c]])
        return pts @ R.T

    layer1 = rotate(coords, +theta_rad / 2)
    layer2 = rotate(coords, -theta_rad / 2)

    # Plot atoms (small dots)
    ax.scatter(layer1[:, 0], layer1[:, 1], s=0.3, c=PALETTE["blue"],
               alpha=0.4, rasterized=True, linewidths=0)
    ax.scatter(layer2[:, 0], layer2[:, 1], s=0.3, c=PALETTE["red"],
               alpha=0.4, rasterized=True, linewidths=0)

    # Mark AA sites (where atoms overlap = moire cell centers)
    # AA sites form a triangular lattice with period lambda_M
    L1 = lambda_M * np.array([1.0, 0.0])
    L2 = lambda_M * np.array([0.5, np.sqrt(3) / 2])
    aa_sites = []
    for i in range(-3, 4):
        for j in range(-3, 4):
            site = i * L1 + j * L2
            if np.linalg.norm(site) < r_clip * 0.9:
                aa_sites.append(site)
    aa_sites = np.array(aa_sites)

    ax.scatter(aa_sites[:, 0], aa_sites[:, 1], s=25,
               c=PALETTE["red"], marker="o", zorder=5,
               edgecolors="none", alpha=0.7, label="AA (cavity)")

    # Draw moire unit cell hexagon at center
    hex_verts = []
    for k in range(7):
        angle = np.pi / 6 + k * np.pi / 3
        hex_verts.append([lambda_M / np.sqrt(3) * np.cos(angle),
                          lambda_M / np.sqrt(3) * np.sin(angle)])
    hex_verts = np.array(hex_verts)
    ax.plot(hex_verts[:, 0], hex_verts[:, 1], color=PALETTE["purple"],
            linewidth=1.2, linestyle="-", label=f"Moire cell", zorder=4)

    # Domain wall network (hexagonal grid connecting AA sites)
    for i in range(len(aa_sites)):
        for j in range(i + 1, len(aa_sites)):
            d = np.linalg.norm(aa_sites[i] - aa_sites[j])
            if abs(d - lambda_M) < lambda_M * 0.1:
                ax.plot([aa_sites[i, 0], aa_sites[j, 0]],
                        [aa_sites[i, 1], aa_sites[j, 1]],
                        color=PALETTE["gray"], linewidth=0.5,
                        alpha=0.5, zorder=2)

    # Scale bar
    bar_x = -r_clip * 0.8
    bar_y = -r_clip * 0.9
    bar_len = 50.0  # Angstroms = 5 nm
    ax.plot([bar_x, bar_x + bar_len], [bar_y, bar_y],
            color="k", linewidth=1.5)
    ax.text(bar_x + bar_len / 2, bar_y - lambda_M * 0.08,
            "5 nm", ha="center", va="top", fontsize=FONT_BODY)

    # Annotations
    ax.text(0.02, 0.98, f"$\\theta = {theta_deg}^\\circ$\n"
            f"$\\lambda_M = {lambda_M_nm:.1f}$ nm\n"
            f"$N = {n_layers}$ layers",
            transform=ax.transAxes, va="top", fontsize=FONT_BODY,
            bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

    ax.set_aspect("equal")
    ax.set_xlim(-r_clip, r_clip)
    ax.set_ylim(-r_clip, r_clip)
    ax.set_xlabel(r"$x$ (\AA)")
    ax.set_ylabel(r"$y$ (\AA)")
    ax.legend(loc="lower right", framealpha=0.8, markerscale=0.8)

    fig.tight_layout()
    _save_figure(fig, output)


# =============================================================================
# 2. Band Structure
# =============================================================================


def plot_band_structure(
    bm_model=None,
    output: str = "band_structure.png",
):
    """Plot band dispersion E(k) along high-symmetry path in moire BZ.

    If bm_model is provided (a BistritzMacDonald solver instance with
    a .bands attribute), uses its output. Otherwise, generates a
    schematic flat-band structure representative of MATBG at magic angle.

    Parameters
    ----------
    bm_model : object, optional
        Bistritzer-MacDonald solver with attributes:
          .k_path: ndarray (N_k,) - path coordinate
          .bands: ndarray (N_k, N_bands) - energies in meV
          .k_labels: list of (position, label) tuples
        If None, generates a schematic band structure.
    output : str
        Output filename.
    """
    _require_matplotlib()
    _apply_style()

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    if bm_model is not None and hasattr(bm_model, "bands"):
        # Use computed band structure
        k_path = bm_model.k_path
        bands = bm_model.bands
        k_labels = getattr(bm_model, "k_labels", None)
    else:
        # Generate schematic MATBG-like band structure
        # Path: Gamma -> K -> M -> Gamma
        N_seg = 100
        k_path = np.linspace(0, 3, 3 * N_seg)

        # Schematic flat bands near E=0 (bandwidth W ~ 5 meV at magic angle)
        W = 5.0  # meV, flat band bandwidth
        E_gap = 30.0  # meV, gap to remote bands

        # Flat bands: cosine dispersion with small bandwidth
        flat_upper = (W / 2) * np.cos(np.pi * k_path / 3)
        flat_lower = -(W / 2) * np.cos(np.pi * k_path / 3 + np.pi / 4)

        # Remote bands: dispersive, starting at +/- E_gap
        remote_upper = E_gap + 20 * (1 - np.cos(2 * np.pi * k_path / 3))
        remote_lower = -E_gap - 20 * (1 - np.cos(2 * np.pi * k_path / 3))

        bands = np.column_stack([remote_lower, flat_lower, flat_upper, remote_upper])
        k_labels = [(0, r"$\Gamma$"), (N_seg, r"$K$"),
                    (2 * N_seg, r"$M$"), (3 * N_seg - 1, r"$\Gamma$")]

    # Plot bands
    n_bands = bands.shape[1]
    for i in range(n_bands):
        # Highlight flat bands (within 10 meV of Fermi level)
        band = bands[:, i]
        is_flat = np.max(np.abs(band)) < 15.0
        color = PALETTE["red"] if is_flat else PALETTE["blue"]
        lw = 1.2 if is_flat else 0.7
        alpha = 1.0 if is_flat else 0.6
        ax.plot(k_path, band, color=color, linewidth=lw, alpha=alpha)

    # Mark bandwidth W
    if bm_model is None:
        W_val = W
        ax.annotate("", xy=(0.3, W / 2), xytext=(0.3, -W / 2),
                    arrowprops=dict(arrowstyle="<->", color=PALETTE["green"],
                                    lw=0.8))
        ax.text(0.5, 0, f"$W = {W_val:.0f}$ meV",
                color=PALETTE["green"], fontsize=FONT_BODY - 1, va="center")

    # Fermi level
    ax.axhline(0, color="k", linewidth=0.4, linestyle="--", alpha=0.5)

    # High-symmetry point labels
    if k_labels is not None:
        tick_positions = [kl[0] if isinstance(kl[0], (int, float))
                          else k_path[kl[0]] for kl in k_labels]
        tick_labels = [kl[1] for kl in k_labels]
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        for pos in tick_positions:
            ax.axvline(pos, color="k", linewidth=0.3, alpha=0.3)

    ax.set_ylabel(r"$E$ (meV)")
    ax.set_xlim(k_path[0], k_path[-1])
    ax.set_ylim(-80, 80)
    ax.set_title("Moire band structure ($\\theta = 1.1^\\circ$)",
                 fontsize=FONT_LABEL)

    fig.tight_layout()
    _save_figure(fig, output)


# =============================================================================
# 3. Phase Diagram (Gap Fold)
# =============================================================================


def plot_phase_diagram(output: str = "phase_diagram.png"):
    """Plot the gap map Delta(alpha) showing saddle-node fold bifurcation.

    The self-consistent gap equation delta = 1 + alpha(ln delta + ell)
    exhibits a fold catastrophe at alpha_c. The physical branch (larger root)
    terminates at the fold edge where it meets the unstable branch.

    Key values (corrections-canon IR, ell = -2.791):
      alpha_c = delta_c = 0.1818
      Working point: lambda* = 0.7628
      Fold-edge gap: Delta_c = 0.109 meV (for Delta_0 = 0.6 meV)

    Parameters
    ----------
    output : str
        Output filename.
    """
    _require_matplotlib()
    _apply_style()

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    # Parameters from qvc.gap_fold
    ell = -2.791  # corrections-canon IR coefficient
    Delta_0 = 0.6  # meV (bare gap = 2J)

    # Solve fold condition: delta_c (1 - ln delta_c - ell) = 1
    # Numerical solution gives delta_c ~ 0.1818
    from functools import reduce
    # Bisection for delta_c
    lo, hi = 1e-12, 1.0
    for _ in range(120):
        mid = 0.5 * (lo + hi)
        val = mid * (1.0 - np.log(mid) - ell) - 1.0
        if val > 0:
            hi = mid
        else:
            lo = mid
    delta_c = 0.5 * (lo + hi)
    alpha_c = delta_c

    # Scan alpha and find branches
    n_pts = 500
    alpha_arr = np.linspace(0.001, alpha_c * 0.999, n_pts)
    phys_branch = np.zeros(n_pts)
    unst_branch = np.zeros(n_pts)

    for idx, alpha in enumerate(alpha_arr):
        # Find roots of: delta - 1 - alpha*(ln(delta) + ell) = 0
        # Scan delta in (0, 1]
        delta_scan = np.linspace(1e-6, 1.0, 5000)
        residual = delta_scan - 1.0 - alpha * (np.log(delta_scan) + ell)
        # Find sign changes
        sign_changes = np.where(np.diff(np.sign(residual)))[0]
        roots = []
        for sc in sign_changes:
            # Linear interpolation
            d1, d2 = delta_scan[sc], delta_scan[sc + 1]
            r1, r2 = residual[sc], residual[sc + 1]
            root = d1 - r1 * (d2 - d1) / (r2 - r1)
            roots.append(root)
        roots = sorted(roots, reverse=True)
        phys_branch[idx] = roots[0] if len(roots) >= 1 else np.nan
        unst_branch[idx] = roots[1] if len(roots) >= 2 else np.nan

    # Convert to physical units
    alpha_normalized = alpha_arr / alpha_c
    Delta_phys = phys_branch * Delta_0
    Delta_unst = unst_branch * Delta_0
    Delta_c_meV = delta_c * Delta_0

    # Working point
    lambda_star = 0.7628
    alpha_star = lambda_star * alpha_c
    # Find delta at working point
    delta_scan = np.linspace(1e-6, 1.0, 10000)
    residual = delta_scan - 1.0 - alpha_star * (np.log(delta_scan) + ell)
    sign_changes = np.where(np.diff(np.sign(residual)))[0]
    if len(sign_changes) > 0:
        sc = sign_changes[0]
        delta_star = delta_scan[sc] - residual[sc] * (delta_scan[sc + 1] - delta_scan[sc]) / (residual[sc + 1] - residual[sc])
        Delta_star = delta_star * Delta_0
    else:
        Delta_star = Delta_c_meV

    # Plot
    ax.plot(alpha_normalized, Delta_phys, color=PALETTE["blue"],
            linewidth=1.2, label="Physical branch")
    ax.plot(alpha_normalized, Delta_unst, color=PALETTE["blue"],
            linewidth=0.8, linestyle="--", alpha=0.6,
            label="Unstable branch")

    # Mark fold point
    ax.plot(1.0, Delta_c_meV, "o", color=PALETTE["red"],
            markersize=5, zorder=5)
    ax.annotate(f"$\\Delta_c = {Delta_c_meV:.3f}$ meV",
                xy=(1.0, Delta_c_meV),
                xytext=(0.7, Delta_c_meV + 0.05),
                fontsize=FONT_BODY - 1,
                arrowprops=dict(arrowstyle="->", color=PALETTE["red"],
                                lw=0.6),
                color=PALETTE["red"])

    # Mark working point
    ax.axvline(lambda_star, color=PALETTE["green"], linewidth=0.6,
               linestyle=":", alpha=0.7)
    ax.plot(lambda_star, Delta_star, "s", color=PALETTE["green"],
            markersize=4, zorder=5)
    ax.text(lambda_star + 0.02, Delta_star + 0.02,
            f"$\\lambda^* = {lambda_star}$",
            color=PALETTE["green"], fontsize=FONT_BODY - 1)

    # Labels
    ax.set_xlabel(r"$\alpha / \alpha_c$")
    ax.set_ylabel(r"$\Delta$ (meV)")
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, Delta_0 * 1.05)
    ax.legend(loc="upper right", framealpha=0.8)
    ax.set_title("Gap fold: saddle-node bifurcation", fontsize=FONT_LABEL)

    fig.tight_layout()
    _save_figure(fig, output)


# =============================================================================
# 4. Hall Staircase
# =============================================================================


def plot_hall_staircase(N: int = 5, output: str = "hall_staircase.png"):
    """Plot sigma_xy vs B/B_0 showing quantized TEVC steps.

    Compares:
      - TEVC staircase: step height e^2/(2Nh) = e^2/10h
      - Standard QHE: step height e^2/h (rescaled for comparison)

    Parameters
    ----------
    N : int
        Number of TEVC layers. Default: 5.
    output : str
        Output filename.
    """
    _require_matplotlib()
    _apply_style()

    fig, ax = plt.subplots(figsize=(FIG_WIDTH, FIG_HEIGHT))

    # Generate field range
    B = np.linspace(0, 8, 2000)

    # TEVC staircase
    step_tevc = 1.0 / (2.0 * N)  # e^2/(2Nh) in units of e^2/h
    base_sigma = float(N)  # N * e^2/h at B=0

    # Sharp staircase
    n_hop = np.floor(B).astype(int)
    sigma_tevc = base_sigma + n_hop * step_tevc

    # Thermal smearing (Fermi function at each step)
    kT_width = 0.03
    sigma_tevc_smooth = base_sigma * np.ones_like(B)
    for n in range(1, int(np.max(B)) + 2):
        fermi = 1.0 / (1.0 + np.exp(-(B - n) / kT_width))
        sigma_tevc_smooth += step_tevc * fermi

    # Standard QHE for comparison (rescaled to same base)
    sigma_qhe = base_sigma * np.ones_like(B)
    step_qhe = 1.0  # e^2/h
    for n in range(1, int(np.max(B)) + 2):
        fermi = 1.0 / (1.0 + np.exp(-(B - n) / kT_width))
        sigma_qhe += step_qhe * fermi

    # Plot
    ax.plot(B, sigma_tevc_smooth, color=PALETTE["blue"], linewidth=1.2,
            label=f"TEVC ($N={N}$): $e^2/{2*N}h$")
    ax.plot(B, sigma_qhe, color=PALETTE["gray"], linewidth=0.8,
            linestyle="--", alpha=0.6, label=r"QHE: $e^2/h$")

    # Annotate step height
    # Arrow showing step height at B=3.5
    y_low = base_sigma + 3 * step_tevc
    y_high = y_low + step_tevc
    ax.annotate("", xy=(3.5, y_high), xytext=(3.5, y_low),
                arrowprops=dict(arrowstyle="<->", color=PALETTE["red"],
                                lw=0.8))
    ax.text(3.7, (y_low + y_high) / 2,
            f"$\\delta\\sigma = e^2/{2*N}h$",
            color=PALETTE["red"], fontsize=FONT_BODY - 1, va="center")

    # Nucleation field markers
    for n in range(1, 8):
        ax.axvline(n, color="k", linewidth=0.2, alpha=0.2)

    ax.set_xlabel(r"$B / B_0$")
    ax.set_ylabel(r"$\sigma_{xy}$ ($e^2/h$)")
    ax.set_xlim(0, 8)
    ax.set_ylim(base_sigma - 0.2, base_sigma + 8 * step_tevc + 0.5)
    ax.legend(loc="upper left", framealpha=0.8)
    ax.set_title("Hall staircase: hopfion nucleation", fontsize=FONT_LABEL)

    fig.tight_layout()
    _save_figure(fig, output)


# =============================================================================
# 5. Mexican Hat Potential
# =============================================================================


def plot_mexican_hat(output: str = "mexican_hat.png"):
    """3D surface plot of Ginzburg-Landau potential V(rho).

    Shows V(rho) = alpha_GL * rho^2 + (beta/2) * rho^4 for:
      - Bare potential: alpha_GL = -Delta_0 (excitonic condensation only)
      - Catalyzon-enhanced: alpha_GL = -Delta_cat (vacuum-assisted, ~30% deeper)

    The deepening of the Mexican hat by the catalyzon mechanism is the
    key visual for the QVC catalysis concept: the vacuum fluctuation
    lowers the barrier without changing the symmetry.

    Parameters
    ----------
    output : str
        Output filename.
    """
    _require_matplotlib()
    _apply_style()

    fig = plt.figure(figsize=(FIG_WIDTH, FIG_WIDTH))
    ax = fig.add_subplot(111, projection="3d")

    # Potential parameters
    Delta_0 = 0.6    # meV, bare excitonic gap
    Delta_cat = 0.78  # meV, catalyzon-enhanced (30% deeper)
    beta = 1.0        # normalized quartic coefficient

    # Grid in order parameter space (rho_x, rho_y)
    N = 100
    rho_max = 1.5
    x = np.linspace(-rho_max, rho_max, N)
    y = np.linspace(-rho_max, rho_max, N)
    X, Y = np.meshgrid(x, y)
    rho_sq = X**2 + Y**2

    # Bare potential: V = -Delta_0 * rho^2 + (beta/2) * rho^4
    V_bare = -Delta_0 * rho_sq + (beta / 2) * rho_sq**2

    # Catalyzon potential: deeper well
    V_cat = -Delta_cat * rho_sq + (beta / 2) * rho_sq**2

    # Plot catalyzon potential (solid, colored)
    surf_cat = ax.plot_surface(X, Y, V_cat, alpha=0.7,
                               cmap="coolwarm", linewidth=0,
                               antialiased=True)

    # Plot bare potential (wireframe, for comparison)
    ax.plot_wireframe(X, Y, V_bare, alpha=0.3, color=PALETTE["gray"],
                      linewidth=0.3, rstride=5, cstride=5)

    # Mark minima
    rho_min_bare = np.sqrt(Delta_0 / beta)
    rho_min_cat = np.sqrt(Delta_cat / beta)
    V_min_bare = -Delta_0**2 / (2 * beta)
    V_min_cat = -Delta_cat**2 / (2 * beta)

    ax.scatter([rho_min_cat], [0], [V_min_cat], color=PALETTE["red"],
               s=20, zorder=5)
    ax.scatter([rho_min_bare], [0], [V_min_bare], color=PALETTE["gray"],
               s=15, zorder=5)

    # Labels
    ax.set_xlabel(r"$\mathrm{Re}\,\Delta$", labelpad=-2)
    ax.set_ylabel(r"$\mathrm{Im}\,\Delta$", labelpad=-2)
    ax.set_zlabel(r"$V(\Delta)$ (meV)", labelpad=-2)
    ax.set_title("Order parameter potential", fontsize=FONT_LABEL, pad=0)

    # Annotation: well depth reduction
    ax.text2D(0.02, 0.92, f"Bare: $V_{{min}} = {V_min_bare:.3f}$ meV",
              transform=ax.transAxes, fontsize=FONT_BODY - 1,
              color=PALETTE["gray"])
    ax.text2D(0.02, 0.85, f"Catalyzon: $V_{{min}} = {V_min_cat:.3f}$ meV",
              transform=ax.transAxes, fontsize=FONT_BODY - 1,
              color=PALETTE["red"])
    depth_reduction = (V_min_cat - V_min_bare) / V_min_bare * 100
    ax.text2D(0.02, 0.78, f"Deepening: {depth_reduction:.0f}%",
              transform=ax.transAxes, fontsize=FONT_BODY - 1,
              color=PALETTE["blue"])

    # View angle
    ax.view_init(elev=25, azim=-60)
    ax.set_zticks([])
    ax.xaxis.set_tick_params(labelsize=6)
    ax.yaxis.set_tick_params(labelsize=6)

    fig.tight_layout()
    _save_figure(fig, output)
