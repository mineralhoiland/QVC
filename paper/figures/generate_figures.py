"""
QVC Publication Figures: 6 PDF figures for arXiv v4 preprint
=============================================================
Generates: gap_fold_collapse, gap_phase_diagram, five_mode_coupling_matrix,
           hall_staircase, symmetry_cascade_energy, catalyzon_dark_modes
Style: Okabe-Ito palette, sans-serif 9pt/7pt, single-column width
Dependencies: numpy, matplotlib (no scipy required)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# === Style Setup ===
oiBlue = '#0072B2'
oiVerm = '#D55E00'
oiOrange = '#E69F00'
oiGreen = '#009E73'
oiSky = '#56B4E9'
oiPink = '#CC79A7'

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 9,
    'axes.titlesize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'axes.linewidth': 0.6,
    'xtick.major.width': 0.6,
    'ytick.major.width': 0.6,
    'lines.linewidth': 1.2,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
})

OUTDIR = '/sessions/busy-eager-gauss/mnt/QVCCursor/tex/QVC_arXiv_v4/figures/'
SINGLE_COL = 3.4  # inches


def bisect(f, a, b, tol=1e-12, maxiter=100):
    """Simple bisection root finder (no scipy needed)."""
    fa, fb = f(a), f(b)
    assert fa * fb < 0, f"No sign change: f({a})={fa}, f({b})={fb}"
    for _ in range(maxiter):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(fc) < tol or (b - a) < tol:
            return c
        if fa * fc < 0:
            b, fb = c, fc
        else:
            a, fa = c, fc
    return 0.5 * (a + b)


# ============================================================
# FIGURE 1: Gap Fold Collapse (Saddle-Node Bifurcation)
# ============================================================
def fig_gap_fold_collapse():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.4))

    ell = 1.0  # representative log-enhancement parameter

    # Parametrize: alpha(delta) = (delta - 1) / (ln(delta) + ell)
    delta = np.linspace(0.38, 2.5, 2000)
    denom = np.log(delta) + ell
    valid = denom > 0.01
    delta_v = delta[valid]
    alpha_v = (delta_v - 1.0) / (np.log(delta_v) + ell)

    # Find fold point: d(alpha)/d(delta) = 0
    # Condition: delta*ln(delta) + delta*ell - delta + 1 = 0
    def fold_cond(d):
        return d * np.log(d) + d * ell - d + 1.0

    delta_c = bisect(fold_cond, 0.5, 2.0)
    alpha_c = (delta_c - 1.0) / (np.log(delta_c) + ell)

    # Split into upper (stable) and lower (unstable) branches at fold
    idx_fold = np.argmax(alpha_v)

    # Upper branch: delta > delta_c
    delta_upper = delta_v[idx_fold:]
    alpha_upper = alpha_v[idx_fold:]

    # Lower branch: delta < delta_c
    delta_lower = delta_v[:idx_fold + 1]
    alpha_lower = alpha_v[:idx_fold + 1]

    ax.plot(alpha_upper, delta_upper, '-', color=oiBlue, linewidth=1.5,
            label=r'$\delta_+(\alpha)$ (stable)')
    ax.plot(alpha_lower, delta_lower, '--', color=oiVerm, linewidth=1.5,
            label=r'$\delta_-(\alpha)$ (unstable)')

    # Mark fold point
    ax.plot(alpha_c, delta_c, 'o', color='k', markersize=6, zorder=5)

    # Annotations
    ax.annotate('Fold (saddle-node)',
                xy=(alpha_c, delta_c),
                xytext=(alpha_c - 0.04, delta_c + 0.5),
                fontsize=7, ha='center',
                arrowprops=dict(arrowstyle='->', color='k', lw=0.8))

    ax.text(0.06, 0.3, r'$\alpha_c = %.4f$ (NOT BKT)' % alpha_c,
            fontsize=7, color=oiVerm,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=oiVerm, alpha=0.8))

    ax.set_xlabel(r'Coupling parameter $\alpha$')
    ax.set_ylabel(r'Reduced gap $\delta = \Delta/\Delta_0$')
    ax.set_xlim(0, 0.25)
    ax.set_ylim(0, 2.5)
    ax.legend(loc='upper left', frameon=False)

    fig.savefig(OUTDIR + 'gap_fold_collapse.pdf', bbox_inches='tight')
    plt.close(fig)
    print("  [1/6] gap_fold_collapse.pdf — alpha_c = %.4f, delta_c = %.4f" % (alpha_c, delta_c))


# ============================================================
# FIGURE 2: Gap Phase Diagram
# ============================================================
def fig_gap_phase_diagram():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.6))

    alpha_c0 = 0.1818

    # Boundary: alpha_c(T) = alpha_c0 * (1 - (T/T*)^2)
    T_ratio = np.linspace(0, 1.0, 300)
    alpha_boundary = alpha_c0 * (1.0 - T_ratio**2)

    # Shading
    ax.fill_between(alpha_boundary, T_ratio, 0, alpha=0.15, color=oiBlue)
    ax.fill_betweenx(T_ratio, alpha_boundary, 0.25, alpha=0.10, color=oiOrange)

    # Boundary curve
    ax.plot(alpha_boundary, T_ratio, '-', color='k', linewidth=1.3)

    # Region labels
    ax.text(0.07, 0.3, 'Gapped\n(condensate)', fontsize=8, color=oiBlue,
            ha='center', fontweight='bold')
    ax.text(0.20, 0.9, 'Normal', fontsize=8, color=oiVerm,
            ha='center', fontweight='bold')

    # Schwinger threshold
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
    ax.text(0.22, 1.03, r'$T/T^* = 1$', fontsize=6, color='gray', ha='right')

    # Operating point
    ax.plot(0.16, 0.17, '*', color=oiGreen, markersize=10, zorder=5,
            markeredgecolor='k', markeredgewidth=0.4)
    ax.annotate('Operating\npoint',
                xy=(0.16, 0.17), xytext=(0.20, 0.35),
                fontsize=6, ha='center', color=oiGreen,
                arrowprops=dict(arrowstyle='->', color=oiGreen, lw=0.8))

    ax.set_xlabel(r'Coupling $\alpha$')
    ax.set_ylabel(r'$T / T^*$')
    ax.set_xlim(0, 0.25)
    ax.set_ylim(0, 1.5)

    # Secondary y-axis with mK
    ax2 = ax.secondary_yaxis('right')
    ax2.set_ylabel(r'$T$ (mK)', fontsize=7)
    T_star_mK = 118.7
    ax2.set_yticks([0, 0.5, 1.0, 1.5])
    ax2.set_yticklabels(['0', '59', '119', '178'])
    ax2.tick_params(labelsize=6)

    fig.savefig(OUTDIR + 'gap_phase_diagram.pdf', bbox_inches='tight')
    plt.close(fig)
    print("  [2/6] gap_phase_diagram.pdf")


# ============================================================
# FIGURE 3: Five-Mode Coupling Matrix
# ============================================================
def fig_five_mode_coupling_matrix():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.0))

    g0 = 0.0982  # meV
    N = 5
    # Democratic coupling: G = g0*(J_N - I_N) visualized with
    # off-diagonal = +g0, diagonal = -g0 per user specification
    G = np.full((N, N), g0)
    np.fill_diagonal(G, -g0)

    cmap = LinearSegmentedColormap.from_list('custom',
        [(0, oiBlue), (0.5, 'white'), (1, oiOrange)])

    vmax = g0 * 1.2
    im = ax.imshow(G, cmap=cmap, vmin=-vmax, vmax=vmax, aspect='equal')

    # Annotate cells
    for i in range(N):
        for j in range(N):
            val = G[i, j]
            color = 'white' if abs(val) > 0.06 else 'black'
            ax.text(j, i, f'{val:+.3f}', ha='center', va='center',
                    fontsize=6.5, color=color, fontweight='bold')

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(r'$G_{ij}$ (meV)', fontsize=7)
    cbar.ax.tick_params(labelsize=6)

    ax.set_xticks(range(N))
    ax.set_xticklabels(range(1, N + 1))
    ax.set_yticks(range(N))
    ax.set_yticklabels(range(1, N + 1))
    ax.set_xlabel('Mode index $j$')
    ax.set_ylabel('Mode index $i$')

    # Restore all spines for matrix plot
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.4)

    # Spectrum annotation below figure
    fig.text(0.5, 0.02,
             r'Spec$(G)$: $\lambda_{\max}=(N{-}1)g_0=0.393$ meV, '
             r'$\lambda_{\min}={-}g_0=-0.098$ meV (mult. 4)',
             fontsize=6.5, ha='center', va='bottom')

    fig.subplots_adjust(bottom=0.12)

    fig.savefig(OUTDIR + 'five_mode_coupling_matrix.pdf', bbox_inches='tight')
    plt.close(fig)
    print("  [3/6] five_mode_coupling_matrix.pdf")


# ============================================================
# FIGURE 4: Hall Staircase
# ============================================================
def fig_hall_staircase():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.6))

    N = 5

    # QVC staircase in units of e^2/(2Nh)
    # sigma_xy = n for n = 0,1,...,5 in units of e^2/(2Nh)
    # Each plateau spans 1 unit of filling parameter
    filling = np.linspace(0, 6, 6000)

    # QVC plateaus: step height = e^2/(2Nh) = e^2/(10h)
    sigma_qvc = np.floor(filling + 0.001)  # integer steps

    # Standard QHE for comparison: step height = e^2/h = 2N times larger
    # In same units of e^2/(2Nh): standard step = 2N = 10
    sigma_std = np.floor(filling + 0.001) * (2 * N)

    # Plot QVC
    ax.step(filling, sigma_qvc, where='post', color=oiBlue, linewidth=1.5,
            label=r'QVC ($N=5$): $\delta\sigma = e^2/10h$', zorder=3)

    # Standard QHE (rescaled to show on same plot: use actual units e^2/(2Nh))
    ax.step(filling, sigma_std, where='post', color=oiGreen, linewidth=1.0,
            linestyle='--', alpha=0.6,
            label=r'Standard QHE: $\delta\sigma = e^2/h$')

    ax.set_xlabel(r'Filling parameter $\nu$')
    ax.set_ylabel(r'$\sigma_{xy}$ $\left[\,e^2/(2Nh)\,\right]$')
    ax.set_xlim(0, 5.5)
    ax.set_ylim(-0.5, 6)

    # Mark fractional step height
    ax.annotate('', xy=(5.2, 2.0), xytext=(5.2, 1.0),
                arrowprops=dict(arrowstyle='<->', color=oiVerm, lw=1.2))
    ax.text(5.35, 1.5, r'$\frac{e^2}{2Nh}$', fontsize=8, color=oiVerm,
            va='center')

    # Annotation box
    ax.text(0.3, 5.2, r'$N=5$ layers $\to$ $1/2N$ quantization',
            fontsize=7, color=oiBlue,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                      edgecolor=oiBlue, alpha=0.7))

    ax.legend(loc='lower right', frameon=False, fontsize=6.5)

    # Thin gridlines at QVC plateau levels
    for n in range(7):
        ax.axhline(n, color='gray', linewidth=0.3, alpha=0.4)

    fig.savefig(OUTDIR + 'hall_staircase.pdf', bbox_inches='tight')
    plt.close(fig)
    print("  [4/6] hall_staircase.pdf")


# ============================================================
# FIGURE 5: Symmetry Cascade Energy Levels
# ============================================================
def fig_symmetry_cascade_energy():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.2))

    # Transition temperatures and labels
    levels = [
        (30.0, r'$U(4){\times}U(4) \to G_1$', 'Valley polarization', r'$P_v$', oiBlue),
        (6.0,  r'$G_1 \to G_2$', r'Kekul\'{e} order', r'$\Delta_{KK\'}$', oiOrange),
        (1.0,  r'$G_2 \to G_{\mathrm{BEC}}$', 'Excitonic condensation', r'$\Delta_{\mathrm{ex}}$', oiGreen),
        (0.12, r'$G_{\mathrm{BEC}} \to G_{\mathrm{QVC}}$', 'ZPF coupling\n(non-thermal)', r'$A_{\mathrm{ZPF}}$', oiPink),
    ]

    # Use log-spaced y positions for visual clarity
    y_positions = [3.5, 2.5, 1.5, 0.5]
    bar_left = 0.2
    bar_right = 0.8

    for i, (T, sym_label, phys_label, op_label, color) in enumerate(levels):
        y = y_positions[i]
        # Horizontal bar
        ax.plot([bar_left, bar_right], [y, y], '-', color=color, linewidth=4.0,
                solid_capstyle='round')
        # Temperature label (right)
        ax.text(bar_right + 0.03, y, f'{T} K', fontsize=7, va='center', color='k')
        # Symmetry breaking label (left)
        ax.text(bar_left - 0.03, y, sym_label, fontsize=6.5, va='center',
                ha='right', color=color)
        # Physics label (above bar)
        ax.text((bar_left + bar_right) / 2, y + 0.15,
                phys_label, fontsize=6, va='bottom', ha='center', color='gray')

    # Downward arrows between levels with order parameter labels
    arrow_x = (bar_left + bar_right) / 2
    for i in range(len(levels) - 1):
        y_top = y_positions[i] - 0.2
        y_bot = y_positions[i + 1] + 0.2
        ax.annotate('', xy=(arrow_x, y_bot),
                    xytext=(arrow_x, y_top),
                    arrowprops=dict(arrowstyle='->', color='gray', lw=1.0))
        # Order parameter on arrow
        ax.text(arrow_x + 0.12, (y_top + y_bot) / 2,
                levels[i + 1][3], fontsize=7, va='center', ha='left',
                color=levels[i + 1][4])

    ax.set_xlim(-0.4, 1.15)
    ax.set_ylim(0, 4.2)

    # Remove all axes (pure diagram)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    ax.set_title('Sequential symmetry-breaking cascade', fontsize=8, pad=8)

    fig.savefig(OUTDIR + 'symmetry_cascade_energy.pdf', bbox_inches='tight')
    plt.close(fig)
    print("  [5/6] symmetry_cascade_energy.pdf")


# ============================================================
# FIGURE 6: Catalyzon Dark Modes
# ============================================================
def fig_catalyzon_dark_modes():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.4))

    N = 5
    g0 = 0.0982  # meV
    lambda_max = (N - 1) * g0   # 0.3928 meV
    lambda_min = -g0             # -0.0982 meV

    # Mode layout: 1 bright (superradiant) + 4 dark (degenerate)
    dark_indices = [1, 2, 4, 5]
    bright_idx = 3

    # Degeneracy lines
    ax.axhline(lambda_min, color='gray', linewidth=0.5, linestyle='--', alpha=0.5,
               xmin=0.05, xmax=0.95)
    ax.axhline(lambda_max, color='gray', linewidth=0.5, linestyle='--', alpha=0.5,
               xmin=0.05, xmax=0.95)

    # Dark modes (open circles)
    for idx in dark_indices:
        ax.plot(idx, lambda_min, 'o', color=oiVerm, markersize=9,
                markerfacecolor='white', markeredgewidth=1.8, zorder=4)

    # Bright mode (filled circle)
    ax.plot(bright_idx, lambda_max, 'o', color=oiBlue, markersize=11,
            markerfacecolor=oiBlue, markeredgewidth=1.5, zorder=4)

    # Annotations
    ax.annotate('Bright (superradiant)',
                xy=(bright_idx + 0.2, lambda_max),
                xytext=(bright_idx + 1.0, lambda_max + 0.06),
                fontsize=7, color=oiBlue, va='center',
                arrowprops=dict(arrowstyle='->', color=oiBlue, lw=0.8))

    ax.annotate(r'Dark (subradiant, mult. $N{-}1=4$)',
                xy=(dark_indices[1] + 0.2, lambda_min),
                xytext=(3.0, lambda_min - 0.10),
                fontsize=7, color=oiVerm, va='top', ha='center',
                arrowprops=dict(arrowstyle='->', color=oiVerm, lw=0.8))

    # Splitting double arrow
    ax.annotate('', xy=(0.6, lambda_max - 0.01), xytext=(0.6, lambda_min + 0.01),
                arrowprops=dict(arrowstyle='<->', color='k', lw=1.0))
    splitting = N * g0
    ax.text(0.25, (lambda_max + lambda_min) / 2,
            r'$\Delta\lambda = Ng_0$' + '\n' + f'= {splitting:.3f} meV',
            fontsize=6, va='center', ha='center')

    # Zero energy reference
    ax.axhline(0, color='gray', linewidth=0.4, linestyle=':')
    ax.text(5.6, 0.01, '0', fontsize=6, color='gray', va='bottom')

    ax.set_xlabel('Mode index')
    ax.set_ylabel('Energy (meV)')
    ax.set_xlim(0, 6)
    ax.set_ylim(-0.25, 0.55)
    ax.set_xticks(range(1, N + 1))

    fig.savefig(OUTDIR + 'catalyzon_dark_modes.pdf', bbox_inches='tight')
    plt.close(fig)
    print("  [6/6] catalyzon_dark_modes.pdf")


# ============================================================
# MAIN
# ============================================================
if __name__ == '__main__':
    print("Generating QVC publication figures...")
    print("  Output directory:", OUTDIR)
    print()
    fig_gap_fold_collapse()
    fig_gap_phase_diagram()
    fig_five_mode_coupling_matrix()
    fig_hall_staircase()
    fig_symmetry_cascade_energy()
    fig_catalyzon_dark_modes()
    print("\nDone. All 6 PDFs written.")
