"""
QVC 3D Figure Generator: Publication-quality 3D matplotlib figures
=================================================================
Generates 4 PDF figures for QVC arXiv v4 preprint:
  1. toroidal_geometry_3d.pdf    — Fat torus Casimir geometry (Sec 4)
  2. mexican_hat_3d_catalyzon.pdf — TQGL sombrero potential (Sec 7)
  3. hopfion_vortex_3d.pdf       — Hopf fibration linked fibers (Sec 8)
  4. catalyzon_formation_3d.pdf  — Catalyzon mode hybridization (Sec 6)

Style: Okabe-Ito palette, sans-serif 9pt labels, 7pt ticks, ~4in wide, PDF output.
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d.proj3d import proj_transform
import matplotlib as mpl
import os

# === Style Setup ===
oiBlue = '#0072B2'
oiVerm = '#D55E00'
oiOrange = '#E69F00'
oiGreen = '#009E73'
oiSky = '#56B4E9'
oiPink = '#CC79A7'
oiBlack = '#000000'

mpl.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 9,
    'axes.labelsize': 9,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.titlesize': 10,
    'figure.dpi': 150,
    'savefig.dpi': 300,
})

OUTDIR = os.path.dirname(os.path.abspath(__file__))


def set_3d_panes_white(ax):
    """Make 3D pane backgrounds white and lighten grid."""
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')
    ax.grid(True, alpha=0.3, linewidth=0.5)


# =========================================================================
# FIGURE 1: Toroidal Geometry 3D (Sec 4)
# =========================================================================
def fig_toroidal_geometry():
    print("  Generating toroidal_geometry_3d.pdf ...")
    R = 50.0  # major radius (nm)
    r = 8.0   # minor radius (nm)

    u = np.linspace(0, 2*np.pi, 60)
    v = np.linspace(0, 2*np.pi, 60)
    U, V = np.meshgrid(u, v)

    X = (R + r * np.cos(V)) * np.cos(U)
    Y = (R + r * np.cos(V)) * np.sin(U)
    Z = r * np.sin(V)

    fig = plt.figure(figsize=(4.2, 3.8))
    ax = fig.add_subplot(111, projection='3d')

    # Surface plot
    ax.plot_surface(X, Y, Z, alpha=0.65, cmap='cool',
                    edgecolor='none', antialiased=True,
                    lightsource=mpl.colors.LightSource(azdeg=45, altdeg=60))

    # Dashed center circle (midline of torus, z=0 plane)
    theta_c = np.linspace(0, 2*np.pi, 100)
    xc = R * np.cos(theta_c)
    yc = R * np.sin(theta_c)
    zc = np.zeros_like(theta_c)
    ax.plot(xc, yc, zc, 'k--', linewidth=1.0, alpha=0.6, label='Torus midline')

    # Annotation arrow for R (major radius)
    ax.plot([0, R*np.cos(np.pi/4)], [0, R*np.sin(np.pi/4)], [0, 0],
            color=oiVerm, linewidth=1.8, zorder=10)
    ax.scatter([0], [0], [0], color=oiBlack, s=15, zorder=11)
    ax.text(R*0.45*np.cos(np.pi/4), R*0.45*np.sin(np.pi/4), 3,
            r'$R = 50$ nm', fontsize=8, color=oiVerm, fontweight='bold')

    # Annotation arrow for r (minor radius)
    u0 = np.pi/4
    # from midline point to outer surface
    mid_x = R * np.cos(u0)
    mid_y = R * np.sin(u0)
    out_x = (R + r) * np.cos(u0)
    out_y = (R + r) * np.sin(u0)
    ax.plot([mid_x, out_x], [mid_y, out_y], [0, 0],
            color=oiGreen, linewidth=1.8, zorder=10)
    ax.text((mid_x + out_x)/2 + 2, (mid_y + out_y)/2 + 2, 2,
            r'$r = 8$ nm', fontsize=8, color=oiGreen, fontweight='bold')

    ax.view_init(elev=25, azim=45)
    set_3d_panes_white(ax)
    ax.set_xlabel('x (nm)', labelpad=5)
    ax.set_ylabel('y (nm)', labelpad=5)
    ax.set_zlabel('z (nm)', labelpad=3)
    ax.set_title('Option A* fat torus geometry', fontsize=10, fontweight='bold', pad=10)

    # Equalize aspect
    max_range = R + r + 5
    ax.set_xlim([-max_range, max_range])
    ax.set_ylim([-max_range, max_range])
    ax.set_zlim([-max_range/3, max_range/3])

    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'toroidal_geometry_3d.pdf'), bbox_inches='tight')
    plt.close()


# =========================================================================
# FIGURE 2: Mexican Hat 3D (Sec 7: TQGL Functional)
# =========================================================================
def fig_mexican_hat():
    print("  Generating mexican_hat_3d_catalyzon.pdf ...")
    alpha_gl = 1.0
    beta_gl = 0.5
    delta_min = np.sqrt(alpha_gl / (2 * beta_gl))  # = 1.0

    # Grid in Cartesian
    x = np.linspace(-1.8, 1.8, 120)
    y = np.linspace(-1.8, 1.8, 120)
    X, Y = np.meshgrid(x, y)
    R2 = X**2 + Y**2
    F = -alpha_gl * R2 + beta_gl * R2**2

    fig = plt.figure(figsize=(4.2, 3.8))
    ax = fig.add_subplot(111, projection='3d')

    # Surface
    surf = ax.plot_surface(X, Y, F, cmap='RdBu_r', alpha=0.75,
                           edgecolor='none', antialiased=True)

    # Minimum ring
    t = np.linspace(0, 2*np.pi, 100)
    xring = delta_min * np.cos(t)
    yring = delta_min * np.sin(t)
    F_min = -alpha_gl * delta_min**2 + beta_gl * delta_min**4
    zring = np.full_like(t, F_min)
    ax.plot(xring, yring, zring, color=oiGreen, linewidth=2.5, zorder=10,
            label=r'$|\Delta^*|$ minimum')

    # Unstable maximum at origin
    ax.scatter([0], [0], [0], color=oiVerm, s=40, zorder=11, edgecolors='k', linewidths=0.5)

    # Annotations
    ax.text(delta_min * np.cos(np.pi/6) + 0.15,
            delta_min * np.sin(np.pi/6) + 0.15,
            F_min - 0.08,
            r'$\Delta^* = \sqrt{\alpha/2\beta}$', fontsize=7.5, color=oiGreen)

    ax.text(-0.6, -0.6, F_min * 0.3,
            r'$U(1)$ manifold', fontsize=7.5, color=oiBlue,
            fontstyle='italic')

    ax.text(0.1, 0.1, 0.08, 'unstable', fontsize=6.5, color=oiVerm)

    ax.view_init(elev=30, azim=-60)
    set_3d_panes_white(ax)
    ax.set_xlabel(r'Re $\Delta$ (meV)', labelpad=5)
    ax.set_ylabel(r'Im $\Delta$ (meV)', labelpad=5)
    ax.set_zlabel(r'$\mathcal{F}$ (arb. units)', labelpad=3)
    ax.set_title('TQGL free energy landscape', fontsize=8, pad=6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'mexican_hat_3d_catalyzon.pdf'), bbox_inches='tight')
    plt.close()


# =========================================================================
# FIGURE 3: Hopfion Vortex 3D (Sec 8: Hopfions)
# =========================================================================
def fig_hopfion_vortex():
    print("  Generating hopfion_vortex_3d.pdf ...")
    # Generate N=5 mutually linked fibers as (p,q)=(1,5) torus knot strands
    # Parametrize on a torus with R_major=2, r_minor=1
    # Fiber k: offset by 2*pi*k/N in toroidal angle
    N = 5
    R_maj = 2.0
    r_min = 0.85
    colors = [oiBlue, oiVerm, oiOrange, oiGreen, oiSky]

    fig = plt.figure(figsize=(4.2, 4.0))
    ax = fig.add_subplot(111, projection='3d')

    t = np.linspace(0, 2*np.pi, 400)

    for k in range(N):
        phi_offset = 2 * np.pi * k / N
        # (1, N) torus knot fiber: winds once around the hole, N times around the tube
        # But for LINKED circles (not a single knot), we want N separate (1,1) curves
        # offset in the toroidal direction — these are the standard Hopf fibers
        # on a Clifford torus.
        #
        # Hopf fibers on a Clifford torus (after stereographic projection to R3):
        # For fiber at angle phi_0 on base S1:
        #   x(t) = (R + r*cos(t)) * cos(t + phi_0)
        #   y(t) = (R + r*cos(t)) * sin(t + phi_0)
        #   z(t) = r * sin(t)
        # Each fiber links every other fiber once.

        x = (R_maj + r_min * np.cos(t)) * np.cos(t + phi_offset)
        y = (R_maj + r_min * np.cos(t)) * np.sin(t + phi_offset)
        z = r_min * np.sin(t)

        ax.plot(x, y, z, color=colors[k], linewidth=2.0, alpha=0.9,
                label=f'Fiber {k+1}')

    ax.view_init(elev=20, azim=30)
    set_3d_panes_white(ax)

    # Annotations
    ax.text(0, 0, -2.0, r'$Q_H = 1/N = 0.2$', fontsize=8, color=oiBlack,
            ha='center', fontweight='bold')
    ax.text(0, 3.2, 1.0, f'$N = {N}$ linked fibers', fontsize=8, color=oiBlack,
            ha='center', fontstyle='italic')

    ax.set_xlabel('x', labelpad=4)
    ax.set_ylabel('y', labelpad=4)
    ax.set_zlabel('z', labelpad=3)
    ax.set_title('Hopf fibration: linked preimage curves', fontsize=10,
                 fontweight='bold', pad=10)

    # Clean up limits
    lim = R_maj + r_min + 0.5
    ax.set_xlim([-lim, lim])
    ax.set_ylim([-lim, lim])
    ax.set_zlim([-1.5, 1.5])

    ax.legend(fontsize=6, loc='upper right', framealpha=0.7, ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'hopfion_vortex_3d.pdf'), bbox_inches='tight')
    plt.close()


# =========================================================================
# FIGURE 4: Catalyzon Formation 3D (Sec 6)
# =========================================================================
def fig_catalyzon_formation():
    print("  Generating catalyzon_formation_3d.pdf ...")
    N = 5
    R_ring = 2.0  # radius of pentagon arrangement
    well_depth = 1.2
    well_width = 0.6  # controls steepness

    # Well positions (pentagon)
    angles_k = [2 * np.pi * k / N for k in range(N)]
    well_pos = [(R_ring * np.cos(a), R_ring * np.sin(a)) for a in angles_k]

    # Create grid
    x = np.linspace(-3.5, 3.5, 150)
    y = np.linspace(-3.5, 3.5, 150)
    X, Y = np.meshgrid(x, y)

    # Energy landscape: sum of inverted Gaussians (wells) + slight coupling saddles
    Z = np.zeros_like(X)
    for (xk, yk) in well_pos:
        dist2 = (X - xk)**2 + (Y - yk)**2
        Z += -well_depth * np.exp(-dist2 / (2 * well_width**2))

    # Add coupling bridges between adjacent wells (shallow saddles)
    for k in range(N):
        k_next = (k + 1) % N
        xm = (well_pos[k][0] + well_pos[k_next][0]) / 2
        ym = (well_pos[k][1] + well_pos[k_next][1]) / 2
        dist2 = (X - xm)**2 + (Y - ym)**2
        Z += -0.3 * np.exp(-dist2 / (2 * 0.5**2))

    # Shift baseline
    Z = Z - Z.min() * 0.1

    fig = plt.figure(figsize=(4.2, 3.8))
    ax = fig.add_subplot(111, projection='3d')

    # Surface with custom coloring: blue for wells, warm for ridges
    norm = plt.Normalize(Z.min(), Z.max())
    ax.plot_surface(X, Y, Z, cmap='coolwarm_r', alpha=0.78,
                    edgecolor='none', antialiased=True, norm=norm)

    # Well labels and bright-mode arrows
    arrow_len = 0.4
    for k, (xk, yk) in enumerate(well_pos):
        z_well = -well_depth + 0.05
        # Mode label
        ax.text(xk * 1.3, yk * 1.3, Z.min() * 0.2,
                f'Mode {k+1}', fontsize=6.5, color=oiBlack, ha='center')
        # Bright mode arrow (all pointing radially outward = in-phase)
        dx = arrow_len * np.cos(angles_k[k])
        dy = arrow_len * np.sin(angles_k[k])
        ax.quiver(xk, yk, Z.min() - 0.1, dx, dy, 0,
                  color=oiVerm, arrow_length_ratio=0.3, linewidth=1.5)

    # Coupling annotation on one bridge
    xm = (well_pos[0][0] + well_pos[1][0]) / 2
    ym = (well_pos[0][1] + well_pos[1][1]) / 2
    ax.text(xm + 0.3, ym + 0.3, Z.min() * 0.55,
            r'$g_0 = 0.098$ meV', fontsize=7, color=oiOrange, fontweight='bold')

    ax.view_init(elev=35, azim=45)
    set_3d_panes_white(ax)
    ax.set_xlabel('x (arb.)', labelpad=4)
    ax.set_ylabel('y (arb.)', labelpad=4)
    ax.set_zlabel('E (meV)', labelpad=3)
    ax.set_title('Catalyzon: 5-mode hybridization landscape', fontsize=8,
                 pad=6)

    plt.tight_layout()
    plt.savefig(os.path.join(OUTDIR, 'catalyzon_formation_3d.pdf'), bbox_inches='tight')
    plt.close()


# =========================================================================
# MAIN
# =========================================================================
if __name__ == '__main__':
    print(f"Output directory: {OUTDIR}")
    print("Generating QVC 3D figures ...")
    fig_toroidal_geometry()
    fig_mexican_hat()
    fig_hopfion_vortex()
    fig_catalyzon_formation()
    print("Done. All 4 PDFs written.")
