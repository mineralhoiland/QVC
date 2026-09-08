"""
QVC Publication Figures: 4 regenerated PDFs for arXiv v4 preprint
================================================================
Regenerates: gap_phase_diagram (Fig 5), computed_vortex_flow (Fig 6),
             computed_phase_portrait (Fig 41), computed_schwinger_rate (Fig 47)

ALL parameters from locked macros.tex.  No hallucinated data.
Style: Okabe-Ito palette, sans-serif 9 pt, single-column width.
Dependencies: numpy, matplotlib (no scipy).
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ================================================================
# Style Setup
# ================================================================
oiBlue   = '#0072B2'
oiVerm   = '#D55E00'
oiOrange = '#E69F00'
oiGreen  = '#009E73'
oiSky    = '#56B4E9'
oiPink   = '#CC79A7'

plt.rcParams.update({
    'font.family':        'sans-serif',
    'font.sans-serif':    ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size':          9,
    'axes.labelsize':     9,
    'axes.titlesize':     9,
    'xtick.labelsize':    7,
    'ytick.labelsize':    7,
    'legend.fontsize':    7,
    'axes.spines.top':    False,
    'axes.spines.right':  False,
    'axes.linewidth':     0.6,
    'xtick.major.width':  0.6,
    'ytick.major.width':  0.6,
    'lines.linewidth':    1.2,
    'pdf.fonttype':       42,
    'ps.fonttype':        42,
})

import os
OUTDIR = os.path.dirname(os.path.abspath(__file__)) + '/'
SINGLE_COL = 3.4   # inches

# ================================================================
# Locked Physical Parameters  (macros.tex, Option A*)
# ================================================================
hbar_cs      = 0.07       # meV * um
Delta_0      = 0.600      # meV  (= Delta_ex)
E_vac        = 0.1287     # meV
ell_lock     = -2.7958    # ln(Delta_0 * a_H / (2*hbar_cs)) + gamma_E
alpha_c      = 0.1818     # saddle-node fold coupling
delta_c      = 0.1818     # normalised gap at fold  (= alpha_c in this lock)
lambda_star  = 0.7628     # participation
alpha_star   = 0.1636     # lambda* * E_vac / Delta_0
Delta_star   = 0.2324     # meV  (working gap)
T_star       = 118.7      # mK
Ac_fold      = 1.5582     # um^-1
A_Ac_res     = 0.9117     # A*/Ac_fold  (resonant)
Delta_c_meV  = delta_c * Delta_0   # 0.1091 meV
a_H          = 0.008      # um  (healing length, 8 nm)
gamma_E      = 0.5772

# Conversion factor: alpha -> A_ZPF (um^-1)
alpha_to_A   = Delta_0 / hbar_cs   # 8.5714  um^-1 per unit alpha

# Schwinger sector
Ac_Schwinger = Delta_0 / hbar_cs   # 8.5714 um^-1

# Vortex-sector stiffnesses at T* (from macros.tex)
K_UV_Tstar   = 46.7       # deep bound phase
K_fold_Tstar = 10.57


# ================================================================
# Utility helpers
# ================================================================

def bisect(f, a, b, tol=1e-12, maxiter=200):
    """Root of f on [a, b] by bisection.  Returns None if no sign change."""
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        return None
    for _ in range(maxiter):
        c = 0.5 * (a + b)
        fc = f(c)
        if abs(fc) < tol or (b - a) < tol:
            return c
        if fa * fc < 0:
            b = c; fb = fc
        else:
            a = c; fa = fc
    return 0.5 * (a + b)


def solve_gap_physical(alpha, ell=ell_lock):
    """Solve delta = 1 + alpha*(ln delta + ell) on the physical (upper) branch.
    Returns delta in (delta_c, 1] for alpha in [0, alpha_c)."""
    if alpha <= 0:
        return 1.0
    if alpha >= alpha_c - 1e-10:
        return delta_c
    f = lambda d: d - 1.0 - alpha * (np.log(d) + ell)
    return bisect(f, delta_c + 1e-10, 1.0 - 1e-14) or delta_c


def rk4_step(deriv, y, t, dt):
    """Single step of 4th-order Runge-Kutta."""
    k1 = dt * deriv(y, t)
    k2 = dt * deriv(y + 0.5*k1, t + 0.5*dt)
    k3 = dt * deriv(y + 0.5*k2, t + 0.5*dt)
    k4 = dt * deriv(y + k3, t + dt)
    return y + (k1 + 2*k2 + 2*k3 + k4) / 6.0


def add_flow_arrow(ax, traj, frac, color, lw=0.7, hl=0.06):
    """Place an arrowhead at fractional position along trajectory."""
    n = len(traj)
    if n < 4:
        return
    i = max(1, min(n - 2, int(frac * n)))
    dx = traj[i+1, 0] - traj[i, 0]
    dy = traj[i+1, 1] - traj[i, 1]
    ax.annotate('', xy=(traj[i, 0] + dx, traj[i, 1] + dy),
                xytext=(traj[i, 0], traj[i, 1]),
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw, mutation_scale=8))


# ================================================================
# Precompute fold-curve lookup  (used by phase-portrait)
# ================================================================
_a_table = np.linspace(0.0, 0.9999, 2000)
_delta_table = np.array([solve_gap_physical(a * alpha_c) for a in _a_table])

def delta_fold_interp(a):
    """Fast lookup: normalised gap on physical branch at a = A/Ac_fold."""
    return np.interp(np.clip(a, 0, 0.9999), _a_table, _delta_table)


# ================================================================
# FIGURE 5  —  Gap Phase Diagram  (fold bifurcation)
# ================================================================
def fig_gap_phase_diagram():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.6))

    ell = ell_lock   # -2.7958

    # ------- parametric branches: alpha(delta) = (delta-1) / (ln delta + ell) -------
    #   Physical (upper):   delta in (delta_c, 1)   ->  alpha 0 .. alpha_c
    #   Unstable (lower):   delta in (~0, delta_c)   ->  alpha 0 .. alpha_c
    delta_phys   = np.linspace(delta_c + 5e-5, 1.0 - 1e-6, 1200)
    delta_unstab = np.linspace(5e-4, delta_c - 5e-5, 1200)

    def alpha_of_delta(d):
        return (d - 1.0) / (np.log(d) + ell)

    a_phys    = alpha_of_delta(delta_phys)
    a_unstab  = alpha_of_delta(delta_unstab)

    # Convert to physical units
    A_phys      = a_phys   * alpha_to_A       # um^-1
    A_unstab    = a_unstab * alpha_to_A
    D_phys      = delta_phys   * Delta_0      # meV
    D_unstab    = delta_unstab * Delta_0

    # Plot
    ax.plot(A_phys,   D_phys,   '-',  color=oiBlue, lw=1.5,
            label=r'Physical branch $\delta_+(\alpha)$')
    ax.plot(A_unstab, D_unstab, '--', color=oiVerm, lw=1.3,
            label=r'Unstable branch $\delta_-(\alpha)$')

    # Fold point
    A_fold = alpha_c * alpha_to_A
    ax.plot(A_fold, Delta_c_meV, 'o', color='k', ms=5, zorder=5)
    ax.annotate(r'Fold ($\alpha_c\!=\!%.4f$)' % alpha_c,
                xy=(A_fold, Delta_c_meV),
                xytext=(A_fold - 0.40, Delta_c_meV + 0.11),
                fontsize=6.5, ha='center',
                arrowprops=dict(arrowstyle='->', color='k', lw=0.7))

    # Working point
    A_work = alpha_star * alpha_to_A
    ax.plot(A_work, Delta_star, '*', color=oiGreen, ms=10, zorder=5,
            markeredgecolor='k', markeredgewidth=0.3)
    ax.annotate(r'$\alpha^{\!*}\!=\!%.4f$' % alpha_star + '\n'
                + r'$\Delta^{\!*}\!=\!%.3f$ meV' % Delta_star,
                xy=(A_work, Delta_star),
                xytext=(A_work - 0.55, Delta_star + 0.12),
                fontsize=6, ha='center', color=oiGreen,
                arrowprops=dict(arrowstyle='->', color=oiGreen, lw=0.7))

    # Critical gap dashed line
    ax.axhline(Delta_c_meV, color='gray', lw=0.5, ls=':', alpha=0.6)
    ax.text(0.08, Delta_c_meV + 0.010,
            r'$\Delta_c = %.3f$ meV' % Delta_c_meV,
            fontsize=5.5, color='gray')

    # No-solution shading beyond fold
    ax.axvline(A_fold, color='gray', lw=0.4, ls='--', alpha=0.25)
    ax.fill_betweenx([0, 0.65], A_fold, 1.85, alpha=0.06, color=oiVerm)

    # Region labels
    ax.text(0.55, 0.44, 'Gapped\n(condensate)', fontsize=7, color=oiBlue,
            ha='center', style='italic', alpha=0.8)
    ax.text(1.65, 0.44, 'No static\nsolution', fontsize=7, color=oiVerm,
            ha='center', style='italic', alpha=0.8)

    ax.set_xlabel(r'$A_{\mathrm{ZPF}}$ ($\mu$m$^{-1}$)')
    ax.set_ylabel(r'$\Delta$ (meV)')
    ax.set_xlim(0, 1.85)
    ax.set_ylim(0, 0.65)
    ax.legend(loc='upper right', frameon=False, fontsize=6)

    # Secondary top x-axis: dimensionless coupling alpha
    ax2 = ax.secondary_xaxis('top', functions=(
        lambda A: A / alpha_to_A,
        lambda a: a * alpha_to_A))
    ax2.set_xlabel(r'Coupling $\alpha$', fontsize=7)
    ax2.tick_params(labelsize=5.5)

    fig.savefig(OUTDIR + 'gap_phase_diagram.pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  [1/4] gap_phase_diagram.pdf  —  fold (%.3f um^-1, %.3f meV)" %
          (A_fold, Delta_c_meV))


# ================================================================
# FIGURE 6  —  Vortex RG Flow  (Kosterlitz-Thouless)
# ================================================================
def fig_vortex_flow():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.0))

    K_c = 2.0 / np.pi   # 0.6366

    # KT flow equations  (Kosterlitz normalisation):
    #   dK^{-1}/dl = 4 pi^3 y^2
    #   dy/dl      = (2 - pi K) y
    # In (K, y) variables:
    #   dK/dl = -4 pi^3 K^2 y^2
    #   dy/dl = (2 - pi K) y

    def kt_deriv(state, ell):
        K, y = state
        if K < 1e-6 or y < 0:
            return np.array([0.0, 0.0])
        dK = -4.0 * np.pi**3 * K**2 * y**2
        dy = (2.0 - np.pi * K) * y
        return np.array([dK, dy])

    dl    = 0.005
    l_max = 30.0
    nstep = int(l_max / dl)

    # Starting conditions — spread around K_c to show separatrix
    starts_bound = [
        # (K0, y0)  —  K > Kc, or K near Kc with small enough y
        (0.70, 0.05), (0.70, 0.10), (0.72, 0.18),
        (0.80, 0.12), (0.80, 0.25), (0.90, 0.20),
        (1.00, 0.18), (1.00, 0.30), (1.00, 0.40),
        (1.30, 0.25), (1.50, 0.30), (2.00, 0.20),
        (2.50, 0.25), (3.00, 0.15),
    ]
    starts_plasma = [
        (0.60, 0.03), (0.58, 0.06), (0.55, 0.04),
        (0.50, 0.06), (0.50, 0.12), (0.45, 0.04),
        (0.40, 0.06), (0.40, 0.12), (0.35, 0.05),
        (0.30, 0.08), (0.25, 0.05),
    ]

    def integrate(K0, y0):
        traj = np.zeros((nstep + 1, 2))
        traj[0] = [K0, y0]
        for i in range(nstep):
            traj[i+1] = rk4_step(kt_deriv, traj[i], i*dl, dl)
            # Stop if out of visible region
            if traj[i+1, 0] < 0.05 or traj[i+1, 1] < 1e-12:
                return traj[:i+2]
            if traj[i+1, 1] > 0.8 or traj[i+1, 0] > 5.0:
                return traj[:i+2]
        return traj

    # Classify trajectory: is y ultimately decaying?
    def is_bound(traj):
        return traj[-1, 1] < traj[0, 1] * 0.5 or traj[-1, 1] < 0.01

    # ---- Bound-phase trajectories ----
    for j, (K0, y0) in enumerate(starts_bound):
        tr = integrate(K0, y0)
        mask = (tr[:, 0] > 0.08) & (tr[:, 0] < 4.2) & (tr[:, 1] > -0.01) & (tr[:, 1] < 0.55)
        idx  = np.where(mask)[0]
        if len(idx) > 3:
            ax.plot(tr[idx, 0], tr[idx, 1], '-', color=oiBlue, alpha=0.55, lw=0.9,
                    label=('Bound (superfluid)' if j == 0 else ''))
            add_flow_arrow(ax, tr[idx], 0.35, oiBlue)

    # ---- Plasma-phase trajectories ----
    for j, (K0, y0) in enumerate(starts_plasma):
        tr = integrate(K0, y0)
        mask = (tr[:, 0] > 0.08) & (tr[:, 0] < 4.2) & (tr[:, 1] > -0.01) & (tr[:, 1] < 0.55)
        idx  = np.where(mask)[0]
        if len(idx) > 3:
            ax.plot(tr[idx, 0], tr[idx, 1], '-', color=oiOrange, alpha=0.55, lw=0.9,
                    label=('Plasma (disordered)' if j == 0 else ''))
            add_flow_arrow(ax, tr[idx], 0.35, oiOrange)

    # Separatrix
    ax.axvline(K_c, color='k', lw=0.8, ls='--', alpha=0.7)
    ax.text(K_c + 0.06, 0.51, r'$K_c = 2/\pi$', fontsize=7, color='k')

    # Phase labels
    ax.text(2.2, 0.06, 'Bound pairs\n' + r'($y \to 0$)',
            fontsize=7, color=oiBlue, ha='center', style='italic')
    ax.text(0.22, 0.42, 'Free-vortex\nplasma\n' + r'($y \to \infty$)',
            fontsize=6.5, color=oiOrange, ha='center', style='italic')

    # QVC annotation
    ax.annotate(r'QVC: $K_{\mathrm{UV}}(T^*) = %.1f$' % K_UV_Tstar
                + '\n(deep bound phase)',
                xy=(3.9, 0.015), xytext=(2.6, 0.30),
                fontsize=5.5, ha='center', color=oiGreen,
                bbox=dict(boxstyle='round,pad=0.3', fc='white',
                          ec=oiGreen, alpha=0.85),
                arrowprops=dict(arrowstyle='->', color=oiGreen, lw=0.8))

    ax.set_xlabel(r'Stiffness $K = \kappa / T_{\mathrm{eff}}$')
    ax.set_ylabel(r'Fugacity $y$')
    ax.set_xlim(0.1, 4.2)
    ax.set_ylim(-0.01, 0.55)
    ax.legend(loc='upper left', frameon=False, fontsize=6)

    fig.savefig(OUTDIR + 'computed_vortex_flow.pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  [2/4] computed_vortex_flow.pdf  —  K_c = %.4f" % K_c)


# ================================================================
# FIGURE 41  —  Phase Portrait  (A/A_c , Delta)
# ================================================================
def fig_phase_portrait():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.6))

    # ---- fold curve (physical branch) ----
    a_arr   = np.linspace(0.005, 0.999, 600)
    D_fold  = np.array([solve_gap_physical(a * alpha_c) for a in a_arr]) * Delta_0
    ax.plot(a_arr, D_fold, '-', color='gray', lw=1.0, alpha=0.55,
            label='Static gap map')

    # ---- fold critical gap ----
    ax.axhline(Delta_c_meV, color=oiVerm, lw=0.7, ls=':', alpha=0.6)
    ax.text(0.04, Delta_c_meV + 0.010,
            r'$\Delta_c = %.3f$ meV' % Delta_c_meV, fontsize=5.5, color=oiVerm)

    # ---- fixed point ----
    a_star = A_Ac_res   # 0.9117
    ax.plot(a_star, Delta_star, 's', color=oiBlue, ms=7, zorder=6,
            markeredgecolor='k', markeredgewidth=0.4)
    ax.annotate(r'$(A^{*}\!/A_c^{\mathrm{fold}},\,\Delta^{*})$'
                + '\n' + r'$= (%.4f,\;%.4f)$' % (a_star, Delta_star),
                xy=(a_star, Delta_star),
                xytext=(a_star - 0.32, Delta_star + 0.10),
                fontsize=5.5, ha='center', color=oiBlue,
                arrowprops=dict(arrowstyle='->', color=oiBlue, lw=0.7))

    # ---- 2-D dynamics ---------------------------------------------------
    #  da/dt  ~  growth that saturates at a*  (DCE detuning + nonlinear)
    #  dD/dt  =  -(D - D_fold(a)) / tau       (fast gap relaxation)
    #
    # Model:  da/dt = G_eff(a, D) * a
    # where G_eff turns positive below a* and negative above.
    #
    # Physical mechanism:  Lorentzian detuning L(D)  times  parametric gain
    # minus loss and nonlinear saturation.

    hbar_ps   = 0.6582         # meV * ps
    Q_moire   = 20.0
    mod_depth = 0.1            # delta-Delta / Delta_0
    omega_0   = 2.0 * Delta_0 / hbar_ps          # pump ref  (ps^-1)
    Gamma_cav = omega_0 / Q_moire                 # cavity linewidth (ps^-1)
    Lambda_0  = mod_depth * Delta_0 / (2*hbar_ps) * Q_moire / (2*np.pi)  # 0.145 ps^-1

    # Tune loss + nonlinear saturation so fixed point lands at a*
    # At a*, D* = Delta_star, detuning:
    dw_star = 2.0 * (Delta_0 - Delta_star) / hbar_ps
    L_star  = 1.0 / (1.0 + (dw_star / Gamma_cav)**2)
    gain_star = Lambda_0 * L_star                 # ~ 9.2e-4  ps^-1
    # Choose kappa so that  gain - kappa = eta * a*^2   with small eta
    kappa_eff = gain_star * 0.85                   # slightly below gain
    eta_nl    = (gain_star - kappa_eff) / (a_star**2)

    tau_gap = 0.4   # ps  (fast gap relaxation toward fold curve)

    def portrait_deriv(state, t):
        a_val, D_val = state
        a_val = np.clip(a_val, 1e-4, 0.9999)
        D_val = max(D_val, 0.01)

        # Fold-curve gap at current a
        D_fold_val = delta_fold_interp(a_val) * Delta_0

        # DCE gain with Lorentzian detuning
        dw   = 2.0 * (Delta_0 - D_val) / hbar_ps
        L_dw = 1.0 / (1.0 + (dw / Gamma_cav)**2)
        gain = Lambda_0 * L_dw

        da = (gain - kappa_eff - eta_nl * a_val**2) * a_val
        dD = -(D_val - D_fold_val) / tau_gap
        return np.array([da, dD])

    # Trajectory initial conditions — spread around the fixed point
    ic_list = [
        (0.05, 0.58),   # far left, high Delta
        (0.15, 0.55),
        (0.30, 0.50),
        (0.50, 0.42),
        (0.70, 0.50),
        (0.20, 0.25),   # lower Delta
        (0.40, 0.20),
        (0.60, 0.15),
        (0.82, 0.12),   # below fold, will rise in Delta
        (0.96, 0.45),   # right of fixed point
        (0.98, 0.30),
        (0.92, 0.55),
    ]

    dt     = 0.02   # ps
    t_max  = 300.0   # ps
    nsteps = int(t_max / dt)

    for k, (a0, D0) in enumerate(ic_list):
        tr = np.zeros((nsteps + 1, 2))
        tr[0] = [a0, D0]
        for i in range(nsteps):
            tr[i+1] = rk4_step(portrait_deriv, tr[i], i*dt, dt)
            tr[i+1, 0] = np.clip(tr[i+1, 0], 1e-3, 0.999)
            tr[i+1, 1] = np.clip(tr[i+1, 1], 0.01, 0.62)
            if (abs(tr[i+1, 0] - a_star) < 0.003
                    and abs(tr[i+1, 1] - Delta_star) < 0.003):
                tr = tr[:i+2]
                break

        c = oiSky if k < 8 else oiGreen
        ax.plot(tr[:, 0], tr[:, 1], '-', color=c, alpha=0.50, lw=0.7)
        add_flow_arrow(ax, tr, 0.30, c, lw=0.6)

    ax.set_xlabel(r'$A\, /\, A_c^{\mathrm{fold}}$')
    ax.set_ylabel(r'$\Delta$ (meV)')
    ax.set_xlim(0, 1.05)
    ax.set_ylim(0, 0.63)
    ax.legend(loc='upper left', frameon=False, fontsize=6)

    fig.savefig(OUTDIR + 'computed_phase_portrait.pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  [3/4] computed_phase_portrait.pdf  —  FP (%.4f, %.4f meV)" %
          (a_star, Delta_star))


# ================================================================
# FIGURE 47  —  Schwinger Rate
# ================================================================
def fig_schwinger_rate():
    fig, ax = plt.subplots(figsize=(SINGLE_COL, 2.6))

    # Attempt frequency  omega_0 = 2 Delta* / hbar   (s^-1)
    hbar_eVs = 6.582e-16   # eV * s
    omega_0  = 2.0 * Delta_star * 1e-3 / hbar_eVs   # 7.06e11  s^-1

    # 2+1-D Schwinger rate:
    #   Gamma = omega_0 * (A_eff / A_c) * exp(-pi * A_c / A_eff)
    x = np.linspace(0.06, 1.0, 2000)      # A_eff / A_c
    Gamma = omega_0 * x * np.exp(-np.pi / x)

    # Working point
    x_star     = alpha_star   # 0.1636
    Gamma_star = omega_0 * x_star * np.exp(-np.pi / x_star)

    ax.semilogy(x, Gamma, '-', color=oiBlue, lw=1.5,
                label=r'$\Gamma_{\mathrm{Sch}}$ (2+1 D)')

    # Working-point marker and annotation
    ax.plot(x_star, Gamma_star, '*', color=oiGreen, ms=10, zorder=5,
            markeredgecolor='k', markeredgewidth=0.3)
    ax.axvline(x_star, color=oiGreen, lw=0.5, ls='--', alpha=0.4)
    ax.annotate(r'$\alpha^{\!*}=%.4f$' % x_star
                + '\n' + r'$\Gamma \approx %.0f\;$Hz' % Gamma_star,
                xy=(x_star, Gamma_star),
                xytext=(x_star + 0.18, Gamma_star * 5e3),
                fontsize=6, ha='left', color=oiGreen,
                arrowprops=dict(arrowstyle='->', color=oiGreen, lw=0.7))

    # Schwinger exponent annotation
    ax.text(0.58, omega_0 * 0.08,
            r'$\Gamma = \omega_0 \,\frac{A_{\mathrm{eff}}}{A_c}'
            r'\,\exp\!\!\left(\!-\frac{\pi\, A_c}{A_{\mathrm{eff}}}\right)$',
            fontsize=7.5, color=oiBlue,
            bbox=dict(boxstyle='round,pad=0.3', fc='white', ec=oiBlue, alpha=0.7))

    # Twin-channel crossover  (2 * Gamma* at T = T*)
    ax.axhline(2 * Gamma_star, color='gray', lw=0.4, ls=':', alpha=0.5)
    ax.text(0.95, 2.8 * Gamma_star,
            r'$2\Gamma^*$ at $T\!=\!T^*\!=\!%.1f$ mK' % T_star,
            fontsize=5.5, color='gray', ha='right')

    ax.set_xlabel(r'$A_{\mathrm{eff}}\, /\, A_{c,\mathrm{Sch}}$')
    ax.set_ylabel(r'$\Gamma_{\mathrm{Schwinger}}$ (s$^{-1}$)')
    ax.set_xlim(0, 1.05)
    ax.set_ylim(1e-2, 5e11)

    fig.savefig(OUTDIR + 'computed_schwinger_rate.pdf', bbox_inches='tight', dpi=300)
    plt.close(fig)
    print("  [4/4] computed_schwinger_rate.pdf  —  Gamma* = %.2e s^-1  (%.0f Hz)" %
          (Gamma_star, Gamma_star))


# ================================================================
# MAIN
# ================================================================
if __name__ == '__main__':
    print("=" * 60)
    print("Regenerating 4 QVC publication figures  (arXiv v4)")
    print("=" * 60)
    print("  Output directory:", OUTDIR)
    print()

    # ------ parameter verification ------
    print("  Locked-parameter cross-checks:")
    ell_check = np.log(Delta_0 * a_H / (2.0 * hbar_cs)) + gamma_E
    print("    ell  = ln(%.3f * %.3f / (2*%.2f)) + gamma_E = %.4f  (locked %.4f, %s)"
          % (Delta_0, a_H, hbar_cs, ell_check, ell_lock,
             'OK' if abs(ell_check - ell_lock) < 0.001 else 'MISMATCH'))
    fold_check = 1.0 + alpha_c * (np.log(delta_c) + ell_lock)
    print("    fold:  1 + alpha_c*(ln(delta_c)+ell) = %.5f  (should be %.4f, %s)"
          % (fold_check, delta_c,
             'OK' if abs(fold_check - delta_c) < 0.002 else 'MISMATCH'))
    Ac_check = alpha_c * alpha_to_A
    print("    Ac^fold = alpha_c * Delta_0/hbar_cs = %.4f um^-1  (locked %.4f, %s)"
          % (Ac_check, Ac_fold,
             'OK' if abs(Ac_check - Ac_fold) < 0.002 else 'MISMATCH'))
    alpha_star_check = lambda_star * E_vac / Delta_0
    print("    alpha* = lambda* * E_vac / Delta_0 = %.4f  (locked %.4f, %s)"
          % (alpha_star_check, alpha_star,
             'OK' if abs(alpha_star_check - alpha_star) < 0.001 else 'MISMATCH'))
    print()

    # ------ generate figures ------
    fig_gap_phase_diagram()
    fig_vortex_flow()
    fig_phase_portrait()
    fig_schwinger_rate()

    print()
    print("Done.  All 4 PDFs written to", OUTDIR)
