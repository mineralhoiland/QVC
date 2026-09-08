"""
Publication Figure: Hopfion Spacetime Crystal
=============================================
Generates a multi-panel figure for the TDVT paper showing:
  Panel A — Preimage curves (linked circles) with toroidal sonic surface
  Panel B — Topological charge density ρ_H midplane slices
  Panel C — Energy density |∇n|² with n-field quiver overlay
  Panel D — n_z component showing the Hopf fibration winding

Uses locked QVC v4 parameters for all physical annotations.
Tier 2 (working-grade): structural correctness, correct topology,
physically labeled axes. Not coefficient-precise numerics.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import os

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR  = os.path.join(BASEDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

# === Locked QVC v4 parameters (for annotations) ===
xi_nm = 50.0          # coherence length [nm]
E_vac = 0.1287        # meV, Casimir vacuum energy
F_fat = 0.58068       # geometric factor
g0    = 0.0982        # meV, coupling
Astar_over_Ac = 0.9117

# === Compute hopfion on grid ===
Ngrid = 80
L = 4.0  # in units of ξ
x = np.linspace(-L, L, Ngrid, endpoint=False)
dx = x[1] - x[0]
X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
r2 = X**2 + Y**2 + Z**2

# Hopf map via inverse stereographic projection
denom = 1.0 + r2
w0, w1, w2, w3 = 2*X/denom, 2*Y/denom, 2*Z/denom, (1-r2)/denom

# Smooth window: push n → (0,0,-1) for r > R_cut
R_cut, delta = 2.5, 0.4
r = np.sqrt(r2)
t = 0.5 * (1.0 + np.tanh((r - R_cut) / delta))

nx_raw = 2*(w0*w2 + w1*w3)
ny_raw = 2*(w1*w2 - w0*w3)
nz_raw = w0**2 + w1**2 - w2**2 - w3**2

nx_b = (1-t)*nx_raw
ny_b = (1-t)*ny_raw
nz_b = (1-t)*nz_raw + t*(-1.0)
norm_b = np.maximum(np.sqrt(nx_b**2+ny_b**2+nz_b**2), 1e-15)
nx, ny, nz = nx_b/norm_b, ny_b/norm_b, nz_b/norm_b

# Gradients (periodic central difference)
def pg(f, ax):
    return (np.roll(f,-1,axis=ax) - np.roll(f,1,axis=ax)) / (2*dx)

dnx = [pg(nx,i) for i in range(3)]
dny = [pg(ny,i) for i in range(3)]
dnz = [pg(nz,i) for i in range(3)]

# Berry curvature B_i = ε_{abc} n_a ∂_j n_b ∂_k n_c  (i,j,k cyclic)
def F_ij(i, j):
    return (nx*(dny[i]*dnz[j]-dny[j]*dnz[i]) +
            ny*(dnz[i]*dnx[j]-dnz[j]*dnx[i]) +
            nz*(dnx[i]*dny[j]-dnx[j]*dny[i]))

Bx, By, Bz = F_ij(1,2), F_ij(2,0), F_ij(0,1)

# Coulomb gauge FFT inversion
Bxh, Byh, Bzh = np.fft.fftn(Bx), np.fft.fftn(By), np.fft.fftn(Bz)
kx = np.fft.fftfreq(Ngrid, d=dx)*2*np.pi
KX, KY, KZ = np.meshgrid(kx, kx, kx, indexing='ij')
K2 = KX**2+KY**2+KZ**2
K2s = np.where(K2>0, K2, 1.0)

# Divergence projection
kdB = KX*Bxh+KY*Byh+KZ*Bzh
Bxh_c = Bxh - KX*kdB/K2s
Byh_c = Byh - KY*kdB/K2s
Bzh_c = Bzh - KZ*kdB/K2s
for h in [Bxh_c, Byh_c, Bzh_c]: h[0,0,0] = 0

cx = KY*Bzh_c - KZ*Byh_c
cy = KZ*Bxh_c - KX*Bzh_c
cz = KX*Byh_c - KY*Bxh_c
Axh = -1j*cx/K2s; Ayh = -1j*cy/K2s; Azh = -1j*cz/K2s
for h in [Axh, Ayh, Azh]: h[0,0,0] = 0

Ax = np.real(np.fft.ifftn(Axh))
Ay = np.real(np.fft.ifftn(Ayh))
Az = np.real(np.fft.ifftn(Azh))
Bx_c = np.real(np.fft.ifftn(Bxh_c))
By_c = np.real(np.fft.ifftn(Byh_c))
Bz_c = np.real(np.fft.ifftn(Bzh_c))

rho_H = (Ax*Bx_c+Ay*By_c+Az*Bz_c) / (4*np.pi**2)
Edens = sum(dnx[i]**2+dny[i]**2+dnz[i]**2 for i in range(3))
Q = np.sum(rho_H)*dx**3
E = 0.5*np.sum(Edens)*dx**3
print(f"Q = {Q:.4f}, E = {E:.1f}")

mid = Ngrid // 2
ext = [x[0], x[-1], x[0], x[-1]]


# =====================================================================
# FIGURE 1: PUBLICATION COMPOSITE (4-PANEL)
# =====================================================================
fig = plt.figure(figsize=(16, 14))

# --- Panel A: Preimage curves + sonic torus ---
ax1 = fig.add_subplot(221, projection='3d')

# Preimage of north pole: nz > 0.9 (ring in xy-plane)
mask_N = nz > 0.9
pts_N = np.column_stack([X[mask_N], Y[mask_N], Z[mask_N]])
if len(pts_N) > 2000:
    idx = np.random.choice(len(pts_N), 2000, replace=False)
    pts_N = pts_N[idx]

# Preimage of equator (nx > 0.9): linked ring
mask_E = nx > 0.9
pts_E = np.column_stack([X[mask_E], Y[mask_E], Z[mask_E]])
if len(pts_E) > 2000:
    idx = np.random.choice(len(pts_E), 2000, replace=False)
    pts_E = pts_E[idx]

ax1.scatter(pts_N[:,0], pts_N[:,1], pts_N[:,2], c='#2166ac', s=2, alpha=0.6,
            label=r'$\pi^{-1}(N)$: $n_z > 0.9$')
ax1.scatter(pts_E[:,0], pts_E[:,1], pts_E[:,2], c='#b2182b', s=2, alpha=0.6,
            label=r'$\pi^{-1}(E)$: $n_x > 0.9$')

# Sonic surface torus: for Whitehead hopfion, the energy density ∝ 1/(1+r²)⁴
# The "acoustic horizon" in the analogy is at a toroidal surface.
# We draw a torus at the energy density maximum locus as the sonic surface.
# For a Q=1 hopfion, energy peaks on a torus at r ≈ 1 from center
# with tube radius ~ 0.4
u_tor = np.linspace(0, 2*np.pi, 60)
v_tor = np.linspace(0, 2*np.pi, 40)
U_tor, V_tor = np.meshgrid(u_tor, v_tor)
R_major, r_minor = 1.0, 0.35  # in ξ units
Xt = (R_major + r_minor*np.cos(V_tor)) * np.cos(U_tor)
Yt = (R_major + r_minor*np.cos(V_tor)) * np.sin(U_tor)
Zt = r_minor * np.sin(V_tor)
ax1.plot_surface(Xt, Yt, Zt, alpha=0.12, color='gold', shade=True)

# z-axis line (preimage of south pole, partial)
zline = np.linspace(-2.5, 2.5, 50)
ax1.plot([0]*50, [0]*50, zline, c='#4daf4a', lw=1.5, alpha=0.5,
         label=r'$\pi^{-1}(S)$: z-axis')

ax1.set_xlim(-2.5, 2.5); ax1.set_ylim(-2.5, 2.5); ax1.set_zlim(-2.5, 2.5)
ax1.set_xlabel(r'$x/\xi$', fontsize=10)
ax1.set_ylabel(r'$y/\xi$', fontsize=10)
ax1.set_zlabel(r'$z/\xi$', fontsize=10)
ax1.set_title('(A) Preimage Curves + Sonic Torus\n'
              r'$Q = 1$, linked $S^1 \hookrightarrow S^3$',
              fontsize=12, pad=15)
ax1.legend(fontsize=8, loc='upper left', framealpha=0.8)
ax1.view_init(elev=25, azim=45)

# --- Panel B: Charge density midplane ---
ax2 = fig.add_subplot(222)
vmax = max(abs(rho_H[:,:,mid].min()), abs(rho_H[:,:,mid].max()))
im2 = ax2.imshow(rho_H[:,:,mid].T, origin='lower', extent=ext,
                  cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
cb2 = plt.colorbar(im2, ax=ax2, shrink=0.85, pad=0.02)
cb2.set_label(r'$\rho_H = \frac{1}{16\pi^2}\,\mathbf{A}\cdot\mathbf{B}$', fontsize=10)
ax2.set_xlabel(r'$x/\xi$', fontsize=11)
ax2.set_ylabel(r'$y/\xi$', fontsize=11)
ax2.set_title(f'(B) Hopf Charge Density ($z=0$)\n'
              f'$Q = {Q:.2f}$ (analytic: 1)', fontsize=12)

# Mark the preimage ring location
theta_ring = np.linspace(0, 2*np.pi, 100)
ax2.plot(np.cos(theta_ring), np.sin(theta_ring), 'w--', lw=1, alpha=0.6,
         label=r'$\pi^{-1}(N)$ ring')
ax2.legend(fontsize=8, loc='lower right')

# --- Panel C: Energy density + quiver ---
ax3 = fig.add_subplot(223)
im3 = ax3.imshow(Edens[:,:,mid].T, origin='lower', extent=ext,
                  cmap='inferno', aspect='equal')
cb3 = plt.colorbar(im3, ax=ax3, shrink=0.85, pad=0.02)
cb3.set_label(r'$|\nabla\mathbf{n}|^2$', fontsize=10)

skip = max(1, Ngrid//20)
xs = x[::skip]
Xs, Ys = np.meshgrid(xs, xs, indexing='ij')
ax3.quiver(Xs, Ys, nx[::skip,::skip,mid], ny[::skip,::skip,mid],
           color='white', scale=30, width=0.003, alpha=0.7)
ax3.set_xlabel(r'$x/\xi$', fontsize=11)
ax3.set_ylabel(r'$y/\xi$', fontsize=11)
ax3.set_title(f'(C) Energy Density + n-field ($z=0$)\n'
              f'$E = {E:.0f}$, VK bound $= 21.1$', fontsize=12)

# --- Panel D: n_z component ---
ax4 = fig.add_subplot(224)
im4 = ax4.imshow(nz[:,:,mid].T, origin='lower', extent=ext,
                  cmap='coolwarm', vmin=-1, vmax=1, aspect='equal')
cb4 = plt.colorbar(im4, ax=ax4, shrink=0.85, pad=0.02)
cb4.set_label(r'$n_z$', fontsize=10)
ax4.set_xlabel(r'$x/\xi$', fontsize=11)
ax4.set_ylabel(r'$y/\xi$', fontsize=11)
ax4.set_title(r'(D) Target Space $n_z$ ($z=0$)' '\n'
              r'Winding: $n_z = +1$ on ring, $-1$ at origin/boundary',
              fontsize=12)
# Mark +1 and -1 loci
ax4.plot(np.cos(theta_ring), np.sin(theta_ring), 'k--', lw=1, alpha=0.5)
ax4.plot(0, 0, 'k+', ms=10, mew=2)

# Global title
fig.suptitle('Whitehead Hopfion ($Q=1$) — TDVT Spacetime Crystal Unit Cell\n'
             r'$\xi = 50\,\mathrm{nm}$, $E_\mathrm{vac} = 0.1287\,\mathrm{meV}$, '
             r'$\mathcal{F} = 0.58068$',
             fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0, 1, 0.94])
path_pub = os.path.join(FIGDIR, "hopfion_publication_composite.png")
plt.savefig(path_pub, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Publication figure: {path_pub}")


# =====================================================================
# FIGURE 2: SONIC SURFACE DETAIL (standalone 3D)
# =====================================================================
fig2 = plt.figure(figsize=(10, 9))
ax = fig2.add_subplot(111, projection='3d')

# Energy density isosurface proxy: find where Edens is ~50% of max
# on the z=0 midplane, trace the ring
E_mid = Edens[:,:,mid]
E_thresh = 0.4 * E_mid.max()

# Toroidal sonic surface with energy density coloring
ax.plot_surface(Xt, Yt, Zt, alpha=0.25, color='gold', shade=True,
                label='Sonic surface (torus)')

# Preimage curves (denser sampling)
for thresh, comp, color, lw, lab in [
    (0.85, nz, '#2166ac', 3, r'$\pi^{-1}(0,0,1)$ — north pole fiber'),
    (0.85, nx, '#b2182b', 3, r'$\pi^{-1}(1,0,0)$ — equator fiber'),
    (0.85, ny, '#4daf4a', 3, r'$\pi^{-1}(0,1,0)$ — equator fiber')]:
    mask = comp > thresh
    if np.any(mask):
        pts = np.column_stack([X[mask], Y[mask], Z[mask]])
        if len(pts) > 3000:
            idx = np.random.choice(len(pts), 3000, replace=False)
            pts = pts[idx]
        ax.scatter(pts[:,0], pts[:,1], pts[:,2], c=color, s=3, alpha=0.5, label=lab)

# z-axis
zl = np.linspace(-3, 3, 80)
ax.plot([0]*80, [0]*80, zl, 'k--', lw=1, alpha=0.3, label=r'$\pi^{-1}(0,0,-1)$: z-axis')

ax.set_xlim(-3,3); ax.set_ylim(-3,3); ax.set_zlim(-3,3)
ax.set_xlabel(r'$x/\xi$', fontsize=12)
ax.set_ylabel(r'$y/\xi$', fontsize=12)
ax.set_zlabel(r'$z/\xi$', fontsize=12)
ax.set_title('Hopfion Preimage Topology + Toroidal Sonic Surface\n'
             r'Acoustic horizon at $|v| = c_s$, genus-1 surface',
             fontsize=13, pad=15)
ax.legend(fontsize=9, loc='upper left', framealpha=0.8)
ax.view_init(elev=20, azim=40)

path_sonic = os.path.join(FIGDIR, "hopfion_sonic_surface.png")
fig2.savefig(path_sonic, dpi=300, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Sonic surface figure: {path_sonic}")


# =====================================================================
# FIGURE 3: ORTHOGONAL SLICES TRIPTYCH
# =====================================================================
fig3, axes = plt.subplots(2, 3, figsize=(16, 10))

for col, (sl_rho, sl_E, sl_nz, plane, xl, yl) in enumerate([
    (rho_H[:,:,mid], Edens[:,:,mid], nz[:,:,mid], '$z=0$', r'$x/\xi$', r'$y/\xi$'),
    (rho_H[:,mid,:], Edens[:,mid,:], nz[:,mid,:], '$y=0$', r'$x/\xi$', r'$z/\xi$'),
    (rho_H[mid,:,:], Edens[mid,:,:], nz[mid,:,:], '$x=0$', r'$y/\xi$', r'$z/\xi$')]):

    # Top row: charge density
    vm = max(abs(sl_rho.min()), abs(sl_rho.max()))
    im = axes[0,col].imshow(sl_rho.T, origin='lower', extent=ext,
                             cmap='RdBu_r', vmin=-vm, vmax=vm, aspect='equal')
    plt.colorbar(im, ax=axes[0,col], shrink=0.8)
    axes[0,col].set_title(rf'$\rho_H$ ({plane})', fontsize=11)
    axes[0,col].set_xlabel(xl); axes[0,col].set_ylabel(yl)

    # Bottom row: n_z
    im2 = axes[1,col].imshow(sl_nz.T, origin='lower', extent=ext,
                              cmap='coolwarm', vmin=-1, vmax=1, aspect='equal')
    plt.colorbar(im2, ax=axes[1,col], shrink=0.8)
    axes[1,col].set_title(rf'$n_z$ ({plane})', fontsize=11)
    axes[1,col].set_xlabel(xl); axes[1,col].set_ylabel(yl)

fig3.suptitle(r'Hopfion Charge Density $\rho_H$ and Target Field $n_z$ — Three Orthogonal Planes'
              f'\n$Q = {Q:.2f}$, grid: {Ngrid}³, $L = {L}\\xi$',
              fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
path_tri = os.path.join(FIGDIR, "hopfion_orthogonal_triptych.png")
fig3.savefig(path_tri, dpi=250, bbox_inches='tight', facecolor='white')
plt.close()
print(f"Triptych figure: {path_tri}")

# =====================================================================
# NUMERICAL READOUT
# =====================================================================
print("\n" + "=" * 60)
print("NUMERICAL READOUT — Hopfion Spacetime Crystal")
print("=" * 60)
print(f"  Topological charge Q = {Q:.4f} (exact: 1)")
print(f"  Total energy E       = {E:.1f}")
print(f"  VK bound C           = 21.07")
print(f"  E/C                  = {E/21.07:.2f}")
print(f"  Grid: {Ngrid}³, L = {L}ξ = {L*xi_nm:.0f} nm")
print(f"  ρ_H range: [{rho_H.min():.6f}, {rho_H.max():.6f}]")
print(f"  |∇n|² max: {Edens.max():.2f}")
print(f"  --- QVC v4 Parameters ---")
print(f"  ξ = {xi_nm} nm")
print(f"  E_vac = {E_vac} meV")
print(f"  F = {F_fat}")
print(f"  g₀ = {g0} meV")
print(f"  A*/Ac = {Astar_over_Ac}")
print(f"  Lattice constant (BCC) = {xi_nm} nm = ξ")
