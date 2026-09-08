"""
Fix script for two issues in generate_structures.py:

1. Twisted bilayer Layer 2: deduplication fails after rotation, giving 7× too many atoms.
   Fix: use fractional-coordinate hashing with 4-decimal rounding.

2. Hopf charge Q = 0: Berry connection A_i = (nx ∂_i ny - ny ∂_i nx)/(1+nz) has
   a gauge singularity at n_z = -1 (the origin for Whitehead hopfion).
   Fix: compute B_i = (1/2)ε_{ijk} F_{jk} from gauge-invariant curvature,
   then invert ∇×A = B via FFT (Coulomb gauge), then Q = (1/4π²) ∫ A·B d³x.
"""
import numpy as np
import os

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR_CIF  = os.path.join(BASEDIR, "cif")
OUTDIR_CUBE = os.path.join(BASEDIR, "cube")
FIGDIR      = os.path.join(BASEDIR, "figures")
os.makedirs(FIGDIR, exist_ok=True)

a_gr  = 2.46    # Å, graphene lattice constant
a_CC  = 1.42    # Å, C-C bond
d_AB  = 3.35    # Å, interlayer spacing


# =========================================================================
# FIX 1: Twisted Bilayer — proper fractional deduplication
# =========================================================================
def write_cif(filename, data_name, a, b, c, alpha, beta, gamma,
              spacegroup, sg_number, atoms, comment=""):
    path = os.path.join(OUTDIR_CIF, filename)
    with open(path, 'w') as f:
        if comment:
            f.write(f"# {comment}\n")
        f.write(f"data_{data_name}\n\n")
        f.write(f"_cell_length_a    {a:.6f}\n")
        f.write(f"_cell_length_b    {b:.6f}\n")
        f.write(f"_cell_length_c    {c:.6f}\n")
        f.write(f"_cell_angle_alpha {alpha:.4f}\n")
        f.write(f"_cell_angle_beta  {beta:.4f}\n")
        f.write(f"_cell_angle_gamma {gamma:.4f}\n\n")
        f.write(f"_symmetry_space_group_name_H-M   '{spacegroup}'\n")
        f.write(f"_symmetry_Int_Tables_number       {sg_number}\n\n")
        f.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n_atom_site_U_iso_or_equiv\n")
        f.write("_atom_site_occupancy\n")
        for label, sym, fx, fy, fz, uiso in atoms:
            f.write(f"{label:6s} {sym:4s} {fx:10.6f} {fy:10.6f} {fz:10.6f} "
                    f"{uiso:.4f} 1.0000\n")
    print(f"  Written: {path}  ({len(atoms)} sites)")
    return path


def gen_twisted_bilayer_fixed(m=2, n=1):
    """Commensurate twisted bilayer with fractional-coordinate deduplication."""
    print(f"\n=== Twisted Bilayer (m={m}, n={n}) — FIXED ===")

    N_cell = m**2 + m*n + n**2
    cos_theta = (m**2 + 4*m*n + n**2) / (2 * N_cell)
    theta = np.arccos(np.clip(cos_theta, -1, 1))
    theta_deg = np.degrees(theta)
    print(f"  θ = {theta_deg:.2f}°, expect {2*N_cell} atoms/layer, {4*N_cell} total")

    # Graphene primitive vectors (Cartesian)
    a1 = a_gr * np.array([1.0, 0.0])
    a2 = a_gr * np.array([0.5, np.sqrt(3)/2])

    # Sublattice basis
    basis_A = np.array([0.0, 0.0])
    basis_B = (a1 + a2) / 3   # (1/3, 2/3) fractional

    # Superlattice vectors
    L1 = m * a1 + n * a2
    L2 = -n * a1 + (m + n) * a2
    inv_L = np.linalg.inv(np.column_stack([L1, L2]))

    L_mag = np.linalg.norm(L1)
    cos_gamma = np.dot(L1, L2) / (np.linalg.norm(L1) * np.linalg.norm(L2))
    gamma_deg = np.degrees(np.arccos(np.clip(cos_gamma, -1, 1)))
    print(f"  |L| = {L_mag:.4f} Å, γ = {gamma_deg:.2f}°")

    # Rotation matrix for layer 2
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])

    # AB-stacking shift (applied to layer 2 before rotation)
    ab_shift = (a1 + a2) / 3

    def collect_layer(use_rotation=False, shift=np.zeros(2)):
        """Generate unique atomic positions for one layer."""
        seen = set()
        positions = []
        # Search range: need enough to cover the supercell after rotation
        R_scan = max(m, n) + 3
        for i in range(-R_scan, R_scan + 1):
            for j in range(-R_scan, R_scan + 1):
                for base in [basis_A, basis_B]:
                    pos = i * a1 + j * a2 + base + shift
                    if use_rotation:
                        pos = R @ pos
                    # Convert to fractional coords of supercell
                    frac = inv_L @ pos
                    fx = frac[0] % 1.0
                    fy = frac[1] % 1.0
                    # Round for deduplication (4 decimal places = 0.0001 tolerance)
                    key = (round(fx, 4), round(fy, 4))
                    # Handle wraparound: 0.9999 → 0.0
                    if key[0] >= 0.9999:
                        key = (0.0, key[1])
                    if key[1] >= 0.9999:
                        key = (key[0], 0.0)
                    if key not in seen:
                        seen.add(key)
                        positions.append((fx, fy))
        return positions

    layer1 = collect_layer(use_rotation=False, shift=np.zeros(2))
    layer2 = collect_layer(use_rotation=True, shift=ab_shift)
    print(f"  Layer 1: {len(layer1)} atoms (expect {2*N_cell})")
    print(f"  Layer 2: {len(layer2)} atoms (expect {2*N_cell})")

    # Build CIF
    c_cell = 20.0
    z1 = 0.5 - d_AB / (2 * c_cell)
    z2 = 0.5 + d_AB / (2 * c_cell)

    atoms = []
    for idx, (fx, fy) in enumerate(layer1, 1):
        atoms.append((f"C{idx}", "C", fx, fy, z1, 0.005))
    for idx, (fx, fy) in enumerate(layer2, len(layer1) + 1):
        atoms.append((f"C{idx}", "C", fx, fy, z2, 0.005))

    return write_cif(
        f"03_twisted_bilayer_m{m}n{n}.cif",
        f"twisted_bilayer_theta_{theta_deg:.1f}",
        np.linalg.norm(L1), np.linalg.norm(L2), c_cell,
        90.0, 90.0, gamma_deg,
        "P 1", 1, atoms,
        f"Commensurate twisted bilayer graphene, theta={theta_deg:.2f} deg, "
        f"(m,n)=({m},{n}), {len(atoms)} atoms"
    ), theta_deg, len(atoms)


# =========================================================================
# FIX 2: Hopfion Charge — Coulomb-gauge FFT inversion
# =========================================================================
def gen_hopfion_cube_fixed(Ngrid=64):
    """
    Whitehead hopfion Q=1 with Fourier-space Hopf invariant computation.

    The gauge-invariant approach:
    1. Compute Berry curvature F_{ij} = ε_{abc} n_a ∂_i n_b ∂_j n_c
    2. Define B_i = (1/2) ε_{ijk} F_{jk}  ("magnetic field" from curvature)
    3. Solve ∇×A = B via FFT in Coulomb gauge: Â(k) = -i k×B̂(k)/|k|²
    4. Q = (1/4π²) ∫ A·B d³x  (Hopf invariant)

    This avoids the gauge singularity that kills the direct A formula.
    """
    print(f"\n=== Hopfion Charge Density — Coulomb Gauge FFT ({Ngrid}³) ===")

    L = 4.0   # box half-width in ξ units (larger than before for convergence)
    x = np.linspace(-L, L, Ngrid, endpoint=False)  # periodic grid for FFT
    dx = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = X**2 + Y**2 + Z**2

    # --- Hopf map: R³ → S³ → S² ---
    denom = 1.0 + r2
    w0 = 2*X / denom
    w1 = 2*Y / denom
    w2 = 2*Z / denom
    w3 = (1.0 - r2) / denom

    # n = π(w): S³ → S² via Hopf
    nx = 2*(w0*w2 + w1*w3)
    ny = 2*(w1*w2 - w0*w3)
    nz = w0**2 + w1**2 - w2**2 - w3**2

    nmag = np.sqrt(nx**2 + ny**2 + nz**2)
    print(f"  |n| check: [{nmag.min():.8f}, {nmag.max():.8f}]")

    # --- Gradients (central differences, periodic) ---
    def pgrad(f, axis, h):
        """Periodic central difference along axis."""
        return (np.roll(f, -1, axis=axis) - np.roll(f, 1, axis=axis)) / (2*h)

    dnx_dx, dnx_dy, dnx_dz = pgrad(nx,0,dx), pgrad(nx,1,dx), pgrad(nx,2,dx)
    dny_dx, dny_dy, dny_dz = pgrad(ny,0,dx), pgrad(ny,1,dx), pgrad(ny,2,dx)
    dnz_dx, dnz_dy, dnz_dz = pgrad(nz,0,dx), pgrad(nz,1,dx), pgrad(nz,2,dx)

    # --- Berry curvature 2-form: F_{ij} = ε_{abc} n_a ∂_i n_b ∂_j n_c ---
    def F_component(i_grads, j_grads):
        """F_{ij} = n_x(∂_i n_y ∂_j n_z - ∂_j n_y ∂_i n_z) + cyclic"""
        di_nx, di_ny, di_nz = i_grads
        dj_nx, dj_ny, dj_nz = j_grads
        return (nx * (di_ny * dj_nz - dj_ny * di_nz) +
                ny * (di_nz * dj_nx - dj_nz * di_nx) +
                nz * (di_nx * dj_ny - dj_nx * di_ny))

    grad_x = (dnx_dx, dny_dx, dnz_dx)
    grad_y = (dnx_dy, dny_dy, dnz_dy)
    grad_z = (dnx_dz, dny_dz, dnz_dz)

    Fyz = F_component(grad_y, grad_z)
    Fzx = F_component(grad_z, grad_x)
    Fxy = F_component(grad_x, grad_y)

    # "Magnetic field" B_i = F_{jk} with (i,j,k) cyclic
    Bx, By, Bz = Fyz, Fzx, Fxy
    print(f"  |B| max: {np.sqrt(Bx**2+By**2+Bz**2).max():.6f}")

    # Check ∇·B = 0 (gauge-invariant consistency)
    divB = pgrad(Bx,0,dx) + pgrad(By,1,dx) + pgrad(Bz,2,dx)
    print(f"  ∇·B max: {np.abs(divB).max():.2e} (should be ~0)")

    # --- FFT inversion: solve ∇×A = B in Coulomb gauge ---
    # Â(k) = -i k × B̂(k) / |k|²
    Bx_hat = np.fft.fftn(Bx)
    By_hat = np.fft.fftn(By)
    Bz_hat = np.fft.fftn(Bz)

    # k-space grid
    kx = np.fft.fftfreq(Ngrid, d=dx) * 2 * np.pi
    ky = np.fft.fftfreq(Ngrid, d=dx) * 2 * np.pi
    kz = np.fft.fftfreq(Ngrid, d=dx) * 2 * np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K2 = KX**2 + KY**2 + KZ**2
    K2[0, 0, 0] = 1.0  # avoid division by zero (set k=0 mode to 0 later)

    # Â = -i (k × B̂) / |k|²
    # (k × B̂)_x = k_y B̂_z - k_z B̂_y, etc.
    cross_x = KY * Bz_hat - KZ * By_hat
    cross_y = KZ * Bx_hat - KX * Bz_hat
    cross_z = KX * By_hat - KY * Bx_hat

    Ax_hat = -1j * cross_x / K2
    Ay_hat = -1j * cross_y / K2
    Az_hat = -1j * cross_z / K2

    # Zero out k=0 mode
    Ax_hat[0,0,0] = 0
    Ay_hat[0,0,0] = 0
    Az_hat[0,0,0] = 0

    Ax = np.real(np.fft.ifftn(Ax_hat))
    Ay = np.real(np.fft.ifftn(Ay_hat))
    Az = np.real(np.fft.ifftn(Az_hat))

    print(f"  |A| max: {np.sqrt(Ax**2+Ay**2+Az**2).max():.6f}")

    # --- Hopf invariant: Q = (1/4π²) ∫ A·B d³x ---
    AdotB = Ax*Bx + Ay*By + Az*Bz
    Q = (1.0 / (4 * np.pi**2)) * np.sum(AdotB) * dx**3
    print(f"  Hopf charge Q = {Q:.6f} (expect ±1)")

    # Hopf density for cube file
    rho_H = (1.0 / (4 * np.pi**2)) * AdotB

    # --- Energy density and VK bound ---
    E_density = (dnx_dx**2 + dnx_dy**2 + dnx_dz**2 +
                 dny_dx**2 + dny_dy**2 + dny_dz**2 +
                 dnz_dx**2 + dnz_dy**2 + dnz_dz**2)
    E_total = 0.5 * np.sum(E_density) * dx**3
    VK_bound = 16 * np.pi**2 * (3**(3/8)) / (2**(7/2))
    print(f"  Energy E = {E_total:.2f}")
    print(f"  VK bound = {VK_bound:.2f}")
    print(f"  E/E_VK   = {E_total/VK_bound:.4f} (≥1 required)")
    print(f"  ρ_H range: [{rho_H.min():.6f}, {rho_H.max():.6f}]")

    # --- Write cube files ---
    origin = np.array([-L, -L, -L])

    def write_cube(fname, data, title1, title2):
        path = os.path.join(OUTDIR_CUBE, fname)
        with open(path, 'w') as f:
            f.write(f"{title1}\n{title2}\n")
            f.write(f"    1 {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
            f.write(f" {Ngrid:4d} {dx:12.6f}  0.000000  0.000000\n")
            f.write(f" {Ngrid:4d}  0.000000 {dx:12.6f}  0.000000\n")
            f.write(f" {Ngrid:4d}  0.000000  0.000000 {dx:12.6f}\n")
            f.write(f"   67  0.000000  0.000000  0.000000  0.000000\n")
            for ix in range(Ngrid):
                for iy in range(Ngrid):
                    vals = data[ix, iy, :]
                    for iz in range(Ngrid):
                        f.write(f" {vals[iz]:13.5E}")
                        if (iz + 1) % 6 == 0:
                            f.write("\n")
                    if Ngrid % 6 != 0:
                        f.write("\n")
        sz = os.path.getsize(path) / 1024
        print(f"  Written: {path} ({sz:.0f} kB)")

    write_cube("hopfion_charge_density.cube", rho_H,
               "Hopfion topological charge density rho_H(r)",
               f"Q={Q:.4f}, Ngrid={Ngrid}, L={L}")
    write_cube("hopfion_energy_density.cube", E_density,
               "Hopfion energy density (1/2)|grad n|^2",
               f"E={E_total:.2f}, VK={VK_bound:.2f}")

    return Q, E_total, VK_bound, rho_H, E_density, nx, ny, nz, x


# =========================================================================
# MATPLOTLIB 3D VISUALIZATIONS
# =========================================================================
def plot_hopfion_3d(nx, ny, nz, x, rho_H, E_density):
    """Publication-quality 3D visualization of hopfion field and topology."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import cm

    Ngrid = len(x)
    dx = x[1] - x[0]

    # --- Figure 1: Hopfion preimage curves (linked circles) ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Find preimage of north pole (nz ≈ 1): should be a circle at |r|=1 in xy-plane
    # Find preimage of a point on the equator: should be a circle linked to the first
    X3, Y3, Z3 = np.meshgrid(x, x, x, indexing='ij')

    # Preimage of north pole: nz > 0.95
    mask_N = nz > 0.95
    if np.any(mask_N):
        ax.scatter(X3[mask_N], Y3[mask_N], Z3[mask_N],
                   c='blue', s=1, alpha=0.3, label=r'$n_z > 0.95$ (north pole)')

    # Preimage of equator point (1,0,0): nx > 0.95
    mask_E = nx > 0.95
    if np.any(mask_E):
        ax.scatter(X3[mask_E], Y3[mask_E], Z3[mask_E],
                   c='red', s=1, alpha=0.3, label=r'$n_x > 0.95$ (equator)')

    # Preimage of another equator point (0,1,0): ny > 0.95
    mask_E2 = ny > 0.95
    if np.any(mask_E2):
        ax.scatter(X3[mask_E2], Y3[mask_E2], Z3[mask_E2],
                   c='green', s=1, alpha=0.3, label=r'$n_y > 0.95$ (equator)')

    ax.set_xlabel('x / ξ', fontsize=12)
    ax.set_ylabel('y / ξ', fontsize=12)
    ax.set_zlabel('z / ξ', fontsize=12)
    ax.set_title('Hopfion Preimage Curves (Linked Circles)\nWhitehead Q = 1',
                 fontsize=14)
    ax.legend(fontsize=10, loc='upper left')
    ax.set_xlim(-3, 3)
    ax.set_ylim(-3, 3)
    ax.set_zlim(-3, 3)

    path1 = os.path.join(FIGDIR, "hopfion_preimage_curves.png")
    plt.savefig(path1, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Figure: {path1}")

    # --- Figure 2: Hopf charge density isosurface (midplane slices) ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    mid = Ngrid // 2
    extent = [x[0], x[-1], x[0], x[-1]]

    # xy-plane (z=0)
    im0 = axes[0].imshow(rho_H[:, :, mid].T, origin='lower', extent=extent,
                          cmap='RdBu_r', aspect='equal')
    axes[0].set_title(f'ρ_H(x,y,z=0)', fontsize=12)
    axes[0].set_xlabel('x / ξ')
    axes[0].set_ylabel('y / ξ')
    plt.colorbar(im0, ax=axes[0], shrink=0.8)

    # xz-plane (y=0)
    im1 = axes[1].imshow(rho_H[:, mid, :].T, origin='lower', extent=extent,
                          cmap='RdBu_r', aspect='equal')
    axes[1].set_title(f'ρ_H(x,y=0,z)', fontsize=12)
    axes[1].set_xlabel('x / ξ')
    axes[1].set_ylabel('z / ξ')
    plt.colorbar(im1, ax=axes[1], shrink=0.8)

    # yz-plane (x=0)
    im2 = axes[2].imshow(rho_H[mid, :, :].T, origin='lower', extent=extent,
                          cmap='RdBu_r', aspect='equal')
    axes[2].set_title(f'ρ_H(x=0,y,z)', fontsize=12)
    axes[2].set_xlabel('y / ξ')
    axes[2].set_ylabel('z / ξ')
    plt.colorbar(im2, ax=axes[2], shrink=0.8)

    fig.suptitle('Hopfion Topological Charge Density Slices', fontsize=14, y=1.02)
    plt.tight_layout()
    path2 = os.path.join(FIGDIR, "hopfion_charge_slices.png")
    plt.savefig(path2, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Figure: {path2}")

    # --- Figure 3: Energy density + n-field arrows (midplane) ---
    fig, ax = plt.subplots(figsize=(8, 7))

    im = ax.imshow(E_density[:, :, mid].T, origin='lower', extent=extent,
                    cmap='inferno', aspect='equal')
    plt.colorbar(im, ax=ax, label=r'$|\nabla \mathbf{n}|^2$')

    # Quiver plot of n-field in xy-plane
    skip = max(1, Ngrid // 20)
    xs = x[::skip]
    Xs, Ys = np.meshgrid(xs, xs, indexing='ij')
    U = nx[::skip, ::skip, mid]
    V = ny[::skip, ::skip, mid]
    C = nz[::skip, ::skip, mid]

    q = ax.quiver(Xs, Ys, U, V, C, cmap='coolwarm', clim=(-1, 1),
                   scale=25, width=0.003, alpha=0.8)

    ax.set_xlabel('x / ξ', fontsize=12)
    ax.set_ylabel('y / ξ', fontsize=12)
    ax.set_title('Energy Density + n-field (z=0 plane)', fontsize=14)

    path3 = os.path.join(FIGDIR, "hopfion_energy_nfield.png")
    plt.savefig(path3, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Figure: {path3}")

    # --- Figure 4: n_z field on three planes ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))

    for ax_i, (data, title) in enumerate([
        (nz[:, :, mid], r'$n_z(x,y,z=0)$'),
        (nz[:, mid, :], r'$n_z(x,y=0,z)$'),
        (nz[mid, :, :], r'$n_z(x=0,y,z)$'),
    ]):
        im = axes[ax_i].imshow(data.T, origin='lower', extent=extent,
                                cmap='coolwarm', vmin=-1, vmax=1, aspect='equal')
        axes[ax_i].set_title(title, fontsize=12)
        plt.colorbar(im, ax=axes[ax_i], shrink=0.8)

    axes[0].set_xlabel('x / ξ'); axes[0].set_ylabel('y / ξ')
    axes[1].set_xlabel('x / ξ'); axes[1].set_ylabel('z / ξ')
    axes[2].set_xlabel('y / ξ'); axes[2].set_ylabel('z / ξ')
    fig.suptitle('Hopfion n_z Component (S² target space)', fontsize=14, y=1.02)
    plt.tight_layout()
    path4 = os.path.join(FIGDIR, "hopfion_nz_slices.png")
    plt.savefig(path4, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"  Figure: {path4}")

    return [path1, path2, path3, path4]


# =========================================================================
# MAIN
# =========================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TEVC Structure Fix — Twisted Bilayer + Hopfion Charge")
    print("=" * 70)

    # Fix twisted bilayers
    gen_twisted_bilayer_fixed(m=2, n=1)
    gen_twisted_bilayer_fixed(m=3, n=2)

    # Fix hopfion computation
    Q, E, VK, rho_H, Edens, nx, ny, nz, grid = gen_hopfion_cube_fixed(Ngrid=64)

    # Generate matplotlib figures
    print("\n=== Generating 3D Visualizations ===")
    figs = plot_hopfion_3d(nx, ny, nz, grid, rho_H, Edens)

    # Summary
    print("\n" + "=" * 70)
    print("FIX SUMMARY")
    print("=" * 70)
    print(f"  Hopf charge Q = {Q:.6f} (analytic: 1)")
    print(f"  Energy E      = {E:.2f}")
    print(f"  VK bound      = {VK:.2f}")
    print(f"  E/E_VK        = {E/VK:.4f}")
    print(f"  Figures: {len(figs)} generated in {FIGDIR}")
