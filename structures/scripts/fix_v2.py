"""
TEVC Fix v2: Corrected twisted bilayer + windowed Hopf charge
=============================================================
Bug 1: Superlattice vectors were L1 = m*a1+n*a2 (wrong).
        Correct: L1 = n*a1+m*a2. Verified by commensurate identity
        R(theta)*(m*a1+n*a2) = n*a1+m*a2 = L1.

Bug 2: Hopf field nz ~ -0.56 at box boundary (r=4) breaks FFT periodicity.
        Fix: smooth window forces n -> (0,0,-1) for r > R_cut, then
        divergence-project B in k-space before Coulomb inversion.
"""
import numpy as np
import os, sys

BASEDIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUTDIR_CIF  = os.path.join(BASEDIR, "cif")
OUTDIR_CUBE = os.path.join(BASEDIR, "cube")
FIGDIR      = os.path.join(BASEDIR, "figures")
for d in [OUTDIR_CIF, OUTDIR_CUBE, FIGDIR]:
    os.makedirs(d, exist_ok=True)

a_gr, a_CC, d_AB = 2.46, 1.42, 3.35

def write_cif(filename, data_name, a, b, c, alpha, beta, gamma,
              sg, sg_num, atoms, comment=""):
    path = os.path.join(OUTDIR_CIF, filename)
    with open(path, 'w') as f:
        if comment: f.write(f"# {comment}\n")
        f.write(f"data_{data_name}\n\n")
        f.write(f"_cell_length_a    {a:.6f}\n_cell_length_b    {b:.6f}\n")
        f.write(f"_cell_length_c    {c:.6f}\n")
        f.write(f"_cell_angle_alpha {alpha:.4f}\n_cell_angle_beta  {beta:.4f}\n")
        f.write(f"_cell_angle_gamma {gamma:.4f}\n\n")
        f.write(f"_symmetry_space_group_name_H-M   '{sg}'\n")
        f.write(f"_symmetry_Int_Tables_number       {sg_num}\n\n")
        f.write("loop_\n_atom_site_label\n_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n")
        f.write("_atom_site_U_iso_or_equiv\n_atom_site_occupancy\n")
        for lab, sym, fx, fy, fz, u in atoms:
            f.write(f"{lab:6s} {sym:4s} {fx:10.6f} {fy:10.6f} {fz:10.6f} {u:.4f} 1.0000\n")
    print(f"  Written: {path}  ({len(atoms)} sites)")
    return path


# ======================================================================
# TWISTED BILAYER — CORRECTED SUPERLATTICE
# ======================================================================
def gen_twisted_bilayer(m, n):
    """
    Commensurate twisted bilayer graphene.

    Superlattice (CORRECTED):
        L1 = n*a1 + m*a2    (layer 1, unrotated)
        L2 = -m*a1 + (n+m)*a2

    Commensurate identity: R(theta) * (m*a1 + n*a2) = n*a1 + m*a2 = L1
    So layer 2 periodicity folds correctly into L1, L2.

    cos(theta) = (m^2 + 4mn + n^2) / [2(m^2 + mn + n^2)]
    N_cell = m^2 + mn + n^2  (unit cells per layer in supercell)
    Total atoms = 4 * N_cell  (2 sublattices × 2 layers)
    """
    print(f"\n=== Twisted Bilayer (m={m}, n={n}) ===")
    N_cell = m**2 + m*n + n**2
    cos_t = (m**2 + 4*m*n + n**2) / (2*N_cell)
    theta = np.arccos(np.clip(cos_t, -1, 1))
    theta_deg = np.degrees(theta)

    a1 = a_gr * np.array([1.0, 0.0])
    a2 = a_gr * np.array([0.5, np.sqrt(3)/2])
    basis_A = np.array([0.0, 0.0])
    basis_B = (a1 + a2) / 3

    # CORRECTED superlattice vectors (n,m swapped from wrong version)
    L1 = n * a1 + m * a2       # was m*a1 + n*a2
    L2 = -m * a1 + (n+m) * a2  # was -n*a1 + (m+n)*a2
    inv_L = np.linalg.inv(np.column_stack([L1, L2]))

    # Verify commensurate identity
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])
    check = R @ (m*a1 + n*a2)
    err = np.linalg.norm(check - L1)
    print(f"  θ = {theta_deg:.2f}°, N_cell = {N_cell}, expect {4*N_cell} atoms")
    print(f"  Commensurate check: |R(m·a1+n·a2) - L1| = {err:.2e}")

    L_mag = np.linalg.norm(L1)
    cos_g = np.dot(L1, L2) / (np.linalg.norm(L1)*np.linalg.norm(L2))
    gamma_deg = np.degrees(np.arccos(np.clip(cos_g, -1, 1)))
    print(f"  |L| = {L_mag:.4f} Å, γ = {gamma_deg:.1f}°")

    ab_shift = (a1 + a2) / 3  # AB-stacking offset

    def collect_atoms(rotate=False, shift=np.zeros(2)):
        seen = set()
        out = []
        R_scan = max(m, n) + 3
        for i in range(-R_scan, R_scan+1):
            for j in range(-R_scan, R_scan+1):
                for base in [basis_A, basis_B]:
                    pos = i*a1 + j*a2 + base + shift
                    if rotate:
                        pos = R @ pos
                    frac = inv_L @ pos
                    fx, fy = frac[0] % 1.0, frac[1] % 1.0
                    key = (round(fx, 4), round(fy, 4))
                    if key[0] >= 0.9999: key = (0.0, key[1])
                    if key[1] >= 0.9999: key = (key[0], 0.0)
                    if key not in seen:
                        seen.add(key)
                        out.append((fx, fy))
        return out

    L1_atoms = collect_atoms(rotate=False)
    L2_atoms = collect_atoms(rotate=True, shift=ab_shift)
    print(f"  Layer 1: {len(L1_atoms)} (expect {2*N_cell})")
    print(f"  Layer 2: {len(L2_atoms)} (expect {2*N_cell})")

    c_cell = 20.0
    z1 = 0.5 - d_AB/(2*c_cell)
    z2 = 0.5 + d_AB/(2*c_cell)
    atoms = []
    for i, (fx,fy) in enumerate(L1_atoms, 1):
        atoms.append((f"C{i}", "C", fx, fy, z1, 0.005))
    for i, (fx,fy) in enumerate(L2_atoms, len(L1_atoms)+1):
        atoms.append((f"C{i}", "C", fx, fy, z2, 0.005))

    return write_cif(
        f"03_twisted_bilayer_m{m}n{n}.cif",
        f"twisted_bilayer_{theta_deg:.1f}deg",
        np.linalg.norm(L1), np.linalg.norm(L2), c_cell,
        90.0, 90.0, gamma_deg, "P 1", 1, atoms,
        f"Commensurate twisted bilayer, theta={theta_deg:.2f}, (m,n)=({m},{n}), {len(atoms)} atoms"
    )


# ======================================================================
# HOPFION — WINDOWED + DIVERGENCE-PROJECTED FFT
# ======================================================================
def gen_hopfion(Ngrid=80):
    """
    Whitehead hopfion Q=1 with smooth boundary window.

    Window: for |r| > R_cut, smoothly rotate n(r) toward (0,0,-1)
    using normalized linear interpolation (preserves |n|=1).
    Then B = curl-of-A is periodic and the FFT inversion converges.

    After windowing, project B divergence-free in k-space:
        B_clean = B - k(k·B)/|k|^2   [removes ∇·B artifacts]
    """
    print(f"\n=== Hopfion Q=1 — Windowed FFT ({Ngrid}³) ===")

    L = 4.0
    x = np.linspace(-L, L, Ngrid, endpoint=False)
    dx = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = X**2 + Y**2 + Z**2
    r = np.sqrt(r2)

    # --- Hopf map ---
    denom = 1.0 + r2
    w0, w1, w2, w3 = 2*X/denom, 2*Y/denom, 2*Z/denom, (1-r2)/denom
    nx_raw = 2*(w0*w2 + w1*w3)
    ny_raw = 2*(w1*w2 - w0*w3)
    nz_raw = w0**2 + w1**2 - w2**2 - w3**2

    # --- Smooth window: blend n toward (0,0,-1) for r > R_cut ---
    R_cut = 2.5  # keep full hopfion structure inside this radius
    delta = 0.4  # transition width
    # t(r): 0 for r < R_cut, 1 for r >> R_cut
    t = 0.5 * (1.0 + np.tanh((r - R_cut) / delta))

    # Linear blend + renormalize to stay on S²
    nx_blend = (1-t) * nx_raw + t * 0.0
    ny_blend = (1-t) * ny_raw + t * 0.0
    nz_blend = (1-t) * nz_raw + t * (-1.0)
    norm_blend = np.sqrt(nx_blend**2 + ny_blend**2 + nz_blend**2)
    norm_blend = np.maximum(norm_blend, 1e-15)
    nx = nx_blend / norm_blend
    ny = ny_blend / norm_blend
    nz = nz_blend / norm_blend

    print(f"  |n| range: [{np.sqrt(nx**2+ny**2+nz**2).min():.8f}, "
          f"{np.sqrt(nx**2+ny**2+nz**2).max():.8f}]")
    print(f"  n_z at boundary: {nz[0,Ngrid//2,Ngrid//2]:.4f} "
          f"(should be ≈ -1.0)")

    # --- Periodic gradients ---
    def pg(f, ax):
        return (np.roll(f,-1,axis=ax) - np.roll(f,1,axis=ax)) / (2*dx)

    dnx = [pg(nx,0), pg(nx,1), pg(nx,2)]
    dny = [pg(ny,0), pg(ny,1), pg(ny,2)]
    dnz = [pg(nz,0), pg(nz,1), pg(nz,2)]

    # --- Berry curvature F_{ij} = ε_{abc} n_a ∂_i n_b ∂_j n_c ---
    def Fcomp(di, dj):
        dix, diy, diz = [dnx[di], dny[di], dnz[di]]  # wrong indexing
        djx, djy, djz = [dnx[dj], dny[dj], dnz[dj]]
        return (nx*(diy*djz - djy*diz) +
                ny*(diz*djx - djz*dix) +
                nz*(dix*djy - djx*diy))

    # Indices: 0=x, 1=y, 2=z
    # F_{ij} where i,j are coordinate indices
    # Need gradients indexed by coordinate
    # ∂_i n_a means: for coordinate i, gradient of component a
    # dnx[i] = ∂_i(n_x), dny[i] = ∂_i(n_y), dnz[i] = ∂_i(n_z)

    def F_ij(i, j):
        """F_{ij} = ε_{abc} n_a (∂_i n_b)(∂_j n_c)"""
        return (nx * (dny[i]*dnz[j] - dny[j]*dnz[i]) +
                ny * (dnz[i]*dnx[j] - dnz[j]*dnx[i]) +
                nz * (dnx[i]*dny[j] - dnx[j]*dny[i]))

    Bx = F_ij(1, 2)  # F_{yz}
    By = F_ij(2, 0)  # F_{zx}
    Bz = F_ij(0, 1)  # F_{xy}

    divB_raw = pg(Bx,0) + pg(By,1) + pg(Bz,2)
    print(f"  ∇·B (raw) max: {np.abs(divB_raw).max():.3e}")

    # --- FFT: divergence-project B, then Coulomb inversion ---
    Bxh = np.fft.fftn(Bx)
    Byh = np.fft.fftn(By)
    Bzh = np.fft.fftn(Bz)

    kx = np.fft.fftfreq(Ngrid, d=dx) * 2*np.pi
    ky = np.fft.fftfreq(Ngrid, d=dx) * 2*np.pi
    kz = np.fft.fftfreq(Ngrid, d=dx) * 2*np.pi
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing='ij')
    K2 = KX**2 + KY**2 + KZ**2
    K2_safe = np.where(K2 > 0, K2, 1.0)

    # Project: B_clean = B - k(k·B)/|k|²
    kdotB = KX*Bxh + KY*Byh + KZ*Bzh
    Bxh_c = Bxh - KX * kdotB / K2_safe
    Byh_c = Byh - KY * kdotB / K2_safe
    Bzh_c = Bzh - KZ * kdotB / K2_safe
    Bxh_c[0,0,0] = 0; Byh_c[0,0,0] = 0; Bzh_c[0,0,0] = 0

    # Coulomb gauge: A_hat = -i (k × B_hat) / |k|²
    cx = KY*Bzh_c - KZ*Byh_c
    cy = KZ*Bxh_c - KX*Bzh_c
    cz = KX*Byh_c - KY*Bxh_c

    Axh = -1j * cx / K2_safe
    Ayh = -1j * cy / K2_safe
    Azh = -1j * cz / K2_safe
    Axh[0,0,0] = 0; Ayh[0,0,0] = 0; Azh[0,0,0] = 0

    Ax = np.real(np.fft.ifftn(Axh))
    Ay = np.real(np.fft.ifftn(Ayh))
    Az = np.real(np.fft.ifftn(Azh))

    # Verify cleaned ∇·B
    Bx_c = np.real(np.fft.ifftn(Bxh_c))
    By_c = np.real(np.fft.ifftn(Byh_c))
    Bz_c = np.real(np.fft.ifftn(Bzh_c))
    divB_c = pg(Bx_c,0) + pg(By_c,1) + pg(Bz_c,2)
    print(f"  ∇·B (projected) max: {np.abs(divB_c).max():.3e}")

    # --- Hopf invariant ---
    AdotB = Ax*Bx_c + Ay*By_c + Az*Bz_c
    Q = np.sum(AdotB) * dx**3 / (4*np.pi**2)
    print(f"  *** Hopf charge Q = {Q:.6f} (expect ±1) ***")

    rho_H = AdotB / (4*np.pi**2)

    # --- Energy ---
    Edens = sum(dnx[i]**2 + dny[i]**2 + dnz[i]**2 for i in range(3))
    E = 0.5 * np.sum(Edens) * dx**3
    VK = 16*np.pi**2 * 3**(3/8) / 2**(7/2)
    print(f"  Energy E = {E:.2f}, VK bound = {VK:.2f}, E/VK = {E/VK:.3f}")

    # --- Write cube files ---
    origin = np.array([-L, -L, -L])
    def write_cube(fname, data, t1, t2):
        path = os.path.join(OUTDIR_CUBE, fname)
        with open(path, 'w') as f:
            f.write(f"{t1}\n{t2}\n")
            f.write(f"    1 {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
            for ax_i in range(3):
                step = [0.0, 0.0, 0.0]; step[ax_i] = dx
                f.write(f" {Ngrid:4d} {step[0]:12.6f} {step[1]:12.6f} {step[2]:12.6f}\n")
            f.write(f"   67  0.000000  0.000000  0.000000  0.000000\n")
            for ix in range(Ngrid):
                for iy in range(Ngrid):
                    for iz in range(Ngrid):
                        f.write(f" {data[ix,iy,iz]:13.5E}")
                        if (iz+1) % 6 == 0: f.write("\n")
                    if Ngrid % 6 != 0: f.write("\n")
        print(f"  Cube: {path} ({os.path.getsize(path)/1024:.0f} kB)")

    write_cube("hopfion_charge_density.cube", rho_H,
               "Hopfion topological charge density rho_H",
               f"Q={Q:.4f} Ngrid={Ngrid} L={L} Rcut={R_cut}")
    write_cube("hopfion_energy_density.cube", Edens,
               "Hopfion energy density |grad n|^2",
               f"E={E:.2f} VK={VK:.2f}")

    return Q, E, VK, rho_H, Edens, nx, ny, nz, x


# ======================================================================
# MATPLOTLIB FIGURES
# ======================================================================
def make_figures(nx, ny, nz, x, rho_H, Edens):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    Ng = len(x); mid = Ng//2
    ext = [x[0], x[-1], x[0], x[-1]]

    # --- Fig 1: Preimage curves ---
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    X3, Y3, Z3 = np.meshgrid(x, x, x, indexing='ij')
    for thresh, comp, color, lab in [
        (0.95, nz, 'blue', r'$n_z > 0.95$'),
        (0.95, nx, 'red',  r'$n_x > 0.95$'),
        (0.95, ny, 'green', r'$n_y > 0.95$')]:
        mask = comp > thresh
        if np.any(mask):
            pts = np.column_stack([X3[mask], Y3[mask], Z3[mask]])
            # Subsample if too many
            if len(pts) > 3000:
                idx = np.random.choice(len(pts), 3000, replace=False)
                pts = pts[idx]
            ax.scatter(pts[:,0], pts[:,1], pts[:,2], c=color, s=1, alpha=0.4, label=lab)
    ax.set_xlabel('x / ξ'); ax.set_ylabel('y / ξ'); ax.set_zlabel('z / ξ')
    ax.set_title('Hopfion Preimage Curves (Q=1)', fontsize=14)
    ax.legend(fontsize=10); ax.set_xlim(-3,3); ax.set_ylim(-3,3); ax.set_zlim(-3,3)
    p1 = os.path.join(FIGDIR, "hopfion_preimage_curves.png")
    plt.savefig(p1, dpi=200, bbox_inches='tight'); plt.close()
    print(f"  Fig: {p1}")

    # --- Fig 2: Charge density slices ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    vmax = max(abs(rho_H.min()), abs(rho_H.max()))
    for ax_i, (sl, tit, xl, yl) in enumerate([
        (rho_H[:,:,mid], r'$\rho_H(x,y,z{=}0)$', 'x / ξ', 'y / ξ'),
        (rho_H[:,mid,:], r'$\rho_H(x,y{=}0,z)$', 'x / ξ', 'z / ξ'),
        (rho_H[mid,:,:], r'$\rho_H(x{=}0,y,z)$', 'y / ξ', 'z / ξ')]):
        im = axes[ax_i].imshow(sl.T, origin='lower', extent=ext,
                                cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='equal')
        axes[ax_i].set_title(tit, fontsize=12)
        axes[ax_i].set_xlabel(xl); axes[ax_i].set_ylabel(yl)
        plt.colorbar(im, ax=axes[ax_i], shrink=0.8)
    fig.suptitle('Hopfion Topological Charge Density', fontsize=14, y=1.02)
    plt.tight_layout()
    p2 = os.path.join(FIGDIR, "hopfion_charge_slices.png")
    plt.savefig(p2, dpi=200, bbox_inches='tight'); plt.close()
    print(f"  Fig: {p2}")

    # --- Fig 3: Energy + quiver ---
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(Edens[:,:,mid].T, origin='lower', extent=ext, cmap='inferno', aspect='equal')
    plt.colorbar(im, ax=ax, label=r'$|\nabla \mathbf{n}|^2$')
    skip = max(1, Ng//20)
    xs = x[::skip]
    Xs, Ys = np.meshgrid(xs, xs, indexing='ij')
    ax.quiver(Xs, Ys, nx[::skip,::skip,mid], ny[::skip,::skip,mid],
              nz[::skip,::skip,mid], cmap='coolwarm', clim=(-1,1),
              scale=25, width=0.003, alpha=0.8)
    ax.set_xlabel('x / ξ'); ax.set_ylabel('y / ξ')
    ax.set_title('Energy Density + n-field (z = 0)', fontsize=14)
    p3 = os.path.join(FIGDIR, "hopfion_energy_nfield.png")
    plt.savefig(p3, dpi=200, bbox_inches='tight'); plt.close()
    print(f"  Fig: {p3}")

    # --- Fig 4: nz slices ---
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    for ax_i, (sl, tit) in enumerate([
        (nz[:,:,mid], r'$n_z(z{=}0)$'),
        (nz[:,mid,:], r'$n_z(y{=}0)$'),
        (nz[mid,:,:], r'$n_z(x{=}0)$')]):
        im = axes[ax_i].imshow(sl.T, origin='lower', extent=ext,
                                cmap='coolwarm', vmin=-1, vmax=1, aspect='equal')
        axes[ax_i].set_title(tit, fontsize=12)
        plt.colorbar(im, ax=axes[ax_i], shrink=0.8)
    fig.suptitle('Hopfion Target Space Component $n_z$', fontsize=14, y=1.02)
    plt.tight_layout()
    p4 = os.path.join(FIGDIR, "hopfion_nz_slices.png")
    plt.savefig(p4, dpi=200, bbox_inches='tight'); plt.close()
    print(f"  Fig: {p4}")

    return [p1, p2, p3, p4]


# ======================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TEVC Fix v2 — Corrected Structures")
    print("=" * 70)

    gen_twisted_bilayer(m=2, n=1)
    gen_twisted_bilayer(m=3, n=2)

    Q, E, VK, rho_H, Edens, nx, ny, nz, grid = gen_hopfion(Ngrid=80)

    print("\n=== Generating Figures ===")
    make_figures(nx, ny, nz, grid, rho_H, Edens)

    print(f"\n{'='*70}")
    print(f"FINAL: Q = {Q:.6f}, E = {E:.2f}, E/VK = {E/VK:.3f}")
    print(f"{'='*70}")
