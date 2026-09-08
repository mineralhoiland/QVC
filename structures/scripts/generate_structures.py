"""
TEVC Crystal Structure Generator for VESTA
===========================================
Generates VESTA-compatible CIF files for TDVT/TEVC platform materials:
  1. Graphene monolayer (P6/mmm)
  2. AB-stacked bilayer graphene
  3. Commensurate twisted bilayer graphene (θ = 21.79°, 28 atoms)
  4. BaTiO₃ tetragonal perovskite (piezoelectric substrate, P4mm)
  5. Hopfion spacetime crystal (BCC arrangement, dummy sites at Q=1 centers)

Also generates:
  - Gaussian cube file of hopfion topological charge density
  - Numerical readout summary

Locked QVC v4 parameters used throughout:
  ξ = 50 nm, Δ* = 0.2324 meV, E_vac = 0.1287 meV,
  g₀ = 0.0982 meV, ℏc_s = 0.07 meV·μm, T* = 118.7 mK

Author: QVC Research
Date: 2026-08-28
"""

import numpy as np
import os

# === Output directory ===
OUTDIR_CIF  = os.path.join(os.path.dirname(os.path.dirname(__file__)), "cif")
OUTDIR_CUBE = os.path.join(os.path.dirname(os.path.dirname(__file__)), "cube")
os.makedirs(OUTDIR_CIF, exist_ok=True)
os.makedirs(OUTDIR_CUBE, exist_ok=True)

# === Physical Constants (SI) ===
hbar    = 1.0546e-34     # J·s
kB      = 1.3806e-23     # J/K
eV      = 1.602e-19      # J/eV
meV     = 1.602e-22      # J/meV

# === Locked QVC v4 Parameters ===
xi_nm       = 50.0          # nm, coherence length
Delta_star  = 0.2324        # meV, locked gap
E_vac       = 0.1287        # meV, vacuum energy (fat torus F=0.58068)
g0          = 0.0982        # meV, coupling
hbar_cs     = 0.07          # meV·μm
T_star      = 118.7e-3      # K (118.7 mK)
F_torus     = 0.58068       # fat-torus Casimir prefactor
lambda_star = 0.7628        # dimensionless coupling
A_over_Ac   = 0.9117        # A*/Ac ratio
eta_ret     = 3.1           # retardation parameter

# === Graphene parameters ===
a_CC  = 1.42    # Å, C-C bond length
a_gr  = 2.46    # Å, lattice constant (= a_CC * sqrt(3))
d_AB  = 3.35    # Å, interlayer spacing (AB stacking)

# =========================================================================
# CIF GENERATION FUNCTIONS
# =========================================================================

def write_cif(filename, data_name, a, b, c, alpha, beta, gamma,
              spacegroup, sg_number, atoms, comment=""):
    """Write a CIF file.

    atoms: list of (label, symbol, frac_x, frac_y, frac_z, U_iso)
    """
    path = os.path.join(OUTDIR_CIF, filename)
    with open(path, 'w') as f:
        f.write(f"# {comment}\n" if comment else "")
        f.write(f"data_{data_name}\n\n")
        f.write(f"_cell_length_a    {a:.6f}\n")
        f.write(f"_cell_length_b    {b:.6f}\n")
        f.write(f"_cell_length_c    {c:.6f}\n")
        f.write(f"_cell_angle_alpha {alpha:.4f}\n")
        f.write(f"_cell_angle_beta  {beta:.4f}\n")
        f.write(f"_cell_angle_gamma {gamma:.4f}\n\n")
        f.write(f"_symmetry_space_group_name_H-M   '{spacegroup}'\n")
        f.write(f"_symmetry_Int_Tables_number       {sg_number}\n\n")
        f.write("loop_\n")
        f.write("_atom_site_label\n")
        f.write("_atom_site_type_symbol\n")
        f.write("_atom_site_fract_x\n")
        f.write("_atom_site_fract_y\n")
        f.write("_atom_site_fract_z\n")
        f.write("_atom_site_U_iso_or_equiv\n")
        f.write("_atom_site_occupancy\n")
        for label, sym, fx, fy, fz, uiso in atoms:
            f.write(f"{label:6s} {sym:4s} {fx:10.6f} {fy:10.6f} {fz:10.6f} {uiso:.4f} 1.0000\n")
        f.write("\n")
    print(f"  Written: {path}  ({len(atoms)} sites)")
    return path


# =========================================================================
# 1. GRAPHENE MONOLAYER
# =========================================================================
def gen_graphene_monolayer():
    """Graphene monolayer with vacuum slab."""
    print("\n=== 1. Graphene Monolayer ===")
    c_vac = 20.0  # Å vacuum
    # Hexagonal cell: C at (0,0,0.5) and (1/3, 2/3, 0.5)
    atoms = [
        ("C1", "C", 0.000000, 0.000000, 0.500000, 0.005),
        ("C2", "C", 1/3,      2/3,      0.500000, 0.005),
    ]
    return write_cif("01_graphene_monolayer.cif", "graphene_monolayer",
                     a_gr, a_gr, c_vac, 90.0, 90.0, 120.0,
                     "P 6/m m m", 191, atoms,
                     "Graphene monolayer, a=2.46 A, C-C=1.42 A")


# =========================================================================
# 2. AB-STACKED BILAYER GRAPHENE
# =========================================================================
def gen_bilayer_AB():
    """AB (Bernal) stacked bilayer graphene."""
    print("\n=== 2. AB-Stacked Bilayer Graphene ===")
    c_cell = 20.0  # Å total cell height
    z1 = 0.5 - d_AB / (2 * c_cell)  # lower layer
    z2 = 0.5 + d_AB / (2 * c_cell)  # upper layer

    # Layer 1 (A sublattice at origin, B at (1/3, 2/3))
    # Layer 2 AB-shifted: A' at (0, 0) [over B1], B' at (2/3, 1/3) [over hollow]
    atoms = [
        ("C1_A", "C", 0.000000, 0.000000, z1, 0.005),
        ("C1_B", "C", 1/3,      2/3,      z1, 0.005),
        ("C2_A", "C", 0.000000, 0.000000, z2, 0.005),  # A' sits over B1
        ("C2_B", "C", 2/3,      1/3,      z2, 0.005),  # B' over hollow
    ]
    return write_cif("02_bilayer_AB.cif", "bilayer_AB",
                     a_gr, a_gr, c_cell, 90.0, 90.0, 120.0,
                     "P 1", 1, atoms,
                     "AB (Bernal) bilayer graphene, d=3.35 A")


# =========================================================================
# 3. COMMENSURATE TWISTED BILAYER GRAPHENE
# =========================================================================
def gen_twisted_bilayer(m=2, n=1):
    """
    Commensurate twisted bilayer graphene using (m,n) construction.

    (m,n)=(2,1) → θ ≈ 21.79°, 28 atoms.  Manageable for VESTA.
    (m,n)=(3,2) → θ ≈ 13.17°, 76 atoms.

    Twist angle: cos(θ) = (m² + 4mn + n²) / [2(m² + mn + n²)]
    Number of atoms per cell: 4(m² + mn + n²)
    """
    print(f"\n=== 3. Twisted Bilayer Graphene (m={m}, n={n}) ===")

    N_cell = m**2 + m*n + n**2
    N_atoms = 4 * N_cell
    cos_theta = (m**2 + 4*m*n + n**2) / (2 * N_cell)
    theta = np.arccos(np.clip(cos_theta, -1, 1))
    theta_deg = np.degrees(theta)
    print(f"  Twist angle: θ = {theta_deg:.2f}°")
    print(f"  Atoms in supercell: {N_atoms}")

    # Graphene primitive vectors (Cartesian, Å)
    a1 = a_gr * np.array([1.0, 0.0])
    a2 = a_gr * np.array([0.5, np.sqrt(3)/2])

    # Superlattice vectors for layer 1 (unrotated convention)
    # L1 = m*a1 + n*a2,  L2 = -n*a1 + (m+n)*a2
    L1 = m * a1 + n * a2
    L2 = -n * a1 + (m + n) * a2
    L_mag = np.linalg.norm(L1)
    print(f"  Superlattice constant: {L_mag:.4f} Å")
    # Angle between L1, L2
    cos_L = np.dot(L1, L2) / (np.linalg.norm(L1) * np.linalg.norm(L2))
    gamma_deg = np.degrees(np.arccos(np.clip(cos_L, -1, 1)))

    # Basis atoms in graphene unit cell (Cartesian)
    basis_A = np.array([0.0, 0.0])
    basis_B = (a1 + a2) / 3  # = (1/3, 2/3) in fractional → Cartesian

    # Generate all atoms in supercell for layer 1 (unrotated, z=0)
    # Scan integer combinations of a1, a2 and keep those inside supercell
    inv_L = np.linalg.inv(np.column_stack([L1, L2]))

    layer1_cart = []
    for i in range(-N_cell, N_cell+1):
        for j in range(-N_cell, N_cell+1):
            for base in [basis_A, basis_B]:
                pos = i * a1 + j * a2 + base
                frac = inv_L @ pos
                # Keep if inside [0,1) x [0,1)
                fx, fy = frac[0] % 1.0, frac[1] % 1.0
                if fx >= 1.0 - 1e-8: fx = 0.0
                if fy >= 1.0 - 1e-8: fy = 0.0
                # Check uniqueness
                pos_check = fx * L1 + fy * L2
                is_dup = False
                for existing in layer1_cart:
                    diff = pos_check - existing
                    if np.linalg.norm(diff) < 0.1:
                        is_dup = True
                        break
                if not is_dup:
                    layer1_cart.append(pos_check)

    print(f"  Layer 1 atoms found: {len(layer1_cart)} (expect {2*N_cell})")

    # Layer 2: rotate by θ about origin, then shift up by d_AB
    R = np.array([[np.cos(theta), -np.sin(theta)],
                  [np.sin(theta),  np.cos(theta)]])

    # For layer 2, also apply AB offset before rotation
    ab_shift = (a1 + a2) / 3  # standard AB shift

    layer2_cart = []
    for i in range(-N_cell-2, N_cell+3):
        for j in range(-N_cell-2, N_cell+3):
            for base in [basis_A, basis_B]:
                pos_orig = i * a1 + j * a2 + base + ab_shift
                pos_rot = R @ pos_orig
                frac = inv_L @ pos_rot
                fx, fy = frac[0] % 1.0, frac[1] % 1.0
                if fx >= 1.0 - 1e-8: fx = 0.0
                if fy >= 1.0 - 1e-8: fy = 0.0
                pos_check = fx * L1 + fy * L2
                is_dup = False
                for existing in layer2_cart:
                    diff = pos_check - existing
                    if np.linalg.norm(diff) < 0.1:
                        is_dup = True
                        break
                if not is_dup:
                    layer2_cart.append(pos_check)

    print(f"  Layer 2 atoms found: {len(layer2_cart)} (expect {2*N_cell})")

    # Build CIF: use superlattice as cell (hexagonal-like)
    c_cell = 20.0  # Å vacuum
    z1 = 0.5 - d_AB / (2 * c_cell)
    z2 = 0.5 + d_AB / (2 * c_cell)

    atoms = []
    idx = 1
    for pos in layer1_cart:
        frac = inv_L @ pos
        fx, fy = frac[0] % 1.0, frac[1] % 1.0
        atoms.append((f"C{idx}", "C", fx, fy, z1, 0.005))
        idx += 1
    for pos in layer2_cart:
        frac = inv_L @ pos
        fx, fy = frac[0] % 1.0, frac[1] % 1.0
        atoms.append((f"C{idx}", "C", fx, fy, z2, 0.005))
        idx += 1

    # Cell parameters
    a_cell = np.linalg.norm(L1)
    b_cell = np.linalg.norm(L2)

    return write_cif(
        f"03_twisted_bilayer_m{m}n{n}.cif",
        f"twisted_bilayer_theta_{theta_deg:.1f}",
        a_cell, b_cell, c_cell, 90.0, 90.0, gamma_deg,
        "P 1", 1, atoms,
        f"Commensurate twisted bilayer graphene, theta={theta_deg:.2f} deg, "
        f"(m,n)=({m},{n}), {len(atoms)} atoms"
    )


# =========================================================================
# 4. BaTiO₃ TETRAGONAL PEROVSKITE (PIEZOELECTRIC SUBSTRATE)
# =========================================================================
def gen_BaTiO3():
    """BaTiO₃ in tetragonal P4mm phase — the piezoelectric substrate."""
    print("\n=== 4. BaTiO₃ Tetragonal Perovskite (P4mm) ===")
    # Tetragonal: a = 3.994 Å, c = 4.038 Å
    # Wyckoff positions for P4mm (#99):
    #   Ba at 1a (0, 0, 0)
    #   Ti at 1b (1/2, 1/2, z_Ti) with z_Ti ≈ 0.5126
    #   O1 at 1b (1/2, 1/2, z_O1) with z_O1 ≈ 0.0230
    #   O2 at 2c (1/2, 0, z_O2) with z_O2 ≈ 0.5147
    atoms = [
        ("Ba1", "Ba", 0.000000, 0.000000, 0.000000, 0.006),
        ("Ti1", "Ti", 0.500000, 0.500000, 0.512600, 0.005),
        ("O1",  "O",  0.500000, 0.500000, 0.023000, 0.008),
        ("O2",  "O",  0.500000, 0.000000, 0.514700, 0.008),
    ]
    # Note: O2 at 2c generates (1/2,0,z) and (0,1/2,z) by symmetry
    return write_cif("04_BaTiO3_P4mm.cif", "BaTiO3_tetragonal",
                     3.994, 3.994, 4.038, 90.0, 90.0, 90.0,
                     "P 4 m m", 99, atoms,
                     "BaTiO3 tetragonal perovskite, piezoelectric substrate for TEVC Platform I")


# =========================================================================
# 5. BaTiO₃ 3×3×1 SUPERCELL (for substrate visualization)
# =========================================================================
def gen_BaTiO3_supercell(nx=3, ny=3, nz=1):
    """BaTiO₃ supercell for substrate visualization."""
    print(f"\n=== 5. BaTiO₃ {nx}×{ny}×{nz} Supercell ===")
    a0, c0 = 3.994, 4.038

    # Primitive cell atoms (fractional in primitive cell)
    prim_atoms = [
        ("Ba", 0.0, 0.0, 0.0),
        ("Ti", 0.5, 0.5, 0.5126),
        ("O",  0.5, 0.5, 0.023),
        ("O",  0.5, 0.0, 0.5147),
        ("O",  0.0, 0.5, 0.5147),
    ]

    atoms = []
    idx = 1
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                for sym, fx, fy, fz in prim_atoms:
                    sfx = (fx + ix) / nx
                    sfy = (fy + iy) / ny
                    sfz = (fz + iz) / nz
                    atoms.append((f"{sym}{idx}", sym, sfx, sfy, sfz, 0.006))
                    idx += 1

    return write_cif(
        f"05_BaTiO3_{nx}x{ny}x{nz}_supercell.cif",
        f"BaTiO3_{nx}x{ny}x{nz}",
        a0*nx, a0*ny, c0*nz, 90.0, 90.0, 90.0,
        "P 1", 1, atoms,
        f"BaTiO3 {nx}x{ny}x{nz} supercell, {len(atoms)} atoms"
    )


# =========================================================================
# 6. HOPFION SPACETIME CRYSTAL — BCC lattice of topological charges
# =========================================================================
def gen_hopfion_crystal():
    """
    Hopfion spacetime crystal: BCC arrangement of Q=1 Faddeev-Niemi hopfions.

    Lattice constant set by coherence length ξ = 50 nm.
    We use a scaled representation (1 Å ↔ 1 nm for VESTA display).
    Nodes: Ho (holmium, Z=67, as visual proxy for hopfion center, blue in VESTA)
    Linking ring: He (helium, Z=2, small white dots tracing the Hopf fiber)

    This is a *schematic crystal* for VESTA rendering — not a real material.
    The physical content is the BCC topology with Q=1 at each site.
    """
    print("\n=== 6. Hopfion Spacetime Crystal (BCC, ξ-scale) ===")

    # Use ξ = 50 as lattice constant (in "display Å" = nm)
    a_hopf = 50.0  # "Å" in CIF ↔ 50 nm physical

    atoms = []
    idx = 1

    # BCC: corners + body center
    # Corner hopfion centers
    bcc_sites = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.5),  # body center
    ]
    for fx, fy, fz in bcc_sites:
        atoms.append((f"Ho{idx}", "Ho", fx, fy, fz, 0.050))
        idx += 1

    # Preimage ring (Hopf fiber): a torus around the body center
    # Major radius R = 0.3a, minor radius r = 0.1a
    # Trace the torus with He atoms
    R_major = 0.30  # in fractional units
    r_minor = 0.10
    n_major = 16    # points around major circle
    n_minor = 8     # points around minor circle
    center = np.array([0.5, 0.5, 0.5])

    for i in range(n_major):
        phi = 2 * np.pi * i / n_major
        for j in range(n_minor):
            psi = 2 * np.pi * j / n_minor
            # Torus parametrization in fractional coords
            x = center[0] + (R_major + r_minor * np.cos(psi)) * np.cos(phi)
            y = center[1] + (R_major + r_minor * np.cos(psi)) * np.sin(phi)
            z = center[2] + r_minor * np.sin(psi)
            # Wrap to [0, 1)
            x, y, z = x % 1.0, y % 1.0, z % 1.0
            atoms.append((f"He{idx}", "He", x, y, z, 0.020))
            idx += 1

    # Second preimage: a linked circle (the other Hopf fiber)
    # This is a circle threaded through the torus hole
    n_link = 24
    R_link = 0.20
    for i in range(n_link):
        phi = 2 * np.pi * i / n_link
        x = center[0] + R_link * np.cos(phi)
        y = center[1]
        z = center[2] + R_link * np.sin(phi)
        x, y, z = x % 1.0, y % 1.0, z % 1.0
        atoms.append((f"Ne{idx}", "Ne", x, y, z, 0.015))
        idx += 1

    return write_cif(
        "06_hopfion_crystal_BCC.cif",
        "hopfion_spacetime_crystal",
        a_hopf, a_hopf, a_hopf, 90.0, 90.0, 90.0,
        "P 1", 1, atoms,
        f"Hopfion spacetime crystal (BCC), a=xi=50nm, Q=1 per site, "
        f"Ho=hopfion center, He=preimage torus, Ne=linked fiber. "
        f"{len(atoms)} sites."
    )


# =========================================================================
# 7. GAUSSIAN CUBE FILE — Hopfion topological charge density
# =========================================================================
def gen_hopfion_cube(Ngrid=60):
    """
    Compute the Whitehead hopfion Q=1 field n(r): R³ → S²,
    then evaluate the topological charge density
        ρ_H = (1/4π²) ε^{ijk} A_i F_{jk}
    and write to Gaussian cube format for VESTA isosurface rendering.

    The Hopf map: stereographic projection R³ → S³ → S² (π)
    For point r = (x,y,z) ∈ R³:
        w = (2r, 1-|r|²) / (1+|r|²) ∈ S³   [inverse stereographic]
        n = π(w) = (2Re(w₁w₂*), 2Im(w₁w₂*), |w₁|²-|w₂|²)
    where w₁ = w[0] + i·w[1], w₂ = w[2] + i·w[3].
    """
    print(f"\n=== 7. Hopfion Topological Charge Density (Cube File, {Ngrid}³) ===")

    # Grid in physical space (centered at origin, extend to ±L)
    L = 3.0   # in units of ξ
    x = np.linspace(-L, L, Ngrid)
    dx = x[1] - x[0]
    X, Y, Z = np.meshgrid(x, x, x, indexing='ij')
    r2 = X**2 + Y**2 + Z**2

    # Inverse stereographic: R³ → S³
    denom = 1.0 + r2
    w0 = 2*X / denom
    w1 = 2*Y / denom
    w2 = 2*Z / denom
    w3 = (1.0 - r2) / denom

    # Hopf map: S³ → S² via π(w) treating (w0+iw1, w2+iw3) as C²
    # n_x = 2 Re(z1 z2*) = 2(w0 w2 + w1 w3)
    # n_y = 2 Im(z1 z2*) = 2(w1 w2 - w0 w3)
    # n_z = |z1|² - |z2|² = w0² + w1² - w2² - w3²
    nx = 2*(w0*w2 + w1*w3)
    ny = 2*(w1*w2 - w0*w3)
    nz = w0**2 + w1**2 - w2**2 - w3**2

    # Verify |n| = 1 (should be exact up to numerics)
    n_norm = np.sqrt(nx**2 + ny**2 + nz**2)
    print(f"  |n| range: [{n_norm.min():.6f}, {n_norm.max():.6f}] (should be 1.0)")

    # Compute topological charge density using the Hopf invariant integrand:
    #   ρ = (1/4π²) ε^{ijk} F_ij A_k
    # where F_ij = (1/2) ε_{abc} n_a ∂_i(n_b) ∂_j(n_c)
    # and A_i = (1/(1+n_z)) (n_x ∂_i n_y - n_y ∂_i n_x)  [Hopf connection]

    # Compute gradients (central differences)
    def grad(f):
        """Gradient of 3D array using central differences."""
        gx = np.zeros_like(f)
        gy = np.zeros_like(f)
        gz = np.zeros_like(f)
        gx[1:-1,:,:] = (f[2:,:,:] - f[:-2,:,:]) / (2*dx)
        gy[:,1:-1,:] = (f[:,2:,:] - f[:,:-2,:]) / (2*dx)
        gz[:,:,1:-1] = (f[:,:,2:] - f[:,:,:-2]) / (2*dx)
        return gx, gy, gz

    dnx = grad(nx)  # (∂nx/∂x, ∂nx/∂y, ∂nx/∂z)
    dny = grad(ny)
    dnz = grad(nz)

    # Berry connection A_i = (1/(1+nz)) (nx ∂_i ny - ny ∂_i nx)
    # Use south-pole gauge where 1+nz > 0 almost everywhere
    denom_gauge = 1.0 + nz + 1e-12  # regularize at south pole
    Ax = (nx * dny[0] - ny * dnx[0]) / denom_gauge
    Ay = (nx * dny[1] - ny * dnx[1]) / denom_gauge
    Az = (nx * dny[2] - ny * dnx[2]) / denom_gauge

    # Berry curvature F_ij = (1/2) ε_{abc} n_a ∂_i n_b ∂_j n_c
    # F_xy = nx(∂x ny ∂y nz - ∂y ny ∂x nz) + ny(∂x nz ∂y nx - ∂y nz ∂x nx)
    #        + nz(∂x nx ∂y ny - ∂y nx ∂x ny)
    # Simplified: F_ij = ε_{abc} n_a (∂_i n_b)(∂_j n_c)  [factor of 1/2 absorbed]
    Fxy = nx*(dny[0]*dnz[1] - dny[1]*dnz[0]) + \
          ny*(dnz[0]*dnx[1] - dnz[1]*dnx[0]) + \
          nz*(dnx[0]*dny[1] - dnx[1]*dny[0])
    Fyz = nx*(dny[1]*dnz[2] - dny[2]*dnz[1]) + \
          ny*(dnz[1]*dnx[2] - dnz[2]*dnx[1]) + \
          nz*(dnx[1]*dny[2] - dnx[2]*dny[1])
    Fzx = nx*(dny[2]*dnz[0] - dny[0]*dnz[2]) + \
          ny*(dnz[2]*dnx[0] - dnz[0]*dnx[2]) + \
          nz*(dnx[2]*dny[0] - dnx[0]*dny[2])

    # Hopf density: ρ = (1/4π²) ε^{ijk} F_{ij} A_k
    # = (1/4π²) [F_{xy} A_z + F_{yz} A_x + F_{zx} A_y
    #            - F_{xy} A_z is already there...
    # Actually: ε^{ijk} F_{ij} A_k = 2(F_{xy}A_z + F_{yz}A_x + F_{zx}A_y)
    rho_H = (1.0 / (4*np.pi**2)) * 2 * (Fxy*Az + Fyz*Ax + Fzx*Ay)

    # Integrate to get Q (should be 1 for Whitehead hopfion)
    Q_numerical = np.sum(rho_H[1:-1, 1:-1, 1:-1]) * dx**3
    print(f"  Hopf charge Q = ∫ρ_H d³x = {Q_numerical:.4f} (expect 1.0)")
    print(f"  ρ_H range: [{rho_H.min():.6f}, {rho_H.max():.6f}]")

    # Also compute energy density |∇n|²
    E_density = sum(dnx[i]**2 + dny[i]**2 + dnz[i]**2 for i in range(3))
    E_total = 0.5 * np.sum(E_density[1:-1,1:-1,1:-1]) * dx**3
    # Vakulenko-Kapitansky bound: E ≥ C|Q|^{3/4} with C = 16π²·(3^{3/8}/2^{7/2})
    # ≈ 163.5 for Q=1
    VK_bound = 16 * np.pi**2 * (3**(3/8)) / (2**(7/2))
    print(f"  Energy E = (1/2)∫|∇n|² = {E_total:.2f}")
    print(f"  VK bound  E ≥ C|Q|^(3/4) = {VK_bound:.2f}")
    print(f"  E/E_VK = {E_total/VK_bound:.4f} (>1 required)")

    # Write Gaussian cube file
    # Convention: atomic units (bohr), but we scale so 1 unit = ξ
    bohr = 0.529177  # Å per bohr (but we use our own length scale)
    cube_path = os.path.join(OUTDIR_CUBE, "hopfion_charge_density.cube")
    with open(cube_path, 'w') as f:
        f.write("Hopfion topological charge density rho_H(r)\n")
        f.write(f"Q={Q_numerical:.4f}, grid={Ngrid}^3, L={L:.1f}xi\n")
        # Number of atoms, origin
        origin = np.array([-L, -L, -L])
        f.write(f"    1 {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        # Voxel steps
        f.write(f" {Ngrid:4d} {dx:12.6f}  0.000000  0.000000\n")
        f.write(f" {Ngrid:4d}  0.000000 {dx:12.6f}  0.000000\n")
        f.write(f" {Ngrid:4d}  0.000000  0.000000 {dx:12.6f}\n")
        # Dummy atom at origin (required by format)
        f.write(f"   67  0.000000  0.000000  0.000000  0.000000\n")
        # Volumetric data (Fortran-style: z varies fastest)
        for ix in range(Ngrid):
            for iy in range(Ngrid):
                vals = rho_H[ix, iy, :]
                for iz in range(Ngrid):
                    f.write(f" {vals[iz]:13.5E}")
                    if (iz + 1) % 6 == 0:
                        f.write("\n")
                if Ngrid % 6 != 0:
                    f.write("\n")

    print(f"  Written: {cube_path}")

    # Also write energy density cube
    ecube_path = os.path.join(OUTDIR_CUBE, "hopfion_energy_density.cube")
    with open(ecube_path, 'w') as f:
        f.write("Hopfion energy density (1/2)|grad n|^2\n")
        f.write(f"E_total={E_total:.4f}, VK_bound={VK_bound:.4f}\n")
        f.write(f"    1 {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        f.write(f" {Ngrid:4d} {dx:12.6f}  0.000000  0.000000\n")
        f.write(f" {Ngrid:4d}  0.000000 {dx:12.6f}  0.000000\n")
        f.write(f" {Ngrid:4d}  0.000000  0.000000 {dx:12.6f}\n")
        f.write(f"   67  0.000000  0.000000  0.000000  0.000000\n")
        for ix in range(Ngrid):
            for iy in range(Ngrid):
                vals = E_density[ix, iy, :]
                for iz in range(Ngrid):
                    f.write(f" {vals[iz]:13.5E}")
                    if (iz + 1) % 6 == 0:
                        f.write("\n")
                if Ngrid % 6 != 0:
                    f.write("\n")

    print(f"  Written: {ecube_path}")

    return Q_numerical, E_total, VK_bound, rho_H, E_density, nx, ny, nz, x


# =========================================================================
# MAIN EXECUTION
# =========================================================================
if __name__ == "__main__":
    print("=" * 70)
    print("TEVC Structure Generator — VESTA-Compatible CIF + Cube Files")
    print("=" * 70)

    # Crystal structures
    gen_graphene_monolayer()
    gen_bilayer_AB()
    gen_twisted_bilayer(m=2, n=1)   # θ = 21.79°
    gen_twisted_bilayer(m=3, n=2)   # θ = 13.17°
    gen_BaTiO3()
    gen_BaTiO3_supercell(3, 3, 2)
    gen_hopfion_crystal()

    # Hopfion field computation + cube files
    Q, E, VK, rho, Edens, nx, ny, nz, grid = gen_hopfion_cube(Ngrid=60)

    # === Numerical Readout Summary ===
    print("\n" + "=" * 70)
    print("NUMERICAL READOUT SUMMARY")
    print("=" * 70)
    print(f"\n--- Locked QVC v4 Parameters ---")
    print(f"  F (fat torus)    = {F_torus}")
    print(f"  E_vac            = {E_vac} meV")
    print(f"  g₀               = {g0} meV")
    print(f"  λ*               = {lambda_star}")
    print(f"  Δ*               = {Delta_star} meV")
    print(f"  T*               = {T_star*1e3:.1f} mK")
    print(f"  ℏc_s             = {hbar_cs} meV·μm")
    print(f"  ξ                = {xi_nm} nm")
    print(f"  η_ret            = {eta_ret}")
    print(f"  A*/A_c           = {A_over_Ac}")

    print(f"\n--- Graphene Lattice ---")
    print(f"  a (lattice)      = {a_gr} Å")
    print(f"  C-C bond         = {a_CC} Å")
    print(f"  d (interlayer)   = {d_AB} Å")

    print(f"\n--- BaTiO₃ Substrate ---")
    print(f"  a                = 3.994 Å")
    print(f"  c                = 4.038 Å")
    print(f"  c/a              = {4.038/3.994:.4f}")
    print(f"  Space group      = P4mm (#99)")
    print(f"  Ti displacement  = {(0.5126-0.5)*4.038:.4f} Å from centrosym")

    print(f"\n--- Twisted Bilayer (m=2,n=1) ---")
    N21 = 4*(4+2+1)
    theta21 = np.degrees(np.arccos((4+8+1)/(2*7)))
    L21 = a_gr * np.sqrt(7)
    print(f"  θ                = {theta21:.2f}°")
    print(f"  Supercell atoms  = {N21}")
    print(f"  Moiré wavelength = {L21:.4f} Å")

    print(f"\n--- Hopfion Topology (Whitehead, Q=1) ---")
    print(f"  Hopf charge Q    = {Q:.4f} (analytic: 1)")
    print(f"  Energy E         = {E:.2f} (dimensionless)")
    print(f"  VK bound         = {VK:.2f}")
    print(f"  E/E_VK           = {E/VK:.4f} (must be ≥ 1)")
    print(f"  BCC lattice const= {xi_nm} nm (= ξ)")

    print(f"\n--- Output Files ---")
    for d in [OUTDIR_CIF, OUTDIR_CUBE]:
        for fn in sorted(os.listdir(d)):
            fpath = os.path.join(d, fn)
            size_kb = os.path.getsize(fpath) / 1024
            print(f"  {fpath}  ({size_kb:.1f} kB)")

    print("\nDone. Open CIF files in VESTA for 3D rendering.")
    print("Open .cube files in VESTA → File → Open → Isosurface for volumetric data.")
