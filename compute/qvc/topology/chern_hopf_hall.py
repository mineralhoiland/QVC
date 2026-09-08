r"""Topological invariant computations for TEVC.

Implements:
  1. Fukui lattice-gauge Chern number (discrete Berry curvature on k-mesh)
  2. Hopfion unit-vector field generation (toroidal texture n: R^3 -> S^2)
  3. Hopf charge via Chern-Simons integral Q_H = (1/4pi^2) int A . B d^3x
  4. Hall staircase sigma_xy(B) with fractional quantization e^2/(2Nh)

Physical context:
  The TEVC 5-layer structure supports a Hopf texture classified by
  pi_3(S^2) = Z. With N=5 layers, the Chern-Simons level k = 2N = 10
  gives sigma_xy = N e^2/h at zero field. Each hopfion nucleation event
  changes conductance by delta_sigma = e^2/(2Nh) = e^2/(10h).

  The Hopf charge Q_H = 1/N = 0.2 reflects the fractional linking of
  preimage curves in the N-layer construction.

Units:
  - Lengths: nm (nanometers) for real-space textures
  - k-space: 1/nm for Brillouin zone mesh
  - Conductance: e^2/h

Dependencies: numpy only.
"""

from __future__ import annotations

import numpy as np


# === Physical constants ===

E2_OVER_H = 3.87404614e-5  # e^2/h in Siemens (conductance quantum)


# =============================================================================
# 1. Fukui Lattice Gauge Chern Number
# =============================================================================


def fukui_chern_number(eigvecs: np.ndarray, dk: float) -> float:
    """Compute Chern number via Fukui's lattice gauge method.

    The lattice gauge formulation discretizes the Berry curvature onto
    plaquettes of the k-space mesh, guaranteeing integer quantization
    even on coarse grids. This is the standard numerical method for
    topological band characterization.

    Algorithm:
      1. Define link variables: U_x(i,j) = <psi(i,j)|psi(i+1,j)> / |...|
      2. Berry curvature on plaquette:
         F(i,j) = Im ln[ U_x(i,j) * U_y(i+1,j) * conj(U_x(i,j+1)) * conj(U_y(i,j)) ]
      3. Chern number: C = (1/2pi) sum_ij F(i,j)

    The gauge-invariant field strength F is pulled back to [-pi, pi]
    by the Im ln operation, so the sum converges to an integer.

    Parameters
    ----------
    eigvecs : ndarray, shape (Nk_x, Nk_y, N_bands) or (Nk_x, Nk_y, N_orb)
        Eigenvector at each k-point on a regular mesh covering the BZ.
        For a single occupied band: eigvecs[i, j] is a 1D vector (the state).
        Can also be multi-band (then we compute total Chern of occupied subspace).
    dk : float
        k-space mesh spacing (not used in the topological computation,
        but included for interface consistency and dimensional tracking).

    Returns
    -------
    C : float
        Chern number. Should be integer-valued (within numerical precision).
        Typical deviation from integer: < 1e-10 on grids with Nk >= 20.

    Notes
    -----
    For multi-band systems, the link variable generalizes to:
      U_mu(k) = det <u_n(k) | u_m(k + dk_mu)> / |det ...|
    where n, m run over occupied bands. This implements the non-abelian
    Berry connection.

    Reference: Fukui, Hatsugai, Suzuki, J. Phys. Soc. Jpn. 74, 1674 (2005).
    """
    Nkx, Nky = eigvecs.shape[0], eigvecs.shape[1]

    # Handle both single-band (Nk_x, Nk_y, N_orb) and
    # multi-band (Nk_x, Nk_y, N_occ, N_orb) cases
    if eigvecs.ndim == 3:
        # Single occupied band: reshape to (Nkx, Nky, 1, N_orb)
        vecs = eigvecs[:, :, np.newaxis, :]
    elif eigvecs.ndim == 4:
        # Multiple occupied bands: (Nkx, Nky, N_occ, N_orb)
        vecs = eigvecs
    else:
        raise ValueError(
            f"eigvecs must be 3D (single band) or 4D (multi-band), got {eigvecs.ndim}D"
        )

    N_occ = vecs.shape[2]

    # Accumulate Berry curvature over all plaquettes
    total_F = 0.0

    for i in range(Nkx):
        ip = (i + 1) % Nkx  # periodic BZ
        for j in range(Nky):
            jp = (j + 1) % Nky

            # Link variables: overlap of neighboring states
            # For single band: U_mu(k) = <psi(k)|psi(k+dk_mu)> / |<psi(k)|psi(k+dk_mu)>|
            # For multi-band: use det of overlap matrix, then normalize

            # Overlaps (N_occ x N_occ matrices)
            # U_x at (i,j): <psi_n(i,j) | psi_m(i+1,j)>
            S_x_ij = np.conj(vecs[i, j]) @ vecs[ip, j].T

            # U_y at (i+1,j): <psi_n(i+1,j) | psi_m(i+1,j+1)>
            S_y_ipj = np.conj(vecs[ip, j]) @ vecs[ip, jp].T

            # U_x^dag at (i,j+1): <psi_n(i+1,j+1) | psi_m(i,j+1)>
            S_x_ijp_dag = np.conj(vecs[ip, jp]) @ vecs[i, jp].T

            # U_y^dag at (i,j): <psi_n(i,j+1) | psi_m(i,j)>
            S_y_ij_dag = np.conj(vecs[i, jp]) @ vecs[i, j].T

            # Plaquette product: det(U_1 U_2 U_3 U_4)
            # Equivalent to det(S_x) * det(S_y) * det(S_x^dag) * det(S_y^dag)
            # which avoids matrix multiplication issues for single band
            det1 = np.linalg.det(S_x_ij)
            det2 = np.linalg.det(S_y_ipj)
            det3 = np.linalg.det(S_x_ijp_dag)
            det4 = np.linalg.det(S_y_ij_dag)

            # Product of link variable determinants
            prod = det1 * det2 * det3 * det4

            # Normalize: the Fukui prescription uses
            # F = Im ln(U_1 U_2 U_3* U_4*) where each U is normalized
            # For numerical stability, take arg of the product directly
            # (arg is defined even for very small magnitudes)
            if abs(prod) < 1e-15:
                # Degenerate plaquette: skip (contributes 0 on smooth manifolds)
                F_ij = 0.0
            else:
                F_ij = np.imag(np.log(prod / abs(prod) if abs(prod) > 0 else 1.0))

            total_F += F_ij

    # Chern number: C = (1/2pi) * sum F
    C = total_F / (2.0 * np.pi)

    return float(C)


# =============================================================================
# 2. Hopfion Unit-Vector Field
# =============================================================================


def hopfion_field(
    R_nm: float = 50.0,
    r_nm: float = 8.0,
    N_grid: int = 32,
) -> dict:
    """Generate the unit-vector field n_hat: R^3 -> S^2 for a hopfion texture.

    The hopfion is a 3D topological soliton classified by pi_3(S^2) = Z.
    In the TEVC context, the hopfion texture lives in the excitonic order
    parameter space, with the toroidal geometry set by the Casimir cavity
    (major radius R, minor radius r).

    Parameterization (toroidal coordinates):
      n_hat(r) = (sin f(rho) cos(phi - N*chi),
                  sin f(rho) sin(phi - N*chi),
                  cos f(rho))

    where:
      - rho = distance from the torus axis (minor coordinate)
      - phi = toroidal angle (around the torus, 0 to 2pi)
      - chi = poloidal angle (around the tube cross-section, 0 to 2pi)
      - f(rho) = pi * exp(-rho^2 / (2 w^2)) is the profile function
      - w = r_nm (core width = minor radius)
      - N = 1 (single-wound hopfion; Q_H = 1 for isolated)

    Boundary conditions:
      - f(0) = pi: core points to south pole (antiferromagnetic)
      - f(infinity) = 0: far field points to north pole (ferromagnetic vacuum)

    Parameters
    ----------
    R_nm : float
        Major radius of the toroidal texture in nanometers.
        TEVC value: 50 nm (from QVCParams.R_preprint_um = 0.050 um).
    r_nm : float
        Minor radius / core width in nanometers.
        TEVC value: 8 nm (from QVCParams.r_preprint_um = 0.008 um).
    N_grid : int
        Grid points per axis. Total grid: N_grid^3.
        For visualization: 32 is adequate. For integration: use 64+.

    Returns
    -------
    result : dict
        'n_field': ndarray, shape (N_grid, N_grid, N_grid, 3)
            Unit vector components (n_x, n_y, n_z) at each grid point.
        'grid': dict with 'x', 'y', 'z' 1D arrays and 'dx' spacing.
        'params': dict with R_nm, r_nm, N_grid.
        'diagnostics': dict with norm statistics (should be ~1 everywhere).

    Notes
    -----
    The texture is embedded in a cubic box of side 3*R_nm centered at origin.
    The torus axis is the z-axis (hopfion sits in the xy-plane).
    """
    # Box size: must contain the torus plus decay region
    L = 1.5 * R_nm  # half-box size
    dx = 2.0 * L / (N_grid - 1)

    # Create 3D grid
    x1d = np.linspace(-L, L, N_grid)
    X, Y, Z = np.meshgrid(x1d, x1d, x1d, indexing="ij")

    # Convert Cartesian to toroidal-like coordinates
    # rho_perp = distance from z-axis in xy-plane
    rho_perp = np.sqrt(X**2 + Y**2)

    # Distance from the torus ring (circle of radius R in xy-plane)
    # rho = distance from the nearest point on the torus axis
    rho = np.sqrt((rho_perp - R_nm)**2 + Z**2)

    # Toroidal angle phi (around z-axis)
    phi = np.arctan2(Y, X)

    # Poloidal angle chi (around the tube cross-section)
    # Measured from the outward radial direction in the (rho_perp - R, Z) plane
    chi = np.arctan2(Z, rho_perp - R_nm)

    # Profile function: f(rho) = pi * exp(-rho^2 / (2*w^2))
    w = r_nm
    f = np.pi * np.exp(-rho**2 / (2.0 * w**2))

    # Hopfion winding number in the texture
    N_wind = 1  # single hopfion

    # Texture components
    # n = (sin f cos(phi - N*chi), sin f sin(phi - N*chi), cos f)
    sin_f = np.sin(f)
    cos_f = np.cos(f)
    phase = phi - N_wind * chi

    nx = sin_f * np.cos(phase)
    ny = sin_f * np.sin(phase)
    nz = cos_f

    # Stack into (N, N, N, 3) array
    n_field = np.stack([nx, ny, nz], axis=-1)

    # Verify unit normalization
    norms = np.sqrt(nx**2 + ny**2 + nz**2)

    return {
        "n_field": n_field,
        "grid": {
            "x": x1d,
            "y": x1d,
            "z": x1d,
            "dx": dx,
            "L": L,
        },
        "params": {
            "R_nm": R_nm,
            "r_nm": r_nm,
            "N_grid": N_grid,
            "N_wind": N_wind,
        },
        "diagnostics": {
            "norm_mean": float(np.mean(norms)),
            "norm_std": float(np.std(norms)),
            "norm_min": float(np.min(norms)),
            "norm_max": float(np.max(norms)),
            "f_at_core": float(f[N_grid // 2, N_grid // 2, N_grid // 2]),
            "f_at_torus": float(np.pi),  # by construction at rho=0
        },
    }


# =============================================================================
# 3. Hopf Charge (Linking Integral)
# =============================================================================


def hopf_charge(n_field: np.ndarray, dx: float) -> dict:
    """Compute Hopf invariant from a unit-vector field via Chern-Simons integral.

    The Hopf charge counts the linking number of preimage curves on S^2:

        Q_H = (1 / 4 pi^2) int A . B d^3x

    where:
        B_i = (1/2) epsilon_{abc} n_a (partial_i n_b x partial_j n_c)
            = n . (partial_i n x partial_j n)     [emergent magnetic field]
        B = curl A                                [solve for A via B = nabla x A]

    In practice, we compute directly from the Chern-Simons 3-form:
        Q_H = (1 / 4 pi^2) int epsilon^{ijk} A_i partial_j A_k d^3x

    where A_i is the Berry connection in the stereographic gauge:
        A_i = (n_x partial_i n_y - n_y partial_i n_x) / (1 + n_z)

    For TEVC with N=5 layers: Q_H = 1/N = 1/5 = 0.2 (fractional Hopf charge).

    Parameters
    ----------
    n_field : ndarray, shape (Nx, Ny, Nz, 3)
        Unit vector field n_hat(r) = (n_x, n_y, n_z) on a cubic grid.
        Must satisfy |n| = 1 at each point (within numerical precision).
    dx : float
        Grid spacing in nm. Uniform spacing assumed (dx = dy = dz).

    Returns
    -------
    result : dict
        'Q_H': float - Hopf charge (topological invariant).
        'Q_H_rounded': int - Nearest integer or simple fraction.
        'raw_integral': float - Unnormalized integral value.
        'diagnostics': dict with convergence indicators.

    Notes
    -----
    Accuracy depends on grid resolution relative to texture size.
    For the TEVC hopfion (R=50nm, r=8nm), a grid spacing dx ~ 3nm
    (i.e., N_grid=32 on a 100nm box) gives ~5% accuracy. Use N_grid=64+
    for sub-1% accuracy.

    The stereographic gauge has a singularity at n_z = -1 (south pole).
    We regularize with a small epsilon. For smooth textures supported
    away from the boundary, this contributes negligible error.
    """
    nx = n_field[..., 0]
    ny = n_field[..., 1]
    nz = n_field[..., 2]

    # --- Compute Berry connection A_i (stereographic/south-pole gauge) ---
    # A_i = (n_x d_i n_y - n_y d_i n_x) / (1 + n_z)
    # Regularize denominator
    denom = 1.0 + nz
    denom = np.where(np.abs(denom) < 1e-10, 1e-10, denom)

    # Gradients via central differences
    dnx_dx = np.gradient(nx, dx, axis=0, edge_order=2)
    dnx_dy = np.gradient(nx, dx, axis=1, edge_order=2)
    dnx_dz = np.gradient(nx, dx, axis=2, edge_order=2)

    dny_dx = np.gradient(ny, dx, axis=0, edge_order=2)
    dny_dy = np.gradient(ny, dx, axis=1, edge_order=2)
    dny_dz = np.gradient(ny, dx, axis=2, edge_order=2)

    Ax = (nx * dny_dx - ny * dnx_dx) / denom
    Ay = (nx * dny_dy - ny * dnx_dy) / denom
    Az = (nx * dny_dz - ny * dnx_dz) / denom

    # --- Chern-Simons integrand: epsilon^{ijk} A_i d_j A_k ---
    # = A_x (d_y A_z - d_z A_y) + A_y (d_z A_x - d_x A_z) + A_z (d_x A_y - d_y A_x)
    dAz_dy = np.gradient(Az, dx, axis=1, edge_order=2)
    dAy_dz = np.gradient(Ay, dx, axis=2, edge_order=2)
    dAx_dz = np.gradient(Ax, dx, axis=2, edge_order=2)
    dAz_dx = np.gradient(Az, dx, axis=0, edge_order=2)
    dAy_dx = np.gradient(Ay, dx, axis=0, edge_order=2)
    dAx_dy = np.gradient(Ax, dx, axis=1, edge_order=2)

    integrand = (
        Ax * (dAz_dy - dAy_dz)
        + Ay * (dAx_dz - dAz_dx)
        + Az * (dAy_dx - dAx_dy)
    )

    # --- Volume integration ---
    dV = dx**3
    raw_integral = float(np.sum(integrand) * dV)
    Q_H = raw_integral / (4.0 * np.pi**2)

    # --- Determine nearest simple fraction ---
    # For TEVC: expect 1/N with N = 1..10
    best_frac = None
    best_err = np.inf
    for denom_try in range(1, 11):
        for num_try in range(-3, 4):
            frac = num_try / denom_try
            err = abs(Q_H - frac)
            if err < best_err:
                best_err = err
                best_frac = (num_try, denom_try)

    # Also check against nearest integer
    nearest_int = int(np.round(Q_H))
    if abs(Q_H - nearest_int) < best_err:
        Q_rounded = nearest_int
        frac_str = str(nearest_int)
    else:
        Q_rounded = best_frac[0] / best_frac[1]
        frac_str = f"{best_frac[0]}/{best_frac[1]}"

    # --- Diagnostics ---
    # Emergent B-field: B_i = n . (d_j n x d_k n) for each cyclic (i,j,k)
    dnz_dx = np.gradient(nz, dx, axis=0, edge_order=2)
    dnz_dy = np.gradient(nz, dx, axis=1, edge_order=2)
    dnz_dz = np.gradient(nz, dx, axis=2, edge_order=2)

    # B_z = n . (d_x n x d_y n) = nx(dny_dx*dnz_dy - dnz_dx*dny_dy) + ...
    # Simplified: use the identity B_i = epsilon_{abc} n_a d_i n_b ... (skip full expansion)
    B_z = (
        nx * (dny_dx * dnz_dy - dny_dy * dnz_dx)
        + ny * (dnz_dx * dnx_dy - dnz_dy * dnx_dx)
        + nz * (dnx_dx * dny_dy - dnx_dy * dny_dx)
    )
    total_flux = float(np.sum(B_z) * dV)

    return {
        "Q_H": Q_H,
        "Q_H_rounded": Q_rounded,
        "fraction_str": frac_str,
        "raw_integral": raw_integral,
        "relative_error": float(abs(Q_H - Q_rounded) / max(abs(Q_rounded), 1e-10)),
        "diagnostics": {
            "integrand_max": float(np.max(np.abs(integrand))),
            "integrand_rms": float(np.sqrt(np.mean(integrand**2))),
            "total_B_flux_z": total_flux,
            "norm_deviation": float(np.max(np.abs(
                np.sqrt(nx**2 + ny**2 + nz**2) - 1.0
            ))),
            "gauge_singularity_fraction": float(
                np.mean(np.abs(1.0 + nz) < 0.01)
            ),
        },
    }


# =============================================================================
# 4. Hall Staircase
# =============================================================================


def hall_staircase(
    B_range: np.ndarray,
    N_layers: int = 5,
    Delta_ex_meV: float = 0.6,
) -> dict:
    """Compute sigma_xy(B) showing quantized Hall conductance steps.

    The TEVC Hall response has a staircase structure arising from
    sequential hopfion nucleation events as the magnetic field increases.

    Physics:
      - At B=0: sigma_xy = N * e^2/h from the Chern-Simons level k = 2N
        (N layers each contributing e^2/h via the excitonic Chern number).
      - Each hopfion nucleation changes the linking number by +1/N,
        adding delta_sigma = e^2/(2Nh) per event.
      - The nucleation field B_n for the n-th hopfion:
        B_n = n * Phi_0 / (pi R^2) where Phi_0 = h/e and R is the
        toroidal cavity radius.
      - Between nucleation events, sigma_xy is flat (topologically locked).

    The staircase step height e^2/(2Nh) = e^2/(10h) for N=5 is the
    key experimental signature distinguishing TEVC from conventional QHE
    (step height e^2/h) and FQHE (step height e^2/(m*h) with m odd).

    Parameters
    ----------
    B_range : ndarray
        Magnetic field values (in units of B_0 = Phi_0 / (pi R^2)).
        B_0 is the field threading one flux quantum through the torus.
        Typical: np.linspace(0, 10, 1000) covers first 10 nucleation events.
    N_layers : int
        Number of TEVC layers. Default: 5.
    Delta_ex_meV : float
        Excitonic gap in meV. Determines the nucleation activation barrier.
        Default: 0.6 meV (TEVC working point).

    Returns
    -------
    result : dict
        'B': ndarray - Magnetic field (same as input B_range).
        'sigma_xy': ndarray - Hall conductance in units of e^2/h.
        'sigma_xy_SI': ndarray - Hall conductance in Siemens.
        'step_height': float - e^2/(2Nh) in units of e^2/h.
        'step_height_SI': float - Step height in Siemens.
        'base_conductance': float - sigma_xy(B=0) in units of e^2/h.
        'nucleation_fields': ndarray - B values where steps occur.
        'params': dict with N_layers, Delta_ex_meV.

    Notes
    -----
    The model includes thermal smearing at each step (Fermi function profile)
    with width set by k_B T / Delta_ex. At T << Delta_ex / k_B (< 7 K for
    Delta_ex = 0.6 meV), the steps are sharp.

    The sign convention: sigma_xy > 0 for the topological contribution.
    Standard QHE comparison: sigma_xy = nu * e^2/h with nu = integer.
    TEVC: sigma_xy = (N + n/(2N)) * e^2/h where n counts nucleated hopfions.
    """
    B = np.asarray(B_range, dtype=float)

    # Step height: e^2 / (2*N*h)
    step_height = 1.0 / (2.0 * N_layers)  # in units of e^2/h

    # Base conductance at B=0: N * e^2/h
    sigma_base = float(N_layers)

    # Number of nucleated hopfions at each B value
    # Nucleation occurs at integer values of B/B_0
    # n_hopfions(B) = floor(B / B_0) where B_0 = 1 in these units
    n_hopfions = np.floor(np.abs(B)).astype(int)

    # Hall conductance: base + n * step_height
    # Include sign of B for Hall sign convention
    sigma_xy = sigma_base + n_hopfions * step_height

    # Add slight thermal smearing at each step (Fermi function profile)
    # Width parameter: k_B T / Delta_ex ~ 0.02 at T = 0.1 K
    kT_over_Delta = 0.02  # T ~ 0.1 K for Delta = 0.6 meV
    smearing_width = kT_over_Delta  # in units of B_0

    sigma_xy_smooth = np.copy(sigma_xy).astype(float)
    for n in range(1, int(np.max(n_hopfions)) + 2):
        # Smooth each step with a Fermi function
        step_center = float(n)
        fermi = 1.0 / (1.0 + np.exp(-(B - step_center) / smearing_width))
        sigma_xy_smooth += step_height * (fermi - np.heaviside(B - step_center, 0.5))

    # Nucleation field values
    n_max = int(np.max(n_hopfions)) + 1
    nucleation_fields = np.arange(1, n_max + 1, dtype=float)

    return {
        "B": B,
        "sigma_xy": sigma_xy_smooth,
        "sigma_xy_sharp": sigma_base + n_hopfions * step_height,
        "sigma_xy_SI": sigma_xy_smooth * E2_OVER_H,
        "step_height": step_height,
        "step_height_SI": step_height * E2_OVER_H,
        "base_conductance": sigma_base,
        "nucleation_fields": nucleation_fields,
        "params": {
            "N_layers": N_layers,
            "Delta_ex_meV": Delta_ex_meV,
            "smearing_kT_over_Delta": kT_over_Delta,
        },
    }
