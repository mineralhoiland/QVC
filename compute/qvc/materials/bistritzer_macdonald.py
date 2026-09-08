"""
qvc.materials.bistritzer_macdonald — Bistritzer-MacDonald continuum model for twisted multilayer graphene.

Implements the continuum Hamiltonian for the K-valley of twisted bilayer/multilayer
graphene (TBG/TLG) as described in Bistritzer & MacDonald, PNAS 108, 12233 (2011).

The BM Hamiltonian couples two Dirac cones (one per layer) via interlayer tunneling
modulated at the moire periodicity. For MATBG at theta ~ 1.1 deg, the resulting
moire minibands become extremely flat (W < 10 meV), enabling correlated phases.

Physical picture:
    - Each layer contributes a massless Dirac cone (sublattice A, B)
    - Rotation by +/-theta/2 shifts each layer's K-point
    - Interlayer tunneling T(r) has moire periodicity, with Fourier components T_j
    - Three q-vectors (120 deg apart) mediate the dominant tunneling processes
    - Plane-wave expansion truncated at n_shells moire reciprocal lattice shells

For N-layer generalization (TEVC: N=5, alternating +/-theta/2):
    Block-tridiagonal Hamiltonian with N Dirac blocks and N-1 tunneling blocks.

Units: Angstroms for lengths/momenta (1/Angstrom), eV for energies.

Key parameters at theta = 1.1 deg:
    a = 2.46 Angstrom (graphene lattice constant)
    hbar*v_F = 5.253 eV*Angstrom (for v_F ~ 0.8e6 m/s)
    w_AA = 0.0797 eV, w_AB = 0.1100 eV
    Moire period lambda_M ~ 128 Angstrom (12.8 nm)
    Expected: bandwidth W < 10 meV at magic angle, Chern number C = +/-1 per valley

References:
    [1] Bistritzer & MacDonald, PNAS 108, 12233 (2011)
    [2] Tarnopolsky, Kruchkov, Vishwanath, PRL 122, 106405 (2019)
    [3] Song, Wang, Bernevig, PRL 123, 036401 (2019)
"""

from __future__ import annotations

import numpy as np
from typing import Tuple, Optional


# === Physical Constants ===

A_GRAPHENE: float = 2.46  # Graphene lattice constant in Angstroms
HBAR_EV_S: float = 6.582119569e-16  # hbar in eV*s


# === Helper Functions ===


def _rotation_matrix_2d(theta_rad: float) -> np.ndarray:
    """2D rotation matrix for angle theta (radians)."""
    c, s = np.cos(theta_rad), np.sin(theta_rad)
    return np.array([[c, -s], [s, c]])


def _generate_g_vectors(
    G1: np.ndarray,
    G2: np.ndarray,
    n_shells: int,
) -> np.ndarray:
    """Generate moire reciprocal lattice vectors within circular cutoff.

    Includes all G = n1*G1 + n2*G2 with |G| <= n_shells * |G1|.
    This gives a well-converged plane-wave basis for the BM model.

    Parameters
    ----------
    G1, G2 : ndarray, shape (2,)
        Moire reciprocal lattice basis vectors (1/Angstrom).
    n_shells : int
        Number of reciprocal lattice shells to include.

    Returns
    -------
    G_list : ndarray, shape (N_G, 2)
        All G-vectors within the cutoff, sorted by magnitude.
    """
    cutoff = n_shells * np.linalg.norm(G1)
    G_list = []

    # Search over a box large enough to contain the circular cutoff
    n_max = n_shells + 1
    for n1 in range(-n_max, n_max + 1):
        for n2 in range(-n_max, n_max + 1):
            G = n1 * G1 + n2 * G2
            if np.linalg.norm(G) <= cutoff + 1e-10:
                G_list.append(G)

    G_list = np.array(G_list)
    # Sort by magnitude (primary) then angle (secondary) for reproducibility
    mags = np.linalg.norm(G_list, axis=1)
    angles = np.arctan2(G_list[:, 1], G_list[:, 0])
    sort_idx = np.lexsort((angles, mags))
    return G_list[sort_idx]


def k_path_high_symmetry(K_M: float, n_points: int = 100) -> Tuple[np.ndarray, np.ndarray, list]:
    """Generate high-symmetry k-path through the moire Brillouin zone.

    Path: Gamma -> K -> M -> Gamma

    The moire BZ is hexagonal with corners at distance K_M from Gamma.
    Convention: K = K_M*(0, 1), M = K_M*(sqrt(3)/4, 3/4).

    Parameters
    ----------
    K_M : float
        Moire BZ corner distance from Gamma (1/Angstrom).
        K_M = 2 * K_D * sin(theta/2) where K_D = 4*pi/(3*a).
    n_points : int
        Total number of k-points along the full path.

    Returns
    -------
    k_path : ndarray, shape (n_points, 2)
        k-point coordinates along the path (1/Angstrom).
    k_dist : ndarray, shape (n_points,)
        Cumulative distance along the path (for plotting x-axis).
    tick_positions : list of (float, str)
        (distance, label) pairs for high-symmetry point markers.
    """
    # High-symmetry points in our convention
    Gamma = np.array([0.0, 0.0])
    K = K_M * np.array([0.0, 1.0])
    M = K_M * np.array([np.sqrt(3) / 4, 3.0 / 4])

    segments = [(Gamma, K, "Gamma", "K"),
                (K, M, "K", "M"),
                (M, Gamma, "M", "Gamma")]

    # Distribute points proportional to segment lengths
    seg_lengths = [np.linalg.norm(b - a) for a, b, _, _ in segments]
    total_length = sum(seg_lengths)
    seg_npoints = [max(2, int(round(n_points * L / total_length))) for L in seg_lengths]
    # Adjust to match total
    seg_npoints[-1] = n_points - sum(seg_npoints[:-1])

    k_path = []
    k_dist = []
    tick_positions = []
    cumulative_dist = 0.0

    for seg_idx, (start, end, label_start, label_end) in enumerate(segments):
        n_seg = seg_npoints[seg_idx]
        # Don't duplicate endpoints between segments
        if seg_idx == 0:
            t_values = np.linspace(0, 1, n_seg)
            tick_positions.append((cumulative_dist, label_start))
        else:
            t_values = np.linspace(0, 1, n_seg + 1)[1:]  # skip first (= previous last)

        for t in t_values:
            k = start + t * (end - start)
            if len(k_path) > 0:
                cumulative_dist += np.linalg.norm(k - k_path[-1])
            k_path.append(k)
            k_dist.append(cumulative_dist)

        if seg_idx < len(segments) - 1:
            tick_positions.append((cumulative_dist, label_end))

    tick_positions.append((cumulative_dist, segments[-1][3]))

    return np.array(k_path), np.array(k_dist), tick_positions


def fukui_chern(eigvecs_mesh: np.ndarray) -> Tuple[float, np.ndarray]:
    """Compute Chern number via Fukui lattice gauge method.

    Uses the gauge-invariant discretization of Berry curvature on a regular mesh.
    The Berry curvature at each plaquette is:
        F_12(k) = Im[ln(U_1(k) * U_2(k+d1) * U_1^*(k+d2) * U_2^*(k))]
    where U_i(k) = <u(k)|u(k+di)> are link variables (overlap of adjacent eigenstates).

    The Chern number C = (1/2pi) * sum_k F_12(k) is guaranteed integer for a gapped band.

    Parameters
    ----------
    eigvecs_mesh : ndarray, shape (Nkx, Nky, dim)
        Eigenvector of the target band on a regular 2D mesh in the BZ.
        The mesh must be periodic (last point = first point is handled internally).

    Returns
    -------
    chern : float
        Chern number (should be close to integer for gapped bands).
    berry_curvature : ndarray, shape (Nkx-1, Nky-1)
        Berry curvature on each plaquette of the mesh.
    """
    Nkx, Nky, dim = eigvecs_mesh.shape

    # Link variables along kx (direction 1) and ky (direction 2)
    # U1[i,j] = <u(i,j)|u(i+1,j)>
    # U2[i,j] = <u(i,j)|u(i,j+1)>
    U1 = np.zeros((Nkx, Nky), dtype=complex)
    U2 = np.zeros((Nkx, Nky), dtype=complex)

    for i in range(Nkx):
        for j in range(Nky):
            ip1 = (i + 1) % Nkx
            jp1 = (j + 1) % Nky
            U1[i, j] = np.vdot(eigvecs_mesh[i, j], eigvecs_mesh[ip1, j])
            U2[i, j] = np.vdot(eigvecs_mesh[i, j], eigvecs_mesh[i, jp1])

    # Normalize link variables (lattice gauge fixing)
    U1 /= np.abs(U1)
    U2 /= np.abs(U2)

    # Berry curvature on each plaquette: F = Im[ln(U1 * U2_shifted * U1_shifted^* * U2^*)]
    berry_curvature = np.zeros((Nkx, Nky))
    for i in range(Nkx):
        for j in range(Nky):
            ip1 = (i + 1) % Nkx
            jp1 = (j + 1) % Nky
            # Plaquette product (counterclockwise)
            plaq = U1[i, j] * U2[ip1, j] * np.conj(U1[i, jp1]) * np.conj(U2[i, j])
            berry_curvature[i, j] = np.imag(np.log(plaq))

    chern = np.sum(berry_curvature) / (2 * np.pi)

    return chern, berry_curvature


# === Main Class ===


class BistritzerMacDonald:
    """Bistritzer-MacDonald continuum model for twisted multilayer graphene.

    Constructs and diagonalizes the moire Hamiltonian for N layers of graphene
    with alternating twist angles +/-theta/2 (relative twist theta between
    adjacent layers).

    The Hamiltonian is block-tridiagonal:
        - Diagonal blocks: Dirac Hamiltonians h_l(k+G) for each layer l
        - Off-diagonal blocks: Interlayer tunneling T connecting adjacent layers

    Attributes
    ----------
    theta_deg : float
        Twist angle in degrees.
    n_layers : int
        Number of graphene layers (2 for TBG, 5 for TEVC).
    hbar_vf : float
        hbar * v_F in eV*Angstrom (Dirac cone velocity scale).
    w_AA : float
        AA-site interlayer tunneling (eV). Typically ~0.080 eV.
    w_AB : float
        AB-site interlayer tunneling (eV). Typically ~0.110 eV.
    K_D : float
        Dirac point distance from Gamma in monolayer BZ (1/Angstrom).
    K_theta : float
        Moire wavevector magnitude = 2*K_D*sin(theta/2) (1/Angstrom).
    G1, G2 : ndarray
        Moire reciprocal lattice vectors (1/Angstrom).
    q_vectors : ndarray, shape (3, 2)
        Three tunneling q-vectors connecting the Dirac cones.
    """

    def __init__(
        self,
        theta_deg: float,
        n_layers: int = 2,
        v_F: float = 1e6,
        w_AA: float = 0.080,
        w_AB: float = 0.110,
        sublattice_mass: float = 0.0,
    ):
        """Initialize BM model parameters.

        Parameters
        ----------
        theta_deg : float
            Twist angle between adjacent layers (degrees).
        n_layers : int
            Number of graphene layers. Default 2 (bilayer).
            For TEVC platform: n_layers=5.
        v_F : float
            Fermi velocity in m/s. Default 1e6 m/s.
            Note: hbar*v_F = 5.253 eV*Angstrom for v_F ~ 0.8e6 m/s
            (renormalized graphene value from literature).
        w_AA : float
            AA-tunneling amplitude in eV. Default 0.080 eV.
        w_AB : float
            AB-tunneling amplitude in eV. Default 0.110 eV.
            The ratio w_AA/w_AB controls particle-hole asymmetry.
            w_AA = 0 is the chiral limit (exact flat bands at magic angle).
        sublattice_mass : float
            Sublattice potential m*sigma_z (eV). Default 0.
            Models aligned hBN substrate that breaks C2 symmetry and gaps
            the flat-band Dirac points, enabling well-defined Chern numbers.
            Typical hBN-induced gap: ~15-30 meV (0.015-0.030 eV).
        """
        self.theta_deg = theta_deg
        self.theta = np.radians(theta_deg)
        self.n_layers = n_layers
        self.v_F = v_F
        self.w_AA = w_AA
        self.w_AB = w_AB
        self.sublattice_mass = sublattice_mass

        # hbar * v_F in eV*Angstrom
        # hbar(eV*s) * v_F(m/s) * 1e10(Angstrom/m)
        self.hbar_vf = HBAR_EV_S * v_F * 1e10  # eV*Angstrom

        # Graphene K-point magnitude: K_D = 4*pi / (3*a)
        self.K_D = 4 * np.pi / (3 * A_GRAPHENE)  # 1/Angstrom

        # Moire wavevector: K_theta = 2*K_D*sin(theta/2)
        # This equals the mBZ corner-to-center distance
        self.K_theta = 2 * self.K_D * np.sin(self.theta / 2)  # 1/Angstrom

        # Moire period: lambda_M = a / (2*sin(theta/2))
        self.lambda_M = A_GRAPHENE / (2 * np.sin(self.theta / 2))  # Angstrom

        # Three q-vectors: tunneling momentum transfers
        # Convention: q1 points along -y, q2 and q3 are 120-degree rotations
        # q1 + q2 + q3 = 0 (triangle in reciprocal space)
        self.q_vectors = self.K_theta * np.array([
            [0.0, -1.0],                         # q1: along -y
            [np.sqrt(3) / 2, 0.5],               # q2: 120 deg from q1
            [-np.sqrt(3) / 2, 0.5],              # q3: 240 deg from q1
        ])

        # Moire reciprocal lattice vectors: G1 = q2 - q1, G2 = q3 - q1
        # |G1| = |G2| = sqrt(3) * K_theta, angle between them = 60 degrees
        self.G1 = self.q_vectors[1] - self.q_vectors[0]
        self.G2 = self.q_vectors[2] - self.q_vectors[0]

        # Dimensionless coupling (alpha parameter from TKV):
        # alpha = w_AB / (hbar*v_F * K_theta)
        # Magic angle at alpha ~ 0.586 (chiral limit)
        self.alpha = self.w_AB / (self.hbar_vf * self.K_theta)

        # Build tunneling matrices T_j (2x2 in sublattice space)
        self._build_tunneling_matrices()

    def _build_tunneling_matrices(self):
        """Construct the three 2x2 tunneling matrices T_1, T_2, T_3.

        T_j = w_AA * sigma_0 + w_AB * (cos(phi_j) sigma_x + sin(phi_j) sigma_y)
        where phi_j = 2*pi*(j-1)/3

        In sublattice basis (A, B):
            T_j = [[w_AA, w_AB * exp(-i*phi_j)],
                   [w_AB * exp(+i*phi_j), w_AA]]

        Each T_j is Hermitian (since w_AA, w_AB are real).
        """
        self.T_matrices = np.zeros((3, 2, 2), dtype=complex)

        for j in range(3):
            phi = 2 * np.pi * j / 3  # phi_1=0, phi_2=2pi/3, phi_3=4pi/3
            self.T_matrices[j] = np.array([
                [self.w_AA, self.w_AB * np.exp(-1j * phi)],
                [self.w_AB * np.exp(+1j * phi), self.w_AA],
            ])

    def _dirac_hamiltonian(self, p: np.ndarray, layer_index: int) -> np.ndarray:
        """Single-layer Dirac Hamiltonian at momentum p.

        h(p) = hbar*v_F * sigma . R(alpha_l) * p

        where alpha_l = +theta/2 for odd layers (1,3,5,...) and
              alpha_l = -theta/2 for even layers (2,4,...).

        The rotation R(alpha_l) accounts for the sublattice orientation
        of the rotated graphene layer.

        Parameters
        ----------
        p : ndarray, shape (2,)
            Momentum vector (kx, ky) in 1/Angstrom.
            For odd layers: p = k + G (measured from mBZ Gamma).
            For even layers: p = k + G + q1 (offset by the moire q-vector).
        layer_index : int
            Layer number (0-indexed): 0, 1, 2, ..., N-1.

        Returns
        -------
        h : ndarray, shape (2, 2), complex
            Dirac Hamiltonian in sublattice space.
            Eigenvalues: +/- hbar*v_F*|R(alpha)*p|.
        """
        # Rotation angle: odd layers (0,2,4,...) get +theta/2, even (1,3,...) get -theta/2
        if layer_index % 2 == 0:
            angle = +self.theta / 2
        else:
            angle = -self.theta / 2

        # Rotate momentum by the layer angle
        R = _rotation_matrix_2d(angle)
        p_rot = R @ p  # rotated momentum components (px', py')

        # Dirac Hamiltonian: h = hbar*v_F * (px' sigma_x + py' sigma_y) + m * sigma_z
        # In matrix form: [[m, px'-i*py'], [px'+i*py', -m]]
        # The mass term (sublattice potential) opens a gap at the Dirac point.
        px, py = p_rot[0], p_rot[1]
        m = self.sublattice_mass
        h = self.hbar_vf * np.array([
            [0.0, px - 1j * py],
            [px + 1j * py, 0.0],
        ])
        # Add sublattice mass (sigma_z term): +m on A, -m on B
        h[0, 0] += m
        h[1, 1] -= m

        return h

    def hamiltonian(self, k: np.ndarray, n_shells: int = 4) -> np.ndarray:
        """Build the full BM Hamiltonian matrix at moire momentum k.

        The Hamiltonian is block-tridiagonal for N layers:
            H = diag(H_1, H_2, ..., H_N) + off-diag(T_12, T_23, ..., T_{N-1,N})

        Each layer block H_l is block-diagonal in G (2x2 Dirac at each G-point).
        Tunneling blocks T_{l,l+1} connect G in layer l to G' in layer l+1
        when G - G' = 0, G1, or G2 (the three tunneling channels).

        Parameters
        ----------
        k : ndarray, shape (2,)
            Crystal momentum in the moire BZ (1/Angstrom).
        n_shells : int
            Number of moire reciprocal lattice shells in the plane-wave expansion.
            Larger = more converged but slower.
            n_shells=4 typically gives well-converged flat bands.

        Returns
        -------
        H : ndarray, shape (dim, dim), complex
            Hermitian Hamiltonian matrix.
            dim = 2 * N_G * n_layers (sublattice * plane waves * layers).
        """
        k = np.asarray(k, dtype=float)

        # Generate G-vectors within cutoff
        G_list = _generate_g_vectors(self.G1, self.G2, n_shells)
        N_G = len(G_list)

        # Total dimension: 2 (sublattice) * N_G (plane waves) * n_layers
        dim = 2 * N_G * self.n_layers
        H = np.zeros((dim, dim), dtype=complex)

        # --- Intralayer blocks (diagonal in layer, diagonal in G) ---
        for layer in range(self.n_layers):
            for ig, G in enumerate(G_list):
                # Momentum entering the Dirac Hamiltonian:
                # Odd layers (0,2,4): p = k + G (cone at Gamma of mBZ)
                # Even layers (1,3): p = k + G + q1 (cone shifted to K-point of mBZ)
                if layer % 2 == 0:
                    p = k + G
                else:
                    p = k + G + self.q_vectors[0]  # offset by q1

                h = self._dirac_hamiltonian(p, layer)

                # Place 2x2 block in the full matrix
                # Index: layer * (2*N_G) + ig * 2 : layer * (2*N_G) + ig * 2 + 2
                row = layer * (2 * N_G) + ig * 2
                H[row:row + 2, row:row + 2] = h

        # --- Interlayer tunneling blocks (off-diagonal in layer) ---
        # Build a lookup: G-vector index from (approximate) G-vector value
        # For fast matching of G - G' = delta conditions
        g_index_map = {}
        for ig, G in enumerate(G_list):
            # Round to avoid floating-point matching issues
            key = (round(G[0] * 1e6), round(G[1] * 1e6))
            g_index_map[key] = ig

        # The three delta vectors for the tunneling connection rule.
        # From momentum conservation: G_layer2 = G_layer1 + (q_j - q_1)
        # Rearranging: G_layer1 - G_layer2 = -(q_j - q_1) = q_1 - q_j
        # delta_1 = 0, delta_2 = q1 - q2 = -G1, delta_3 = q1 - q3 = -G2
        deltas = [np.array([0.0, 0.0]), -self.G1, -self.G2]

        for layer in range(self.n_layers - 1):
            # Tunneling from layer to layer+1
            layer_top = layer      # "layer 1" in the bilayer pair
            layer_bot = layer + 1  # "layer 2" in the bilayer pair

            for ig_top, G_top in enumerate(G_list):
                for j, delta in enumerate(deltas):
                    # Connection: G_top - G_bot = delta_j
                    # So G_bot = G_top - delta_j
                    G_bot_target = G_top - delta
                    key = (round(G_bot_target[0] * 1e6), round(G_bot_target[1] * 1e6))

                    if key in g_index_map:
                        ig_bot = g_index_map[key]

                        # T_j matrix element: connects top layer at G_top to bot at G_bot
                        T_j = self.T_matrices[j]

                        # Position in full Hamiltonian:
                        # Top layer (layer_top): row = layer_top * 2*N_G + ig_top * 2
                        # Bot layer (layer_bot): col = layer_bot * 2*N_G + ig_bot * 2
                        row = layer_top * (2 * N_G) + ig_top * 2
                        col = layer_bot * (2 * N_G) + ig_bot * 2

                        # Upper-right block: H[top, bot] += T_j^dagger = T_j (Hermitian)
                        H[row:row + 2, col:col + 2] += T_j

                        # Lower-left block: H[bot, top] += T_j (Hermitian conjugate)
                        H[col:col + 2, row:row + 2] += T_j.conj().T

        # Verify Hermiticity (sanity check in debug)
        assert np.allclose(H, H.conj().T, atol=1e-12), "Hamiltonian is not Hermitian!"

        return H

    def band_structure(
        self,
        k_path: np.ndarray,
        n_bands: int = 10,
        n_shells: int = 4,
    ) -> np.ndarray:
        """Compute eigenvalues along a k-path through the moire BZ.

        Parameters
        ----------
        k_path : ndarray, shape (N_k, 2)
            Array of k-points along the desired path (1/Angstrom).
        n_bands : int
            Number of bands to return (centered around E=0).
            Returns the n_bands/2 bands above and below the charge neutrality point.
        n_shells : int
            Plane-wave truncation (see hamiltonian method).

        Returns
        -------
        energies : ndarray, shape (N_k, n_bands)
            Eigenvalues in eV along the path, sorted by energy.
            Centered around E=0 (charge neutrality).
        """
        k_path = np.asarray(k_path)
        N_k = len(k_path)

        # Compute all eigenvalues at first k-point to determine total dimension
        H0 = self.hamiltonian(k_path[0], n_shells=n_shells)
        dim = H0.shape[0]

        # Index range: take n_bands centered around the middle
        mid = dim // 2
        idx_lo = mid - n_bands // 2
        idx_hi = idx_lo + n_bands

        energies = np.zeros((N_k, n_bands))

        for ik, k in enumerate(k_path):
            if ik == 0:
                H = H0
            else:
                H = self.hamiltonian(k, n_shells=n_shells)

            # Hermitian eigenvalue decomposition (returns sorted eigenvalues)
            evals = np.linalg.eigh(H)[0]
            energies[ik] = evals[idx_lo:idx_hi]

        return energies

    def berry_curvature(
        self,
        k_mesh_size: int = 50,
        band_index: int = 0,
        n_shells: int = 4,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Compute Berry curvature on a k-mesh via Fukui lattice gauge method.

        Discretizes the moire BZ on a regular mesh and computes the
        gauge-invariant Berry curvature for a specified band.

        Parameters
        ----------
        k_mesh_size : int
            Number of mesh points along each direction (total = k_mesh_size^2).
        band_index : int
            Band index relative to charge neutrality:
            0 = highest valence band, -1 = next below, +1 = lowest conduction, etc.
        n_shells : int
            Plane-wave truncation.

        Returns
        -------
        kx_mesh : ndarray, shape (k_mesh_size, k_mesh_size)
            kx coordinates of mesh points.
        ky_mesh : ndarray, shape (k_mesh_size, k_mesh_size)
            ky coordinates of mesh points.
        omega : ndarray, shape (k_mesh_size, k_mesh_size)
            Berry curvature Omega(k) at each plaquette.
        """
        # Span the moire BZ using reciprocal lattice vectors
        # The BZ is spanned by fractional coordinates (f1, f2) in [0, 1)
        # k = f1*G1 + f2*G2, with the BZ centered at Gamma
        # Use offset to center: f in [-0.5, 0.5)
        N = k_mesh_size
        f_vals = np.linspace(-0.5, 0.5, N, endpoint=False)

        # Determine band index in the full spectrum
        H_test = self.hamiltonian(np.array([0.0, 0.0]), n_shells=n_shells)
        dim = H_test.shape[0]
        mid = dim // 2

        # band_index: 0 = highest valence, 1 = lowest conduction, -1 = next valence
        if band_index >= 0:
            target_band = mid - 1 - band_index  # valence bands count down from mid-1
        else:
            target_band = mid + (-band_index - 1)  # conduction bands count up from mid

        # Compute eigenvectors on the full mesh
        eigvecs_mesh = np.zeros((N, N, dim), dtype=complex)

        for i, f1 in enumerate(f_vals):
            for j, f2 in enumerate(f_vals):
                k = f1 * self.G1 + f2 * self.G2
                H = self.hamiltonian(k, n_shells=n_shells)
                evals, evecs = np.linalg.eigh(H)
                eigvecs_mesh[i, j] = evecs[:, target_band]

        # Apply Fukui method
        chern, omega = fukui_chern(eigvecs_mesh)

        # Build coordinate meshes for plotting
        kx_mesh = np.zeros((N, N))
        ky_mesh = np.zeros((N, N))
        for i, f1 in enumerate(f_vals):
            for j, f2 in enumerate(f_vals):
                k = f1 * self.G1 + f2 * self.G2
                kx_mesh[i, j] = k[0]
                ky_mesh[i, j] = k[1]

        return kx_mesh, ky_mesh, omega

    def chern_number(
        self,
        k_mesh_size: int = 50,
        band_index: int = 0,
        n_shells: int = 4,
    ) -> float:
        """Compute the Chern number for a specified band.

        Uses the Fukui lattice gauge method on a discretized BZ mesh.
        For MATBG flat bands: expected C = +/-1 per valley per spin.

        Parameters
        ----------
        k_mesh_size : int
            Mesh density. Larger = more accurate but slower.
            k_mesh_size >= 30 typically sufficient for convergence.
        band_index : int
            Band index (0 = highest valence band).
        n_shells : int
            Plane-wave truncation.

        Returns
        -------
        C : float
            Chern number (should be close to integer for gapped bands).
        """
        _, _, omega = self.berry_curvature(k_mesh_size, band_index, n_shells)
        C = np.sum(omega) / (2 * np.pi)
        return C

    def flat_band_width(self, n_shells: int = 4, n_k: int = 50) -> float:
        """Compute the bandwidth W of the flat bands near charge neutrality.

        Samples the moire BZ on a regular mesh and finds the bandwidth
        (max - min energy) of the bands closest to E=0.

        For MATBG at the magic angle: W < 10 meV.
        Away from magic angle: W increases approximately as |alpha - alpha_magic|.

        Parameters
        ----------
        n_shells : int
            Plane-wave truncation.
        n_k : int
            Number of k-points along each BZ direction (total mesh = n_k^2).

        Returns
        -------
        W : float
            Flat-band bandwidth in eV.
        """
        # Sample BZ on regular mesh
        f_vals = np.linspace(-0.5, 0.5, n_k, endpoint=False)

        # Get dimension and identify flat-band indices
        H_test = self.hamiltonian(np.array([0.0, 0.0]), n_shells=n_shells)
        dim = H_test.shape[0]
        mid = dim // 2

        # For bilayer: 2 flat bands (indices mid-1 and mid in sorted eigenvalues)
        # For N-layer: approximately N flat bands centered around mid
        # The flat bands are the N bands closest to E=0 in the moire spectrum
        n_flat = self.n_layers

        idx_lo = mid - n_flat // 2
        idx_hi = mid + (n_flat + 1) // 2

        # Collect flat-band energies across BZ
        flat_band_evals = []

        for f1 in f_vals:
            for f2 in f_vals:
                k = f1 * self.G1 + f2 * self.G2
                H = self.hamiltonian(k, n_shells=n_shells)
                evals = np.linalg.eigh(H)[0]
                flat_band_evals.append(evals[idx_lo:idx_hi])

        flat_band_evals = np.array(flat_band_evals)

        # Bandwidth of the narrowest band(s) closest to E=0
        # Find the band with smallest bandwidth
        n_bands_collected = idx_hi - idx_lo
        bandwidths = np.zeros(n_bands_collected)
        for ib in range(n_bands_collected):
            bandwidths[ib] = np.max(flat_band_evals[:, ib]) - np.min(flat_band_evals[:, ib])

        # Return the minimum bandwidth (the flattest band)
        # This is the relevant scale for correlation physics
        return float(np.min(bandwidths))

    def density_of_states(
        self,
        energies: np.ndarray,
        n_k: int = 100,
        broadening: float = 0.001,
        n_shells: int = 4,
    ) -> np.ndarray:
        """Compute density of states via Lorentzian broadening on a k-mesh.

        DOS(E) = (1/N_k) * sum_{k,n} delta(E - E_n(k))
               ~ (1/N_k) * sum_{k,n} (eta/pi) / ((E - E_n(k))^2 + eta^2)

        Parameters
        ----------
        energies : ndarray, shape (N_E,)
            Energy values at which to evaluate the DOS (eV).
        n_k : int
            Number of k-points per BZ direction (total mesh = n_k^2).
        broadening : float
            Lorentzian half-width eta in eV. Typical: 0.001 eV (1 meV).
        n_shells : int
            Plane-wave truncation.

        Returns
        -------
        dos : ndarray, shape (N_E,)
            Density of states at each energy value (states/eV/cell).
        """
        energies = np.asarray(energies)
        N_E = len(energies)
        dos = np.zeros(N_E)

        f_vals = np.linspace(-0.5, 0.5, n_k, endpoint=False)
        N_k_total = n_k * n_k
        eta = broadening

        for f1 in f_vals:
            for f2 in f_vals:
                k = f1 * self.G1 + f2 * self.G2
                H = self.hamiltonian(k, n_shells=n_shells)
                evals = np.linalg.eigh(H)[0]

                # Add Lorentzian contribution from each eigenvalue
                for E_n in evals:
                    dos += (eta / np.pi) / ((energies - E_n) ** 2 + eta ** 2)

        dos /= N_k_total

        return dos

    def get_moire_bz_info(self) -> dict:
        """Return key geometric information about the moire Brillouin zone.

        Useful for setting up calculations and verifying parameters.

        Returns
        -------
        info : dict
            Dictionary with moire geometry parameters.
        """
        return {
            "theta_deg": self.theta_deg,
            "n_layers": self.n_layers,
            "lambda_M_Angstrom": self.lambda_M,
            "lambda_M_nm": self.lambda_M / 10,
            "K_theta_inv_Angstrom": self.K_theta,
            "G_magnitude_inv_Angstrom": np.linalg.norm(self.G1),
            "hbar_vf_eV_Angstrom": self.hbar_vf,
            "alpha_BM": self.alpha,
            "alpha_magic_chiral": 0.586,  # TKV value for w_AA=0
            "K_point_mBZ": -self.q_vectors[0],  # K-point of moire BZ
        }


# === Convenience constructor ===


def matbg_model(theta_deg: float = 1.1, **kwargs) -> BistritzerMacDonald:
    """Create a BM model for magic-angle twisted bilayer graphene.

    Uses literature parameters:
        hbar*v_F = 5.253 eV*Angstrom (v_F = 7.98e5 m/s)
        w_AA = 0.0797 eV
        w_AB = 0.1100 eV

    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees. Default: 1.1 (first magic angle).
    **kwargs
        Additional keyword arguments passed to BistritzerMacDonald.

    Returns
    -------
    model : BistritzerMacDonald
        Configured BM model instance.
    """
    # v_F chosen to reproduce hbar*v_F = 5.253 eV*Angstrom
    v_F = 5.253 / (HBAR_EV_S * 1e10)  # m/s (approximately 7.98e5)
    defaults = dict(v_F=v_F, w_AA=0.0797, w_AB=0.1100)
    defaults.update(kwargs)
    return BistritzerMacDonald(theta_deg, **defaults)


def tevc_model(theta_deg: float = 1.1, n_layers: int = 5, **kwargs) -> BistritzerMacDonald:
    """Create a BM model for the TEVC platform (N=5 layer stack).

    Parameters
    ----------
    theta_deg : float
        Twist angle between adjacent layers (degrees).
    n_layers : int
        Number of layers. Default: 5 (TEVC).
    **kwargs
        Additional keyword arguments passed to BistritzerMacDonald.

    Returns
    -------
    model : BistritzerMacDonald
        Configured BM model for TEVC.
    """
    v_F = 5.253 / (HBAR_EV_S * 1e10)
    defaults = dict(v_F=v_F, w_AA=0.0797, w_AB=0.1100, n_layers=n_layers)
    defaults.update(kwargs)
    return BistritzerMacDonald(theta_deg, **defaults)


# === Validation ===


if __name__ == "__main__":
    print("=" * 70)
    print("Bistritzer-MacDonald Model: Validation at theta = 1.1 deg")
    print("=" * 70)

    # Create model with literature parameters
    model = matbg_model(theta_deg=1.1)
    info = model.get_moire_bz_info()

    print(f"\n--- Moire Geometry ---")
    print(f"  Twist angle:        {info['theta_deg']:.2f} deg")
    print(f"  Moire period:       {info['lambda_M_nm']:.2f} nm")
    print(f"  K_theta:            {info['K_theta_inv_Angstrom']:.5f} 1/Angstrom")
    print(f"  |G_moire|:          {info['G_magnitude_inv_Angstrom']:.5f} 1/Angstrom")
    print(f"  hbar*v_F:           {info['hbar_vf_eV_Angstrom']:.3f} eV*Angstrom")
    print(f"  alpha (BM):         {info['alpha_BM']:.4f}")
    print(f"  alpha (magic, chiral): {info['alpha_magic_chiral']:.4f}")

    # Check Hamiltonian dimension and Hermiticity
    print(f"\n--- Hamiltonian Structure (n_shells=4) ---")
    H_test = model.hamiltonian(np.array([0.0, 0.0]), n_shells=4)
    print(f"  Matrix dimension:   {H_test.shape[0]} x {H_test.shape[1]}")
    print(f"  Hermitian check:    {np.allclose(H_test, H_test.conj().T)}")

    G_list = _generate_g_vectors(model.G1, model.G2, n_shells=4)
    print(f"  Plane waves (N_G):  {len(G_list)}")
    print(f"  States per layer:   {2 * len(G_list)}")

    # Compute flat-band width
    print(f"\n--- Flat-Band Width ---")
    print(f"  Computing on 30x30 k-mesh with n_shells=3...")
    W = model.flat_band_width(n_shells=3, n_k=30)
    print(f"  Bandwidth W:        {W * 1000:.2f} meV")
    print(f"  (Expected: W ~ 10 meV near magic angle)")
    print(f"  (alpha = {model.alpha:.4f}, magic = 0.586 in chiral limit)")

    # Verify exact flat bands in chiral limit at magic angle
    print(f"\n--- Chiral Limit Validation (alpha = 0.586) ---")
    K_D_val = 4 * np.pi / (3 * A_GRAPHENE)
    K_theta_magic = model.w_AB / (model.hbar_vf * 0.586)
    theta_magic = np.degrees(2 * np.arcsin(K_theta_magic / (2 * K_D_val)))
    v_F_lit = 5.253 / (HBAR_EV_S * 1e10)
    model_chiral = BistritzerMacDonald(theta_magic, v_F=v_F_lit, w_AA=0.0, w_AB=0.110)
    W_chiral = model_chiral.flat_band_width(n_shells=4, n_k=20)
    print(f"  Magic angle theta:  {theta_magic:.3f} deg")
    print(f"  W (chiral, magic):  {W_chiral * 1000:.4f} meV  (should be ~0)")

    # Band structure along high-symmetry path
    print(f"\n--- Band Structure (Gamma -> K -> M -> Gamma) ---")
    k_path, k_dist, ticks = k_path_high_symmetry(model.K_theta, n_points=60)
    bands = model.band_structure(k_path, n_bands=6, n_shells=3)
    print(f"  Energy range of flat bands: [{bands[:, 2].min()*1000:.2f}, "
          f"{bands[:, 3].max()*1000:.2f}] meV")
    print(f"  Gap above flat bands:       {(bands[:, 4].min() - bands[:, 3].max())*1000:.2f} meV")
    print(f"  Gap below flat bands:       {(bands[:, 2].min() - bands[:, 1].max())*1000:.2f} meV")

    # Chern number with C2-breaking sublattice mass (hBN substrate)
    # Without the mass, flat bands touch at K and Gamma (C2-protected),
    # making individual Chern numbers ill-defined. With a sublattice mass
    # (modeling aligned hBN), the gap opens and C = +/-1 per band.
    print(f"\n--- Topology (Chern number with hBN mass = 0.015 eV) ---")
    model_hbn = matbg_model(theta_deg=1.1, sublattice_mass=0.015)
    # Verify the gap opens
    H_K = model_hbn.hamiltonian(-model_hbn.q_vectors[0], n_shells=3)
    evals_K = np.linalg.eigh(H_K)[0]
    mid_K = len(evals_K) // 2
    print(f"  Gap at K-point:     {(evals_K[mid_K] - evals_K[mid_K-1])*1000:.2f} meV")
    C = model_hbn.chern_number(k_mesh_size=25, band_index=0, n_shells=3)
    print(f"  Chern number C:     {C:.3f}")
    print(f"  (Expected: |C| = 1 per valley for gapped flat band)")

    print(f"\n{'=' * 70}")
    print(f"Validation complete.")
    print(f"{'=' * 70}")
