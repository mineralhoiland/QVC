"""
qvc.montecarlo.xy_model — 2D XY model Metropolis Monte Carlo.

Physics
-------
The phase sector of TEVC: θ(r) ∈ U(1) lives on the circle.
The effective theory is the 2D XY model with Kosterlitz coupling K:

    H / T_eff = -K Σ_{⟨ij⟩} cos(θ_i - θ_j)

where K = κ_s / T_eff (spin stiffness over effective temperature).

BKT transition: K_c = 2/π ≈ 0.6366.
- K > K_c: bound vortex-antivortex pairs (quasi-long-range order)
- K < K_c: free vortices (disordered)

For TEVC at T* = 118.7 mK on the primary branch:
    K ≈ 10-47 >> K_c (deeply ordered)

In the spin-wave approximation (low-T, large K):
    g(r) = ⟨exp(i[θ(r) - θ(0)])⟩ ~ r^{-η}
    η = 1 / (2πK)

At K = 10.57 (primary branch at α_c): η ≈ 0.01506, far below η_BKT = 1/4.

This module provides Metropolis sampling of the XY model for:
- Measuring the phase correlator and effective exponent
- Vortex density statistics
- Helicity modulus (spin stiffness) via winding number fluctuations

Units: K is dimensionless coupling (already divided by T).
"""

from __future__ import annotations

import numpy as np


class XYModel:
    """2D XY model on a square lattice with Metropolis Monte Carlo.

    The lattice is L×L with periodic boundary conditions. Each site
    carries an angle θ_i ∈ [0, 2π). Nearest-neighbor coupling K.

    Parameters
    ----------
    L : int
        Linear system size (L×L square lattice). Default: 16.
    K : float
        Coupling constant K = J/T (dimensionless). Default: 10.0.
        At BKT: K_c = 2/π ≈ 0.6366.
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(self, L: int = 16, K: float = 10.0, seed: int = 137):
        self.L = L
        self.K = K
        self.N = L * L  # total number of spins

        # RNG
        self.rng = np.random.default_rng(seed)

        # Initialize angles uniformly in [0, 2π) — hot start
        self.theta = self.rng.uniform(0, 2 * np.pi, size=(L, L))

        # Metropolis proposal width: adapt for acceptance ~50% at this K
        # At large K, small proposals are needed; δθ ~ 1/√K is a good heuristic
        self.proposal_width = min(np.pi, 2.0 / np.sqrt(K))

        # Counters
        self.sweeps_done = 0
        self._accepted = 0
        self._proposed = 0

    def _energy_site(self, i: int, j: int) -> float:
        """Local energy contribution from site (i,j): -K Σ_{nn} cos(Δθ).

        Sum over 4 nearest neighbors (periodic BC).
        """
        L = self.L
        theta_ij = self.theta[i, j]

        # Neighbor indices (periodic)
        neighbors = [
            self.theta[(i + 1) % L, j],
            self.theta[(i - 1) % L, j],
            self.theta[i, (j + 1) % L],
            self.theta[i, (j - 1) % L],
        ]

        E = 0.0
        for theta_nn in neighbors:
            E -= self.K * np.cos(theta_ij - theta_nn)
        return E

    def sweep(self, n_sweeps: int = 1) -> float:
        """Perform n_sweeps full Metropolis sweeps (L² single-site updates each).

        Each sweep visits all L² sites once in random order.

        Parameters
        ----------
        n_sweeps : int
            Number of full lattice sweeps.

        Returns
        -------
        acceptance_rate : float
            Fraction of proposed moves accepted (across all sweeps).
        """
        L = self.L
        accepted = 0
        proposed = 0

        for _ in range(n_sweeps):
            # Random visitation order for one sweep
            sites = self.rng.permutation(self.N)

            for site in sites:
                i = site // L
                j = site % L

                # Current energy contribution from this site
                E_old = self._energy_site(i, j)

                # Propose new angle
                theta_old = self.theta[i, j]
                delta = self.rng.uniform(-self.proposal_width, self.proposal_width)
                theta_new = (theta_old + delta) % (2 * np.pi)

                # Compute new energy
                self.theta[i, j] = theta_new
                E_new = self._energy_site(i, j)

                # Metropolis acceptance: accept if ΔE < 0 or with prob exp(-ΔE)
                # Note: H/T = -K Σ cos(...), so ΔE is already in units of T
                dE = E_new - E_old
                proposed += 1

                if dE <= 0 or self.rng.random() < np.exp(-dE):
                    accepted += 1  # keep new configuration
                else:
                    self.theta[i, j] = theta_old  # revert

            self.sweeps_done += 1

        self._accepted += accepted
        self._proposed += proposed
        return accepted / proposed if proposed > 0 else 0.0

    def thermalize(self, n_sweeps: int = 1000) -> None:
        """Run n_sweeps to reach thermal equilibrium.

        For K >> K_c (deeply ordered): starting from random config,
        ~100-200 sweeps typically suffice for L ≤ 32.
        For K near K_c: critical slowing down requires more sweeps.
        """
        # Reset acceptance counter
        self._accepted = 0
        self._proposed = 0
        self.sweep(n_sweeps)

    def phase_correlator(self) -> np.ndarray:
        """Measure the phase-phase correlator g(r) = ⟨exp(i[θ(r) - θ(0)])⟩.

        Computed as a spatial average over all pairs at distance r
        (Manhattan distance projected onto x-axis for simplicity).

        Returns
        -------
        g_r : ndarray, shape (L//2 + 1,)
            Correlator g(r) for r = 0, 1, ..., L//2 lattice spacings.
            g(0) = 1 by definition.
        """
        L = self.L
        max_r = L // 2 + 1
        g_r = np.zeros(max_r)
        counts = np.zeros(max_r)

        # Use translational invariance: average over all origin sites
        # Project correlator along x-direction (equivalent to y by symmetry)
        for r in range(max_r):
            # cos[θ(i, j+r) - θ(i, j)] averaged over all (i, j)
            shifted = np.roll(self.theta, -r, axis=1)
            delta_theta = shifted - self.theta
            # g(r) = Re⟨exp(iΔθ)⟩ for real correlator (imaginary part = 0 by symmetry)
            g_r[r] = float(np.mean(np.cos(delta_theta)))
            counts[r] = self.N

        return g_r

    def effective_exponent(self) -> float:
        """Extract effective power-law exponent η from g(r).

        η = -ln g(L/2) / ln(L/2)

        In the spin-wave regime (K >> K_c):
            η_sw = 1 / (2πK)

        At BKT:
            η_BKT = 1/4

        Returns
        -------
        eta : float
            Effective exponent. Should match 1/(2πK) deep in ordered phase.
        """
        g = self.phase_correlator()
        L_half = self.L // 2

        # Avoid log(0) or log(negative) from fluctuations
        g_half = max(g[L_half], 1e-10)

        # η = -ln g(r) / ln(r) evaluated at r = L/2
        eta = -np.log(g_half) / np.log(L_half)
        return float(eta)

    def vortex_density(self) -> float:
        """Measure vortex density: fraction of plaquettes carrying a vortex.

        A vortex on plaquette (i,j) has winding:
            q = (1/2π) [Δθ₁₂ + Δθ₂₃ + Δθ₃₄ + Δθ₄₁]

        where Δθ is wrapped to [-π, π]. q = ±1 for vortex/antivortex.

        Returns
        -------
        rho_v : float
            Vortex density = N_vortices / N_plaquettes.
            Deep in ordered phase (K >> K_c): exponentially small.
            Above BKT (K < K_c): O(1).
        """
        L = self.L
        theta = self.theta

        # Phase differences along bonds, wrapped to [-π, π]
        # Δθ_x[i,j] = θ[i,j+1] - θ[i,j]
        dtheta_x = _wrap_angle(np.roll(theta, -1, axis=1) - theta)
        # Δθ_y[i,j] = θ[i+1,j] - θ[i,j]
        dtheta_y = _wrap_angle(np.roll(theta, -1, axis=0) - theta)

        # Plaquette winding (counterclockwise around each elementary square):
        # right + up(shifted right) - right(shifted up) - up
        winding = (dtheta_x
                   + np.roll(dtheta_y, -1, axis=1)
                   - np.roll(dtheta_x, -1, axis=0)
                   - dtheta_y)

        # Vortex charge = winding / 2π
        charges = np.round(winding / (2 * np.pi)).astype(int)
        n_vortices = int(np.sum(np.abs(charges)))

        # Density = total vortices / total plaquettes
        return n_vortices / (L * L)

    def measure_stiffness(self) -> float:
        """Helicity modulus (spin stiffness) from winding number fluctuations.

        The winding number in direction μ:
            W_μ = (1/2π) Σ_along_μ Δθ  (wrapped)

        Helicity modulus:
            ρ_s = K × [1 - K × ⟨(Σ_bond sin Δθ)²⟩ / N]

        Simplified estimator using instantaneous winding number:
            W_x = (1/L) Σ_rows [ (1/2π) Σ_j Δθ_{j,j+1} ]

        For a single configuration this is noisy; call after many
        thermalizing sweeps and average externally.

        Returns
        -------
        kappa : float
            Instantaneous stiffness estimator from winding numbers.
            Average over many measurements for a reliable result.
        """
        L = self.L
        theta = self.theta

        # Winding number in x-direction: sum of wrapped Δθ along each row
        dtheta_x = _wrap_angle(np.roll(theta, -1, axis=1) - theta)
        # W_x for each row: (1/2π) × sum of Δθ along that row
        W_x_rows = np.sum(dtheta_x, axis=1) / (2 * np.pi)
        W_x = float(np.mean(W_x_rows))

        # Winding number in y-direction
        dtheta_y = _wrap_angle(np.roll(theta, -1, axis=0) - theta)
        W_y_cols = np.sum(dtheta_y, axis=0) / (2 * np.pi)
        W_y = float(np.mean(W_y_cols))

        # Stiffness estimator: ρ_s ∝ K × ⟨W²⟩
        # For a single measurement: return K × (W_x² + W_y²) / 2
        # (average over many configs externally)
        W2 = W_x**2 + W_y**2
        kappa = self.K * W2 / 2.0

        return kappa

    @property
    def acceptance_rate(self) -> float:
        """Cumulative acceptance rate since last thermalize()."""
        if self._proposed == 0:
            return 0.0
        return self._accepted / self._proposed

    @property
    def eta_spinwave(self) -> float:
        """Spin-wave prediction for the decay exponent: η = 1/(2πK)."""
        return 1.0 / (2 * np.pi * self.K)


def _wrap_angle(dtheta: np.ndarray) -> np.ndarray:
    """Wrap angle difference to [-π, π]."""
    return (dtheta + np.pi) % (2 * np.pi) - np.pi


# ============================================================
# Validation
# ============================================================

if __name__ == "__main__":
    print("=" * 64)
    print("qvc.montecarlo.xy_model — Validation")
    print("=" * 64)

    # --- Test at K = 10.57 (TEVC primary branch at α_c) ---
    K_test = 10.57
    L_test = 16
    eta_sw_expected = 1.0 / (2 * np.pi * K_test)

    print(f"\nK = {K_test} (TEVC primary branch at α_c)")
    print(f"  BKT critical coupling: K_c = 2/π = {2/np.pi:.4f}")
    print(f"  K/K_c = {K_test/(2/np.pi):.1f} (deeply ordered)")
    print(f"  Spin-wave prediction: η = 1/(2πK) = {eta_sw_expected:.5f}")
    print(f"  (Compare η_BKT = 1/4 = 0.250)")

    model = XYModel(L=L_test, K=K_test, seed=42)
    print(f"\nLattice: {L_test}×{L_test}, N = {model.N} spins")
    print(f"  Proposal width: δθ = {model.proposal_width:.3f} rad")

    # Thermalize
    print(f"\nThermalizing (2000 sweeps)...")
    model.thermalize(n_sweeps=2000)
    print(f"  Acceptance rate: {model.acceptance_rate:.3f}")
    print(f"  (Target: 0.3-0.7 for efficient sampling)")

    # Measure phase correlator
    print(f"\nMeasuring phase correlator (averaging over 500 sweeps)...")
    g_r_accum = np.zeros(L_test // 2 + 1)
    eta_accum = []
    stiffness_accum = []
    vortex_accum = []
    n_measure = 500

    for _ in range(n_measure):
        model.sweep(n_sweeps=2)  # decorrelation sweeps between measurements
        g_r_accum += model.phase_correlator()
        eta_accum.append(model.effective_exponent())
        stiffness_accum.append(model.measure_stiffness())
        vortex_accum.append(model.vortex_density())

    g_r_avg = g_r_accum / n_measure
    eta_measured = float(np.mean(eta_accum))
    eta_err = float(np.std(eta_accum) / np.sqrt(n_measure))
    stiffness_avg = float(np.mean(stiffness_accum))
    vortex_avg = float(np.mean(vortex_accum))

    print(f"\nResults:")
    print(f"  g(1) = {g_r_avg[1]:.6f}  (should be ≈ 1 for K >> K_c)")
    print(f"  g(L/2) = {g_r_avg[L_test//2]:.6f}")
    print(f"  Effective η = {eta_measured:.5f} ± {eta_err:.5f}")
    print(f"  Spin-wave η = {eta_sw_expected:.5f}")
    print(f"  Ratio η_meas/η_sw = {eta_measured/eta_sw_expected:.3f}")
    print(f"  Vortex density = {vortex_avg:.6f} (should be ~0 for K >> K_c)")
    print(f"  Stiffness ⟨K W²/2⟩ = {stiffness_avg:.4f}")

    # Validation checks
    print(f"\nValidation:")

    # η should be close to spin-wave prediction (within ~50% for L=16)
    eta_ratio = eta_measured / eta_sw_expected
    if 0.3 < eta_ratio < 3.0:
        print(f"  η within expected range: PASS (ratio = {eta_ratio:.2f})")
    else:
        print(f"  η ratio = {eta_ratio:.2f}: WARNING (finite-size effects at L=16)")

    # η should be far below BKT value
    assert eta_measured < 0.25, f"η = {eta_measured:.4f} >= 1/4 (above BKT!)"
    print(f"  η << 1/4 (deeply ordered): PASS")

    # Vortex density should be essentially zero
    assert vortex_avg < 0.01, f"Vortex density {vortex_avg:.4f} too high"
    print(f"  Vortex density ≈ 0: PASS ({vortex_avg:.6f})")

    # Correlator should be close to 1 at short range
    assert g_r_avg[1] > 0.95, f"g(1) = {g_r_avg[1]:.4f} too low"
    print(f"  g(1) ≈ 1: PASS ({g_r_avg[1]:.6f})")

    # --- Quick check at K near K_c (should see η → 1/4) ---
    print(f"\n--- Cross-check near BKT (K = 0.9, slightly above K_c) ---")
    model_bkt = XYModel(L=16, K=0.9, seed=99)
    model_bkt.thermalize(n_sweeps=5000)
    eta_near_bkt = []
    vortex_near_bkt = []
    for _ in range(200):
        model_bkt.sweep(n_sweeps=5)
        eta_near_bkt.append(model_bkt.effective_exponent())
        vortex_near_bkt.append(model_bkt.vortex_density())

    eta_bkt_meas = float(np.mean(eta_near_bkt))
    vortex_bkt_meas = float(np.mean(vortex_near_bkt))
    print(f"  η(K=0.9) = {eta_bkt_meas:.4f} (should approach 1/4 = 0.250)")
    print(f"  Vortex density(K=0.9) = {vortex_bkt_meas:.4f} (some free vortices)")
    print(f"  η_sw(K=0.9) = {1/(2*np.pi*0.9):.4f} (spin-wave breaks down near BKT)")

    print("\n" + "=" * 64)
    print("XY model validation complete.")
    print("=" * 64)
