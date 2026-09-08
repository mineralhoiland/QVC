"""B1: G_ij from Wigner–Seitz Gaussian and ring-Fourier/Wannier overlaps.

Status
------
This is a **computed kernel model**. It is not torsion (H-TEGR1), not a
linearized TQGL theorem, and it does not insert G = g0 (J − I) by hand.

Two constructions share the Macdonald K0 vacuum kernel
K(r,r') = K0(|r−r'| / ξ_K) with ξ_K = R (torus major):

1. Existing ``run_tqgl_ws_G`` (Gaussian blobs on the moiré WS hexagon),
   re-run at a finer grid.
2. NEW ring modes: angular harmonics n = 0 … N−1 of the hopfion envelope
   on the locked torus, Wannier-transformed to N localized packets.
   G_ij = ∬ Φ_i*(r) K(r,r') Φ_j(r') d²r d²r'  (i ≠ j),  G_ii = 0.

Democratic Spec(G) is locked *algebra* at whatever g0 is used:
    λ_max = (N−1) g0,  λ_min = −g0 (multiplicity N−1).
For N=5 that is (4 g0, −g0 × 4). Compare the computed spectrum to that
shape at the *computed* g0_eff (mean off-diagonal). Do not replace g0_eff
by the BRIDGE identification g0★ = λ★ E_vac (meV); they are different
objects (kernel units vs meV).

Optional diagnostic: GPE quadratic (kinetic) overlap ⟨Φ_i| −(ħ²/2m*)∇² |Φ_j⟩
with the same G_ii=0 convention. That kernel is *not* Macdonald K0.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.coupling_matrix import democratic_coupling_matrix, democratic_spectrum  # noqa: E402
from qvc.gij_tqgl_ws import (  # noqa: E402
    kernel_matrix_K0,
    overlap_G,
    run_tqgl_ws_G,
)
from qvc.tqgl.gpe import hopfion_initial, kinetic_coeff_mev_um2, make_grid  # noqa: E402
from qvc.tqgl.locks import FrozenTQGL  # noqa: E402
from qvc_l1.locks import L1Locks, build_l1_locks  # noqa: E402


# ---------------------------------------------------------------------------
# Spectrum comparison (never builds G from J−I except as a reference)
# ---------------------------------------------------------------------------


def compare_to_democratic(G: np.ndarray, g0_eff: float) -> dict[str, Any]:
    """Shape comparison of a computed Hermitian G to democratic Spec at g0_eff."""
    n = int(G.shape[0])
    eigs = np.sort(np.linalg.eigvalsh(np.real(0.5 * (G + G.conj().T))))
    G_dem = democratic_coupling_matrix(n, g0_eff)
    spec = democratic_spectrum(n, g0_eff)
    fro = float(np.linalg.norm(np.real(G) - G_dem, ord="fro"))
    rel_fro = fro / (float(np.linalg.norm(G_dem, ord="fro")) + 1e-18)
    off = np.real(G)[np.triu_indices(n, k=1)]
    return {
        "eigenvalues": eigs.tolist(),
        "lambda_max": float(eigs[-1]),
        "lambda_min": float(eigs[0]),
        "lambda_min_mult_numeric": int(np.sum(np.isclose(eigs, eigs[0], atol=1e-8, rtol=0))),
        "democratic_lambda_max": float(spec["lambda_max"]),
        "democratic_lambda_min": float(spec["lambda_min"]),
        "democratic_lambda_min_mult": int(spec["lambda_min_mult"]),
        "rel_err_lambda_max": abs(float(eigs[-1]) - spec["lambda_max"])
        / (abs(spec["lambda_max"]) + 1e-18),
        "rel_err_lambda_min": abs(float(eigs[0]) - spec["lambda_min"])
        / (abs(spec["lambda_min"]) + 1e-18),
        "frobenius_to_democratic": fro,
        "rel_frobenius_to_democratic": rel_fro,
        "offdiag_cv": float(np.std(off) / (abs(g0_eff) + 1e-18)),
        "g0_eff_mean_offdiag": float(g0_eff),
    }


def _g0_eff_from_G(G: np.ndarray) -> float:
    n = G.shape[0]
    off = np.real(G)[np.triu_indices(n, k=1)]
    return float(np.mean(off)) if off.size else 0.0


def _acceptance(G: np.ndarray, cmp: dict[str, Any], *, rel_fro_cut: float = 0.15) -> dict[str, Any]:
    eigs = np.asarray(cmp["eigenvalues"], dtype=float)
    return {
        "hermitian": bool(np.allclose(G, G.conj().T, atol=1e-10)),
        "diag_zero": bool(np.allclose(np.diag(G), 0.0, atol=1e-12)),
        "trace_zero": abs(float(np.trace(G))) < 1e-10,
        "lambda_min_lt_0": float(eigs[0]) < 0.0,
        "near_democratic": bool(cmp["rel_frobenius_to_democratic"] < rel_fro_cut),
        "not_hand_inserted_J_minus_I": True,
    }


# ---------------------------------------------------------------------------
# Complex overlap (K0) and GPE kinetic overlap
# ---------------------------------------------------------------------------


def overlap_G_complex(
    Phi: np.ndarray,
    K: np.ndarray,
    dA: float,
    *,
    zero_diagonal: bool = True,
) -> np.ndarray:
    """G_ij = Σ_{ab} Φ_i*(r_a) K_ab Φ_j(r_b) dA². Hermitian; optional G_ii=0.

    Phi: (N, M), K: (M, M) real symmetric.
    """
    gk = K @ Phi.T  # (M, N)
    g = (Phi.conj() @ gk) * (dA * dA)
    g = 0.5 * (g + g.conj().T)
    if zero_diagonal:
        np.fill_diagonal(g, 0.0)
    return g


def kinetic_overlap_G(
    Phi_2d: np.ndarray,
    k2: np.ndarray,
    dA: float,
    ke_coef: float,
    *,
    zero_diagonal: bool = True,
) -> np.ndarray:
    """G_ij = ⟨Φ_i| ke_coef (−∇²) |Φ_j⟩, then G_ii=0.

    Phi_2d: (N, Ny, Nx). Units: meV if Φ is L2-normalised on the grid.
    This is the GPE *quadratic* kernel, not Macdonald K0.
    """
    n = Phi_2d.shape[0]
    g = np.zeros((n, n), dtype=np.complex128)
    laps = np.empty_like(Phi_2d)
    for j in range(n):
        laps[j] = np.fft.ifft2(k2 * np.fft.fft2(Phi_2d[j]))
    for i in range(n):
        for j in range(n):
            g[i, j] = ke_coef * np.sum(np.conj(Phi_2d[i]) * laps[j]) * dA
    g = 0.5 * (g + g.conj().T)
    if zero_diagonal:
        np.fill_diagonal(g, 0.0)
    return g


# ---------------------------------------------------------------------------
# Ring Fourier harmonics and Wannier packets from hopfion_initial
# ---------------------------------------------------------------------------


def hopfion_envelope(grid, lock: FrozenTQGL, *, n0: float = 1.0) -> np.ndarray:
    """Torus envelope: ring Gaussian × |hopfion_initial| (noise off).

    hopfion_initial is a filled disk with a 10% ring dip; multiplying by
    ``gauss_ring`` confines the modes to the locked torus. Amplitude only:
    hopfion ≠ 2D phase vortex ≠ CS fluxon.
    """
    psi = hopfion_initial(
        grid, lock, n0=n0, core_depletion=0.10, noise_amp=0.0, seed=0
    )
    return (grid.gauss_ring * np.abs(psi)).astype(float)


def angular_harmonic_fields(
    envelope: np.ndarray,
    theta: np.ndarray,
    n_modes: int,
) -> np.ndarray:
    """Φ_n(r) = env(ρ) e^{i n θ} / ‖·‖, n = 0 … n_modes−1. Shape (N, Ny, Nx)."""
    fields = np.empty((n_modes,) + envelope.shape, dtype=np.complex128)
    for n in range(n_modes):
        fields[n] = envelope * np.exp(1j * n * theta)
    return fields


def _real_peak_gauge(fields: np.ndarray, dA: float) -> np.ndarray:
    """Remove a global phase per mode so the peak is real (real-orbital gauge)."""
    out = np.empty_like(fields)
    for i, f in enumerate(fields):
        peak = np.unravel_index(int(np.argmax(np.abs(f))), f.shape)
        phase = float(np.angle(f[peak]))
        out[i] = f * np.exp(-1j * phase)
    return _l2_normalize_stack(out, dA)


def wannier_from_harmonics(
    harmonics: np.ndarray,
    theta: np.ndarray,
    *,
    dA: float,
) -> np.ndarray:
    """Localized ring packets: w_j = N^{-1/2} Σ_n Φ_n exp(−i n θ_j).

    θ_j = 2π j / N. Built *from* angular harmonics n=0..N−1; not inserted
    as Gaussians by hand, and not G = g0(J−I). Real-peak gauge makes G
    essentially real for a real kernel.
    """
    n_modes = harmonics.shape[0]
    sites = np.empty_like(harmonics)
    for j in range(n_modes):
        th_j = 2.0 * math.pi * j / n_modes
        acc = np.zeros(harmonics.shape[1:], dtype=np.complex128)
        for n in range(n_modes):
            acc += harmonics[n] * np.exp(-1j * n * th_j)
        sites[j] = acc / math.sqrt(n_modes)
    return _real_peak_gauge(sites, dA)


def _l2_normalize_stack(fields: np.ndarray, dA: float) -> np.ndarray:
    out = np.empty_like(fields)
    for i, f in enumerate(fields):
        nrm = math.sqrt(float(np.sum(np.abs(f) ** 2) * dA))
        if nrm <= 0.0:
            raise RuntimeError("mode normalisation vanished")
        out[i] = f / nrm
    return out


def ring_mask(gauss_ring: np.ndarray, *, thresh: float = 0.02) -> np.ndarray:
    return gauss_ring >= thresh


def ring_theta_profiles(
    fields: np.ndarray,
    grid,
    *,
    n_theta: int = 360,
) -> dict[str, Any]:
    """Interpolate |Φ| onto the major circle ρ = R."""
    ny, nx = grid.X.shape
    # nearest-neighbour sample on the circle
    theta = np.linspace(-math.pi, math.pi, n_theta, endpoint=False)
    xs = grid.R_um * np.cos(theta)
    ys = grid.R_um * np.sin(theta)
    # map to indices
    ix = np.clip(np.rint((xs - grid.x[0]) / grid.dx).astype(int), 0, nx - 1)
    iy = np.clip(np.rint((ys - grid.y[0]) / grid.dx).astype(int), 0, ny - 1)
    amp = np.empty((fields.shape[0], n_theta), dtype=float)
    real = np.empty((fields.shape[0], n_theta), dtype=float)
    for i, f in enumerate(fields):
        z = f[iy, ix]
        amp[i] = np.abs(z)
        real[i] = np.real(z)
    return {"theta": theta, "amp": amp, "real": real, "R_um": grid.R_um}


def _sample_on_mask(fields: np.ndarray, mask: np.ndarray) -> np.ndarray:
    """(N, n_pts) complex samples."""
    return np.stack([f[mask] for f in fields], axis=0)


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------


def run_ws_gaussian_G(
    locks: L1Locks | None = None,
    *,
    n_grid: int = 80,
) -> dict[str, Any]:
    """Re-run ``run_tqgl_ws_G`` at a finer grid. G is computed, not J−I."""
    if locks is None:
        locks = build_l1_locks()
    raw = run_tqgl_ws_G(n_grid=n_grid)
    g0_eff = float(raw["g0_eff_mean_offdiag"])
    G = np.asarray(raw["G"], dtype=float)
    cmp = compare_to_democratic(G, g0_eff)
    return {
        "status": "PROGRAM_COMPUTED_kernel_model_not_torsion_not_TQGL_theorem",
        "construction": "WS Gaussian × Macdonald K0 on moiré hexagon (run_tqgl_ws_G)",
        "n_grid": int(raw["n_grid"]),
        "n_pts": int(raw["n_pts"]),
        "lambda_M_nm": raw["lambda_M_nm"],
        "R_hex_nm": raw["R_hex_nm"],
        "sigma_nm": raw["sigma_nm"],
        "xi_kernel_nm": raw["xi_kernel_nm"],
        "G": G,
        "g0_eff_kernel": g0_eff,
        "g0_star_meV_BRIDGE": locks.g0_star_meV,
        "g0_eff_is_not_g0_star": True,
        "units": "Macdonald K0 overlap (dimensionless kernel), not meV",
        **cmp,
        "acceptance": _acceptance(G, cmp),
        "hermitian": bool(np.allclose(G, G.T)),
        "notes": (
            "Finer-grid WS overlap. Democratic form is an OUTPUT of long-range K, "
            "not an input. g0_eff (kernel) is not identified with g0★ = λ★ E_vac."
        ),
    }


def run_ring_fourier_G(
    locks: L1Locks | None = None,
    *,
    n_grid: int = 96,
    L_um: float = 0.200,
    mask_thresh: float = 0.02,
    n_theta_profile: int = 360,
) -> dict[str, Any]:
    """Ring Fourier n=0..N−1 → Wannier packets; K0 overlap; G_ii=0."""
    if locks is None:
        locks = build_l1_locks()
    lock = locks.tqgl
    n = locks.N_layers
    grid = make_grid(lock, n_grid=n_grid, L_um=L_um)
    dA = grid.dx * grid.dx
    env = hopfion_envelope(grid, lock)
    harm = _l2_normalize_stack(angular_harmonic_fields(env, grid.theta, n), dA)
    wannier = wannier_from_harmonics(harm, grid.theta, dA=dA)

    mask = ring_mask(grid.gauss_ring, thresh=mask_thresh)
    pts = np.stack([grid.X[mask], grid.Y[mask]], axis=1)
    xi = lock.R_um  # same Macdonald range as WS (torus major)
    K = kernel_matrix_K0(pts, xi)

    Phi_w = _sample_on_mask(wannier, mask)
    G = overlap_G_complex(Phi_w, K, dA, zero_diagonal=True)
    G_use = np.real(0.5 * (G + G.conj().T))
    imag_max = float(np.max(np.abs(np.imag(G))))

    g0_eff = _g0_eff_from_G(G_use)
    cmp = compare_to_democratic(G_use, g0_eff)

    # Fourier-basis kernel (rotationally invariant K is nearly diagonal).
    # Reported as a diagnostic; NOT zeroed, because that would erase Spec(K).
    Phi_f = _sample_on_mask(harm, mask)
    G_fourier = overlap_G_complex(Phi_f, K, dA, zero_diagonal=False)
    fourier_eigs = np.sort(np.linalg.eigvalsh(np.real(0.5 * (G_fourier + G_fourier.conj().T))))

    # GPE quadratic kernel on the same Wannier packets (diagnostic).
    ke = kinetic_coeff_mev_um2(lock)
    G_ke = kinetic_overlap_G(wannier, grid.K2, dA, ke, zero_diagonal=True)
    G_ke_real = np.real(G_ke)
    g0_ke = _g0_eff_from_G(G_ke_real)
    cmp_ke = compare_to_democratic(G_ke_real, g0_ke)

    profiles_w = ring_theta_profiles(wannier, grid, n_theta=n_theta_profile)
    profiles_f = ring_theta_profiles(harm, grid, n_theta=n_theta_profile)

    # 2D |w_0| slice for heatmap of a single packet
    return {
        "status": "PROGRAM_COMPUTED_kernel_model_not_torsion_not_TQGL_theorem",
        "construction": (
            "hopfion_initial envelope × angular harmonics n=0..N-1, "
            "Wannier-localized on the locked torus, overlapped with Macdonald K0"
        ),
        "n_grid": n_grid,
        "L_um": L_um,
        "n_pts": int(mask.sum()),
        "xi_kernel_nm": xi * 1e3,
        "n_modes": n,
        "harmonics": list(range(n)),
        "G": G_use,
        "G_imag_max": imag_max,
        "g0_eff_kernel": g0_eff,
        "g0_star_meV_BRIDGE": locks.g0_star_meV,
        "g0_eff_is_not_g0_star": True,
        "units": "Macdonald K0 overlap (dimensionless kernel), not meV",
        **cmp,
        "acceptance": _acceptance(G_use, cmp, rel_fro_cut=0.15),
        "fourier_kernel_eigenvalues": fourier_eigs.tolist(),
        "fourier_kernel_note": (
            "K0 is radial, so the angular-harmonic basis nearly diagonalises it. "
            "Coupling G uses the Wannier transform of those harmonics, then G_ii=0."
        ),
        "gpe_quadratic": {
            "status": "diagnostic_kinetic_overlap_not_K0",
            "units": "meV (GPE ke_coef × Laplacian overlap, G_ii=0)",
            "ke_coef_meV_um2": ke,
            "G": G_ke_real,
            "g0_eff_meV": g0_ke,
            **cmp_ke,
            "acceptance": _acceptance(G_ke_real, cmp_ke, rel_fro_cut=0.15),
        },
        "profiles": {
            "theta": profiles_w["theta"],
            "wannier_amp": profiles_w["amp"],
            "fourier_amp": profiles_f["amp"],
            "fourier_real": profiles_f["real"],
        },
        "wannier_abs_0": np.abs(wannier[0]),
        "fourier_abs_stack": np.abs(harm),
        "x_um": grid.x,
        "y_um": grid.y,
        "envelope": env,
        "notes": (
            "Computed kernel model. Circulant ring geometry has nn and nnn chords, "
            "so Spec(G) need not collapse to democratic unless K is long-range. "
            "That is a result, not a failure of the lock (the lock is algebra)."
        ),
    }


def run_b1(
    locks: L1Locks | None = None,
    *,
    ws_n_grid: int = 80,
    ring_n_grid: int = 96,
) -> dict[str, Any]:
    """Full B1: finer WS Gaussian + ring Fourier/Wannier."""
    if locks is None:
        locks = build_l1_locks()
    ws = run_ws_gaussian_G(locks, n_grid=ws_n_grid)
    ring = run_ring_fourier_G(locks, n_grid=ring_n_grid)
    return {
        "status": "PROGRAM_COMPUTED_kernel_model_not_torsion_not_TQGL_theorem",
        "locks": {
            "r_nm": locks.r_nm,
            "R_nm": locks.R_nm,
            "Delta0_meV": locks.Delta0_meV,
            "E_vac_meV": locks.E_vac_meV,
            "lambda_star": locks.lambda_star,
            "g0_star_meV": locks.g0_star_meV,
            "N_layers": locks.N_layers,
            "democratic_algebra": {
                "lambda_max": f"{locks.N_layers - 1}*g0",
                "lambda_min": "-g0",
                "lambda_min_mult": locks.dem_lambda_min_mult,
                "inserted_J_minus_I": False,
            },
        },
        "ws_gaussian": ws,
        "ring_fourier": ring,
        "computed_vs_ansatz": {
            "computed": [
                "WS Gaussian × K0 overlap on finer grid",
                "ring Wannier packets from hopfion envelope × harmonics n=0..4",
                "Hermitian G with G_ii=0",
                "Spec(G) and rel Frobenius vs democratic at g0_eff",
                "optional GPE kinetic overlap (labeled diagnostic)",
            ],
            "locked_algebra": [
                "Democratic Spec: λ_max=(N-1)g0, λ_min=−g0 (mult N−1)",
                "g0★ = λ★ E_vac (BRIDGE, meV)",
            ],
            "not": [
                "torsion / H-TEGR1",
                "linearized TQGL theorem",
                "hand-inserted G=g0(J−I)",
                "identification g0_eff = g0★",
            ],
        },
    }


# Re-export for tests that want the WS primitive.
__all__ = [
    "compare_to_democratic",
    "overlap_G",
    "overlap_G_complex",
    "run_b1",
    "run_ring_fourier_G",
    "run_tqgl_ws_G",
    "run_ws_gaussian_G",
]
