"""S3 — Quasiparticle / catalyzon spectral solve.

Two channels, compared, never forced equal:

1. **Democratic lock** (L1, phenomenological until H-TEGR1 is computed):
   G = g₀(J−I), Spec(G) = {(N−1)g₀ ×1, −g₀ ×(N−1)}, Helmert dark basis,
   Δ_cat = Δ_ex − g₀. Parameter dependence g₀(θ, N) is *open* in the paper;
   it is exposed here as a sweep, not a new lock.
2. **Wigner–Seitz overlap** G^(H)_ij = ∬ Φ_i K₀(|r−r'|/ξ_K) Φ_j (Phase Two
   WS2, :mod:`qvc.gij_tqgl_ws`) — a computed microscopic model. We report its
   spectrum and the Frobenius distance to g₀(J−I) as a diagnostic.

Honesty flags carried on every output: five-mode catalyzon is phenomenological;
H-FN1 is not accepted on the truncated TQGL (one discrete U(1) Goldstone);
H-TOR1 (torsion → G_ij) is parked and never feeds this module.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.coupling_matrix import (
    catalyzon_dark_basis,
    democratic_coupling_matrix,
    democratic_spectrum,
)
from qvc.gij_tqgl_ws import run_tqgl_ws_G
from qvc.workbench.locks import Layer, WorkbenchLocks, get_locks, stamp

# Paper C §Five-Mode table (boundary modes at ω₀ = 2Δ/ℏ; residual U(1)_fiber × C6v)
FIVE_MODES: list[dict[str, str]] = [
    {"origin": "LA phonon (radial)", "boundary_mode": "Bogoliubov amplitude", "sym": "A1", "color": "#d62728"},
    {"origin": "TA1 phonon", "boundary_mode": "Josephson plasma", "sym": "E1", "color": "#1f77b4"},
    {"origin": "TA2 phonon", "boundary_mode": "Flat-band plasmon", "sym": "E1", "color": "#17becf"},
    {"origin": "Torsion wave", "boundary_mode": "Hopfion fiber oscillation", "sym": "A2", "color": "#9467bd"},
    {"origin": "Fiber magnon", "boundary_mode": "Out-of-phase Bogoliubov", "sym": "B1", "color": "#2ca02c"},
]

HONESTY = {
    "five_mode_status": "phenomenological (Paper C §Five-Mode); democratic G is a lock, not a derivation",
    "H_FN1": "NOT accepted on truncated TQGL: one discrete U(1) Goldstone, not five",
    "H_TEGR1": "open — microscopic derivation of democratic G from teleparallel overlaps",
    "H_TOR1": "parked — torsion is never fed into G_ij here",
    "WS_overlap": "computed microscopic model (PROGRAM_COMPUTED), compared to the lock, not identified with it",
}


def democratic_panel(N: int | None = None, g0_meV: float | None = None, *, Delta_ex_meV: float | None = None) -> dict[str, Any]:
    """Democratic G, exact spectrum, numeric spectrum, Helmert dark basis, Δ_cat."""
    L = get_locks()
    N = L.N if N is None else int(N)
    g0 = L.g0_star_meV if g0_meV is None else float(g0_meV)
    Delta_ex = L.Delta0_meV if Delta_ex_meV is None else float(Delta_ex_meV)
    G = democratic_coupling_matrix(N, g0)
    spec = democratic_spectrum(N, g0)
    eigs, vecs = np.linalg.eigh(G)
    B = catalyzon_dark_basis(N)  # (N-1, N)
    bright = np.ones(N) / math.sqrt(N)
    # dark subspace projector check
    dark_proj = B.T @ B
    resid = float(np.linalg.norm(G @ dark_proj - (-g0) * dark_proj))
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, f"Democratic G = g₀(J−I), N={N}, g₀={g0:.4f} meV — phenomenological lock"),
        "N": N,
        "g0_meV": g0,
        "G": G,
        "eigenvalues_numeric": eigs,
        "eigenvectors": vecs,
        "lambda_max_meV": spec["lambda_max"],
        "lambda_min_meV": spec["lambda_min"],
        "lambda_min_multiplicity": spec["lambda_min_mult"],
        "bright_mode": bright,
        "helmert_dark_basis": B,
        "dark_subspace_residual": resid,
        "Delta_ex_meV": Delta_ex,
        "Delta_cat_meV": Delta_ex - g0,
        "reduction_pct": 100.0 * g0 / Delta_ex,
        "omega0_ps_inv": 2.0 * Delta_ex / L.tqgl.hbar_meV_ps,
        "trace_zero": abs(float(np.trace(G))) < 1e-12,
        "honesty": HONESTY,
    }


def spec_sweep(
    *,
    g0_range_meV: np.ndarray | None = None,
    N_values: tuple[int, ...] = (3, 4, 5, 6, 7),
) -> dict[str, Any]:
    """Live Spec(G) vs g₀ and N (open parameter dependence — sweep, not lock)."""
    L = get_locks()
    if g0_range_meV is None:
        g0_range_meV = np.linspace(0.0, 2.0 * L.g0_star_meV, 41)
    rows = {}
    for N in N_values:
        rows[N] = {
            "lambda_max": (N - 1) * g0_range_meV,
            "lambda_min": -g0_range_meV,
            "Delta_cat": L.Delta0_meV - g0_range_meV,
        }
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, "Spec(G) vs g₀ — g₀(θ,N) dependence is open in the paper; this is a sweep"),
        "g0_meV": g0_range_meV,
        "by_N": rows,
        "g0_star_meV": L.g0_star_meV,
        "E_vac_meV": L.E_vac_meV,
        "lambda_star": L.lambda_star,
        "note": "g₀* = λ* E_vac ties the catalyzon channel to the same participation as the fold; it is not added on top of Δ*",
    }


def ws_overlap_comparison(
    *,
    xi_over_lambda_values: tuple[float, ...] = (0.5, 1.0, 2.0, 3.85, 8.0),
    n_grid: int = 40,
    N: int | None = None,
) -> dict[str, Any]:
    """Wigner–Seitz G^(H) for several kernel ranges vs democratic g₀(J−I)."""
    L = get_locks()
    N = L.N if N is None else int(N)
    rows = []
    for xi in xi_over_lambda_values:
        r = run_tqgl_ws_G(N=N, n_grid=n_grid, xi_over_lambda=xi)
        eig = np.asarray(r["eigenvalues"])
        g0e = r["g0_eff_mean_offdiag"] or 1e-300
        eig_n = eig / g0e  # dimensionless: democratic → {−1 ×(N−1), N−1 ×1}
        mult_rel = int(np.sum(np.abs(eig_n - eig_n[0]) < 0.05 * max(abs(eig_n[0]), 1e-12)))
        rows.append(
            {
                "xi_over_lambda_M": xi,
                "xi_nm": r["xi_kernel_nm"],
                "eigenvalues": eig,
                "eigenvalues_over_g0eff": eig_n,
                "g0_eff_meV_units": r["g0_eff_mean_offdiag"],
                "lambda_min": r["lambda_min"],
                "lambda_max": r["lambda_max"],
                "lambda_min_over_g0eff": float(eig_n[0]),
                "lambda_max_over_g0eff": float(eig_n[-1]),
                "lambda_min_mult_numeric": r["lambda_min_mult_numeric"],
                "lambda_min_mult_rel5pct": mult_rel,
                "rel_frobenius_to_democratic": r["rel_frobenius_to_democratic"],
                "offdiag_cv": r["offdiag_cv"],
                "near_democratic": r["acceptance"]["near_democratic"],
                "G": np.asarray(r["G"]),
            }
        )
    default = run_tqgl_ws_G(N=N, n_grid=n_grid)  # ξ_K = R_preprint
    d_eig = np.asarray(default["eigenvalues"]) / (default["g0_eff_mean_offdiag"] or 1e-300)
    return {
        "default_eigenvalues_over_g0eff": d_eig,
        "layer": Layer.L1.value,
        "caption": stamp(
            Layer.L1,
            "WS overlap G^(H) = ∬Φ_i K₀(|r−r'|/ξ_K)Φ_j — computed model, compared to democratic lock (not forced)",
        ),
        "N": N,
        "rows": rows,
        "default_xi_nm": default["xi_kernel_nm"],
        "default_rel_frobenius": default["rel_frobenius_to_democratic"],
        "default_eigenvalues": np.asarray(default["eigenvalues"]),
        "default_G": np.asarray(default["G"]),
        "status": default["status"],
        "honesty": HONESTY,
    }


# --------------------------------------------------------------------------
# 3D render helpers (five mode rings / layers, dark-multiplet isosurface)
# --------------------------------------------------------------------------


def mode_ring_geometry(
    *,
    R_nm: float | None = None,
    r_nm: float | None = None,
    layer_spacing_nm: float = 0.34,
    n_pts: int = 200,
    amplitudes: np.ndarray | None = None,
) -> dict[str, Any]:
    """Five colored rings on the torus, one per layer (z = i·d⊥), radius modulated by |amplitude|.

    ``amplitudes`` defaults to the bright mode (equal weights); pass a Helmert
    dark vector to see a sum-zero (dark) pattern.
    """
    L = get_locks()
    R = L.R_nm if R_nm is None else R_nm
    r = L.r_nm if r_nm is None else r_nm
    N = len(FIVE_MODES)
    amp = np.ones(N) / math.sqrt(N) if amplitudes is None else np.asarray(amplitudes, float)
    phi = np.linspace(0, 2 * np.pi, n_pts)
    rings = []
    z0 = -0.5 * (N - 1) * layer_spacing_nm * 10.0  # exaggerate spacing ×10 for visibility
    for i, m in enumerate(FIVE_MODES):
        rad = R + 0.6 * r * amp[i] * np.cos(3 * phi)  # 3-fold modulation = C6v sub-harmonic sketch
        rings.append(
            {
                "name": m["boundary_mode"],
                "origin": m["origin"],
                "sym": m["sym"],
                "color": m["color"],
                "amplitude": float(amp[i]),
                "x": rad * np.cos(phi),
                "y": rad * np.sin(phi),
                "z": np.full_like(phi, z0 + i * layer_spacing_nm * 10.0),
            }
        )
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, "Five phenomenological catalyzon modes as layer rings (amplitude = eigenvector weight)"),
        "rings": rings,
        "R_nm": R,
        "r_nm": r,
        "note": "layer spacing exaggerated ×10; ring modulation is a symmetry sketch, not a computed mode shape",
    }


def dark_multiplet_field(
    *,
    n_grid: int = 48,
    dark_index: int = 0,
    N: int | None = None,
    sigma_over_lambda: float = 0.18,
    centre_over_Rhex: float = 0.45,
) -> dict[str, Any]:
    """Ψ_dark(r) = Σ_i B[dark_index, i] Φ_i(r) on the WS hexagon (isosurface / heatmap).

    Uses the same Gaussian mode profiles as the WS-overlap model so the
    dark-multiplet picture and G^(H) share one geometry.
    """
    from qvc.gij_tqgl_ws import gaussian_modes, mode_centres_polygon, ws_hexagon_mask
    from qvc.params import PARAMS

    L = get_locks()
    N = L.N if N is None else int(N)
    lam = PARAMS.lambda_M_um * 1e3  # nm
    R_hex = lam / math.sqrt(3.0)
    lim = R_hex * 1.05
    axis = np.linspace(-lim, lim, n_grid)
    X, Y = np.meshgrid(axis, axis, indexing="xy")
    dA = (axis[1] - axis[0]) ** 2
    mask = ws_hexagon_mask(X, Y, R_hex)
    centres = mode_centres_polygon(N, centre_over_Rhex * R_hex)
    Phi = gaussian_modes(X, Y, centres, sigma_over_lambda * lam, mask, dA)
    B = catalyzon_dark_basis(N)
    coeff = B[dark_index % (N - 1)]
    psi = np.zeros(X.shape)
    psi[mask] = coeff @ Phi
    bright = np.zeros(X.shape)
    bright[mask] = (np.ones(N) / math.sqrt(N)) @ Phi
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, f"Dark (sum-zero) multiplet #{dark_index+1} on the moiré WS cell; bright mode for contrast"),
        "x_nm": axis,
        "y_nm": axis,
        "psi_dark": psi,
        "psi_bright": bright,
        "mask": mask,
        "centres_nm": centres,
        "coefficients": coeff,
        "lambda_M_nm": lam,
    }


def catalyzon_report() -> dict[str, Any]:
    """Compact JSON-able S3 summary for the acceptance report."""
    d = democratic_panel()
    ws = ws_overlap_comparison(n_grid=32)
    return {
        "democratic": {
            "N": d["N"],
            "g0_meV": d["g0_meV"],
            "lambda_max_meV": d["lambda_max_meV"],
            "lambda_min_meV": d["lambda_min_meV"],
            "lambda_min_multiplicity": d["lambda_min_multiplicity"],
            "numeric_eigenvalues": d["eigenvalues_numeric"].tolist(),
            "Delta_cat_meV": d["Delta_cat_meV"],
            "reduction_pct": d["reduction_pct"],
            "dark_subspace_residual": d["dark_subspace_residual"],
        },
        "ws_overlap": [
            {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in row.items() if k != "G"}
            for row in ws["rows"]
        ],
        "ws_default_rel_frobenius": ws["default_rel_frobenius"],
        "honesty": HONESTY,
    }
