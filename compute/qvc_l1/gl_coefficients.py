"""Track A GL coefficients (matched, not Hubbard).

- α = gap-map coupling α = λ E_vac / Δ₀
- α_GL = Landau mass with |α_GL| = Δ₀  (α_GL = −Δ₀)
- κ_GL = Δ₀ a_H²  with a_H = r (Option A*)
- λ = λ★
- β from Mexican-hat / GPE n_ref on hopfion_initial
- κ_S: Derrick estimate on the locked torus (ansatz)

Do not import holographic γ_QM ≈ 3.73 meV·nm² as κ_GL.
Do not use cartoon α0 = −0.5, α_QVC = −0.35.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.tqgl.gpe import hopfion_initial, kinetic_coeff_mev_um2, make_grid  # noqa: E402
from qvc.tqgl.live_gap import mean_density_weighted  # noqa: E402

from qvc_l1.locks import HOLOGRAPHIC_GAMMA_QM_MEV_NM2, L1Locks

# GPE hopfion IC used by the truncated solver (read the field; do not cartoon).
N0_GPE = 1.0
CORE_DEPLETION = 0.10
GPE_L_UM = 0.200
GPE_N_GRID = 128


def kappa_GL_mev_um2(locks: L1Locks) -> float:
    """κ_GL = Δ₀ a_H² = ħ²/(2m*)  [meV·µm²]. Status: derived given a_H = r."""
    return float(kinetic_coeff_mev_um2(locks.tqgl))


def alpha_gap_map(locks: L1Locks, *, lam: float | None = None) -> dict[str, Any]:
    """Dimensionless gap-map coupling α = λ E_vac / Δ₀. Status: derived."""
    lam = locks.lambda_star if lam is None else lam
    alpha = float(lam * locks.E_vac_meV / locks.Delta0_meV)
    return {
        "alpha": alpha,
        "lambda": float(lam),
        "formula": "α = λ E_vac / Δ₀",
        "status": "derived",
        "alpha_over_ac": alpha / locks.alpha_c if locks.alpha_c else float("nan"),
        "note": "Gap-map coupling, not the Landau mass α_GL.",
    }


def hopfion_n_ref(locks: L1Locks, *, n0: float = N0_GPE, noise_amp: float = 0.0) -> dict[str, Any]:
    """Ring-weighted ⟨|ψ|²⟩ on hopfion_initial. Status: computed."""
    grid = make_grid(locks.tqgl, n_grid=GPE_N_GRID, L_um=GPE_L_UM)
    psi = hopfion_initial(
        grid, locks.tqgl, n0=n0, core_depletion=CORE_DEPLETION, noise_amp=noise_amp, seed=0
    )
    n_ref = float(mean_density_weighted(psi, grid.weight))
    n_mean = float(np.mean(np.abs(psi) ** 2))
    return {
        "n_ref": n_ref,
        "n_mean_box": n_mean,
        "n0": float(n0),
        "core_depletion": CORE_DEPLETION,
        "noise_amp": float(noise_amp),
        "L_um": GPE_L_UM,
        "n_grid": GPE_N_GRID,
        "status": "computed",
        "formula": "n_ref = ⟨|ψ|²⟩_{ring} on hopfion_initial",
        "note": (
            "Same n_ref construction as the truncated GPE/DCE solver. "
            "Cartoon α0=−0.5, α_QVC=−0.35 is not used."
        ),
    }


def beta_from_n_ref(locks: L1Locks, n_ref: float) -> dict[str, Any]:
    """Mexican-hat quartic: ρ_eq² = |α_GL|/β = n_ref ⇒ β = Δ₀ / n_ref.

    Status: matched (GL minimum sits on the GPE hopfion density).
    """
    if n_ref <= 0.0:
        raise ValueError("n_ref must be positive")
    alpha_GL = -locks.Delta0_meV
    beta = float(locks.Delta0_meV / n_ref)
    return {
        "alpha_GL_meV": alpha_GL,
        "alpha_GL_abs_meV": locks.Delta0_meV,
        "beta_meV": beta,
        "n_ref": float(n_ref),
        "rho_eq2": float(n_ref),
        "formula": "β = |α_GL| / n_ref = Δ₀ / n_ref",
        "status": "matched",
        "gpe_cubic_quintic_beta_at_n0": float(locks.Delta0_meV / N0_GPE),
        "note": (
            "GL Mexican hat uses this β so the vacuum manifold sits at the "
            "hopfion n_ref. The truncated GPE EOM uses cubic-quintic "
            "α_nl=−Δ₀, β_nl=Δ₀/n0 with n0=1 (dimensionless), a related but "
            "not identical matching."
        ),
    }


def kappa_S_derrick(locks: L1Locks, *, kappa_GL: float | None = None) -> dict[str, Any]:
    """Crude Faddeev–Skyrme Derrick balance on the locked torus tube.

    3D hopfion: gradient ~ λ, Skyrme ~ λ^{-1}. Stationarity at the locked
    tube scale a_H gives κ_S ~ κ_GL a_H² = Δ₀ a_H⁴. The torus length 2πR
    cancels. Status: ansatz.
    """
    k_gl = kappa_GL_mev_um2(locks) if kappa_GL is None else float(kappa_GL)
    a_H = locks.a_H_um
    kappa_S = float(k_gl * a_H * a_H)  # = Δ₀ a_H⁴
    kappa_S_R = float(k_gl * locks.R_um * locks.R_um)
    return {
        "kappa_S_meV_um4": kappa_S,
        "kappa_S_meV_nm4": kappa_S * 1e12,  # µm⁴ → nm⁴ is 1e12
        "formula": "κ_S ≈ κ_GL a_H² = Δ₀ a_H⁴  (Derrick balance at tube scale)",
        "status": "ansatz",
        "ansatz_label": "Faddeev–Skyrme 3D Derrick on locked tube a_H; 2πR cancels",
        "alternate_major_radius_kappa_S_meV_um4": kappa_S_R,
        "alternate_note": (
            "κ_S^(R) ≈ κ_GL R² uses the major radius instead of a_H. "
            "Labeled alternate, not the tube-scale ansatz."
        ),
        "not_holographic_gamma_QM": True,
        "gamma_QM_meV_nm2_not_imported": HOLOGRAPHIC_GAMMA_QM_MEV_NM2,
        "note": (
            "Crude scaling, not a calibrated Skyrme coefficient. "
            "Do not import holographic γ_QM≈3.73 meV·nm² as κ_GL."
        ),
    }


def gl_coefficients(locks: L1Locks) -> dict[str, Any]:
    """Full Track A GL card."""
    k_gl = kappa_GL_mev_um2(locks)
    nref = hopfion_n_ref(locks)
    beta = beta_from_n_ref(locks, nref["n_ref"])
    alpha_star = alpha_gap_map(locks, lam=locks.lambda_star)
    alpha_fold = alpha_gap_map(locks, lam=locks.lambda_fold)
    ks = kappa_S_derrick(locks, kappa_GL=k_gl)
    return {
        "alpha_gap_map_star": alpha_star,
        "alpha_gap_map_fold": alpha_fold,
        "alpha_GL": {
            "value_meV": beta["alpha_GL_meV"],
            "abs_meV": beta["alpha_GL_abs_meV"],
            "formula": "α_GL = −Δ₀  (|α_GL|=Δ₀)",
            "status": "matched",
            "note": "Landau mass, not the dimensionless gap-map α.",
        },
        "kappa_GL": {
            "value_meV_um2": k_gl,
            "formula": "κ_GL = Δ₀ a_H² with a_H = r",
            "status": "derived",
            "a_H_um": locks.a_H_um,
            "not_gamma_QM": True,
        },
        "lambda_star": {
            "value": locks.lambda_star,
            "status": "locked",
            "rule": "λ★ = 0.90 λ_fold",
        },
        "lambda_fold": {
            "value": locks.lambda_fold,
            "status": "derived",
            "formula": "λ_fold = α_c Δ₀ / E_vac",
        },
        "n_ref": nref,
        "beta": beta,
        "kappa_S": ks,
        "cartoon_not_used": {
            "alpha0": -0.5,
            "alpha_QVC": -0.35,
            "used": False,
        },
    }
