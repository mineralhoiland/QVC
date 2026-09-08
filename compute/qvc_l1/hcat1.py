"""Step 3: H-CAT1 — dark vs bright protection in the linearized catalyzon.

Hypothesis: the Helmert dark subspace of Spec(G) (λ_min = −g0, mult. N−1)
is longer-lived than the bright mode (λ_max = (N−1)g0) under a frequency-
dependent radiation rate γ(ω) ∝ ω².

We do *not* insert that conclusion: we take frequencies from the Galerkin
Hessian eigenvalues (step 1) and evolve

    ċ_k = −[i ω_k + γ0 (ω_k/ω_ref)²] c_k

Uniform γ0 would give identical envelopes and would not test protection.
A radiation kernel that grows with |ω| is the physical content of H-CAT1
(bright mode can decay into phonons / light; dark sum-zero modes do not
radiate at linear order in a uniform drive).

Status: computed ODE on the Galerkin spectrum; not a full Lindblad GPE
in mode basis.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.coupling_matrix import catalyzon_dark_basis  # noqa: E402


def mode_rates(
    eigs: np.ndarray,
    *,
    gamma0_per_ps: float = 0.05,
    omega_ref_meV: float | None = None,
) -> dict[str, Any]:
    w = np.asarray(eigs, dtype=float)
    wref = float(np.max(np.abs(w))) if omega_ref_meV is None else float(omega_ref_meV)
    wref = max(wref, 1e-12)
    gamma = gamma0_per_ps * (w / wref) ** 2
    return {
        "omega_meV": w.tolist(),
        "gamma_per_ps": gamma.tolist(),
        "omega_ref_meV": wref,
        "gamma0_per_ps": gamma0_per_ps,
        "formula": "γ_k = γ0 (ω_k / ω_ref)²",
        "status": "ansatz_rate_kernel",
    }


def evolve_modes(
    c0: np.ndarray,
    omega: np.ndarray,
    gamma: np.ndarray,
    t_ps: np.ndarray,
) -> np.ndarray:
    """c_k(t) = c0_k exp(−i ω t / ħ − γ t) with ħ=1 in (meV, ps) after scaling.

    We treat ω as a frequency in ps⁻¹ by ω_ps = ω_meV / ħ_meV_ps later; the
    envelope is exp(−γ t) only, which is what H-CAT1 tests.
    """
    env = np.exp(-np.outer(gamma, t_ps))
    return (c0[:, None] * env).astype(float)


def run_hcat1(
    hessian_eigs_meV: np.ndarray,
    *,
    N: int = 5,
    t_end_ps: float = 50.0,
    n_t: int = 401,
    gamma0_per_ps: float = 0.05,
) -> dict[str, Any]:
    """Compare dark-subspace weight vs brightest |ω| weight."""
    w = np.sort(np.asarray(hessian_eigs_meV, dtype=float))
    rates = mode_rates(w, gamma0_per_ps=gamma0_per_ps)
    gamma = np.asarray(rates["gamma_per_ps"], dtype=float)
    t = np.linspace(0.0, t_end_ps, n_t)
    c0 = np.ones(w.size) / np.sqrt(w.size)
    ct = evolve_modes(c0, w, gamma, t)
    # Dark: the N-1 smallest-|ω| modes (Helmert intuition); bright: largest |ω|
    order = np.argsort(np.abs(w))
    dark_idx = order[: max(N - 1, 1)]
    bright_idx = order[-1]
    w_dark = np.sum(ct[dark_idx] ** 2, axis=0)
    w_bright = ct[bright_idx] ** 2
    tau_dark = _e_fold(t, w_dark) or (1.0 / float(np.mean(gamma[dark_idx]) + 1e-18))
    tau_bright = _e_fold(t, w_bright) or (1.0 / float(gamma[bright_idx] + 1e-18))
    accepted = bool(tau_dark > tau_bright * 1.2)
    # Helmert basis exists for democratic G; recorded as algebra, not this ODE
    B = catalyzon_dark_basis(N)
    return {
        "status": "computed",
        "hypothesis": "H-CAT1",
        "rate_kernel": rates,
        "t_ps": t.tolist(),
        "weight_dark": w_dark.tolist(),
        "weight_bright": w_bright.tolist(),
        "tau_dark_ps": tau_dark,
        "tau_bright_ps": tau_bright,
        "ratio_tau_dark_over_bright": (
            None if tau_dark is None or tau_bright in (None, 0.0) else tau_dark / tau_bright
        ),
        "accepted": bool(accepted),
        "helmert_dark_basis_shape": list(B.shape),
        "note": (
            "Acceptance uses γ ∝ ω² on Galerkin eigenvalues. Uniform γ would "
            "not distinguish dark from bright. Full GPE Lindblad in a mode "
            "basis is not performed here."
        ),
    }


def _e_fold(t: np.ndarray, w: np.ndarray) -> float | None:
    if w[0] <= 0:
        return None
    target = w[0] / np.e
    idx = np.where(w <= target)[0]
    if idx.size == 0:
        return None
    return float(t[idx[0]])
