"""Step 5: catalyzon spectral function A_cat(ω).

    A_cat(ω) = Σ_k |c_k|²  (Γ_k/π) / [(ω − ω_k)² + Γ_k²]

ω_k from Galerkin Hessian (step 1). Γ_k from the H-CAT1 radiation kernel
or a single Γ_cat. Resonance remains ω0 = 2Δ/ℏ at the two-particle
threshold of the phonon section; this A_cat is the *mode* spectral
density, not a tunneling measurement.

Status: computed from the linearized spectrum; Γ is phenomenological.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc_l1.hcat1 import mode_rates


def lorentzian(omega: np.ndarray, omega0: float, gamma: float) -> np.ndarray:
    g = max(float(gamma), 1e-12)
    return (g / np.pi) / ((omega - omega0) ** 2 + g * g)


def run_spectral_cat(
    hessian_eigs_meV: np.ndarray,
    *,
    Delta0_meV: float,
    g0_star_meV: float,
    n_omega: int = 801,
) -> dict[str, Any]:
    w = np.sort(np.asarray(hessian_eigs_meV, dtype=float))
    rates = mode_rates(w, gamma0_per_ps=0.05)
    # Convert γ (1/ps) to meV: Γ = ħ γ with ħ ≈ 0.658 meV·ps
    hbar = 0.6582119569
    Gamma = hbar * np.asarray(rates["gamma_per_ps"], dtype=float)
    om = np.linspace(float(w[0]) - 0.05, float(w[-1]) + 0.05, n_omega)
    A = np.zeros_like(om)
    ck2 = np.ones_like(w) / w.size
    for wk, gk, p in zip(w, Gamma, ck2):
        A += p * lorentzian(om, float(wk), float(gk))
    omega0 = 2.0 * Delta0_meV  # two-particle threshold in meV (ħ=1 units of energy)
    Delta_cat = Delta0_meV - g0_star_meV
    return {
        "status": "computed",
        "omega_meV": om.tolist(),
        "A_cat": A.tolist(),
        "poles_meV": w.tolist(),
        "Gamma_meV": Gamma.tolist(),
        "omega0_two_particle_meV": omega0,
        "Delta_cat_meV": Delta_cat,
        "peak_omega_meV": float(om[int(np.argmax(A))]),
        "note": (
            "A_cat is the Galerkin-mode Lorentzian mix. The catalyzon gap "
            "Δ_cat = Δ0 − g0★ remains the bridge identification. ω0 = 2Δ "
            "is the phonon two-particle threshold, a distinct energy."
        ),
    }
