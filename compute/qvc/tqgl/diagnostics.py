"""Optional Schwinger twin-channel and Berry GOE/GUE diagnostics.

These are *not* part of the GPE/DCE EOM. Berry Monte Carlo is a random-matrix
probe of level repulsion; it is separate from the locked democratic Spec(G)
where λ_min = −g₀ (multiplicity N−1).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from qvc.coupling_matrix import democratic_spectrum
from qvc.gij_tqgl_ws import run_tqgl_ws_G
from qvc.tqgl.locks import FrozenTQGL


def schwinger_twin_channel(lock: FrozenTQGL) -> dict[str, Any]:
    """Schwinger + thermal rates with Option A* E_vac and democratic |λ_min|.

    E_pair = 2 |λ_min| with |λ_min| = g₀ = λ★ E_vac (same participation as
    the fold operating point). The withdrawn 0.172 meV is not used.
    A_c,Sch = Δ₀ / (ℏ c_s) is a *different* scale from the fold A_c.
    """
    g0 = lock.lambda_star * lock.E_vac_meV
    spec = democratic_spectrum(lock.N_layers, g0)
    lam_min = float(spec["lambda_min"])
    E_pair = 2.0 * abs(lam_min)

    A_c_sch = lock.Delta0_meV / lock.hbar_cs_meV_um
    A_eff = lock.A_star_um_inv
    ratio = A_eff / A_c_sch if A_c_sch else float("nan")
    omega0 = 2.0 * lock.Delta_star_meV / lock.hbar_meV_ps
    sch_exp = np.pi / ratio if ratio > 0.0 else float("inf")
    Gamma_sch = float(omega0 * ratio * np.exp(-sch_exp)) if np.isfinite(sch_exp) else 0.0
    T_star = float(E_pair / (lock.kB_meV_per_K * sch_exp)) if np.isfinite(sch_exp) else float("nan")

    ratio_sweep = np.linspace(0.05, 0.95, 200)
    Gamma_sweep = omega0 * ratio_sweep * np.exp(-np.pi / ratio_sweep)
    T_arr = np.linspace(0.05, 0.80, 200)
    Gamma_th = omega0 * np.exp(-E_pair / (lock.kB_meV_per_K * T_arr))

    ws = run_tqgl_ws_G(n_grid=28)
    return {
        "E_vac_meV": lock.E_vac_meV,
        "g0_democratic_meV": g0,
        "lambda_min_democratic_meV": lam_min,
        "lambda_min_is_minus_g0": True,
        "E_pair_meV": E_pair,
        "E_pair_not_0p172_unless_rederived": abs(E_pair - 0.172) > 1e-6,
        "A_c_schwinger_um_inv": A_c_sch,
        "A_c_fold_um_inv": lock.A_c_fold_um_inv,
        "A_eff_um_inv": A_eff,
        "A_eff_over_Ac_schwinger": ratio,
        "A_star_over_Ac_fold": lock.A_star_um_inv / lock.A_c_fold_um_inv,
        "omega0_ps_inv": omega0,
        "schwinger_exponent": float(sch_exp),
        "Gamma_sch_ps_inv": Gamma_sch,
        "T_star_K": T_star,
        "ratio_sweep": ratio_sweep,
        "Gamma_sch_sweep": Gamma_sweep,
        "T_K": T_arr,
        "Gamma_thermal": Gamma_th,
        "Gamma_total": Gamma_sch + Gamma_th,
        "ws_G": {
            "lambda_min": ws["lambda_min"],
            "lambda_max": ws["lambda_max"],
            "g0_eff_mean_offdiag": ws["g0_eff_mean_offdiag"],
            "near_democratic": ws["acceptance"]["near_democratic"],
            "units": "overlap kernel, not meV",
            "note": "WS G is a shape check; meV |λ_min| uses democratic g0=λ★ E_vac.",
        },
        "note": (
            "Twin-channel rates use democratic |λ_min|=g0=λ★ E_vac. "
            "Fold A_c and Schwinger A_c=Δ₀/(ℏ c_s) are different objects."
        ),
    }


def berry_goe_gue(
    *,
    n_samples: int = 2000,
    n_modes: int = 5,
    g_bar: float = 0.05,
    seed: int = 42,
) -> dict[str, Any]:
    """GOE vs GUE Monte Carlo on zero-diagonal random matrices.

    Not the democratic G. Enhancement R is reported as computed.
    """
    rng = np.random.default_rng(seed)
    sigma = g_bar / np.sqrt(2.0)
    lmin_goe = np.empty(n_samples)
    lmin_gue = np.empty(n_samples)
    for i in range(n_samples):
        raw = rng.normal(0.0, sigma, (n_modes, n_modes))
        G_goe = 0.5 * (raw + raw.T) * np.sqrt(2.0)
        np.fill_diagonal(G_goe, 0.0)
        lmin_goe[i] = float(np.min(np.linalg.eigvalsh(G_goe)))

        mag = np.abs(raw)
        phi = rng.uniform(0.0, 2.0 * np.pi, (n_modes, n_modes))
        phi = phi - phi.T
        G_gue = mag * np.exp(1j * phi)
        G_gue = 0.5 * (G_gue + G_gue.conj().T)
        np.fill_diagonal(G_gue, 0.0)
        lmin_gue[i] = float(np.min(np.linalg.eigvalsh(G_gue)))

    abs_goe = np.abs(lmin_goe)
    abs_gue = np.abs(lmin_gue)
    R_mean = float(np.mean(abs_gue) / (np.mean(abs_goe) + 1e-18))
    R_med = float(np.median(abs_gue) / (np.median(abs_goe) + 1e-18))
    return {
        "n_samples": n_samples,
        "n_modes": n_modes,
        "g_bar_meV": g_bar,
        "lambda_min_GOE": lmin_goe,
        "lambda_min_GUE": lmin_gue,
        "mean_abs_GOE": float(np.mean(abs_goe)),
        "mean_abs_GUE": float(np.mean(abs_gue)),
        "R_mean": R_mean,
        "R_median": R_med,
        "separate_from_democratic_G": True,
        "not_forced_to_1p43": True,
        "note": (
            "R = ⟨|λ_min|⟩_GUE / ⟨|λ_min|⟩_GOE on random zero-diagonal ensembles. "
            "Democratic lock λ_min=−g0 is a different object and is not sampled here."
        ),
    }
