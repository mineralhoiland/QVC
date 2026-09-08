r"""Lindblad open-system TQGL dynamics (Phase Two, Workstream 5).

Models the gap-mode quantum state under dissipation:
  dρ/dt = -i[H_eff, ρ] + Σ_k γ_k (L_k ρ L_k† - ½{L_k†L_k, ρ})

Truncated to a finite Fock space (n_max ~ 20 quanta of the gap mode).

Physical scales (Option A*):
  ω₀ = 2Δ*/ℏ ~ 0.708 ps⁻¹  (gap oscillation)
  γ_loss ~ 1/30 ps⁻¹ = 0.033 ps⁻¹  (cavity photon loss)
  γ_deph ~ 1/50 ps⁻¹ = 0.020 ps⁻¹  (pure dephasing)

Key question: τ_formation ≪ τ_decoherence → gap survives.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any

import numpy as np
from scipy.linalg import expm


@dataclass
class LindbladConfig:
    """Configuration for Lindblad gap-mode simulation."""

    # Hilbert space truncation
    n_max: int = 20

    # Hamiltonian parameters (meV)
    Delta_star_meV: float = 0.2324  # Option A* steady gap
    g0_meV: float = 0.0982  # catalyzon coupling = anharmonicity scale
    hbar_meV_ps: float = 0.6582  # ℏ in meV·ps

    # Dissipation rates (ps⁻¹)
    gamma_loss: float = 0.033  # 1/τ_cav, cavity loss ~ 30 ps
    gamma_deph: float = 0.020  # 1/T₂*, pure dephasing ~ 50 ps
    gamma_thermal: float = 0.001  # thermal excitation at 20 mK (negligible)

    # DCE pump
    pump_amplitude_meV: float = 0.05  # parametric driving
    pump_detuning: float = 0.0  # resonant

    # Time stepping
    dt_ps: float = 0.1
    t_end_ps: float = 200.0

    @property
    def omega0(self) -> float:
        """Gap oscillation frequency in ps⁻¹."""
        return 2.0 * self.Delta_star_meV / self.hbar_meV_ps

    @property
    def n_steps(self) -> int:
        return int(self.t_end_ps / self.dt_ps)


def _annihilation(n: int) -> np.ndarray:
    """Annihilation operator â in Fock basis (n×n)."""
    a = np.zeros((n, n), dtype=complex)
    for i in range(n - 1):
        a[i, i + 1] = math.sqrt(i + 1)
    return a


def _number(n: int) -> np.ndarray:
    """Number operator n̂ = â†â."""
    return np.diag(np.arange(n, dtype=float))


def build_hamiltonian(cfg: LindbladConfig) -> np.ndarray:
    """Effective gap-mode Hamiltonian in meV.

    H = ω₀ n̂ + U n̂(n̂−1)/2 + V_pump (â² + â†²) / 2
    where U = g₀ (Kerr anharmonicity from catalyzon coupling).
    """
    n = cfg.n_max
    nhat = _number(n)
    a = _annihilation(n)
    adag = a.conj().T

    H = cfg.hbar_meV_ps * cfg.omega0 * nhat
    H += 0.5 * cfg.g0_meV * nhat @ (nhat - np.eye(n))

    # Parametric pump (DCE source term)
    H += 0.5 * cfg.pump_amplitude_meV * (a @ a + adag @ adag)

    return H


def build_lindblad_operators(cfg: LindbladConfig) -> list[tuple[float, np.ndarray]]:
    """Return list of (γ_k, L_k) pairs."""
    n = cfg.n_max
    a = _annihilation(n)
    nhat = _number(n)

    ops = []
    # Cavity loss: L₁ = â
    if cfg.gamma_loss > 0:
        ops.append((cfg.gamma_loss, a))

    # Pure dephasing: L₂ = n̂
    if cfg.gamma_deph > 0:
        ops.append((cfg.gamma_deph, nhat))

    # Thermal excitation: L₃ = â† (weak at 20 mK)
    if cfg.gamma_thermal > 0:
        ops.append((cfg.gamma_thermal, a.conj().T))

    return ops


def lindblad_rhs(
    rho: np.ndarray,
    H: np.ndarray,
    ops: list[tuple[float, np.ndarray]],
    hbar: float,
) -> np.ndarray:
    """Compute dρ/dt from Lindblad master equation."""
    # Unitary part: -i/ℏ [H, ρ]
    drho = -1j / hbar * (H @ rho - rho @ H)

    # Dissipative part: Σ γ_k (L ρ L† - ½{L†L, ρ})
    for gamma, L in ops:
        Ldag = L.conj().T
        LdL = Ldag @ L
        drho += gamma * (L @ rho @ Ldag - 0.5 * (LdL @ rho + rho @ LdL))

    return drho


def rk4_step(
    rho: np.ndarray,
    H: np.ndarray,
    ops: list[tuple[float, np.ndarray]],
    hbar: float,
    dt: float,
) -> np.ndarray:
    """Single RK4 step for density matrix evolution."""
    k1 = lindblad_rhs(rho, H, ops, hbar)
    k2 = lindblad_rhs(rho + 0.5 * dt * k1, H, ops, hbar)
    k3 = lindblad_rhs(rho + 0.5 * dt * k2, H, ops, hbar)
    k4 = lindblad_rhs(rho + dt * k3, H, ops, hbar)
    return rho + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)


def coherent_state(n_max: int, alpha: complex) -> np.ndarray:
    """Coherent state |α⟩ as density matrix."""
    psi = np.zeros(n_max, dtype=complex)
    for n in range(n_max):
        psi[n] = (alpha**n / math.sqrt(math.factorial(n))) * math.exp(
            -0.5 * abs(alpha) ** 2
        )
    psi /= np.linalg.norm(psi)
    return np.outer(psi, psi.conj())


def thermal_state(n_max: int, n_bar: float) -> np.ndarray:
    """Thermal density matrix with mean occupation n_bar."""
    probs = np.array(
        [(n_bar**k / (1 + n_bar) ** (k + 1)) for k in range(n_max)]
    )
    probs /= probs.sum()
    return np.diag(probs.astype(complex))


def gap_from_rho(rho: np.ndarray, Delta0: float, n_ref: float) -> float:
    """Extract effective gap Δ from ⟨n̂⟩.

    Model: Δ_eff = Δ₀ × (⟨n̂⟩ / n_ref)^{1/2} for condensate interpretation,
    or Δ_eff = Δ* × (1 - γ_loss/γ_pump) at steady state.

    Simplest: Δ ~ Δ₀ √(⟨n̂⟩ / n_ref) where n_ref = (Δ*/Δ₀)².
    """
    n_hat = _number(rho.shape[0])
    n_mean = max(float(np.real(np.trace(n_hat @ rho))), 0.0)
    if n_ref <= 0:
        return 0.0
    return Delta0 * math.sqrt(n_mean / n_ref)


def run_lindblad(cfg: LindbladConfig | None = None) -> dict[str, Any]:
    """Run Lindblad dynamics and return timeseries.

    Returns dict with:
      t_ps: time array
      n_mean: ⟨n̂⟩ vs time
      Delta_meV: effective gap vs time
      purity: Tr(ρ²) vs time
      entropy: von Neumann entropy vs time
      steady_state: final density matrix properties
    """
    if cfg is None:
        cfg = LindbladConfig()

    H = build_hamiltonian(cfg)
    ops = build_lindblad_operators(cfg)

    # Initial state: vacuum (no gap mode quanta)
    rho = np.zeros((cfg.n_max, cfg.n_max), dtype=complex)
    rho[0, 0] = 1.0

    n_hat = _number(cfg.n_max)
    n_ref = (cfg.Delta_star_meV / 0.600) ** 2  # reference occupation for Δ*

    # Storage
    n_steps = cfg.n_steps
    t_arr = np.linspace(0, cfg.t_end_ps, n_steps + 1)
    n_mean = np.zeros(n_steps + 1)
    Delta_arr = np.zeros(n_steps + 1)
    purity = np.zeros(n_steps + 1)
    entropy = np.zeros(n_steps + 1)

    def _vn_entropy(r: np.ndarray) -> float:
        evals = np.real(np.linalg.eigvalsh(r))
        evals = evals[evals > 1e-15]
        return -float(np.sum(evals * np.log(evals)))

    # Record IC
    n_mean[0] = float(np.real(np.trace(n_hat @ rho)))
    Delta_arr[0] = gap_from_rho(rho, 0.600, n_ref)
    purity[0] = float(np.real(np.trace(rho @ rho)))
    entropy[0] = _vn_entropy(rho)

    # Evolve
    for step in range(n_steps):
        rho = rk4_step(rho, H, ops, cfg.hbar_meV_ps, cfg.dt_ps)
        # Enforce Hermiticity and trace
        rho = 0.5 * (rho + rho.conj().T)
        tr = np.trace(rho)
        if abs(tr) > 1e-10:
            rho /= tr

        n_mean[step + 1] = float(np.real(np.trace(n_hat @ rho)))
        Delta_arr[step + 1] = gap_from_rho(rho, 0.600, n_ref)
        purity[step + 1] = float(np.real(np.trace(rho @ rho)))
        entropy[step + 1] = _vn_entropy(rho)

    # Steady-state analysis
    ss_n = n_mean[-1]
    ss_Delta = Delta_arr[-1]
    ss_purity = purity[-1]

    # Formation timescale: time to reach 90% of steady-state ⟨n⟩
    if ss_n > 0.01:
        threshold = 0.9 * ss_n
        idx = np.argmax(n_mean >= threshold)
        tau_form = float(t_arr[idx]) if idx > 0 else float("inf")
    else:
        tau_form = float("inf")

    return {
        "config": {
            "n_max": cfg.n_max,
            "Delta_star_meV": cfg.Delta_star_meV,
            "gamma_loss_ps_inv": cfg.gamma_loss,
            "gamma_deph_ps_inv": cfg.gamma_deph,
            "pump_meV": cfg.pump_amplitude_meV,
            "t_end_ps": cfg.t_end_ps,
        },
        "t_ps": t_arr.tolist(),
        "n_mean": n_mean.tolist(),
        "Delta_meV": Delta_arr.tolist(),
        "purity": purity.tolist(),
        "entropy": entropy.tolist(),
        "steady_state": {
            "n_mean": ss_n,
            "Delta_meV": ss_Delta,
            "Delta_over_Deltastar": ss_Delta / cfg.Delta_star_meV if cfg.Delta_star_meV > 0 else 0,
            "purity": ss_purity,
            "tau_formation_ps": tau_form,
        },
        "timescale_comparison": {
            "tau_formation_ps": tau_form,
            "tau_cavity_ps": 1.0 / cfg.gamma_loss if cfg.gamma_loss > 0 else float("inf"),
            "tau_dephasing_ps": 1.0 / cfg.gamma_deph if cfg.gamma_deph > 0 else float("inf"),
            "separation_ratio": (1.0 / cfg.gamma_loss) / max(tau_form, 0.01),
            "gap_survives": tau_form < (1.0 / cfg.gamma_loss),
        },
    }


def run_coherent_comparison(cfg: LindbladConfig | None = None) -> dict[str, Any]:
    """Run BOTH coherent (γ=0) and Lindblad dynamics for direct comparison."""
    if cfg is None:
        cfg = LindbladConfig()

    # Lindblad run
    result_open = run_lindblad(cfg)

    # Coherent run (no dissipation)
    cfg_closed = LindbladConfig(
        n_max=cfg.n_max,
        Delta_star_meV=cfg.Delta_star_meV,
        g0_meV=cfg.g0_meV,
        hbar_meV_ps=cfg.hbar_meV_ps,
        gamma_loss=0.0,
        gamma_deph=0.0,
        gamma_thermal=0.0,
        pump_amplitude_meV=cfg.pump_amplitude_meV,
        dt_ps=cfg.dt_ps,
        t_end_ps=cfg.t_end_ps,
    )
    result_closed = run_lindblad(cfg_closed)

    return {
        "open_system": result_open,
        "closed_system": result_closed,
        "comparison": {
            "Delta_open_final_meV": result_open["steady_state"]["Delta_meV"],
            "Delta_closed_final_meV": result_closed["steady_state"]["Delta_meV"],
            "ratio": (
                result_open["steady_state"]["Delta_meV"]
                / max(result_closed["steady_state"]["Delta_meV"], 1e-10)
            ),
            "purity_drop": 1.0 - result_open["steady_state"]["purity"],
        },
    }


def robustness_sweep(
    gamma_loss_range: np.ndarray | None = None,
) -> dict[str, Any]:
    """Sweep γ_loss to find critical decoherence rate where gap dies."""
    if gamma_loss_range is None:
        gamma_loss_range = np.logspace(-3, 0, 20)  # 0.001 to 1.0 ps⁻¹

    results = []
    for gl in gamma_loss_range:
        cfg = LindbladConfig(gamma_loss=float(gl), t_end_ps=150.0, dt_ps=0.2)
        res = run_lindblad(cfg)
        results.append({
            "gamma_loss": float(gl),
            "tau_cav_ps": 1.0 / float(gl),
            "Delta_ss_meV": res["steady_state"]["Delta_meV"],
            "n_ss": res["steady_state"]["n_mean"],
        })

    # Find critical gamma where Δ drops below 50% of Δ*
    Delta_star = 0.2324
    gammas = [r["gamma_loss"] for r in results]
    deltas = [r["Delta_ss_meV"] for r in results]
    gamma_crit = None
    for i in range(len(deltas) - 1):
        if deltas[i] >= 0.5 * Delta_star > deltas[i + 1]:
            # Linear interpolation
            frac = (0.5 * Delta_star - deltas[i + 1]) / (deltas[i] - deltas[i + 1])
            gamma_crit = gammas[i + 1] - frac * (gammas[i + 1] - gammas[i])
            break

    return {
        "sweep": results,
        "gamma_critical_ps_inv": gamma_crit,
        "tau_critical_ps": 1.0 / gamma_crit if gamma_crit else None,
        "interpretation": (
            f"Gap survives for τ_cav > {1.0/gamma_crit:.1f} ps"
            if gamma_crit
            else "Gap robust across full sweep range"
        ),
    }
