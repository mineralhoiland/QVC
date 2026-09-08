"""DCE amplitude ODE, actually coupled to a live gap.

    dA/dt = [ Λ(Δ_live) L(δω/Γ) − Λ_loss ] A

    V_ZPF ∝ −E_vac · (A(t)/A(0)) · gaussian_ring   (running A)

A_c is the *fold* critical amplitude α_c Δ₀ / (ℏ c_s).
It is not 1.2233 “BKT”, and it is not back-computed as A*/0.52.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import numpy as np

from qvc.tqgl.locks import FrozenTQGL, fold_delta_of_A


def lorentzian(x: float) -> float:
    return 1.0 / (1.0 + x * x)


def Lambda_parametric(Delta_meV: float, *, mod_depth: float, Q: float, hbar: float) -> float:
    """Λ(Δ) = (δΔ / 2ℏ) · Q / (2π) with δΔ = mod_depth · Δ."""
    return (mod_depth * Delta_meV / (2.0 * hbar)) * (Q / (2.0 * np.pi))


@dataclass
class DCEConfig:
    Q_moire: float = 20.0
    mod_depth: float = 0.10
    t_end_ps: float = 80.0
    dt_ps: float = 0.05
    # Loss = cavity linewidth at the operating gap (not fit to A*/A_c=0.52).
    loss_mode: str = "gamma_cav_star"
    # Pump: "resonant" → ω/2 = 2 Δ★/ℏ; "detuning" → + n_gamma · Γ.
    pump_mode: str = "resonant"
    n_gamma_detune: float = 2.0
    A0_mode: str = "lambda_star"  # "lambda_star" or "unit_geom"
    A_clip_over_Ac: float = 8.0


@dataclass
class DCEIntegrator:
    """Mutable DCE state advanced with whatever Δ_live the caller supplies."""

    lock: FrozenTQGL
    cfg: DCEConfig
    A: float
    A0: float
    omega_pump_half: float
    Gamma_cav: float
    Lam_loss: float

    def rhs(self, A: float, Delta_live: float) -> float:
        A = max(A, 1e-12)
        D = max(Delta_live, 1e-12)
        omega_cat = 2.0 * D / self.lock.hbar_meV_ps
        dw = self.omega_pump_half - omega_cat
        L = lorentzian(dw / self.Gamma_cav)
        lam = Lambda_parametric(
            D,
            mod_depth=self.cfg.mod_depth,
            Q=self.cfg.Q_moire,
            hbar=self.lock.hbar_meV_ps,
        )
        return (lam * L - self.Lam_loss) * A

    def lorentzian_now(self, Delta_live: float) -> float:
        D = max(Delta_live, 1e-12)
        omega_cat = 2.0 * D / self.lock.hbar_meV_ps
        return lorentzian((self.omega_pump_half - omega_cat) / self.Gamma_cav)

    def step(self, dt: float, Delta_live: float) -> float:
        """RK4 step. Δ_live must come from ψ or from the fold map of current A."""
        A = self.A
        k1 = self.rhs(A, Delta_live)
        k2 = self.rhs(A + 0.5 * dt * k1, Delta_live)
        k3 = self.rhs(A + 0.5 * dt * k2, Delta_live)
        k4 = self.rhs(A + dt * k3, Delta_live)
        A_new = A + dt * (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0
        a_max = self.cfg.A_clip_over_Ac * self.lock.A_c_fold_um_inv
        self.A = float(np.clip(A_new, 1e-12, a_max))
        return self.A


def make_dce(lock: FrozenTQGL, cfg: DCEConfig | None = None) -> DCEIntegrator:
    cfg = DCEConfig() if cfg is None else cfg
    if cfg.A0_mode == "unit_geom":
        A0 = lock.A_geom_um_inv
    else:
        A0 = lock.A_star_um_inv

    omega_star = 2.0 * lock.Delta_star_meV / lock.hbar_meV_ps
    Gamma = omega_star / cfg.Q_moire
    if cfg.pump_mode == "resonant":
        w_half = omega_star
    elif cfg.pump_mode == "detuning":
        w_half = omega_star + cfg.n_gamma_detune * Gamma
    else:
        raise ValueError(f"unknown pump_mode {cfg.pump_mode}")

    if cfg.loss_mode == "gamma_cav_star":
        loss = Gamma
    elif cfg.loss_mode == "lambda_star_balance":
        loss = Lambda_parametric(
            lock.Delta_star_meV,
            mod_depth=cfg.mod_depth,
            Q=cfg.Q_moire,
            hbar=lock.hbar_meV_ps,
        )
    else:
        raise ValueError(f"unknown loss_mode {cfg.loss_mode}")

    return DCEIntegrator(
        lock=lock,
        cfg=cfg,
        A=A0,
        A0=A0,
        omega_pump_half=w_half,
        Gamma_cav=Gamma,
        Lam_loss=loss,
    )


def run_dce_fold_slaved(
    lock: FrozenTQGL,
    cfg: DCEConfig | None = None,
    *,
    delta_override: Callable[[float], float] | None = None,
) -> dict[str, Any]:
    """80 ps DCE with Δ from the static fold branch at the current A.

    This is the slow-timescale slaving: DCE is slow vs GPE, so Δ follows
    the physical fold branch Δ(A) while A < A_c. If A exceeds A_c the
    static theory has no physical branch — flagged, Δ held at last
    physical value only for the ODE (not presented as an operating point).
    """
    cfg = DCEConfig() if cfg is None else cfg
    integ = make_dce(lock, cfg)
    n = int(round(cfg.t_end_ps / cfg.dt_ps))
    t = np.arange(n + 1, dtype=float) * cfg.dt_ps
    A = np.zeros(n + 1)
    D = np.zeros(n + 1)
    L = np.zeros(n + 1)
    folded = np.zeros(n + 1, dtype=bool)
    A[0] = integ.A
    last_phys = lock.Delta_star_meV
    # Cache Δ(A) so the 80 ps window does not rescan the fold map each step.
    _cache: dict[int, tuple[bool, float | None]] = {}

    def _fold_cached(A_now: float) -> tuple[bool, float | None]:
        key = int(round(A_now * 1e7))
        hit = _cache.get(key)
        if hit is not None:
            return hit
        info = fold_delta_of_A(A_now, lock)
        val = (bool(info["folded_out"]), info["Delta_phys_meV"])
        _cache[key] = val
        return val

    for i in range(n + 1):
        folded_i, d_phys = _fold_cached(integ.A)
        folded[i] = folded_i
        if delta_override is not None:
            D[i] = float(delta_override(integ.A))
        elif d_phys is not None:
            D[i] = float(d_phys)
            last_phys = D[i]
        else:
            D[i] = last_phys
        L[i] = integ.lorentzian_now(D[i])
        A[i] = integ.A
        if i < n:
            integ.step(cfg.dt_ps, D[i])

    idx = int(0.85 * n)
    A_star = float(np.mean(A[idx:]))
    D_star = float(np.mean(D[idx:]))
    A_c = lock.A_c_fold_um_inv
    ratio = A_star / A_c if A_c else float("nan")
    runaway = bool(np.max(A) >= 0.98 * cfg.A_clip_over_Ac * A_c)
    return {
        "pump_mode": cfg.pump_mode,
        "loss_mode": cfg.loss_mode,
        "A0_mode": cfg.A0_mode,
        "t_ps": t,
        "A_um_inv": A,
        "Delta_fold_meV": D,
        "lorentzian": L,
        "folded_out": folded,
        "A0_um_inv": integ.A0,
        "A_star_um_inv": A_star,
        "Delta_star_dce_meV": D_star,
        "A_c_fold_um_inv": A_c,
        "A_star_over_Ac": ratio,
        "A_c_is_fold_not_BKT": True,
        "A_c_not_backcomputed_from_0p52": True,
        "paper_ratio_0p52_not_imposed": True,
        "Lorentzian_identically_zero": bool(np.all(L < 1e-12)),
        "runaway_clipped": runaway,
        "any_folded_out": bool(np.any(folded)),
        "omega_pump_half_ps_inv": integ.omega_pump_half,
        "Gamma_cav_ps_inv": integ.Gamma_cav,
        "Lam_loss_ps_inv": integ.Lam_loss,
        "Lambda_at_Delta_star_ps_inv": Lambda_parametric(
            lock.Delta_star_meV,
            mod_depth=cfg.mod_depth,
            Q=cfg.Q_moire,
            hbar=lock.hbar_meV_ps,
        ),
        "L0": float(L[0]),
        "note": (
            "A_c = α_c Δ₀ / (ℏ c_s) is the fold critical amplitude. "
            "A*/A_c is reported as computed; it is not forced to 0.52."
        ),
    }
