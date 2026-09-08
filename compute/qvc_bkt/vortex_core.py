"""C2: 2D GP vortex core energy and core cutoff a_c.

Relaxes a winding-1 Gross–Pitaevskii vortex (NOT a hopfion, NOT a CS fluxon)
in the truncated cubic-quintic GPE used by QVCCompute. Reports E_core/κ and
a_c. The healing-length identification a_c~ξ is labeled when used as a cutoff
in the IR logarithm.

    E_v = π κ ln(R / a_c) + E_core
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.integrate import solve_bvp

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc_bkt.locks import BKTLocks, build_bkt_locks
from qvc_bkt.stiffness import kappa_GL_mev_um2, n0_healing_disk_um2


def _xi_um(kappa_GL: float, Delta0: float) -> float:
    """Amplitude healing length ξ = √(κ_GL / Δ₀) = a_H for the GPE match."""
    return math.sqrt(kappa_GL / Delta0)


def _vortex_ode(kappa_GL: float, alpha_nl: float, beta_nl: float):
    def fun(r: np.ndarray, y: np.ndarray) -> np.ndarray:
        f, p = y
        r_safe = np.maximum(r, 1e-18)
        fpp = -p / r_safe + f / (r_safe**2) + (alpha_nl * f**3 + beta_nl * f**5) / kappa_GL
        return np.vstack((p, fpp))

    return fun


def relax_vortex_radial(
    *,
    kappa_GL: float,
    Delta0_meV: float,
    n_bulk: float = 1.0,
    R_over_xi: float = 40.0,
    n_grid: int = 400,
) -> dict[str, Any]:
    """Stationary winding-1 radial GP vortex, μ=0 at the cubic-quintic vacuum.

    GPE: [-κ_GL ∇² + α n + β n²] ψ with α=−Δ₀, β=Δ₀/n_bulk so the bulk
    density n_bulk is the potential minimum and μ_bulk=0.
    Status: computed.
    """
    alpha_nl = -Delta0_meV
    beta_nl = Delta0_meV / n_bulk
    xi = _xi_um(kappa_GL, Delta0_meV)
    R = R_over_xi * xi
    r_min = max(R / n_grid, 1e-6 * xi)
    r = np.linspace(r_min, R, n_grid)
    amp = math.sqrt(n_bulk)
    f_guess = amp * np.tanh(r / xi)

    fun = _vortex_ode(kappa_GL, alpha_nl, beta_nl)

    def bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
        # Regular vortex: f(r)~c r as r→0 ⇒ f(r_min) = r_min f'(r_min).
        return np.array([ya[0] - r_min * ya[1], yb[0] - amp])

    y0 = np.vstack((f_guess, np.gradient(f_guess, r)))
    sol = solve_bvp(fun, bc, r, y0, max_nodes=max(2 * n_grid, 800), tol=1e-6)
    if not sol.success:
        # Fallback: damped imaginary-time on the same mesh.
        f = f_guess.copy()
        dr = r[1] - r[0]
        dtau = 0.15 * dr * dr / max(kappa_GL, 1e-18)
        for _ in range(8000):
            fp = np.gradient(f, r)
            fpp = np.gradient(fp, r)
            rhs = kappa_GL * (fpp + fp / r - f / (r * r)) - alpha_nl * f**3 - beta_nl * f**5
            f = np.maximum(f + dtau * rhs, 0.0)
            f[0] = 0.0
            f[-1] = amp
        r_sol = r
        f_sol = f
        bvp_ok = False
        bvp_message = "fallback_imag_time"
    else:
        r_sol = r
        f_sol = sol.sol(r)[0]
        bvp_ok = True
        bvp_message = str(sol.message)

    n = np.maximum(f_sol, 0.0) ** 2
    fp = np.gradient(f_sol, r_sol)

    # Energy densities matching E = ∫ [κ_GL |∇ψ|² + (α/2) n² + (β/3) n³] dA
    e_grad = kappa_GL * (fp**2 + n / np.maximum(r_sol, 1e-18) ** 2)
    e_pot = 0.5 * alpha_nl * n**2 + (beta_nl / 3.0) * n**3
    e_pot_bulk = 0.5 * alpha_nl * (n_bulk**2) + (beta_nl / 3.0) * (n_bulk**3)
    integrand = 2.0 * math.pi * r_sol * (e_grad + e_pot - e_pot_bulk)
    trapz = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    E_gpe = float(trapz(integrand, r_sol))  # meV·µm²  (dimensionless |ψ|)

    # Half-density core radius (computed).
    target = 0.5 * n_bulk
    above = np.where(n >= target)[0]
    if above.size:
        a_half = float(r_sol[above[0]])
    else:
        a_half = float("nan")

    # 1-1/e recovery radius
    target_e = n_bulk * (1.0 - math.exp(-1.0))
    rec = np.where(n >= target_e)[0]
    a_rec = float(r_sol[rec[0]]) if rec.size else float("nan")

    return {
        "r_um": r_sol,
        "f": f_sol,
        "n": n,
        "xi_um": float(xi),
        "R_um": float(R),
        "E_v_gpe_meV_um2": E_gpe,
        "a_half_um": a_half,
        "a_recovery_um": a_rec,
        "n_bulk": float(n_bulk),
        "bvp_success": bool(bvp_ok),
        "bvp_message": bvp_message,
        "alpha_nl": float(alpha_nl),
        "beta_nl": float(beta_nl),
        "status": "computed",
        "defect": "2d_phase_vortex_winding_1",
        "not_hopfion": True,
        "not_cs_fluxon": True,
    }


def extract_core(
    radial: dict[str, Any],
    *,
    kappa_meV: float,
    n0_um2: float,
    a_c_um: float | None = None,
) -> dict[str, Any]:
    """Convert GPE vortex energy to E_v = πκ ln(R/a_c) + E_core.

    Physical energy E_true = n0 * E_gpe. If a_c is omitted, the healing
    length ξ is used and that cutoff is labeled ansatz.
    """
    E_true = n0_um2 * float(radial["E_v_gpe_meV_um2"])
    R = float(radial["R_um"])
    xi = float(radial["xi_um"])
    ac_ansatz = xi
    ac_computed = float(radial["a_half_um"])
    ac_used = float(a_c_um) if a_c_um is not None else ac_ansatz
    ac_label = "computed_half_density" if a_c_um is not None else "healing_length_ansatz"
    log_arg = max(R / max(ac_used, 1e-18), 1.0000001)
    E_log = math.pi * kappa_meV * math.log(log_arg)
    E_core = E_true - E_log
    ratio = E_core / kappa_meV if kappa_meV else float("nan")
    # Also report with computed a_half as cutoff.
    log_half = math.log(max(R / max(ac_computed, 1e-18), 1.0000001))
    E_core_half = E_true - math.pi * kappa_meV * log_half
    return {
        "E_v_meV": float(E_true),
        "E_log_meV": float(E_log),
        "E_core_meV": float(E_core),
        "E_core_over_kappa": float(ratio),
        "kappa_meV": float(kappa_meV),
        "a_c_used_um": float(ac_used),
        "a_c_used_label": ac_label,
        "a_c_healing_um": float(xi),
        "a_c_half_um": ac_computed,
        "a_c_recovery_um": float(radial["a_recovery_um"]),
        "E_core_over_kappa_using_a_half": float(E_core_half / kappa_meV if kappa_meV else float("nan")),
        "R_um": R,
        "ln_R_over_ac": float(math.log(log_arg)),
        "n0_um2": float(n0_um2),
        "status": {
            "profile": "computed",
            "a_c_healing": "ansatz" if a_c_um is None else "override",
            "a_c_half": "computed",
            "E_core": "computed_given_cutoff",
        },
    }


def run_c2(
    locks: BKTLocks | None = None,
    *,
    kappa_uv_meV: float | None = None,
    n_grid: int = 400,
    R_over_xi: float = 40.0,
) -> dict[str, Any]:
    """C2 driver. Uses UV conventional+geom κ unless overridden (from C1)."""
    locks = locks or build_bkt_locks()
    k_gl = kappa_GL_mev_um2(locks.tqgl)
    dens = n0_healing_disk_um2(locks.r_um)
    n0 = float(dens["n0_um2"])
    radial = relax_vortex_radial(
        kappa_GL=k_gl,
        Delta0_meV=locks.Delta0_meV,
        n_bulk=1.0,
        R_over_xi=R_over_xi,
        n_grid=n_grid,
    )
    if kappa_uv_meV is None:
        # Conventional UV only for the log coefficient of this GP vortex
        # (the geometric floor is a band property, not in the GP functional).
        kappa_gp = 2.0 * k_gl * n0 * 1.0  # ⟨ρ²⟩=1 at UV bulk
        kappa_source = "GP_conventional_two_sector_UV"
    else:
        kappa_gp = float(kappa_uv_meV)
        kappa_source = "caller"
    core = extract_core(radial, kappa_meV=kappa_gp, n0_um2=n0)
    # Fit κ from the IR tail: E_phase(r>r0) ~ πκ ln, using the GP integrand.
    r = radial["r_um"]
    n = radial["n"]
    xi = radial["xi_um"]
    mask = r > 8.0 * xi
    if np.count_nonzero(mask) > 8:
        e_phase_gpe = k_gl * (n / np.maximum(r, 1e-18) ** 2)
        # cumulative from outside: not needed; slope of remaining log
        kappa_tail = 2.0 * n0 * k_gl * float(np.mean(n[mask]))
    else:
        kappa_tail = kappa_gp

    return {
        "task": "C2",
        "defect": "2d_phase_vortex_winding_1",
        "kappa_for_log_meV": float(kappa_gp),
        "kappa_source": kappa_source,
        "kappa_IR_tail_meV": float(kappa_tail),
        "radial": {
            k: radial[k]
            for k in (
                "xi_um",
                "R_um",
                "E_v_gpe_meV_um2",
                "a_half_um",
                "a_recovery_um",
                "bvp_success",
                "bvp_message",
                "status",
                "defect",
                "not_hopfion",
                "not_cs_fluxon",
            )
        },
        "profile": {"r_um": radial["r_um"], "n": radial["n"], "f": radial["f"]},
        "core": core,
        "closures": {
            "EOM": "truncated GPE (no CS / Berry / Skyrme)",
            "a_c_healing": "ansatz unless replaced by computed half-density radius",
            "n0": dens,
            "hopfion": "not used as the KT vortex",
        },
        "E_core_over_kappa": core["E_core_over_kappa"],
        "a_c_um": core["a_c_used_um"],
        "a_c_label": core["a_c_used_label"],
    }
