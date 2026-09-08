"""Self-consistent gap ↔ vacuum loop under locked rigor.

Gap equation (static / Matsubara ω_n=0 truncation):
    Δ = Δ₀ + ℏ c_s A_ZPF [ ln(Δ a_H / (2 ℏ c_s)) + γ_E ]
with Δ₀ = 2J.

When a_H = ℏ c_s / Δ (running healing), the log argument is exactly 1/2,
so
    Δ = Δ₀ + ℏ c_s A_ZPF ℓ_tree,   ℓ_tree = ln(1/2)+γ_E.

The locked *fold map* uses a renormalized IR coefficient ℓ (Corrections
ℓ≈−2.791) and retains ln δ:
    δ = 1 + α (ln δ + ℓ),   α ∝ A_ZPF.

A_ZPF is taken from the locked mode-enhancement construction
(``qvc.mode_enhancement.A_zpf_eff``) or from a fixed geometric vacuum
energy via α = E_vac / Δ₀ (dimensional bridge — provisional).

All meV outputs inherit vacuum-freeze provisional status until a lock
is promoted to ``frozen``.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any

from qvc.gap_fold import (
    ELL_CORRECTIONS,
    fold_condition,
    solve_fold_branches,
    tree_level_ell,
)
from qvc.mode_enhancement import A_zpf_eff
from qvc.params import PARAMS
from qvc.retardation import eta_ret
from qvc.vacuum_freeze import (
    FrozenVacuumSpec,
    bessel_mode_sum,
    evaluate_vacuum,
    vacuum_self_energy_corrected,
)


@dataclass(frozen=True)
class SelfConsistencyResult:
    geometry_label: str
    status: str
    Delta0_meV: float
    Delta_star_meV: float | None
    delta_star: float | None
    gap_reduction_pct: float | None
    A_zpf_um_inv: float | None
    alpha: float | None
    alpha_over_ac: float | None
    E_vac_tierC_meV: float
    x: float
    F: float
    ell: float
    eta_ret: float
    n_modes: int
    converged: bool
    iterations: int
    notes: str


def gap_from_A_algebraic(
    A_zpf: float,
    *,
    Delta0: float | None = None,
    hbar_cs: float | None = None,
) -> float:
    """Δ = Δ₀ + ℏ c_s A ℓ_tree with a_H = ℏ c_s/Δ (log arg = 1/2)."""
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    return Delta0 + hbar_cs * A_zpf * tree_level_ell()


def alpha_from_A(
    A_zpf: float,
    *,
    Delta0: float | None = None,
    hbar_cs: float | None = None,
) -> float:
    """α = (ℏ c_s A_ZPF) / Δ₀ — amplitude in fold-map units."""
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    return (hbar_cs * A_zpf) / Delta0


def iterate_running_healing_loop(
    R_um: float,
    *,
    hbar_cs: float | None = None,
    Delta0: float | None = None,
    n_modes: int | None = None,
    ell: float = ELL_CORRECTIONS,
    max_iter: int = 80,
    tol: float = 1e-10,
) -> dict[str, Any]:
    """Fixed-point: Δ → a_H(Δ) → A_ZPF(Δ) → fold branch δ(α).

    Uses Corrections ℓ for the fold map (not tree-level) so the
    self-consistency table matches the locked saddle-node programme.
    """
    hbar_cs = PARAMS.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    n_modes = PARAMS.n_modes_default if n_modes is None else n_modes
    delta_c, alpha_c = fold_condition(ell)

    Delta = Delta0
    history: list[dict[str, float]] = []
    converged = False
    for it in range(1, max_iter + 1):
        A = A_zpf_eff(Delta, hbar_cs, R_um, n_modes)
        alpha = alpha_from_A(A, Delta0=Delta0, hbar_cs=hbar_cs)
        br = solve_fold_branches(alpha, ell)
        if br["phys"] is None:
            history.append(
                {
                    "iter": float(it),
                    "Delta": Delta,
                    "A": A,
                    "alpha": alpha,
                    "delta": float("nan"),
                }
            )
            return {
                "converged": False,
                "iterations": it,
                "Delta_star": None,
                "delta_star": None,
                "A_zpf": A,
                "alpha": alpha,
                "alpha_c": alpha_c,
                "delta_c": delta_c,
                "folded_out": True,
                "history": history,
            }
        delta = float(br["phys"])
        Delta_new = delta * Delta0
        history.append(
            {
                "iter": float(it),
                "Delta": Delta_new,
                "A": A,
                "alpha": alpha,
                "delta": delta,
            }
        )
        if abs(Delta_new - Delta) < tol * max(1.0, abs(Delta0)):
            Delta = Delta_new
            converged = True
            break
        # Damped update for stability near fold
        Delta = 0.5 * Delta + 0.5 * Delta_new

    A = A_zpf_eff(Delta, hbar_cs, R_um, n_modes)
    alpha = alpha_from_A(A, Delta0=Delta0, hbar_cs=hbar_cs)
    return {
        "converged": converged,
        "iterations": len(history),
        "Delta_star": Delta,
        "delta_star": Delta / Delta0,
        "A_zpf": A,
        "alpha": alpha,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "alpha_over_ac": alpha / alpha_c if alpha_c else float("nan"),
        "folded_out": False,
        "history": history,
    }


def iterate_fixed_geometry_loop(
    E_vac_meV: float,
    *,
    Delta0: float | None = None,
    n_modes: int | None = None,
    ell: float = ELL_CORRECTIONS,
    enhance: bool = True,
) -> dict[str, Any]:
    """One-shot fold evaluation with fixed (r,R) vacuum energy.

    Bridge (locked dimensional map):
        α_bare = E_vac / Δ₀
        α_eff  = √n · E_vac / Δ₀   (independent-mode quadrature)

    Fold-in requires α ≤ α_c(ℓ). Default reports enhanced bridge; set
    ``enhance=False`` for the bare single-mode diagnostic.
    """
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    n_modes = PARAMS.n_modes_default if n_modes is None else n_modes
    delta_c, alpha_c = fold_condition(ell)
    factor = math.sqrt(n_modes) if enhance else 1.0
    alpha = factor * E_vac_meV / Delta0
    bridge = (
        "alpha = sqrt(n) * E_vac / Delta0"
        if enhance
        else "alpha = E_vac / Delta0 (bare)"
    )
    br = solve_fold_branches(alpha, ell)
    if br["phys"] is None:
        return {
            "converged": False,
            "iterations": 1,
            "Delta_star": None,
            "delta_star": None,
            "A_zpf": alpha * Delta0 / (PARAMS.hbar_cs_meV_um),
            "alpha": alpha,
            "alpha_c": alpha_c,
            "delta_c": delta_c,
            "alpha_over_ac": alpha / alpha_c if alpha_c else float("nan"),
            "folded_out": True,
            "bridge": bridge,
            "n_modes": n_modes,
            "enhance": enhance,
        }
    delta = float(br["phys"])
    return {
        "converged": True,
        "iterations": 1,
        "Delta_star": delta * Delta0,
        "delta_star": delta,
        "A_zpf": alpha * Delta0 / PARAMS.hbar_cs_meV_um,
        "alpha": alpha,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "alpha_over_ac": alpha / alpha_c,
        "folded_out": False,
        "bridge": bridge,
        "n_modes": n_modes,
        "enhance": enhance,
    }


def fold_in_E_vac_bound(
    *,
    Delta0: float | None = None,
    n_modes: int | None = None,
    ell: float = ELL_CORRECTIONS,
) -> dict[str, Any]:
    """Maximum E_vac that still admits a physical fold branch.

    α ≤ α_c  ⇒  E_vac ≤ α_c Δ₀ / √n   (enhanced)
               E_vac ≤ α_c Δ₀         (bare)
    """
    Delta0 = PARAMS.Delta0_meV if Delta0 is None else Delta0
    n_modes = PARAMS.n_modes_default if n_modes is None else n_modes
    delta_c, alpha_c = fold_condition(ell)
    return {
        "ell": ell,
        "alpha_c": alpha_c,
        "delta_c": delta_c,
        "Delta0_meV": Delta0,
        "n_modes": n_modes,
        "E_vac_max_bare_meV": alpha_c * Delta0,
        "E_vac_max_enhanced_meV": alpha_c * Delta0 / math.sqrt(n_modes),
        "notes": (
            "Package Tier-C on fat torus exceeds bare bound → folded out "
            "unless bridge or (ℏc_s,r,R,n) changes."
        ),
    }


def bridge_scan_row(E_vac_meV: float, geometry_label: str) -> dict[str, Any]:
    """Compare bare vs √n bridges for one E_vac."""
    bounds = fold_in_E_vac_bound()
    bare = iterate_fixed_geometry_loop(E_vac_meV, enhance=False)
    enh = iterate_fixed_geometry_loop(E_vac_meV, enhance=True)
    return {
        "geometry_label": geometry_label,
        "E_vac_tierC_meV": E_vac_meV,
        "E_vac_max_bare_meV": bounds["E_vac_max_bare_meV"],
        "E_vac_max_enhanced_meV": bounds["E_vac_max_enhanced_meV"],
        "bare": {
            "alpha_over_ac": bare["alpha_over_ac"],
            "Delta_star_meV": bare.get("Delta_star"),
            "gap_reduction_pct": (
                None
                if bare.get("Delta_star") is None
                else 100.0
                * (PARAMS.Delta0_meV - bare["Delta_star"])
                / PARAMS.Delta0_meV
            ),
            "folded_out": bare["folded_out"],
        },
        "enhanced_sqrt_n": {
            "n_modes": PARAMS.n_modes_default,
            "alpha_over_ac": enh["alpha_over_ac"],
            "Delta_star_meV": enh.get("Delta_star"),
            "gap_reduction_pct": (
                None
                if enh.get("Delta_star") is None
                else 100.0
                * (PARAMS.Delta0_meV - enh["Delta_star"])
                / PARAMS.Delta0_meV
            ),
            "folded_out": enh["folded_out"],
        },
    }


def evaluate_geometry_self_consistency(
    spec: FrozenVacuumSpec,
    *,
    n_modes: int | None = None,
    ell: float = ELL_CORRECTIONS,
    xi_um: float | None = None,
    mode: str = "fixed_geometry",
) -> SelfConsistencyResult:
    """Build one self-consistency row for a vacuum geometry candidate.

    mode:
      - ``fixed_geometry`` (default): α = √n E_vac/Δ₀ from frozen (r,R)
      - ``running_healing``: iterate A_zpf_eff(Δ) (usually exponentially small)
    """
    n_modes = PARAMS.n_modes_default if n_modes is None else n_modes
    xi_um = PARAMS.lambda_M_um if xi_um is None else xi_um
    ev = evaluate_vacuum(spec)
    if mode == "running_healing":
        loop = iterate_running_healing_loop(
            spec.R,
            hbar_cs=spec.hbar_cs,
            n_modes=n_modes,
            ell=ell,
        )
    else:
        loop = iterate_fixed_geometry_loop(
            ev["E_corrected_self_energy_meV"],
            n_modes=n_modes,
            ell=ell,
        )
    Delta0 = PARAMS.Delta0_meV
    Delta_star = loop.get("Delta_star")
    delta_star = loop.get("delta_star")
    gap_pct = None
    if Delta_star is not None:
        gap_pct = 100.0 * (Delta0 - Delta_star) / Delta0

    notes = (
        f"PROVISIONAL self-consistency under locked fold (ℓ={ell}), mode={mode}. "
        f"E_vac Tier-C on fixed (r,R); n={n_modes}. "
        "Do not quote meV as frozen until vacuum lock promotion."
    )
    return SelfConsistencyResult(
        geometry_label=spec.geometry_label or "unlabeled",
        status=spec.status,
        Delta0_meV=Delta0,
        Delta_star_meV=Delta_star,
        delta_star=delta_star,
        gap_reduction_pct=gap_pct,
        A_zpf_um_inv=loop.get("A_zpf"),
        alpha=loop.get("alpha"),
        alpha_over_ac=loop.get("alpha_over_ac"),
        E_vac_tierC_meV=ev["E_corrected_self_energy_meV"],
        x=ev["x"],
        F=ev["F"],
        ell=ell,
        eta_ret=eta_ret(Delta0, xi_um, spec.hbar_cs),
        n_modes=n_modes,
        converged=bool(loop.get("converged")),
        iterations=int(loop.get("iterations", 0)),
        notes=notes,
    )


def self_consistency_table(
    *,
    ell: float = ELL_CORRECTIONS,
    n_modes: int | None = None,
    mode: str = "fixed_geometry",
) -> list[dict[str, Any]]:
    """Compare preprint fat-torus vs Corrections moiré under identical locks."""
    p = PARAMS
    rows: list[SelfConsistencyResult] = []
    for label, r, R in [
        ("preprint_fat_torus", p.r_preprint_um, p.R_preprint_um),
        ("corrections_moire", p.r_corr_um, p.R_corr_um),
    ]:
        spec = FrozenVacuumSpec(
            formula_id="corrected_self_energy_2d",
            r=r,
            R=R,
            hbar_cs=p.hbar_cs_meV_um,
            geometry_label=label,
            status="provisional",
            notes="Self-consistency suite row",
        )
        rows.append(
            evaluate_geometry_self_consistency(
                spec, n_modes=n_modes, ell=ell, mode=mode
            )
        )
    return [asdict(r) for r in rows]


def catalyzon_dressed_gap(g0: float, Delta_ex: float | None = None) -> dict[str, float]:
    """Separate 52%-class mechanism: Δ_cat = Δ_ex + λ_min = Δ_ex - g₀."""
    Delta_ex = PARAMS.Delta0_meV if Delta_ex is None else Delta_ex
    Delta_cat = Delta_ex - g0
    return {
        "Delta_ex_meV": Delta_ex,
        "g0_meV": g0,
        "Delta_cat_meV": Delta_cat,
        "reduction_pct": 100.0 * g0 / Delta_ex,
        "mechanism": "Spec(G) lambda_min = -g0 (not vacuum freeze)",
    }
