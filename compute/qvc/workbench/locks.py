"""S0 — Freeze the experiment.

Option A* locks are *loaded* from the existing solver chain
(:func:`qvc.tqgl.locks.build_frozen_tqgl` → bridge → fold → gap_reduction),
never re-fitted here. Withdrawn numbers live in :data:`DEPRECATED` and any
attempt to display one raises.

Three epistemic layers from the QVC v4 preprint and Hopfion Paper C:

L1  laboratory locks (fold, E_vac, Spec(G), Q_H=1/N, Hall step)  — solver truth
L2  analogue geometry (acoustic metric, torsion proxy, Goldstones) — diagnosed
L3  spacetime ontology                                              — schematic only

Unit conventions used by every viewport
---------------------------------------
atomic crystal   Å            (pymatgen / CIF / LAMMPS)
mesoscale fields nm           (hopfion, torus, Casimir envelope)
GPE / DCE        µm, ps, meV  (existing ``qvc.tqgl`` units, ℏ = 0.6582 meV·ps)
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from typing import Any

from qvc.coupling_matrix import democratic_spectrum
from qvc.lens_chern_simons import anyonic_exchange_angle, fractional_hall_charge
from qvc.self_consistency import catalyzon_dressed_gap
from qvc.tqgl.locks import FrozenTQGL, build_frozen_tqgl

E2_OVER_H_SIEMENS = 3.87404614e-5  # conductance quantum e²/h


class Layer(str, Enum):
    """Epistemic layer stamped on every viewport."""

    L1 = "L1"
    L2 = "L2"
    L3 = "L3"

    @property
    def title(self) -> str:
        return {
            Layer.L1: "L1 · laboratory lock (QVC / TEVC)",
            Layer.L2: "L2 · analogue geometry (TDVT) — diagnosed, not Einstein",
            Layer.L3: "L3 · spacetime ontology — schematic, not locked",
        }[self]


def stamp(layer: Layer, text: str) -> str:
    """Return a viewport caption prefixed with its layer stamp."""
    return f"[{layer.value}] {text}"


# --------------------------------------------------------------------------
# Reject list (docs/PHASE_TWO_PLAN.md “DEPRECATED (do NOT use)”)
# --------------------------------------------------------------------------

DEPRECATED: dict[str, dict[str, Any]] = {
    "E_Cas_0p244_meV": {
        "value": 0.244,
        "unit": "meV",
        "why": "withdrawn preprint Casimir number; Option A* gives E_vac = 0.1287 meV",
    },
    "gap_reduction_92pct": {
        "value": 92.0,
        "unit": "%",
        "why": "δ = 0.032/0.6 lies below the fold δ_c ≈ 0.182; max static reduction is 81.8%",
    },
    "Delta_1ps_0p032_meV": {
        "value": 0.032,
        "unit": "meV",
        "why": "scheduled (not ψ-derived) live gap from the withdrawn solver",
    },
    "Delta_QVC_0p389_meV": {
        "value": 0.389,
        "unit": "meV",
        "why": "35% static loop at E_Cas = 0.244 meV",
    },
    "T_star_285_290_mK": {
        "value": 285.0,
        "unit": "mK",
        "why": "superseded by T* = 118.7 mK from E_pair = 0.1963 meV",
    },
    "A_ZPF_xi_0p296": {
        "value": 0.296,
        "unit": "dimensionless",
        "why": "A_ZPF·ξ is not the coupling constant; bridge is α = λ E_vac/Δ₀",
    },
    "lambda_min_minus_4g0": {
        "value": -4.0,
        "unit": "g₀",
        "why": "democratic G = g₀(J−I) has λ_min = −g₀ (multiplicity N−1); −4g₀/−5g₀ are false",
    },
    "lambda_min_minus_5g0": {
        "value": -5.0,
        "unit": "g₀",
        "why": "see lambda_min_minus_4g0",
    },
}


class DeprecatedNumberError(ValueError):
    """Raised when a withdrawn preprint number is about to be displayed."""


def assert_not_deprecated(value: float, *, unit: str, rel_tol: float = 1e-3) -> None:
    """Raise if ``value`` (in ``unit``) matches a withdrawn folklore number."""
    for key, rec in DEPRECATED.items():
        if rec["unit"] != unit:
            continue
        ref = float(rec["value"])
        if math.isclose(float(value), ref, rel_tol=rel_tol, abs_tol=1e-12):
            raise DeprecatedNumberError(
                f"{value} {unit} matches withdrawn number '{key}': {rec['why']}"
            )


# --------------------------------------------------------------------------
# Frozen workbench locks
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class WorkbenchLocks:
    """Option A* numbers consumed by every viewport (all L1)."""

    tqgl: FrozenTQGL
    # catalyzon (democratic channel, parallel to the fold — never stacked)
    g0_star_meV: float
    lambda_max_meV: float
    lambda_min_meV: float
    lambda_min_multiplicity: int
    Delta_cat_meV: float
    catalyzon_reduction_pct: float
    # lens-space topology
    Q_H: float
    theta_anyon: float
    hall_step_e2_over_h: float
    hall_step_siemens: float
    hall_base_e2_over_h: float
    chern_simons_level_k: int
    # Schwinger diagnostic (documented lock, not recomputed here)
    T_star_mK: float = 118.7
    E_pair_meV: float = 0.1963
    # resonance
    omega0_star_ps_inv: float = field(default=0.0)
    nu0_star_THz: float = field(default=0.0)

    # ---- convenience -----------------------------------------------------
    @property
    def N(self) -> int:
        return self.tqgl.N_layers

    @property
    def Delta0_meV(self) -> float:
        return self.tqgl.Delta0_meV

    @property
    def Delta_star_meV(self) -> float:
        return self.tqgl.Delta_star_meV

    @property
    def E_vac_meV(self) -> float:
        return self.tqgl.E_vac_meV

    @property
    def lambda_star(self) -> float:
        return self.tqgl.lambda_star

    @property
    def r_nm(self) -> float:
        return self.tqgl.r_nm

    @property
    def R_nm(self) -> float:
        return self.tqgl.R_nm

    def table(self) -> list[tuple[str, str, str]]:
        """(symbol, value, layer) rows for the lock panel."""
        t = self.tqgl
        rows = [
            ("(r, R)", f"({t.r_nm:.0f}, {t.R_nm:.0f}) nm", "L1"),
            ("ℏc_s", f"{t.hbar_cs_meV_um:.3f} meV·µm (provisional)", "L1"),
            ("Δ₀ = 2J", f"{t.Delta0_meV:.3f} meV", "L1"),
            ("F(2πr/R)", f"{t.F:.5f}", "L1"),
            ("E_vac", f"{t.E_vac_meV:.4f} meV", "L1"),
            ("ℓ", f"{t.ell:.4f}", "L1"),
            ("α_c = δ_c", f"{t.alpha_c:.4f}", "L1"),
            ("Δ_c (fold)", f"{t.Delta_c_meV:.3f} meV ({t.max_static_reduction_pct:.1f}% max static)", "L1"),
            ("λ_fold", f"{t.lambda_fold:.4f}", "L1"),
            ("λ*", f"{t.lambda_star:.4f}", "L1"),
            ("Δ*", f"{t.Delta_star_meV:.4f} meV ({t.gap_reduction_star_pct:.1f}% static)", "L1"),
            ("g₀* = λ*E_vac", f"{self.g0_star_meV:.4f} meV", "L1"),
            ("Spec(G)", f"λ_max = {self.lambda_max_meV:.4f}, λ_min = {self.lambda_min_meV:.4f} meV (×{self.lambda_min_multiplicity})", "L1"),
            ("Δ_cat = Δ_ex − g₀*", f"{self.Delta_cat_meV:.4f} meV ({self.catalyzon_reduction_pct:.1f}%, parallel channel)", "L1"),
            ("Q_H = 1/N", f"{self.Q_H:.3f}", "L1"),
            ("θ_anyon = π/N", f"{self.theta_anyon:.4f} rad", "L1"),
            ("δσ_xy = e²/(2Nh)", f"{self.hall_step_e2_over_h:.3f} e²/h = {self.hall_step_siemens:.3e} S", "L1"),
            ("k_CS = 2N", f"{self.chern_simons_level_k}", "L1"),
            ("A*/A_c^fold", "reported by DCE run (not imposed)", "L1"),
            ("T*", f"{self.T_star_mK:.1f} mK (Schwinger diagnostic)", "L1"),
            ("ω₀ = 2Δ*/ℏ", f"{self.omega0_star_ps_inv:.3f} ps⁻¹ (ν₀ = {self.nu0_star_THz:.3f} THz)", "L1"),
        ]
        return rows

    def as_dict(self) -> dict[str, Any]:
        d = {k: v for k, v in self.__dict__.items() if k != "tqgl"}
        d["tqgl"] = self.tqgl.as_dict()
        return d


@lru_cache(maxsize=1)
def get_locks() -> WorkbenchLocks:
    """Build (once) the frozen workbench locks from the existing solver chain."""
    t = build_frozen_tqgl()
    N = t.N_layers
    g0 = t.lambda_star * t.E_vac_meV
    spec = democratic_spectrum(N, g0)
    cat = catalyzon_dressed_gap(g0, t.Delta0_meV)
    Q_H = fractional_hall_charge(N)
    step = 1.0 / (2.0 * N)
    omega0 = 2.0 * t.Delta_star_meV / t.hbar_meV_ps  # ps⁻¹
    locks = WorkbenchLocks(
        tqgl=t,
        g0_star_meV=g0,
        lambda_max_meV=spec["lambda_max"],
        lambda_min_meV=spec["lambda_min"],
        lambda_min_multiplicity=spec["lambda_min_mult"],
        Delta_cat_meV=cat["Delta_cat_meV"],
        catalyzon_reduction_pct=cat["reduction_pct"],
        Q_H=Q_H,
        theta_anyon=anyonic_exchange_angle(N),
        hall_step_e2_over_h=step,
        hall_step_siemens=step * E2_OVER_H_SIEMENS,
        hall_base_e2_over_h=float(N),
        chern_simons_level_k=2 * N,
        omega0_star_ps_inv=omega0,
        nu0_star_THz=omega0 / (2.0 * math.pi),  # ps⁻¹ / 2π = THz
    )
    # Guard: the frozen numbers must not coincide with withdrawn folklore.
    assert_not_deprecated(t.E_vac_meV, unit="meV")
    assert_not_deprecated(t.Delta_star_meV, unit="meV")
    assert_not_deprecated(t.gap_reduction_star_pct, unit="%")
    assert_not_deprecated(locks.lambda_min_meV / g0, unit="g₀")
    return locks


def check_against_phase_two() -> dict[str, Any]:
    """S9 acceptance: compare frozen locks to docs/PHASE_TWO_PLAN.md table."""
    L = get_locks()
    t = L.tqgl
    expected = {
        "F": (0.58068, 1e-4),
        "E_vac_meV": (0.1287, 2e-3),
        "ell": (-2.7958, 2e-3),
        "alpha_c": (0.1818, 2e-3),
        "lambda_fold": (0.8475, 2e-3),
        "lambda_star": (0.7628, 2e-3),
        "Delta_star_meV": (0.2324, 2e-3),
        "Delta_c_meV": (0.109, 1e-2),
        "g0_star_meV": (0.0982, 2e-3),
        "lambda_max_meV": (0.3928, 2e-3),
        "lambda_min_meV": (-0.0982, 2e-3),
        "Q_H": (0.2, 1e-9),
        "hall_step_siemens": (3.874e-6, 1e-3),
    }
    got = {
        "F": t.F,
        "E_vac_meV": t.E_vac_meV,
        "ell": t.ell,
        "alpha_c": t.alpha_c,
        "lambda_fold": t.lambda_fold,
        "lambda_star": t.lambda_star,
        "Delta_star_meV": t.Delta_star_meV,
        "Delta_c_meV": t.Delta_c_meV,
        "g0_star_meV": L.g0_star_meV,
        "lambda_max_meV": L.lambda_max_meV,
        "lambda_min_meV": L.lambda_min_meV,
        "Q_H": L.Q_H,
        "hall_step_siemens": L.hall_step_siemens,
    }
    rows = {}
    ok = True
    for k, (ref, rtol) in expected.items():
        val = got[k]
        passed = math.isclose(val, ref, rel_tol=rtol, abs_tol=1e-12)
        ok &= passed
        rows[k] = {"expected": ref, "got": val, "pass": passed}
    rows["lambda_min_multiplicity"] = {
        "expected": L.N - 1,
        "got": L.lambda_min_multiplicity,
        "pass": L.lambda_min_multiplicity == L.N - 1,
    }
    ok &= rows["lambda_min_multiplicity"]["pass"]
    return {"all_pass": bool(ok), "rows": rows}


UNITS = {
    "atomic": "Å",
    "mesoscale": "nm",
    "gpe_length": "µm",
    "gpe_time": "ps",
    "energy": "meV",
    "hbar": "0.6582119569 meV·ps",
}
