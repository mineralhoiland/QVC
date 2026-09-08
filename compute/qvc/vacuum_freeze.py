"""Vacuum spectral problem freeze API — single spectral problem lock.

PROVEN vs PROVISIONAL
---------------------
Proven (Tier C form):
    E_vac = (1 / 4π²) (ħ c_s / r) F(2π r / R)
    F(x) = Σ_{n=1}^∞ K_0(n x)

Audit-failed (old Casimir preprint formula):
    E_old = (1 / 4π) (ħ c_s / r) (R/r) F
    → does NOT reproduce the claimed 0.244 meV on preprint (r,R).

Geometry and which formula to quote in papers remain OPEN until
``freeze(spec)`` writes a lockfile under ``qvc/locks/``.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from scipy.special import k0

from qvc.params import PARAMS

LOCKS_DIR = Path(__file__).resolve().parent / "locks"



@dataclass(frozen=True)
class FrozenVacuumSpec:
    """One concrete vacuum spectral problem (formula + geometry + scale).

    Attributes
    ----------
    formula_id :
        Identifier of the energy formula. Use
        ``corrected_self_energy_2d`` for the locked form, or
        ``casimir_preprint_old`` for the audit-failed formula.
    r, R :
        Geometric radii in µm (tube/healing radius, torus major radius).
    hbar_cs :
        ħ c_s in meV·µm.
    F_definition :
        How F(x) is defined; default Bessel K0 mode sum.
    notes :
        Free-text provenance / provisional markings.
    """

    formula_id: str
    r: float
    R: float
    hbar_cs: float
    F_definition: str = "bessel_K0_sum"
    notes: str = ""
    status: str = "provisional"
    geometry_label: str = ""
    nmax: int = 500


def bessel_mode_sum(x: float, nmax: int = 500) -> float:
    """F(x) = Σ_{n=1}^{nmax} K_0(n x)."""
    if x <= 0:
        raise ValueError("x must be positive")
    return float(sum(k0(n * x) for n in range(1, nmax + 1)))


def casimir_preprint_formula(hbar_cs: float, r: float, R: float, F: float) -> float:
    """OLD (audit-failed) formula: (1/4π)(ħc_s/r)(R/r) F."""
    return (1.0 / (4.0 * math.pi)) * (hbar_cs / r) * (R / r) * F


def vacuum_self_energy_corrected(hbar_cs: float, r: float, F: float) -> float:
    """CORRECTED phonon vacuum self-energy (2D phase-space measure):
    E_vac = (1 / 4π²) (ħ c_s / r) F(2π r / R).
    """
    return (1.0 / (4.0 * math.pi**2)) * (hbar_cs / r) * F


def jacobi_theta3_alternating(alpha_c: float, nmax: int = 200) -> float:
    """θ₃(π, e^{-α}) = Σ (-1)^n exp(-α n²) — NOT equal to F≈0.596."""
    s = 1.0
    for n in range(1, nmax + 1):
        s += 2.0 * ((-1) ** n) * math.exp(-alpha_c * n * n)
    return s


def candidate_geometries() -> dict[str, dict[str, float]]:
    """Two candidate geometries; do not mix silently."""
    p = PARAMS
    return {
        "preprint_fat_torus": {
            "r_um": p.r_preprint_um,
            "R_um": p.R_preprint_um,
            "r_nm": p.r_preprint_um * 1e3,
            "R_nm": p.R_preprint_um * 1e3,
        },
        "corrections_moire": {
            "r_um": p.r_corr_um,
            "R_um": p.R_corr_um,
            "r_nm": p.r_corr_um * 1e3,
            "R_nm": p.R_corr_um * 1e3,
        },
    }


def evaluate_vacuum(spec: FrozenVacuumSpec) -> dict[str, Any]:
    """Compute F, old Casimir, and corrected self-energy for a spec."""
    x = 2.0 * math.pi * spec.r / spec.R
    F = bessel_mode_sum(x, nmax=spec.nmax)
    E_old = casimir_preprint_formula(spec.hbar_cs, spec.r, spec.R, F)
    E_corr = vacuum_self_energy_corrected(spec.hbar_cs, spec.r, F)
    return {
        "x": x,
        "F": F,
        "E_old_casimir_meV": E_old,
        "E_corrected_self_energy_meV": E_corr,
        "claimed_preprint_meV": 0.244,
        "old_reproduces_0.244": abs(E_old - 0.244) < 0.05,
    }


def assert_preprint_0p244_not_reproduced(
    r: float | None = None,
    R: float | None = None,
    hbar_cs: float | None = None,
    tol: float = 0.05,
) -> dict[str, Any]:
    """Explicitly FAIL the false claim that old formula → 0.244 meV.

    Raises AssertionError if the old formula *does* yield ~0.244
    (unexpected); otherwise returns audit details with status FAIL_CLAIM.
    """
    p = PARAMS
    r = p.r_preprint_um if r is None else r
    R = p.R_preprint_um if R is None else R
    hbar_cs = p.hbar_cs_meV_um if hbar_cs is None else hbar_cs
    x = 2.0 * math.pi * r / R
    F = bessel_mode_sum(x)
    E_old = casimir_preprint_formula(hbar_cs, r, R, F)
    claim_ok = abs(E_old - 0.244) <= tol
    result = {
        "r_um": r,
        "R_um": R,
        "hbar_cs": hbar_cs,
        "x": x,
        "F": F,
        "E_old_meV": E_old,
        "claimed_meV": 0.244,
        "claim_reproduced": claim_ok,
        "audit": (
            "FAIL_CLAIM: preprint 0.244 meV is NOT reproduced by "
            "old Casimir formula on stated preprint geometry"
            if not claim_ok
            else "UNEXPECTED: old formula unexpectedly matches 0.244"
        ),
    }
    if claim_ok:
        raise AssertionError(
            "Old Casimir unexpectedly reproduces 0.244 meV; audit assumption broken"
        )
    return result


def provisional_corrections_self_energy_spec() -> FrozenVacuumSpec:
    """Default provisional lock: Corrections geometry + corrected formula."""
    p = PARAMS
    return FrozenVacuumSpec(
        formula_id="corrected_self_energy_2d",
        r=p.r_corr_um,
        R=p.R_corr_um,
        hbar_cs=p.hbar_cs_meV_um,
        F_definition="bessel_K0_sum",
        geometry_label="corrections_moire",
        status="provisional",
        notes=(
            "PROVISIONAL — not for quoting as a final meV in papers. "
            "Uses Corrections (r=0.335 nm, R=13.4 nm) with "
            "E_vac=(1/4π²)(ħc_s/r)F(2πr/R). Geometry and ħc_s still open."
        ),
    )


def freeze(
    spec: FrozenVacuumSpec,
    *,
    lock_name: str | None = None,
    locks_dir: Path | None = None,
) -> Path:
    """Write a JSON lockfile under ``qvc/locks/`` for the chosen geometry.

    Returns the path to the written lockfile.
    """
    out_dir = LOCKS_DIR if locks_dir is None else Path(locks_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if lock_name is None:
        lock_name = f"{spec.status}_{spec.geometry_label or spec.formula_id}.json"
    if not lock_name.endswith(".json"):
        lock_name += ".json"

    eval_result = evaluate_vacuum(spec)
    # Re-assert audit on preprint geometry whenever freezing
    audit = assert_preprint_0p244_not_reproduced()

    payload = {
        "lock_schema": "qvc.FrozenVacuumSpec/v1",
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "spec": asdict(spec),
        "evaluation": eval_result,
        "preprint_0p244_audit": audit,
        "warning": (
            "If status is 'provisional', do NOT quote E_vac meV as a "
            "final paper result. See docs/VACUUM_FREEZE_SPEC.md."
        ),
    }
    path = out_dir / lock_name
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def ensure_default_provisional_lock() -> Path:
    """Write the default provisional Corrections self-energy lock if absent/refresh."""
    return freeze(
        provisional_corrections_self_energy_spec(),
        lock_name="provisional_corrections_self_energy.json",
    )
