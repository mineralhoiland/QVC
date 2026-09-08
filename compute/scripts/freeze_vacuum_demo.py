#!/usr/bin/env python3
"""Demo vacuum evaluation without claiming a final meV until frozen.

Compares preprint vs Corrections geometries and old vs corrected formulas.
Writes/refreshes the provisional Corrections self-energy lockfile.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from qvc.params import PARAMS
from qvc.vacuum_freeze import (
    FrozenVacuumSpec,
    assert_preprint_0p244_not_reproduced,
    candidate_geometries,
    ensure_default_provisional_lock,
    evaluate_vacuum,
    freeze,
    provisional_corrections_self_energy_spec,
)


def main() -> int:
    print("=== QVC Vacuum Freeze Demo (PROVISIONAL) ===\n")
    print("Candidate geometries:")
    print(json.dumps(candidate_geometries(), indent=2))

    print("\n--- Audit: preprint 0.244 meV claim ---")
    audit = assert_preprint_0p244_not_reproduced()
    print(json.dumps(audit, indent=2))

    p = PARAMS
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
            notes="Demo evaluation only — not a paper meV claim.",
        )
        ev = evaluate_vacuum(spec)
        print(f"\n--- {label} ---")
        print(json.dumps(ev, indent=2))

    path = ensure_default_provisional_lock()
    print(f"\nProvisional lock written: {path}")
    print(
        "\nWARNING: Do not quote E_vac meV in papers until "
        "docs/VACUUM_FREEZE_SPEC.md open decisions are closed."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
