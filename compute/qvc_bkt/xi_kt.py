"""C4: continuous Δ→0 diagnostic from the vortex sector.

If unbinding occurs, construct ξ_KT from the KT flow and the H-VRG1 target

    Δ ∼ ħ c_s / ξ_KT

and compare to the locked finite Δ_c. This is a numerical implication of
the vortex sector, not a replacement of the fold theorem.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_bkt.locks import BKTLocks, build_bkt_locks


def xi_from_flow_row(
    row: dict[str, Any],
    *,
    a_c_um: float,
    ell_ir: float,
    y_cut: float = 0.4,
) -> dict[str, Any]:
    """ξ_KT = a_c exp(ℓ*) with ℓ* the plasma scale, else > a_c exp(ℓ_IR)."""
    if row.get("ell_star") is not None:
        ell_star = float(row["ell_star"])
        xi = a_c_um * math.exp(ell_star)
        kind = "unbound_ell_star"
        finite = True
    elif int(row.get("plasma", 0)) == 1:
        # Bare K already below K_c: plasma on microscopic scales.
        xi = a_c_um
        ell_star = 0.0
        kind = "bare_plasma"
        finite = True
    else:
        ell_star = ell_ir
        xi = a_c_um * math.exp(ell_ir)
        kind = "bound_IR_lower_bound"
        finite = False
    return {
        "xi_KT_um": float(xi),
        "ell_star": None if not finite and kind == "bound_IR_lower_bound" else ell_star,
        "kind": kind,
        "finite": bool(finite),
        "y_cut": float(y_cut),
    }


def delta_from_xi(xi_um: float, hbar_cs_meV_um: float) -> float:
    """H-VRG1 target Δ = ħ c_s / ξ. Status: conjecture as a gap identification."""
    if xi_um <= 0.0:
        return float("nan")
    return float(hbar_cs_meV_um / xi_um)


def build_vortex_gap_curve(
    scan: dict[str, Any],
    *,
    locks: BKTLocks,
    a_c_um: float,
    ell_ir: float = 8.0,
) -> dict[str, Any]:
    rows_out = []
    unbound_any = False
    for rec in scan["rows"]:
        x = xi_from_flow_row(rec, a_c_um=a_c_um, ell_ir=ell_ir)
        d_v = delta_from_xi(x["xi_KT_um"], locks.hbar_cs_meV_um) if x["finite"] else float("nan")
        if x["finite"]:
            unbound_any = True
        rows_out.append(
            {
                "alpha": rec["alpha"],
                "alpha_over_ac": rec["alpha_over_ac"],
                "t_A": rec["t_A"],
                "delta_fold": rec["delta"],
                "Delta_fold_meV": rec["Delta_meV"],
                "K_bare": rec["K_bare"],
                "phase": rec["phase"],
                "xi_KT_um": x["xi_KT_um"],
                "xi_kind": x["kind"],
                "Delta_vortex_meV": d_v,
                "finite_unbinding": x["finite"],
            }
        )
    return {
        "unbinding_occurs": bool(unbound_any),
        "a_c_um": float(a_c_um),
        "ell_ir": float(ell_ir),
        "Delta_c_locked_meV": float(locks.Delta_c_meV),
        "hbar_cs_meV_um": float(locks.hbar_cs_meV_um),
        "rows": rows_out,
        "identification": {
            "formula": "Δ ∼ ħ c_s / ξ_KT",
            "status": "conjecture",
            "note": (
                "Vortex-sector implication only. Does not replace the locked "
                "fold theorem, which retains finite Δ_c in the amplitude map."
            ),
        },
    }


def run_c4(
    c3: dict[str, Any],
    c2: dict[str, Any],
    locks: BKTLocks | None = None,
) -> dict[str, Any]:
    """C4 driver. Primary at T*; counterfactual where KT-before-fold if needed."""
    locks = locks or build_bkt_locks()
    a_c = float(c2["a_c_um"])
    ell_ir = float(c3["scan_Tstar_primary"]["ell_ir"])
    primary = build_vortex_gap_curve(c3["scan_Tstar_primary"], locks=locks, a_c_um=a_c, ell_ir=ell_ir)
    conv = build_vortex_gap_curve(c3["scan_Tstar_conv_only"], locks=locks, a_c_um=a_c, ell_ir=ell_ir)

    # Counterfactual: first T on the grid where conv-only hits KT before fold.
    counterfactual = None
    for row in c3["T_eff_grid"]:
        if row["conv_only_kt_before_fold"]:
            T = float(row["T_K"])
            from qvc_bkt.kt_flow import competing_scale_scan

            # rebuild scan at that T from C3 primary's conv rows
            # (C3 already has T_eff_grid flags; reconstruct curve via stored scan if T==T*)
            if abs(T - locks.T_star_K) < 1e-12:
                curve = conv
            else:
                # Use C3 scan_Tstar_conv_only only at T*; for other T we skip a full
                # re-integrate here and record the flag. Full curve is built in suite
                # when requested. Store T and the competing-scale flag.
                curve = None
            counterfactual = {
                "T_K": T,
                "closure": "conv_only healing-disk, no geometric floor",
                "kt_before_fold": True,
                "curve_at_Tstar": curve is not None,
            }
            break

    d_fold = locks.Delta_c_meV
    d_vort_vals = [
        r["Delta_vortex_meV"]
        for r in primary["rows"]
        if r["finite_unbinding"] and np.isfinite(r["Delta_vortex_meV"])
    ]
    comparison = {
        "Delta_c_locked_meV": d_fold,
        "Delta_vortex_min_meV": float(np.nanmin(d_vort_vals)) if d_vort_vals else None,
        "Delta_vortex_at_last_plasma_meV": d_vort_vals[-1] if d_vort_vals else None,
        "unbinding_at_Tstar_primary": primary["unbinding_occurs"],
        "unbinding_at_Tstar_conv_only": conv["unbinding_occurs"],
        "note": (
            "If unbinding does not occur, ξ_KT is IR-limited and the vortex "
            "sector does not generate a continuous Δ→0 replacing Δ_c."
        ),
    }
    return {
        "task": "C4",
        "primary_Tstar": primary,
        "conv_only_Tstar": conv,
        "counterfactual": counterfactual,
        "vs_locked_Delta_c": comparison,
        "status": {
            "xi_KT": "computed_if_unbound_else_IR_bound",
            "Delta_identification": "conjecture",
            "fold_theorem": "not_replaced",
        },
    }
