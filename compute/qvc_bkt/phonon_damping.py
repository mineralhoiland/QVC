"""C6: phonon damping / AHNS overdamped criterion (H-VRG2).

Estimate vortex mobility scales and test ħ/τ_v ≳ πκ. Gao–Khalaf-type
phonon numbers already in the repo are used where present; otherwise the
estimate is labeled. Report jump vs rounded crossover.
"""

from __future__ import annotations

import math
from typing import Any

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.retardation import eta_ret  # noqa: E402

from qvc_bkt.locks import BKTLocks, build_bkt_locks
from qvc_bkt.stiffness import kappa_GL_mev_um2


def gpe_phonon_scale(locks: BKTLocks) -> dict[str, Any]:
    """Phonon drive already used by the truncated GPE/DCE solver.

    g_ph = 0.1 Δ₀, ω₀ = 2Δ₀/ħ. Status: computed from existing solver
    parameters, not a derived vortex mobility.
    """
    g_ph = 0.10 * locks.Delta0_meV
    omega_ps = 2.0 * locks.Delta0_meV / locks.hbar_meV_ps
    # Identify ħ/τ_v with the phonon potential amplitude (labeled).
    hbar_over_tau = g_ph
    tau_ps = locks.hbar_meV_ps / hbar_over_tau
    return {
        "kind": "GPE_phonon_potential",
        "status": "labeled_estimate",
        "g_ph_meV": float(g_ph),
        "omega0_ps_inv": float(omega_ps),
        "hbar_over_tau_meV": float(hbar_over_tau),
        "tau_v_ps": float(tau_ps),
        "formula": "ħ/τ_v ∼ g_ph = 0.1 Δ₀  (existing GPE V_ph amplitude)",
        "note": (
            "This is the phonon scale already in qvc.tqgl.coupled, not a "
            "Bardeen–Stephen mobility from a spectral function."
        ),
    }


def eliashberg_gao_khalaf(locks: BKTLocks, T_K: float, *, lambda_ep: float = 1.0) -> dict[str, Any]:
    """Gao–Khalaf-type electron–phonon scattering OOM: ħ/τ ∼ 2π λ_ep kT.

    λ_ep=1 is a MATBG-typical labeled coupling, not computed from a
    TEVC Eliashberg function.
    """
    kT = T_K * locks.kB_meV_per_K
    hbar_over_tau = 2.0 * math.pi * lambda_ep * kT
    tau_ps = locks.hbar_meV_ps / max(hbar_over_tau, 1e-18)
    return {
        "kind": "Eliashberg_GaoKhalaf_OOM",
        "status": "ansatz",
        "lambda_ep": float(lambda_ep),
        "T_K": float(T_K),
        "hbar_over_tau_meV": float(hbar_over_tau),
        "tau_v_ps": float(tau_ps),
        "formula": "ħ/τ_v ∼ 2π λ_ep k_B T  with λ_ep=1 labeled MATBG-typical",
        "note": (
            "Gao–Khalaf geometric e-ph enhancement is not integrated over the "
            "moiré BZ here; λ_ep is an O(1) placeholder."
        ),
    }


def retardation_scale(locks: BKTLocks) -> dict[str, Any]:
    """Locked retardation η_ret = 2 Δ₀ ξ / (ħ c_s) with ξ = a_H = r."""
    eta = eta_ret(locks.Delta0_meV, locks.r_um, locks.hbar_cs_meV_um)
    return {
        "kind": "retardation_eta",
        "status": "locked_algebra",
        "eta_ret": float(eta),
        "xi_um": float(locks.r_um),
        "note": "Causal phonon propagator identity. Not by itself a vortex damping rate.",
    }


def ahns_criterion(hbar_over_tau_meV: float, kappa_meV: float) -> dict[str, Any]:
    """Overdamped if ħ/τ_v ≳ πκ (H-VRG2)."""
    pi_k = math.pi * kappa_meV
    over = hbar_over_tau_meV >= pi_k
    ratio = hbar_over_tau_meV / max(pi_k, 1e-18)
    return {
        "pi_kappa_meV": float(pi_k),
        "hbar_over_tau_meV": float(hbar_over_tau_meV),
        "ratio": float(ratio),
        "overdamped": bool(over),
        "verdict": "rounded_crossover" if over else "sharp_jump_survives",
        "criterion": "ħ/τ_v ≳ πκ",
        "status": "computed_given_estimates",
    }


def scan_ahns(
    rows_kappa: list[dict[str, Any]],
    hbar_over_tau_meV: float,
) -> dict[str, Any]:
    n_over = 0
    n = 0
    first_over_alpha = None
    recs = []
    for rec in rows_kappa:
        n += 1
        ah = ahns_criterion(hbar_over_tau_meV, float(rec["kappa_meV"]))
        if ah["overdamped"]:
            n_over += 1
            if first_over_alpha is None:
                first_over_alpha = rec["alpha"]
        recs.append(
            {
                "alpha_over_ac": rec["alpha_over_ac"],
                "kappa_meV": rec["kappa_meV"],
                "ratio": ah["ratio"],
                "overdamped": ah["overdamped"],
                "verdict": ah["verdict"],
            }
        )
    return {
        "n": n,
        "n_overdamped": n_over,
        "fraction_overdamped": n_over / max(n, 1),
        "first_overdamped_alpha_over_ac": (
            None if first_over_alpha is None else recs[next(i for i, r in enumerate(recs) if r["overdamped"])]["alpha_over_ac"]
        ),
        "rows": recs,
    }


def run_c6(
    c1: dict[str, Any],
    c3: dict[str, Any],
    locks: BKTLocks | None = None,
) -> dict[str, Any]:
    locks = locks or build_bkt_locks()
    gpe = gpe_phonon_scale(locks)
    elia_star = eliashberg_gao_khalaf(locks, locks.T_star_K)
    elia_1K = eliashberg_gao_khalaf(locks, 1.0)
    eta = retardation_scale(locks)
    k_gl = kappa_GL_mev_um2(locks.tqgl)

    rows = c1["primary"]["rows"]
    scan_gpe = scan_ahns(rows, gpe["hbar_over_tau_meV"])
    scan_el = scan_ahns(rows, elia_star["hbar_over_tau_meV"])

    k_uv = c1["primary"]["kappa_uv_meV"]
    k_fold = c1["primary"]["kappa_fold_meV"]
    uv_gpe = ahns_criterion(gpe["hbar_over_tau_meV"], k_uv)
    fold_gpe = ahns_criterion(gpe["hbar_over_tau_meV"], k_fold)
    uv_el = ahns_criterion(elia_star["hbar_over_tau_meV"], k_uv)
    fold_el = ahns_criterion(elia_star["hbar_over_tau_meV"], k_fold)

    # Would-be KT surface: if C3 finds a crossing, evaluate damping there.
    kt_surface = None
    sc = c3["scan_Tstar_primary"]
    if sc["alpha_kt_bare"] is not None:
        rec = min(rows, key=lambda r: abs(r["alpha"] - sc["alpha_kt_bare"]))
        kt_surface = {
            "alpha_over_ac": rec["alpha_over_ac"],
            "kappa_meV": rec["kappa_meV"],
            "GPE": ahns_criterion(gpe["hbar_over_tau_meV"], rec["kappa_meV"]),
            "Eliashberg_Tstar": ahns_criterion(elia_star["hbar_over_tau_meV"], rec["kappa_meV"]),
        }

    # H-VRG2 reading
    if fold_el["overdamped"] or fold_gpe["overdamped"]:
        h_vrg2 = (
            "Near the fold, labeled phonon scales satisfy ħ/τ_v ≳ πκ: "
            "a would-be universal jump is rounded into a crossover (H-VRG2). "
            "This does not by itself establish H-VRG1."
        )
        jump = "rounded_crossover_near_fold"
    elif uv_el["overdamped"] or uv_gpe["overdamped"]:
        h_vrg2 = "Overdamped already at the UV stiffness under this estimate."
        jump = "rounded_everywhere"
    else:
        h_vrg2 = (
            "Underdamped at T* for the primary κ(α): a sharp stiffness jump "
            "could survive if a KT surface were hit. H-VRG2 is not forced."
        )
        jump = "sharp_jump_possible"

    return {
        "task": "C6",
        "hypothesis": "H-VRG2",
        "gpe_phonon": gpe,
        "eliashberg_Tstar": elia_star,
        "eliashberg_1K": elia_1K,
        "retardation": eta,
        "anchors": {
            "UV_GPE": uv_gpe,
            "fold_GPE": fold_gpe,
            "UV_Eliashberg_Tstar": uv_el,
            "fold_Eliashberg_Tstar": fold_el,
        },
        "scan_GPE": {k: scan_gpe[k] for k in ("n", "n_overdamped", "fraction_overdamped", "first_overdamped_alpha_over_ac")},
        "scan_Eliashberg_Tstar": {
            k: scan_el[k]
            for k in ("n", "n_overdamped", "fraction_overdamped", "first_overdamped_alpha_over_ac")
        },
        "scan_GPE_rows": scan_gpe["rows"],
        "scan_Eliashberg_rows": scan_el["rows"],
        "kt_surface": kt_surface,
        "jump_vs_crossover": jump,
        "h_vrg2_reading": h_vrg2,
        "kappa_GL_meV_um2": k_gl,
        "closures": {
            "g_ph": "labeled_from_existing_GPE_solver",
            "lambda_ep": "ansatz O(1)",
            "AHNS_algebra": "derived given ħ/τ_v and κ",
        },
    }
