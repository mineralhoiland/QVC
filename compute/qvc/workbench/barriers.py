"""S4 — Catalyzed barriers (all L1, phenomenological rate estimates labelled).

1. Mexican hat  V(ρ) = a ρ² + (b/2) ρ⁴  with the *same* cubic-quintic
   coefficients the GPE uses (a = −Δ, b = Δ₀/n₀). Track A: a = −Δ₀ (bare)
   versus catalyzon-tracked mass a = −Δ_cat, and (parallel channel) the static
   fold a = −Δ*. Well depth B = a²/2b, hence B_cat/B₀ = (Δ_cat/Δ₀)².
2. Fold map δ(α) with physical / unstable branches and δ_c = α_c = 0.1818,
   the λ* operating point, and the fold-out unit-participation point.
3. Arrhenius / WKB *sketch* of Γ_QVC/Γ_bare at T ≪ 7 K (order of magnitude,
   not a theorem; no bounce integrator exists in-repo).
4. A_ZPF slider: U_eff(ρ; A) reshaped through the fold branch Δ(A) at
   **fixed total energy** — the vacuum term is traceless (Tr G = 0), so the
   ground-state energy is pinned and only the landscape shape changes
   (catalytic, not energetic).
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.gap_fold import fold_curve
from qvc.tqgl.locks import fold_delta_of_A
from qvc.workbench.locks import Layer, get_locks, stamp

M_E_MEV_PS2_PER_NM2 = 5.685630e-3  # electron mass in meV·ps²/nm²


# --------------------------------------------------------------------------
# 1. Mexican hat
# --------------------------------------------------------------------------


def mexican_hat(
    Delta_meV: float,
    *,
    Delta0_meV: float | None = None,
    n0: float = 1.0,
    rho_max: float | None = None,
    n_pts: int = 301,
) -> dict[str, Any]:
    """V(ρ) = −Δ ρ² + (Δ₀/2n₀) ρ⁴  (GPE cubic-quintic coefficients)."""
    L = get_locks()
    Delta0 = L.Delta0_meV if Delta0_meV is None else Delta0_meV
    a = -float(Delta_meV)
    b = Delta0 / n0
    rho_min = math.sqrt(-a / b) if a < 0 else 0.0
    depth = a * a / (2.0 * b) if a < 0 else 0.0
    rmax = (1.6 * max(rho_min, math.sqrt(n0))) if rho_max is None else rho_max
    rho = np.linspace(0.0, rmax, n_pts)
    V = a * rho**2 + 0.5 * b * rho**4
    return {
        "rho": rho,
        "V_meV": V,
        "a_meV": a,
        "b_meV": b,
        "rho_min": rho_min,
        "n_min": rho_min**2,
        "well_depth_meV": depth,
        "amplitude_stiffness_meV": -4.0 * a if a < 0 else 2.0 * a,  # V''(ρ_min) = −4a
        "Delta_meV": float(Delta_meV),
    }


def mexican_hat_tracks() -> dict[str, Any]:
    """Track A (bare Δ₀) vs catalyzon-tracked (Δ_cat) vs static fold (Δ*) — parallel channels."""
    L = get_locks()
    bare = mexican_hat(L.Delta0_meV)
    cat = mexican_hat(L.Delta_cat_meV)
    fold = mexican_hat(L.Delta_star_meV)
    ratio_cat = cat["well_depth_meV"] / bare["well_depth_meV"]
    ratio_fold = fold["well_depth_meV"] / bare["well_depth_meV"]
    return {
        "layer": Layer.L1.value,
        "caption": stamp(
            Layer.L1,
            "Mexican hat V(ρ)=−Δρ²+(Δ₀/2n₀)ρ⁴ — bare Δ₀, catalyzon Δ_cat, static-fold Δ* are parallel channels (never stacked)",
        ),
        "tracks": {
            "bare_Delta0": bare,
            "catalyzon_Delta_cat": cat,
            "static_fold_Delta_star": fold,
        },
        "B_cat_over_B0": ratio_cat,
        "B_cat_over_B0_expected": (L.Delta_cat_meV / L.Delta0_meV) ** 2,
        "B_fold_over_B0": ratio_fold,
        "B_fold_over_B0_expected": (L.Delta_star_meV / L.Delta0_meV) ** 2,
        "identity": "B_cat/B_0 = (Δ_cat/Δ_0)^2  (exact for the quartic well)",
    }


def mexican_hat_surface(Delta_meV: float, *, n_r: int = 60, n_phi: int = 90, **kw) -> dict[str, Any]:
    """(ρ cos φ, ρ sin φ, V(ρ)) surface for the 3D viewport."""
    mh = mexican_hat(Delta_meV, n_pts=n_r, **kw)
    phi = np.linspace(0, 2 * np.pi, n_phi)
    R, PHI = np.meshgrid(mh["rho"], phi, indexing="ij")
    Vg = np.broadcast_to(mh["V_meV"][:, None], R.shape)
    return {"X": R * np.cos(PHI), "Y": R * np.sin(PHI), "Z": np.array(Vg), **{k: v for k, v in mh.items() if k not in ("rho", "V_meV")}}


# --------------------------------------------------------------------------
# 2. Fold map
# --------------------------------------------------------------------------


def fold_map(n: int = 160) -> dict[str, Any]:
    L = get_locks()
    t = L.tqgl
    fc = fold_curve(ell=t.ell, n=n)
    a = np.asarray(fc["alpha_over_ac"], float)
    phys = np.array([np.nan if v is None else v for v in fc["delta_phys"]], float)
    unst = np.array([np.nan if v is None else v for v in fc["delta_unst"]], float)
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, "Fold map δ = 1 + α(ln δ + ℓ): physical / unstable branches meet at δ_c = α_c = 0.1818"),
        "alpha_over_ac": a,
        "delta_phys": phys,
        "delta_unst": unst,
        "alpha_c": t.alpha_c,
        "delta_c": t.delta_c,
        "ell": t.ell,
        "operating_point": {"alpha_over_ac": t.alpha_star_over_ac, "delta": t.delta_star, "label": "λ* = 0.90 λ_fold"},
        "unit_lambda": {"alpha_over_ac": (t.E_vac_meV / t.Delta0_meV) / t.alpha_c, "label": "λ = 1 (fold-out, no physical branch)"},
        "Delta0_meV": t.Delta0_meV,
        "Delta_meV_phys": phys * t.Delta0_meV,
    }


# --------------------------------------------------------------------------
# 3. Arrhenius / WKB sketch
# --------------------------------------------------------------------------


def rate_catalysis_sketch(
    *,
    T_K: np.ndarray | None = None,
    m_star_over_me: float = 0.5,
    barrier_width_nm: float | None = None,
) -> dict[str, Any]:
    """Γ_QVC/Γ_bare from gap substitution Δ_ex → Δ_cat (order of magnitude).

    Thermal (Arrhenius):  Γ ∝ exp(−Δ/k_BT)  ⇒  ratio = exp(g₀*/k_BT).
    Tunnelling (WKB, rectangular barrier of height Δ and width w = a_H):
        S/ħ ≈ 2 w √(2 m* Δ)/ħ  ⇒  ratio = exp[2w√(2m*)(√Δ_ex − √Δ_cat)/ħ].
    Labelled sketch: no bounce/instanton integrator exists in-repo.
    """
    L = get_locks()
    t = L.tqgl
    T = np.linspace(0.1, 7.0, 140) if T_K is None else np.asarray(T_K, float)
    kT = t.kB_meV_per_K * T
    ratio_arr = np.exp(L.g0_star_meV / kT)
    ratio_fold = np.exp((t.Delta0_meV - t.Delta_star_meV) / kT)
    w = t.r_nm if barrier_width_nm is None else barrier_width_nm
    m = m_star_over_me * M_E_MEV_PS2_PER_NM2
    S_bare = 2.0 * w * math.sqrt(2.0 * m * t.Delta0_meV) / t.hbar_meV_ps
    S_cat = 2.0 * w * math.sqrt(2.0 * m * L.Delta_cat_meV) / t.hbar_meV_ps
    S_fold = 2.0 * w * math.sqrt(2.0 * m * t.Delta_star_meV) / t.hbar_meV_ps
    return {
        "layer": Layer.L1.value,
        "caption": stamp(
            Layer.L1,
            "Rate catalysis SKETCH: Arrhenius exp(g₀*/k_BT) and WKB exponent shift — order of magnitude, not a theorem",
        ),
        "T_K": T,
        "arrhenius_ratio_catalyzon": ratio_arr,
        "arrhenius_ratio_fold_channel": ratio_fold,
        "T_7K_marker": 7.0,
        "wkb": {
            "m_star_over_me": m_star_over_me,
            "width_nm": w,
            "S_over_hbar_bare": S_bare,
            "S_over_hbar_catalyzon": S_cat,
            "S_over_hbar_fold": S_fold,
            "ratio_catalyzon": math.exp(S_bare - S_cat),
            "ratio_fold": math.exp(S_bare - S_fold),
        },
        "status": "phenomenological sketch (no bounce/instanton solver in-repo)",
    }


# --------------------------------------------------------------------------
# 4. A_ZPF slider at fixed total energy
# --------------------------------------------------------------------------


def u_eff_of_A(
    A_um_inv: float,
    *,
    n0: float = 1.0,
    n_pts: int = 301,
) -> dict[str, Any]:
    """U_eff(ρ; A) with the mass tracked by the static fold branch Δ(A).

    The vacuum term is traceless, so the condensate energy at the well minimum
    is held fixed: we subtract the constant c(A) = U_A(ρ_min(A)) − U_0(ρ_min(0)).
    What the slider changes is the *shape* (ρ_min, curvature), not the energy.
    """
    L = get_locks()
    t = L.tqgl
    info = fold_delta_of_A(float(A_um_inv), t)
    folded = bool(info["folded_out"])
    Delta = t.Delta0_meV if A_um_inv <= 0 else (info["Delta_phys_meV"] if not folded else None)
    base = mexican_hat(t.Delta0_meV, n0=n0, n_pts=n_pts)
    if Delta is None:
        # past the fold: no physical branch — return the bare landscape flagged
        return {
            **base,
            "A_um_inv": float(A_um_inv),
            "A_over_Ac": float(A_um_inv) / t.A_c_fold_um_inv,
            "alpha": info["alpha"],
            "folded_out": True,
            "Delta_track_meV": None,
            "V_meV_fixed_energy": base["V_meV"],
            "energy_offset_meV": 0.0,
            "note": "A > A_c: static theory has no physical branch (fold-out) — landscape not defined",
        }
    cur = mexican_hat(Delta, n0=n0, n_pts=n_pts, rho_max=float(base["rho"][-1]))
    offset = cur["V_meV"].min() - base["V_meV"].min()
    return {
        **cur,
        "A_um_inv": float(A_um_inv),
        "A_over_Ac": float(A_um_inv) / t.A_c_fold_um_inv,
        "alpha": info["alpha"],
        "alpha_over_ac": info["alpha_over_ac"],
        "folded_out": False,
        "Delta_track_meV": Delta,
        "V_meV_fixed_energy": cur["V_meV"] - offset,
        "energy_offset_meV": float(offset),
        "bare_V_meV": base["V_meV"],
        "well_depth_ratio": cur["well_depth_meV"] / base["well_depth_meV"],
        "note": "total condensate energy pinned (Tr G = 0): only the landscape shape moves — catalytic, not energetic",
    }


def a_zpf_sweep(n_A: int = 25, n0: float = 1.0) -> dict[str, Any]:
    L = get_locks()
    t = L.tqgl
    A = np.linspace(0.0, 0.999 * t.A_c_fold_um_inv, n_A)
    rows = [u_eff_of_A(float(a), n0=n0, n_pts=121) for a in A]
    return {
        "layer": Layer.L1.value,
        "caption": stamp(Layer.L1, "A_ZPF reshapes U_eff along the fold branch at fixed total energy (A_c^fold = α_cΔ₀/ℏc_s)"),
        "A_um_inv": A,
        "A_over_Ac": A / t.A_c_fold_um_inv,
        "Delta_track_meV": np.array([r["Delta_track_meV"] if r["Delta_track_meV"] is not None else np.nan for r in rows]),
        "rho_min": np.array([r["rho_min"] for r in rows]),
        "well_depth_meV": np.array([r["well_depth_meV"] for r in rows]),
        "stiffness_meV": np.array([r["amplitude_stiffness_meV"] for r in rows]),
        "A_star_um_inv": t.A_star_um_inv,
        "A_c_fold_um_inv": t.A_c_fold_um_inv,
        "A_star_over_Ac": t.A_star_um_inv / t.A_c_fold_um_inv,
        "rows": rows,
    }


def barriers_report() -> dict[str, Any]:
    mh = mexican_hat_tracks()
    fm = fold_map(n=40)
    rk = rate_catalysis_sketch(T_K=np.array([0.1, 0.5, 1.0, 4.0, 7.0]))
    sw = a_zpf_sweep(n_A=6)
    return {
        "mexican_hat": {
            "B_cat_over_B0": mh["B_cat_over_B0"],
            "B_cat_over_B0_expected": mh["B_cat_over_B0_expected"],
            "B_fold_over_B0": mh["B_fold_over_B0"],
            "well_depth_bare_meV": mh["tracks"]["bare_Delta0"]["well_depth_meV"],
            "well_depth_cat_meV": mh["tracks"]["catalyzon_Delta_cat"]["well_depth_meV"],
            "well_depth_fold_meV": mh["tracks"]["static_fold_Delta_star"]["well_depth_meV"],
        },
        "fold": {"alpha_c": fm["alpha_c"], "delta_c": fm["delta_c"], "operating_point": fm["operating_point"]},
        "rates": {
            "T_K": rk["T_K"].tolist(),
            "arrhenius_ratio_catalyzon": rk["arrhenius_ratio_catalyzon"].tolist(),
            "wkb": rk["wkb"],
            "status": rk["status"],
        },
        "a_zpf": {
            "A_star_over_Ac": sw["A_star_over_Ac"],
            "Delta_track_meV": sw["Delta_track_meV"].tolist(),
            "energy_pinned": all(abs(r.get("V_meV_fixed_energy", np.zeros(1)).min() - sw["rows"][0]["V_meV"].min()) < 1e-9 for r in sw["rows"] if not r["folded_out"]),
        },
    }
