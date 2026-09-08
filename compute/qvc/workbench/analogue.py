"""S7 — L2 analogue geometry (optional, tagged panel). Simulated as geometry, not Einstein.

* Unruh–Visser acoustic metric from the superfluid velocity v = (ħ/m*)∇θ on a
  2D slice; sonic surface |v| = c_s; surface gravity κ and T_H = ħκ/2πk_B as a
  **diagnostic** against the locked T* = 118.7 mK (H-HAWK1), never as a lock.
* Torsion proxy T_axial ∝ J_Hopf ∝ ρ_H (Hopf density). H-TOR1 (torsion → G_ij)
  is parked: the toggle only annotates; nothing is fed into Spec(G).
* Goldstone dispersion cartoon from Paper C: 2 TA (helicity ±1) + 1 LA + 2
  fibre magnons. H-EIN1 diagnostic: Frobenius overlap of the TA strain tensor
  with a transverse-traceless (helicity ±2) polarisation is zero — TA phonons
  are not gravitons.

Nothing in this module is promoted to L1 or L3.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.workbench.locks import Layer, get_locks, stamp

KB_MEV_PER_K = 0.08617333262
HBAR_MEV_PS = 0.6582119569


def hbar_over_mstar_um2_ps() -> float:
    """ħ/m* from the GPE kinetic coefficient ħ²/2m* = Δ₀ r² (meV·µm²)."""
    L = get_locks()
    ke = L.Delta0_meV * L.tqgl.r_um**2
    return 2.0 * ke / HBAR_MEV_PS  # µm²/ps


def sound_speed_um_ps() -> float:
    L = get_locks()
    return L.tqgl.hbar_cs_meV_um / HBAR_MEV_PS  # µm/ps


def acoustic_metric_slice(
    *,
    n_grid: int = 96,
    L_um: float = 0.12,
    Q_H: float | None = None,
    drain_v0_over_cs: float = 0.0,
    drain_R_um: float | None = None,
) -> dict[str, Any]:
    """g_μν (Unruh–Visser) on the z=0 plane for v = (ħ/m*) Q_H φ̂/ρ + optional radial drain.

    ds² = −(c_s² − v²)dt² − 2 v·dx dt + dx².  The sonic surface is |v| = c_s.
    With the pure fractional vortex (Q_H = 1/N) the sonic radius ρ_h =
    (ħ/m*)Q_H/c_s sits far inside the healing core — reported honestly.
    ``drain_v0_over_cs`` adds v_r = −v₀ (R/ρ) (draining bathtub) so a horizon
    can be placed outside the core for the κ / T_H diagnostic.
    """
    L = get_locks()
    Q_H = L.Q_H if Q_H is None else Q_H
    cs = sound_speed_um_ps()
    hm = hbar_over_mstar_um2_ps()
    Rd = L.tqgl.R_um if drain_R_um is None else drain_R_um
    x = np.linspace(-L_um, L_um, n_grid)
    X, Y = np.meshgrid(x, x, indexing="xy")
    rho = np.sqrt(X * X + Y * Y)
    rho_safe = np.maximum(rho, 0.25 * (x[1] - x[0]))
    v_phi = hm * Q_H / rho_safe
    v_r = -drain_v0_over_cs * cs * (Rd / rho_safe)
    vx = v_r * X / rho_safe - v_phi * Y / rho_safe
    vy = v_r * Y / rho_safe + v_phi * X / rho_safe
    v2 = vx * vx + vy * vy
    g_tt = -(cs * cs - v2)
    g_tx, g_ty = -vx, -vy
    mach = np.sqrt(v2) / cs
    # sonic radius (analytic, radial part dominates the drain)
    rho_h_vortex = hm * Q_H / cs
    if drain_v0_over_cs > 0:
        # |v|² = v_r² + v_φ² = c_s² → solve for ρ numerically on a radial line
        r_line = np.linspace(0.2 * (x[1] - x[0]), L_um, 4000)
        vr_l = drain_v0_over_cs * cs * Rd / r_line
        vp_l = hm * Q_H / r_line
        f = np.sqrt(vr_l**2 + vp_l**2) - cs
        idx = np.where(np.diff(np.sign(f)) != 0)[0]
        if idx.size:
            i = idx[-1]
            rho_h = float(r_line[i] - f[i] * (r_line[i + 1] - r_line[i]) / (f[i + 1] - f[i]))
            # κ = d(|v_n| − c)/dn at horizon (radial normal): d|v_r|/dρ
            kappa = abs(float(np.gradient(vr_l, r_line)[i]))  # ps⁻¹
        else:
            rho_h, kappa = float("nan"), float("nan")
    else:
        rho_h = rho_h_vortex
        kappa = hm * Q_H / rho_h**2 if rho_h > 0 else float("nan")  # d v_φ/dρ magnitude
    T_H_mK = 1e3 * HBAR_MEV_PS * kappa / (2.0 * math.pi * KB_MEV_PER_K) if np.isfinite(kappa) else float("nan")
    return {
        "layer": Layer.L2.value,
        "caption": stamp(
            Layer.L2,
            f"Unruh–Visser acoustic metric on z=0; sonic radius ρ_h = {rho_h*1e3:.3g} nm "
            f"(core a_H = {L.r_nm:.0f} nm); T_H = {T_H_mK:.3g} mK vs locked T* = {L.T_star_mK} mK — DIAGNOSTIC (H-HAWK1)",
        ),
        "x_um": x,
        "g_tt": g_tt,
        "g_tx": g_tx,
        "g_ty": g_ty,
        "mach": mach,
        "vx": vx,
        "vy": vy,
        "c_s_um_ps": cs,
        "hbar_over_mstar_um2_ps": hm,
        "rho_h_um": rho_h,
        "rho_h_vortex_only_um": rho_h_vortex,
        "horizon_inside_core": bool(rho_h < L.tqgl.r_um),
        "kappa_ps_inv": kappa,
        "T_H_mK": T_H_mK,
        "T_star_mK": L.T_star_mK,
        "drain_v0_over_cs": drain_v0_over_cs,
        "note": "Analogue kinematics only (T_H = κ/2π). Not Einstein gravity; T* is a Schwinger-channel lock, not T_H.",
    }


def torsion_proxy(rho_H: np.ndarray, *, h_tor1_toggle: bool = False) -> dict[str, Any]:
    """Axial torsion proxy T ∝ J_Hopf ∝ ρ_H. H-TOR1 is parked; the toggle only annotates."""
    scale = float(np.max(np.abs(rho_H))) or 1.0
    return {
        "layer": Layer.L2.value,
        "caption": stamp(Layer.L2, "Axial torsion proxy T_axial ∝ J_Hopf (normalised Hopf density). H-TOR1 parked: never fed into G_ij"),
        "T_axial_normalised": rho_H / scale,
        "h_tor1_toggle": h_tor1_toggle,
        "feeds_Gij": False,
        "note": (
            "H-TOR1 toggle is ON: this remains a parked experiment; Spec(G) shown elsewhere is unchanged."
            if h_tor1_toggle
            else "H-TOR1 parked (toggle off)."
        ),
    }


def goldstone_cartoon(
    *,
    c_LA_m_s: float = 21400.0,
    c_TA_m_s: float = 13600.0,
    a_nm: float = 50.0,
    torsion_zeta: float = 0.08,
    magnon_gap_meV: float = 0.05,
    c_magnon_m_s: float = 8000.0,
    n_k: int = 200,
) -> dict[str, Any]:
    """Five Goldstone branches vs k up to the BZ edge π/a (Paper C Table; cartoon)."""
    k = np.linspace(0.0, math.pi / a_nm, n_k)  # nm⁻¹
    conv = HBAR_MEV_PS * 1e-3  # ħ·(m/s)·(nm⁻¹) → meV : ħ[meV ps]·v[m/s]=v[nm/ps]·1e-3 ... v[m/s]=1e-3 nm/ps
    w_LA = conv * c_LA_m_s * k
    w_TA = conv * c_TA_m_s * k * (1.0 - torsion_zeta * (k * a_nm / math.pi) ** 2)
    w_m1 = np.sqrt(magnon_gap_meV**2 + (conv * c_magnon_m_s * k) ** 2)
    w_m2 = np.sqrt((1.3 * magnon_gap_meV) ** 2 + (conv * c_magnon_m_s * k) ** 2)
    return {
        "layer": Layer.L2.value,
        "caption": stamp(Layer.L2, "Goldstone cartoon: 2 TA (helicity ±1, degenerate at Γ) + 1 LA + 2 fibre magnons — not gravitons (H-EIN1)"),
        "k_nm_inv": k,
        "branches": {
            "LA": w_LA,
            "TA1": w_TA,
            "TA2": w_TA.copy(),
            "magnon_1 (fiber)": w_m1,
            "magnon_2 (fiber)": w_m2,
        },
        "params": {"c_LA_m_s": c_LA_m_s, "c_TA_m_s": c_TA_m_s, "a_nm": a_nm, "torsion_zeta": torsion_zeta},
        "note": "TA roton-like correction O((ka)²) from torsion is a sketch; speeds from phonon_folding defaults",
    }


def h_ein1_overlap(k_dir=(0.0, 0.0, 1.0), e_dir=(1.0, 0.0, 0.0)) -> dict[str, Any]:
    """Frobenius overlap of a TA strain tensor with helicity-±2 TT polarisations (→ 0)."""
    k = np.asarray(k_dir, float)
    k /= np.linalg.norm(k)
    e = np.asarray(e_dir, float)
    e -= k * (e @ k)
    e /= np.linalg.norm(e)
    strain = 0.5 * (np.outer(k, e) + np.outer(e, k))  # helicity ±1 (vector) mode
    # TT basis in the plane ⊥ k
    f = np.cross(k, e)
    e_plus = np.outer(e, e) - np.outer(f, f)
    e_cross = np.outer(e, f) + np.outer(f, e)
    ov_plus = float(np.sum(strain * e_plus)) / (np.linalg.norm(strain) * np.linalg.norm(e_plus))
    ov_cross = float(np.sum(strain * e_cross)) / (np.linalg.norm(strain) * np.linalg.norm(e_cross))
    return {
        "layer": Layer.L2.value,
        "caption": stamp(Layer.L2, f"H-EIN1: TA strain vs TT (+,×) overlap = ({ov_plus:.2e}, {ov_cross:.2e}) — TA is helicity ±1, not a graviton"),
        "overlap_plus": ov_plus,
        "overlap_cross": ov_cross,
        "strain_tensor": strain,
        "is_graviton_like": bool(max(abs(ov_plus), abs(ov_cross)) > 1e-8),
    }


L3_SCHEMATIC_NOTE = stamp(
    Layer.L3,
    "Spacetime overlay is a SCHEMATIC: no Sakharov G, no hopfion-phonon gravitons, no E8 taxonomy, meV locks are not cosmology.",
)


def analogue_report(rho_H: np.ndarray | None = None) -> dict[str, Any]:
    am = acoustic_metric_slice()
    am_d = acoustic_metric_slice(drain_v0_over_cs=1.2)
    ein = h_ein1_overlap()
    out = {
        "acoustic_vortex_only": {
            "rho_h_nm": am["rho_h_um"] * 1e3,
            "horizon_inside_core": am["horizon_inside_core"],
            "T_H_mK": am["T_H_mK"],
            "T_star_mK": am["T_star_mK"],
        },
        "acoustic_with_drain_1p2cs": {"rho_h_nm": am_d["rho_h_um"] * 1e3, "T_H_mK": am_d["T_H_mK"], "kappa_ps_inv": am_d["kappa_ps_inv"]},
        "h_ein1": {"overlap_plus": ein["overlap_plus"], "overlap_cross": ein["overlap_cross"], "is_graviton_like": ein["is_graviton_like"]},
        "h_tor1": "parked",
        "L3": L3_SCHEMATIC_NOTE,
    }
    if rho_H is not None:
        tp = torsion_proxy(rho_H)
        out["torsion_proxy_max"] = float(np.max(np.abs(tp["T_axial_normalised"])))
    return out
