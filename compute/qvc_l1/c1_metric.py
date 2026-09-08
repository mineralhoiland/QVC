"""C1 stiffness: two-sector identity on the physical fold branch.

    κ(α) = 2 κ_GL n0 [δ_+(α)]² + κ_geom

κ_GL = Δ₀ a_H² is derived given Option A*. n0 = 1/(π a_H²) remains ansatz.
κ_geom is the Brillouin-zone quantum metric of a Qi–Wu–Zhang C=1 two-band
toy (Peotta–Törmä geometric weight), not a MATBG continuum integral.

Competing-scale test: K(α)=κ(α)/T_eff versus K_c=2/π along δ_+ up to the
fold. A KT hit before α_c would be a C1-level sufficient condition for
H-VRG1; absence of a hit is not a disproof of the vortex-RG implication.

Status: computed toy-band estimator + derived two-sector algebra.
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc_bkt.kt_flow import K_C_BOSONIC  # noqa: E402
from qvc_bkt.stiffness import (  # noqa: E402
    K_of_kappa,
    kappa_GL_mev_um2,
    kappa_conv_meV,
    n0_healing_disk_um2,
    physical_branch_delta,
    twist_stiffness_uniform,
)
from qvc_l1.cs_eom import k_c_anyon
from qvc_l1.locks import L1Locks, build_l1_locks


def qi_wu_zhang_metric_kappa(
    *,
    Delta_meV: float,
    lambda_M_um: float,
    n_k: int = 48,
) -> dict[str, Any]:
    """2-band QWZ Chern insulator; κ_geom from the quantum metric of the lower band."""
    m = -1.0
    kx = np.linspace(-math.pi, math.pi, n_k, endpoint=False)
    ky = np.linspace(-math.pi, math.pi, n_k, endpoint=False)
    KX, KY = np.meshgrid(kx, ky, indexing="xy")
    dx = np.sin(KX)
    dy = np.sin(KY)
    dz = m + np.cos(KX) + np.cos(KY)
    d = np.sqrt(dx * dx + dy * dy + dz * dz) + 1e-18
    nx, ny, nz = dx / d, dy / d, dz / d

    def deriv(arr, axis):
        return np.gradient(arr, kx[1] - kx[0], axis=axis)

    ngx = (deriv(nx, 1), deriv(ny, 1), deriv(nz, 1))
    ngy = (deriv(nx, 0), deriv(ny, 0), deriv(nz, 0))
    gxx = 0.5 * (ngx[0] ** 2 + ngx[1] ** 2 + ngx[2] ** 2)
    gyy = 0.5 * (ngy[0] ** 2 + ngy[1] ** 2 + ngy[2] ** 2)
    gxy = 0.5 * (ngx[0] * ngy[0] + ngx[1] * ngy[1] + ngx[2] * ngy[2])
    omega = 0.5 * (
        nx * (deriv(ny, 1) * deriv(nz, 0) - deriv(nz, 1) * deriv(ny, 0))
        + ny * (deriv(nz, 1) * deriv(nx, 0) - deriv(nx, 1) * deriv(nz, 0))
        + nz * (deriv(nx, 1) * deriv(ny, 0) - deriv(ny, 1) * deriv(nx, 0))
    )
    dk2 = (kx[1] - kx[0]) * (ky[1] - ky[0])
    chern = float(np.sum(omega) * dk2 / (2.0 * math.pi))
    tr_g = float(np.mean(gxx + gyy))
    kappa_geom = (Delta_meV / (2.0 * math.pi)) * tr_g
    area_cell = 0.5 * math.sqrt(3.0) * lambda_M_um * lambda_M_um
    return {
        "status": "computed_toy_band",
        "model": "Qi-Wu-Zhang m=-1",
        "chern_numerical": chern,
        "chern_target": 1.0,
        "mean_tr_g": tr_g,
        "kappa_geom_meV": kappa_geom,
        "Delta_meV": Delta_meV,
        "lambda_M_um": lambda_M_um,
        "A_cell_um2": area_cell,
        "n_k": n_k,
        "not_MATBG_continuum": True,
        "note": (
            "Geometric stiffness from a C=1 two-band toy. Replaces the frozen "
            "floor Δ0/2π only as an estimator; the healing-disk n0 remains ansatz."
        ),
    }


def _scan_branch(
    locks: L1Locks,
    *,
    kappa_geom_uv: float,
    geom: str,
    n_alpha: int,
) -> dict[str, Any]:
    """κ(α), K(α) on δ_+. geom: frozen_toy | tracks_gap | floor | zero."""
    k_gl = kappa_GL_mev_um2(locks.tqgl)
    dens = n0_healing_disk_um2(locks.a_H_um)
    n0 = float(dens["n0_um2"])
    floor = locks.Delta0_meV / (2.0 * math.pi)
    T_star = locks.T_star_K
    kB = locks.kB_meV_per_K
    alphas = np.linspace(1e-6, 0.995 * locks.alpha_c, n_alpha)
    rows: list[dict[str, Any]] = []
    for a in alphas:
        d = physical_branch_delta(float(a), locks.ell)
        if d is None:
            continue
        rho2 = d * d
        kc = kappa_conv_meV(kappa_GL=k_gl, n0_um2=n0, rho2=rho2)
        if geom == "tracks_gap":
            kg = kappa_geom_uv * d
        elif geom == "floor":
            kg = floor
        elif geom == "zero":
            kg = 0.0
        else:
            kg = kappa_geom_uv
        kappa = kc + kg
        rec = {
            "alpha": float(a),
            "alpha_over_ac": float(a / locks.alpha_c),
            "delta": float(d),
            "kappa_conv_meV": float(kc),
            "kappa_geom_meV": float(kg),
            "kappa_meV": float(kappa),
            "K_Tstar": K_of_kappa(kappa, T_star, kB),
            "K_T1K": K_of_kappa(kappa, 1.0, kB),
        }
        rows.append(rec)

    def kappa_at(delta: float) -> float:
        kc = kappa_conv_meV(kappa_GL=k_gl, n0_um2=n0, rho2=delta * delta)
        if geom == "tracks_gap":
            kg = kappa_geom_uv * delta
        elif geom == "floor":
            kg = floor
        elif geom == "zero":
            kg = 0.0
        else:
            kg = kappa_geom_uv
        return kc + kg

    k_uv = kappa_at(1.0)
    k_fold = kappa_at(locks.delta_c)
    K_uv = K_of_kappa(k_uv, T_star, kB)
    K_fold = K_of_kappa(k_fold, T_star, kB)
    Kc = K_C_BOSONIC
    hit_before_fold = bool(rows) and (min(r["K_Tstar"] for r in rows) <= Kc)
    # K decreases toward the fold if geom is not a huge frozen floor
    return {
        "geom": geom,
        "kappa_uv_meV": float(k_uv),
        "kappa_fold_meV": float(k_fold),
        "K_uv_Tstar": float(K_uv),
        "K_fold_Tstar": float(K_fold),
        "K_c_bosonic": Kc,
        "K_c_anyon": k_c_anyon(locks.N_layers),
        "KT_hit_before_fold_at_Tstar": hit_before_fold,
        "K_always_above_Kc_Tstar": bool(rows) and min(r["K_Tstar"] for r in rows) > Kc,
        "T_star_K": T_star,
        "kT_star_meV": locks.kT_star_meV,
        "n_rows": len(rows),
        "alpha_over_ac": [r["alpha_over_ac"] for r in rows],
        "delta": [r["delta"] for r in rows],
        "kappa_meV": [r["kappa_meV"] for r in rows],
        "K_Tstar": [r["K_Tstar"] for r in rows],
        "K_T1K": [r["K_T1K"] for r in rows],
    }


def run_c1_metric(locks: L1Locks | None = None, *, n_alpha: int = 48) -> dict[str, Any]:
    locks = locks or build_l1_locks()
    toy = qi_wu_zhang_metric_kappa(Delta_meV=locks.Delta0_meV, lambda_M_um=0.013)
    n0 = n0_healing_disk_um2(locks.a_H_um)
    k_gl = kappa_GL_mev_um2(locks.tqgl)
    k_geom_toy = float(toy["kappa_geom_meV"])
    k_geom_floor = locks.Delta0_meV / (2.0 * math.pi)
    k_conv_uv = 2.0 * k_gl * n0["n0_um2"] * 1.0
    twist = twist_stiffness_uniform(kappa_GL=k_gl, n0_um2=n0["n0_um2"], rho2=1.0)
    frozen = _scan_branch(
        locks, kappa_geom_uv=k_geom_toy, geom="frozen_toy", n_alpha=n_alpha
    )
    tracks = _scan_branch(
        locks, kappa_geom_uv=k_geom_toy, geom="tracks_gap", n_alpha=n_alpha
    )
    conv = _scan_branch(locks, kappa_geom_uv=0.0, geom="zero", n_alpha=n_alpha)
    return {
        "status": "computed",
        "identity": "κ = 2 κ_GL n0 [δ_+]² + κ_geom",
        "kappa_GL_meV_um2": k_gl,
        "n0_healing_disk": n0,
        "quantum_metric_toy": toy,
        "previous_frozen_floor_meV": k_geom_floor,
        "kappa_conv_UV_meV": k_conv_uv,
        "kappa_UV_with_toy_geom_meV": k_conv_uv + k_geom_toy,
        "kappa_UV_with_floor_meV": k_conv_uv + k_geom_floor,
        "twist_check": twist,
        "branch_frozen_toy": frozen,
        "branch_tracks_gap": tracks,
        "branch_conv_only": conv,
        "note": (
            "C1 identity is derived given n0 and κ_geom. The QWZ toy is not a "
            "TEVC Brillouin-zone theorem. Competing-scale: K(α)/T* vs 2/π on "
            "the physical branch. A miss at T* does not kill H-VRG1 as an "
            "implication; it means this C1 estimator does not by itself force "
            "a KT hit before the fold."
        ),
    }
