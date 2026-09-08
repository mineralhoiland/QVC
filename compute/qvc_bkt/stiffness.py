"""C1: microscopic-leaning phase stiffness κ(α) from the two-sector identity.

    κ ≡ ρ_s = 2 κ_GL ⟨ρ²⟩ + κ_geom

This is *not* the old diagnostic κ∝Δ² with an inserted T_BKT^bare.
κ_GL is matched to the truncated GPE kinetic coefficient. ⟨ρ²⟩ comes from
the locked physical fold branch (and optionally from a live GPE field).
κ_geom is a labeled Peotta–Törmä-type estimator, not a derived TEVC
quantum-metric integral.

Defect vocabulary: a hopfion is a Casimir cavity (Q_H=1/N), not a 2D
phase vortex. Using a hopfion |ψ|² as ⟨ρ²⟩ reads the amplitude sector only.
"""

from __future__ import annotations

import math
from typing import Any, Literal

import numpy as np

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.gap_fold import solve_fold_branches  # noqa: E402
from qvc.tqgl.gpe import (  # noqa: E402
    hopfion_initial,
    kinetic_coeff_mev_um2,
    make_grid,
)
from qvc.tqgl.live_gap import mean_density_weighted  # noqa: E402
from qvc.tqgl.locks import FrozenTQGL  # noqa: E402

from qvc_bkt.locks import BKTLocks, build_bkt_locks

GeomKind = Literal["frozen_band", "tracks_gap", "zero"]
DensityKind = Literal["healing_disk", "moire_cell"]


def kappa_GL_mev_um2(lock: FrozenTQGL) -> float:
    """GL gradient weight κ_GL = ħ²/(2m*) = Δ₀ a_H²  [meV·µm²].

    Status: derived given the GPE healing identification a_H = r (Option A*).
    """
    return float(kinetic_coeff_mev_um2(lock))


def n0_healing_disk_um2(a_H_um: float) -> dict[str, Any]:
    """One quantum per π a_H². Ansatz for converting dimensionless GPE density."""
    n0 = 1.0 / (math.pi * a_H_um * a_H_um)
    return {
        "n0_um2": n0,
        "kind": "healing_disk",
        "status": "ansatz",
        "formula": "n0 = 1/(π a_H²)",
        "note": (
            "Converts dimensionless GPE |ψ|²~1 into a 2D density so that "
            "κ = 2 n0 κ_GL ⟨ρ²⟩ has energy units. With this choice, "
            "κ_conv(α=0) = 2Δ₀/π independent of a_H."
        ),
    }


def n0_moire_cell_um2(lambda_M_um: float) -> dict[str, Any]:
    """One pair per hexagonal moiré cell. Alternate density ansatz."""
    area = 0.5 * math.sqrt(3.0) * lambda_M_um * lambda_M_um
    n0 = 1.0 / area
    return {
        "n0_um2": n0,
        "kind": "moire_cell",
        "status": "ansatz",
        "formula": "n0 = 1 / (√3/2 λ_M²)",
        "A_cell_um2": area,
        "note": "Moiré hexagonal cell as the TEVC density unit. Labeled ansatz.",
    }


def kappa_conv_meV(
    *,
    kappa_GL: float,
    n0_um2: float,
    rho2: float,
) -> float:
    """Conventional piece 2 κ_GL n0 ⟨ρ²⟩  [meV].

    Matching κ_GL |∇Ψ|² to (ρ_s/2)|∇θ|² with Ψ=√n0 ρ e^{iθ} gives
    κ = 2 n0 κ_GL ⟨ρ²⟩. Status: derived given the density ansatz for n0.
    """
    return float(2.0 * kappa_GL * n0_um2 * rho2)


def kappa_geom_meV(
    locks: BKTLocks,
    *,
    delta: float,
    kind: GeomKind = "frozen_band",
    C_eff: float = 1.0,
) -> dict[str, Any]:
    """Peotta–Törmä-type geometric weight.

    Frozen-band (default): κ_geom = (C_eff / 2π) Δ₀, independent of the
    running excitonic amplitude. This is the interesting competing-scale
    floor: a Chern-band quantum-metric contribution need not vanish with δ.

    tracks_gap: κ_geom = (C_eff / 2π) Δ(α) — still ansatz, tracks Δ.

    C_eff=1 is a Chern-band OOM. C_eff=Q_H=1/N is a labeled Hopf-fraction
    variant; it does *not* identify a hopfion with a 2D phase vortex.
    """
    pre = C_eff / (2.0 * math.pi)
    if kind == "zero":
        kg = 0.0
        formula = "κ_geom = 0"
    elif kind == "tracks_gap":
        kg = pre * locks.Delta0_meV * delta
        formula = "κ_geom = (C_eff/2π) Δ(α)"
    else:
        kg = pre * locks.Delta0_meV
        formula = "κ_geom = (C_eff/2π) Δ₀  (frozen band geometry)"
    return {
        "kappa_geom_meV": float(kg),
        "kind": kind,
        "C_eff": float(C_eff),
        "status": "ansatz",
        "formula": formula,
        "note": (
            "Peotta–Törmä 2015 geometric superfluid weight, mapped to a "
            "neutral excitonic stiffness. Not a computed Brillouin-zone "
            "integral of the TEVC quantum metric."
        ),
    }


def physical_branch_delta(alpha: float, ell: float) -> float | None:
    br = solve_fold_branches(float(alpha), ell)
    phys = br["phys"]
    return None if phys is None else float(phys)


def T_to_meV(T_K: float, kB: float) -> float:
    return float(T_K * kB)


def K_of_kappa(kappa_meV: float, T_K: float, kB: float) -> float:
    """Kosterlitz K = κ / T_eff  (T_eff in energy units)."""
    t = T_to_meV(T_K, kB)
    if t <= 0.0:
        return float("inf")
    return float(kappa_meV / t)


def twist_stiffness_uniform(
    *,
    kappa_GL: float,
    n0_um2: float,
    rho2: float,
    n_grid: int = 48,
    L_um: float = 0.200,
    q_um_inv: float | None = None,
) -> dict[str, Any]:
    """Current–current (phase-twist) stiffness on a uniform condensate.

    For ψ = √⟨ρ²⟩ exp(i q x) with q a periodic grid mode, the diamagnetic
    piece is exact and the paramagnetic correlator vanishes. Status: computed
    (consistency check of the two-sector algebra on a discrete grid).
    """
    dx = L_um / n_grid
    # q must be a Fourier mode or the field is not periodic.
    q = 2.0 * math.pi / L_um if q_um_inv is None else float(q_um_inv)
    x = np.linspace(-0.5 * L_um, 0.5 * L_um - dx, n_grid)
    X, _Y = np.meshgrid(x, x, indexing="xy")
    amp = math.sqrt(max(rho2, 0.0))
    psi0 = np.full((n_grid, n_grid), amp, dtype=np.complex128)
    psiq = psi0 * np.exp(1j * q * X)
    kx = 2.0 * math.pi * np.fft.fftfreq(n_grid, d=dx)
    KX, KY = np.meshgrid(kx, kx, indexing="xy")
    K2 = KX * KX + KY * KY
    area = L_um * L_um
    ntot = float(n_grid * n_grid)

    def e_kin(psi: np.ndarray) -> float:
        pk = np.fft.fft2(psi)
        # numpy 2D Parseval: ∑_r |x|² = N⁻² ∑_k |FFT x|²
        # ∫ |∇ψ|² dA = dx² ∑_r |∇ψ|² = dx² N⁻² ∑_k k² |FFT ψ|²
        grad2 = float(np.sum(np.abs(pk) ** 2 * K2) / ntot * dx * dx)
        return kappa_GL * grad2

    e0 = e_kin(psi0)
    e1 = e_kin(psiq)
    dE_gpe = e1 - e0  # meV·µm²
    dE = n0_um2 * dE_gpe  # meV
    kappa_twist = 2.0 * dE / (q * q * area)
    kappa_formula = kappa_conv_meV(kappa_GL=kappa_GL, n0_um2=n0_um2, rho2=rho2)
    return {
        "kappa_twist_meV": float(kappa_twist),
        "kappa_formula_meV": float(kappa_formula),
        "rel_err": float(abs(kappa_twist - kappa_formula) / max(abs(kappa_formula), 1e-30)),
        "q_um_inv": float(q),
        "n_grid": int(n_grid),
        "status": "computed",
        "note": (
            "Static phase twist on a uniform field, q=2π/L (periodic). "
            "Equals 2 n0 κ_GL ⟨ρ²⟩ up to discretisation; paramagnetic χ_JJ "
            "vanishes for this state."
        ),
    }


def hopfion_rho2_weighted(locks: BKTLocks, *, n_grid: int = 64) -> dict[str, Any]:
    """Ring-weighted ⟨|ψ|²⟩ of the truncated-GPE hopfion initial data.

    Status: computed amplitude-sector density. The hopfion is not a 2D
    phase vortex; this is only ⟨ρ²⟩ for the cavity condensate.
    """
    grid = make_grid(locks.tqgl, n_grid=n_grid, L_um=0.200)
    psi = hopfion_initial(grid, locks.tqgl, n0=1.0, noise_amp=0.0, seed=0)
    rho2 = mean_density_weighted(psi, grid.weight)
    return {
        "rho2_hopfion_ring": float(rho2),
        "n_grid": int(n_grid),
        "status": "computed",
        "defect": "hopfion_cavity_amplitude_not_2d_vortex",
        "Q_H": float(locks.Q_H),
        "note": (
            "Live-gap / hopfion |ψ|² is blind to 2D phase-vortex unbinding. "
            "Used only as an amplitude-sector ⟨ρ²⟩ estimator at the UV IC."
        ),
    }


def scan_physical_branch(
    locks: BKTLocks | None = None,
    *,
    n_alpha: int = 80,
    density: DensityKind = "healing_disk",
    geom: GeomKind = "frozen_band",
    C_eff: float = 1.0,
    T_grid_K: tuple[float, ...] | None = None,
) -> dict[str, Any]:
    """κ(α), K(α)=κ/T_eff on the locked physical branch. Frozen α (Option A*)."""
    locks = locks or build_bkt_locks()
    k_gl = kappa_GL_mev_um2(locks.tqgl)
    if density == "moire_cell":
        dens = n0_moire_cell_um2(locks.lambda_M_um)
    else:
        dens = n0_healing_disk_um2(locks.r_um)
    n0 = float(dens["n0_um2"])
    T_star = locks.T_star_K
    if T_grid_K is None:
        T_grid_K = (T_star, 0.020, 0.050, 0.250, 0.500, 1.000)
    T_grid_K = tuple(float(t) for t in T_grid_K)

    alphas = np.linspace(1e-6, 0.995 * locks.alpha_c, n_alpha)
    rows: list[dict[str, Any]] = []
    for a in alphas:
        d = physical_branch_delta(float(a), locks.ell)
        if d is None:
            continue
        rho2 = d * d  # live-gap map: Δ/Δ₀ = √(n/n_ref) ⇒ ⟨ρ²⟩/n_ref = δ²
        kc = kappa_conv_meV(kappa_GL=k_gl, n0_um2=n0, rho2=rho2)
        kg_info = kappa_geom_meV(locks, delta=d, kind=geom, C_eff=C_eff)
        kg = float(kg_info["kappa_geom_meV"])
        kappa = kc + kg
        tA = (locks.alpha_c - a) / locks.alpha_c
        rec: dict[str, Any] = {
            "alpha": float(a),
            "alpha_over_ac": float(a / locks.alpha_c),
            "t_A": float(tA),
            "delta": float(d),
            "Delta_meV": float(d * locks.Delta0_meV),
            "rho2": float(rho2),
            "kappa_conv_meV": float(kc),
            "kappa_geom_meV": float(kg),
            "kappa_meV": float(kappa),
        }
        for T in T_grid_K:
            rec[f"K_T{T:.4g}"] = K_of_kappa(kappa, T, locks.kB_meV_per_K)
        rec["K_Tstar"] = K_of_kappa(kappa, T_star, locks.kB_meV_per_K)
        rows.append(rec)

    # Analytic UV and fold anchors
    d0, dc = 1.0, locks.delta_c
    kg0 = float(kappa_geom_meV(locks, delta=d0, kind=geom, C_eff=C_eff)["kappa_geom_meV"])
    kgc = float(kappa_geom_meV(locks, delta=dc, kind=geom, C_eff=C_eff)["kappa_geom_meV"])
    k_uv = kappa_conv_meV(kappa_GL=k_gl, n0_um2=n0, rho2=1.0) + kg0
    k_fold = kappa_conv_meV(kappa_GL=k_gl, n0_um2=n0, rho2=dc * dc) + kgc

    twist = twist_stiffness_uniform(kappa_GL=k_gl, n0_um2=n0, rho2=1.0)
    hopf = hopfion_rho2_weighted(locks)

    closures = {
        "two_sector": {
            "formula": "κ ≡ ρ_s = 2 κ_GL ⟨ρ²⟩ + κ_geom",
            "status": "derived",
            "convention": "Kosterlitz K=κ/T_eff so K_c=2/π",
        },
        "kappa_GL": {
            "value_meV_um2": k_gl,
            "formula": "κ_GL = Δ₀ a_H² with a_H=r (Option A* frozen tube)",
            "status": "derived",
        },
        "rho2": {
            "formula": "⟨ρ²⟩ = [δ_+(α)]² from locked gap map / live-gap map",
            "status": "computed",
            "source": "qvc.gap_fold.solve_fold_branches physical branch",
        },
        "n0": dens,
        "kappa_geom": kappa_geom_meV(locks, delta=1.0, kind=geom, C_eff=C_eff),
        "not_used": {
            "T_BKT_bare": "NOT an input. Old diagnostic begged the question.",
            "kappa_propto_Delta2_only": (
                "Conventional piece does scale as δ², but κ_geom need not, "
                "and K is κ/T_eff from computed-or-ansatz κ, not from T_BKT^bare."
            ),
        },
        "option_A_star": "cavity frozen; α does not run with KT ℓ",
        "CS": "truncated",
    }

    return {
        "task": "C1",
        "locks_T_star_K": T_star,
        "kT_star_meV": locks.kT_star_meV,
        "density_kind": density,
        "geom_kind": geom,
        "C_eff": float(C_eff),
        "kappa_GL_meV_um2": k_gl,
        "n0_um2": n0,
        "kappa_uv_meV": float(k_uv),
        "kappa_fold_meV": float(k_fold),
        "K_uv_Tstar": K_of_kappa(k_uv, T_star, locks.kB_meV_per_K),
        "K_fold_Tstar": K_of_kappa(k_fold, T_star, locks.kB_meV_per_K),
        "K_c_bosonic": locks.k_c_bosonic,
        "T_grid_K": list(T_grid_K),
        "rows": rows,
        "twist_uniform": twist,
        "hopfion_rho2": hopf,
        "closures": closures,
        "status": {
            "kappa_conv": "derived_given_n0_ansatz",
            "kappa_geom": "ansatz",
            "rho2_branch": "computed",
            "twist_check": "computed",
            "theorem": False,
        },
    }


def run_c1(
    locks: BKTLocks | None = None,
    *,
    n_alpha: int = 80,
) -> dict[str, Any]:
    """C1 driver: healing-disk + frozen geometric floor, plus labeled variants."""
    locks = locks or build_bkt_locks()
    primary = scan_physical_branch(
        locks, n_alpha=n_alpha, density="healing_disk", geom="frozen_band", C_eff=1.0
    )
    variants = {
        "conv_only": scan_physical_branch(
            locks, n_alpha=n_alpha, density="healing_disk", geom="zero", C_eff=1.0
        ),
        "geom_QH": scan_physical_branch(
            locks,
            n_alpha=n_alpha,
            density="healing_disk",
            geom="frozen_band",
            C_eff=locks.Q_H,
        ),
        "moire_plus_geom": scan_physical_branch(
            locks, n_alpha=n_alpha, density="moire_cell", geom="frozen_band", C_eff=1.0
        ),
    }
    # compact variant table (anchors only)
    compact = {}
    for name, sc in variants.items():
        compact[name] = {
            "kappa_uv_meV": sc["kappa_uv_meV"],
            "kappa_fold_meV": sc["kappa_fold_meV"],
            "K_uv_Tstar": sc["K_uv_Tstar"],
            "K_fold_Tstar": sc["K_fold_Tstar"],
            "density_kind": sc["density_kind"],
            "geom_kind": sc["geom_kind"],
            "C_eff": sc["C_eff"],
            "rows": sc["rows"],
        }
    return {
        "task": "C1",
        "primary": primary,
        "variants": compact,
        "primary_closure": (
            "healing-disk n0=1/(π a_H²) + Peotta–Törmä frozen κ_geom=(1/2π)Δ₀. "
            "C1 is partly ansatz (n0, κ_geom): not a BKT theorem."
        ),
    }
