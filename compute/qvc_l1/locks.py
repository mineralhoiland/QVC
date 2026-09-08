"""Option A* locks for the L1 G_ij / Lindblad programs.

BRIDGE identification (algebra, not a kernel number):
    g0★ = λ★ E_vac ≈ 0.0982 meV.

The overlap mean off-diagonal g0_eff is a *computed* kernel number in
Macdonald-K0 units. Do not force g0_eff = g0★.

Democratic Spec(G) is locked algebra: λ_max = 4 g0, λ_min = −g0 (mult 4)
for N=5. Computed G must not insert J−I by hand.

Never quote 0.244 / 0.389 / 0.032 / 92% as physics.
Hopfion ≠ 2D phase vortex ≠ CS fluxon.
CS / Berry / Skyrme stay truncated in the GPE.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.coupling_matrix import democratic_spectrum  # noqa: E402
from qvc.tqgl.locks import FrozenTQGL, build_frozen_tqgl  # noqa: E402

Status = Literal["locked", "derived", "computed", "ansatz", "truncated"]

FOLKLORE_E_CAS_MEV = 0.244
FOLKLORE_DELTA_QVC_MEV = 0.389
FOLKLORE_DELTA_1PS_MEV = 0.032
FOLKLORE_PCT = 92.0

# Published coherent 1 ps truncated GPE+DCE baseline (n_grid=128).
COHERENT_DELTA_LIVE_1PS_MEV = 0.6305
COHERENT_G1_1PS = 0.8832
LIVE_GAP_1PS_MEV = COHERENT_DELTA_LIVE_1PS_MEV
G1_1PS = COHERENT_G1_1PS

# Track A dictionary guards (not physics results).
RETIRED_G0_MEV = 0.21
RETIRED_KAPPA_TOR = 0.25
CARTOON_ALPHA0 = -0.5
CARTOON_ALPHA_QVC = -0.35
CARTOON_BARRIER_PCT = 42.0
HOLOGRAPHIC_GAMMA_QM_MEV_NM2 = 3.73
ETA_RET_ALTERNATE = 3.1
XI_ALTERNATE_UM = 0.181
TSTAR_WITHDRAWN_K = 0.285
R_WS_WITHDRAWN_UM = 0.91
BERRY_R_MEAN_RATIO = 0.908


@dataclass(frozen=True)
class L1Locks:
    """Immutable Option A* numbers consumed by B1 and B2."""

    tqgl: FrozenTQGL
    r_nm: float
    R_nm: float
    hbar_cs_meV_um: float
    Delta0_meV: float
    E_vac_meV: float
    lambda_star: float
    g0_star_meV: float
    N_layers: int
    Q_H: float
    dem_lambda_max_over_g0: float
    dem_lambda_min_over_g0: float
    dem_lambda_min_mult: int
    cs_status: str = "truncated"
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.pop("tqgl", None)
        d["tqgl_keys"] = sorted(self.tqgl.as_dict().keys())
        return d

    @property
    def r_um(self) -> float:
        return 1e-3 * self.r_nm

    @property
    def R_um(self) -> float:
        return 1e-3 * self.R_nm

    @property
    def a_H_um(self) -> float:
        """Option A* GL/GPE healing freeze: a_H = r."""
        return self.r_um

    @property
    def F(self) -> float:
        return float(self.tqgl.F)

    @property
    def ell(self) -> float:
        return float(self.tqgl.ell)

    @property
    def alpha_c(self) -> float:
        return float(self.tqgl.alpha_c)

    @property
    def delta_c(self) -> float:
        return float(self.tqgl.delta_c)

    @property
    def Delta_c_meV(self) -> float:
        return float(self.tqgl.Delta_c_meV)

    @property
    def lambda_fold(self) -> float:
        return float(self.tqgl.lambda_fold)

    @property
    def Delta_star_meV(self) -> float:
        return float(self.tqgl.Delta_star_meV)

    @property
    def kB_meV_per_K(self) -> float:
        return float(self.tqgl.kB_meV_per_K)

    @property
    def hbar_meV_ps(self) -> float:
        return float(self.tqgl.hbar_meV_ps)

    @property
    def g0_fold_meV(self) -> float:
        return float(self.lambda_fold * self.E_vac_meV)

    @property
    def lambda_min_star_meV(self) -> float:
        return float(-self.g0_star_meV)

    @property
    def kT_star_meV(self) -> float:
        return self.T_star_K * self.kB_meV_per_K

    @property
    def T_star_mK(self) -> float:
        return 1e3 * self.T_star_K

    def twin_channel_light(self) -> dict[str, float]:
        """E_pair and T* from democratic g0★. Does not run WS G_ij / Berry."""
        import math

        g0 = self.g0_star_meV
        e_pair = 2.0 * abs(g0)
        a_c_sch = self.Delta0_meV / self.hbar_cs_meV_um
        ratio = self.tqgl.A_star_um_inv / a_c_sch if a_c_sch else float("nan")
        sch_exp = math.pi / ratio if ratio > 0.0 else float("inf")
        t_star = e_pair / (self.kB_meV_per_K * sch_exp) if math.isfinite(sch_exp) else float("nan")
        return {
            "E_pair_meV": float(e_pair),
            "T_star_K": float(t_star),
            "schwinger_exponent": float(sch_exp),
        }

    @property
    def E_pair_meV(self) -> float:
        return self.twin_channel_light()["E_pair_meV"]

    @property
    def T_star_K(self) -> float:
        return self.twin_channel_light()["T_star_K"]


def build_l1_locks(*, lambda_star_frac: float = 0.90) -> L1Locks:
    """Assemble Option A* locks. Democratic spectrum is algebra, not a GPE term."""
    tqgl = build_frozen_tqgl(lambda_star_frac=lambda_star_frac)
    n = int(tqgl.N_layers)
    g0_star = float(tqgl.lambda_star * tqgl.E_vac_meV)
    spec = democratic_spectrum(n, 1.0)
    notes = (
        "Fat torus (r,R)=(8,50) nm, ħ c_s=0.07 meV·µm, Δ0=0.600 meV.",
        "g0★ = λ★ E_vac is BRIDGE algebra. Overlap g0_eff is a kernel number.",
        "Democratic Spec(G): λ_max=(N-1)g0, λ_min=−g0 (mult N−1). Not inserted into G_ij.",
        "Hopfion ≠ 2D phase vortex ≠ CS fluxon.",
        "CS / Berry / Skyrme truncated in the GPE.",
        "Never quote 0.244/0.389/0.032/92% as physics.",
    )
    return L1Locks(
        tqgl=tqgl,
        r_nm=float(tqgl.r_nm),
        R_nm=float(tqgl.R_nm),
        hbar_cs_meV_um=float(tqgl.hbar_cs_meV_um),
        Delta0_meV=float(tqgl.Delta0_meV),
        E_vac_meV=float(tqgl.E_vac_meV),
        lambda_star=float(tqgl.lambda_star),
        g0_star_meV=g0_star,
        N_layers=n,
        Q_H=float(tqgl.Q_H),
        dem_lambda_max_over_g0=float(spec["lambda_max"]),
        dem_lambda_min_over_g0=float(spec["lambda_min"]),
        dem_lambda_min_mult=int(spec["lambda_min_mult"]),
        cs_status="truncated",
        notes=notes,
    )


def folklore_values() -> tuple[float, float, float, float]:
    return (
        FOLKLORE_E_CAS_MEV,
        FOLKLORE_DELTA_QVC_MEV,
        FOLKLORE_DELTA_1PS_MEV,
        FOLKLORE_PCT,
    )
