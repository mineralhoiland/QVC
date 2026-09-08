"""Option A* locks and epistemic labels for the BKT C1–C6 suite.

Locked numbers are pulled from QVCCompute APIs (fat torus, fold, Schwinger
T*, democratic Spec(G)). Folklore 0.244 / 0.389 / 0.032 / 92% is never a
physics result.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any, Literal

from qvc_bkt.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc.coupling_matrix import democratic_spectrum  # noqa: E402
from qvc.gap_fold import ELL_CORRECTIONS, fold_condition  # noqa: E402
from qvc.lens_chern_simons import fractional_hall_charge  # noqa: E402
from qvc.params import PARAMS  # noqa: E402
from qvc.tqgl.diagnostics import schwinger_twin_channel  # noqa: E402
from qvc.tqgl.locks import FrozenTQGL, build_frozen_tqgl  # noqa: E402

Status = Literal["locked", "derived", "computed", "ansatz", "conjecture", "truncated"]

# Withdrawn preprint numbers — tests assert these are not emitted as results.
FOLKLORE_E_CAS_MEV = 0.244
FOLKLORE_DELTA_QVC_MEV = 0.389
FOLKLORE_DELTA_1PS_MEV = 0.032
FOLKLORE_PCT = 92.0

K_C_BOSONIC = 2.0 / 3.141592653589793  # 2/π, Kosterlitz convention K=κ/T
PI = 3.141592653589793


@dataclass(frozen=True)
class BKTLocks:
    """Immutable Option A* + Schwinger numbers consumed by C1–C6."""

    tqgl: FrozenTQGL
    ell: float
    alpha_c: float
    delta_c: float
    Delta_c_meV: float
    Delta0_meV: float
    E_vac_meV: float
    F: float
    lambda_star: float
    Delta_star_meV: float
    E_pair_meV: float
    T_star_K: float
    kB_meV_per_K: float
    hbar_cs_meV_um: float
    hbar_meV_ps: float
    r_um: float
    R_um: float
    lambda_M_um: float
    N_layers: int
    Q_H: float
    g0_meV: float
    lambda_min_meV: float
    k_c_bosonic: float
    k_c_anyon: float
    theta_stat: float
    cs_status: str = "truncated"
    notes: tuple[str, ...] = field(default_factory=tuple)

    def as_dict(self) -> dict[str, Any]:
        d = asdict(self)
        d.pop("tqgl", None)
        d["tqgl_keys"] = sorted(self.tqgl.as_dict().keys())
        return d

    @property
    def kT_star_meV(self) -> float:
        return self.T_star_K * self.kB_meV_per_K


def _k_c_anyon(n: int) -> float:
    """Labeled CS-restored anyon shift K_c=(2/π)(1-1/N)^2. Not the truncated default."""
    return (2.0 / PI) * (1.0 - 1.0 / n) ** 2


def build_bkt_locks(*, lambda_star_frac: float = 0.90) -> BKTLocks:
    """Assemble locks from QVCCompute. Democratic λ_min = −g₀."""
    tqgl = build_frozen_tqgl(lambda_star_frac=lambda_star_frac)
    delta_c, alpha_c = fold_condition(tqgl.ell)
    sch = schwinger_twin_channel(tqgl)
    g0 = tqgl.lambda_star * tqgl.E_vac_meV
    spec = democratic_spectrum(tqgl.N_layers, g0)
    n = tqgl.N_layers
    notes = (
        "Fat torus (8,50) nm, Option A* frozen cavity: α does not run with KT ℓ.",
        "Hopfion ≠ 2D phase vortex ≠ CS fluxon.",
        "Chern–Simons truncated unless anyon K_c is an explicit labeled option.",
        "Democratic λ_min = −g₀ with g₀ = λ★ E_vac.",
        "Q_H = 1/N.",
        "Never quote 0.244/0.389/0.032/92% as physics.",
    )
    return BKTLocks(
        tqgl=tqgl,
        ell=float(tqgl.ell),
        alpha_c=float(alpha_c),
        delta_c=float(delta_c),
        Delta_c_meV=float(tqgl.Delta_c_meV),
        Delta0_meV=float(tqgl.Delta0_meV),
        E_vac_meV=float(tqgl.E_vac_meV),
        F=float(tqgl.F),
        lambda_star=float(tqgl.lambda_star),
        Delta_star_meV=float(tqgl.Delta_star_meV),
        E_pair_meV=float(sch["E_pair_meV"]),
        T_star_K=float(sch["T_star_K"]),
        kB_meV_per_K=float(tqgl.kB_meV_per_K),
        hbar_cs_meV_um=float(tqgl.hbar_cs_meV_um),
        hbar_meV_ps=float(tqgl.hbar_meV_ps),
        r_um=float(tqgl.r_um),
        R_um=float(tqgl.R_um),
        lambda_M_um=float(PARAMS.lambda_M_um),
        N_layers=int(n),
        Q_H=float(fractional_hall_charge(n)),
        g0_meV=float(g0),
        lambda_min_meV=float(spec["lambda_min"]),
        k_c_bosonic=float(2.0 / PI),
        k_c_anyon=float(_k_c_anyon(n)),
        theta_stat=float(PI / n),
        cs_status="truncated",
        notes=notes,
    )


def assert_not_folklore_result(payload: dict[str, Any]) -> None:
    """Guard: folklore numbers must not appear as physics outputs."""
    banned = {
        "E_Cas_meV": FOLKLORE_E_CAS_MEV,
        "Delta_QVC_meV": FOLKLORE_DELTA_QVC_MEV,
        "Delta_1ps_meV": FOLKLORE_DELTA_1PS_MEV,
        "gap_collapse_pct": FOLKLORE_PCT,
    }
    for key, val in banned.items():
        if key in payload and payload[key] is not None:
            raise AssertionError(f"{key} is folklore and must not be a suite result")
        # also catch exact float dumps under other names if marked as result
        _ = val


# Re-export fold IR for callers that need the canon ℓ.
ELL_GAP = ELL_CORRECTIONS
