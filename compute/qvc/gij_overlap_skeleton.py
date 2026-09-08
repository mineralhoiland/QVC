"""P1 starter: G_ij overlap skeleton vs democratic Spec lock.

Does NOT claim a microscopic derivation yet. Builds:
  1. Democratic G = g0 (J - I) and locked spectrum (verification)
  2. Toy WS-cell overlap quadrature → Hermitian 5×5
  3. Acceptance tests for λ_min < 0 and distance to democratic Spec

Run:
  .venv/bin/python -m qvc.gij_overlap_skeleton
"""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from qvc.params import PARAMS


def democratic_G(N: int, g0: float) -> np.ndarray:
    """G = g0 (J - I): all-to-all equal off-diagonal couplings."""
    J = np.ones((N, N), dtype=float)
    return g0 * (J - np.eye(N))


def locked_spec(N: int, g0: float) -> dict[str, Any]:
    """Closed-form Spec(G): λ_max=(N-1)g0, λ_min=-g0 with mult. N-1."""
    return {
        "lambda_max": (N - 1) * g0,
        "lambda_min": -g0,
        "lambda_min_mult": N - 1,
        "trace": 0.0,
    }


def numeric_spec(G: np.ndarray) -> dict[str, Any]:
    w = np.linalg.eigvalsh(G)
    w = np.sort(w)
    # multiplicity of near-min
    lam_min = float(w[0])
    mult = int(np.sum(np.isclose(w, lam_min, rtol=0, atol=1e-10)))
    return {
        "eigenvalues": w.tolist(),
        "lambda_max": float(w[-1]),
        "lambda_min": lam_min,
        "lambda_min_mult": mult,
        "trace": float(np.trace(G)),
    }


def toy_ws_overlap_G(
    N: int = 5,
    *,
    g0_scale: float = 0.21,
    seed: int = 0,
    offdiag_jitter: float = 0.05,
) -> np.ndarray:
    """Placeholder WS-cell overlap: near-democratic with small jitter.

    Replace with TQGL+ZPF correlator quadrature when P1 derivation lands.
    """
    rng = np.random.default_rng(seed)
    G = democratic_G(N, g0_scale)
    # Symmetric real jitter on off-diagonals (keeps Hermiticity)
    noise = rng.normal(0.0, offdiag_jitter * g0_scale, size=(N, N))
    noise = 0.5 * (noise + noise.T)
    np.fill_diagonal(noise, 0.0)
    return G + noise


def compare_to_democratic(G: np.ndarray, g0: float) -> dict[str, Any]:
    N = G.shape[0]
    locked = locked_spec(N, g0)
    num = numeric_spec(G)
    return {
        "locked": locked,
        "numeric": num,
        "lambda_min_negative": num["lambda_min"] < 0,
        "rel_err_lambda_max": abs(num["lambda_max"] - locked["lambda_max"])
        / abs(locked["lambda_max"]),
        "rel_err_lambda_min": abs(num["lambda_min"] - locked["lambda_min"])
        / abs(locked["lambda_min"]),
        "frobenius_to_democratic": float(
            np.linalg.norm(G - democratic_G(N, g0), ord="fro")
        ),
        "acceptance": {
            "lambda_min_lt_0": num["lambda_min"] < 0,
            "trace_near_0": abs(num["trace"]) < 1e-9 * max(1.0, abs(g0) * N),
        },
    }


def run_p1_skeleton() -> dict[str, Any]:
    N = PARAMS.N_layers
    g0 = 0.21  # meV working estimate; not frozen from vacuum
    G_exact = democratic_G(N, g0)
    exact = compare_to_democratic(G_exact, g0)
    G_toy = toy_ws_overlap_G(N, g0_scale=g0)
    toy = compare_to_democratic(G_toy, g0)
    return {
        "status": "skeleton_not_a_derivation",
        "N": N,
        "g0_meV": g0,
        "democratic_exact": exact,
        "toy_ws_jitter": toy,
        "next_derivation_steps": [
            "Discretize TQGL order-parameter modes on WS cell",
            "Build ZPF correlator ⟨φ(r)φ(r′)⟩ from Tier-C torus modes",
            "G_ij = ∫ ψ_i*(r) K_ZPF(r,r′) ψ_j(r′) d²r d²r′",
            "Accept only if λ_min < 0 without inserting g0(J-I) by hand",
        ],
    }


if __name__ == "__main__":
    import json

    print(json.dumps(run_p1_skeleton(), indent=2))
