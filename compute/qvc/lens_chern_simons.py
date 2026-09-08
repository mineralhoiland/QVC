"""Lens-space Chern–Simons / linking fractionalization — Tier C lock.

CLOSED-FORM LOCK
----------------
CS[A_k] = k²/N, lk=1/N, Q_H=1/N, θ=π/N
H₂(L(N,1))=0 — 2-cycle proof path is INVALID.
"""

from __future__ import annotations

import math


def chern_simons_flat_u1(k: int, N: int) -> float:
    """CS[A_k] = k² / N for flat U(1) on L(N,1)."""
    if N <= 0:
        raise ValueError("N must be positive")
    return (k * k) / float(N)


def linking_form_self(N: int) -> float:
    """ℓk([γ],[γ]) = 1/N mod 1 on H₁ = ℤ_N."""
    if N <= 0:
        raise ValueError("N must be positive")
    return 1.0 / float(N)


def fractional_hall_charge(N: int) -> float:
    """Q_H = 1/N — physical identification used in QVC."""
    return 1.0 / float(N)


def anyonic_exchange_angle(N: int) -> float:
    """θ = π / N from the linking form (abelian anyon program prediction)."""
    return math.pi / float(N)


def h2_vanishes_note(N: int = 5) -> dict[str, str]:
    """Documentary lock: H₂(L(N,1))=0 invalidates 2-cycle proofs."""
    return {
        "space": f"L({N},1)",
        "H1": f"Z_{N}",
        "H2": "0",
        "status": "PROVEN (standard algebraic topology)",
        "consequence": (
            "Preprint 2-cycle fractionalization proof path is INVALID; "
            "use CS/linking on H1 instead."
        ),
    }
