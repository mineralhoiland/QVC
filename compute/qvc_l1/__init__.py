"""Layer-1 programs: Track A parameter dictionary, G_ij overlaps, Lindblad GPE.

Package name is ``qvc_l1`` so it does not shadow QVCCompute's ``qvc``.
Independent of Hubbard, CS-in-EOM, H-TOR1 vertex, and holography.

Submodules
----------
- ``dictionary``: Track A L1 parameter card (lengths, GL, catalyzon, Hall)
- ``gij_ring``: WS Gaussian K0 overlap (finer grid) and ring Fourier/Wannier G_ij
- ``lindblad``: number-damping / dephasing on the truncated split-step GPE
- ``suite``: B1/B2 driver writing summary.json + figures
"""

from qvc_l1.bootstrap import ensure_qvccompute

ensure_qvccompute()

from qvc_l1.locks import L1Locks, build_l1_locks  # noqa: E402

__all__ = [
    "L1Locks",
    "build_l1_locks",
    "ensure_qvccompute",
]
