"""Make QVCCompute's ``qvc`` importable without shadowing this package.

QVCCursor and QVCCompute both cannot own the top-level name ``qvc``.
This suite lives as ``qvc_bkt`` and imports solvers from QVCCompute.
"""

from __future__ import annotations

import sys
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_DEFAULTS = (
    _HERE.parent,  # sibling of qvc_bkt (this compute/ tree)
    Path("/Users/mineralhoiland/Code/QVCCompute"),
    Path(__file__).resolve().parents[2] / "QVCCompute",
    Path(__file__).resolve().parents[1].parent / "QVCCompute",
)


def ensure_qvccompute(path: str | Path | None = None) -> Path:
    """Insert QVCCompute on ``sys.path`` (first) and return the root used."""
    candidates: list[Path] = []
    if path is not None:
        candidates.append(Path(path))
    candidates.extend(_DEFAULTS)
    for root in candidates:
        if (root / "qvc" / "gap_fold.py").is_file():
            s = str(root)
            if s in sys.path:
                sys.path.remove(s)
            sys.path.insert(0, s)
            return root
    raise ImportError(
        "QVCCompute not found. Expected qvc/gap_fold.py under "
        + ", ".join(str(p) for p in candidates)
    )


ensure_qvccompute()
