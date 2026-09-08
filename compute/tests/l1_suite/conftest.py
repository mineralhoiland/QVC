"""Ensure compute/ is on sys.path; qvc is the sibling package in this snapshot."""

from __future__ import annotations

import sys
from pathlib import Path

COMPUTE = Path(__file__).resolve().parents[2]
if str(COMPUTE) not in sys.path:
    sys.path.insert(0, str(COMPUTE))

from qvc_l1.bootstrap import ensure_qvccompute  # noqa: E402

ensure_qvccompute()
