"""Put the compute/ root on sys.path so qvc, qvc_l1, and qvc_bkt import."""

from __future__ import annotations

import sys
from pathlib import Path

COMPUTE = Path(__file__).resolve().parents[1]
if str(COMPUTE) not in sys.path:
    sys.path.insert(0, str(COMPUTE))
