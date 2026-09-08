#!/usr/bin/env python3
"""Regenerate vortex-RG PGFPlots tables for QVC_arXiv_v4.

Run from anywhere:
  PYTHONPATH=/Users/mineralhoiland/Code/QVCCompute python tex/QVC_arXiv_v4/vortex_rg_run.py
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
QVCCOMPUTE = Path("/Users/mineralhoiland/Code/QVCCompute")
if str(QVCCOMPUTE) not in sys.path:
    sys.path.insert(0, str(QVCCOMPUTE))

from qvc.vortex_rg import run_diagnostic  # noqa: E402


def main() -> None:
    data = HERE / "data"
    summary = run_diagnostic(data, do_mc=True)
    print("wrote", data)
    print("path_scan_kt_before_fold", summary["path_scan_kt_before_fold"])
    print("coupled_from_star_winner", summary["coupled_from_star_winner"])
    print("coupled_from_uv_winner", summary["coupled_from_uv_winner"])
    print("T_eff_fold_win_K", summary["T_eff_fold_win_K"])
    print("alpha_kt_bare_over_ac", summary["alpha_kt_bare_over_ac"])
    print("C5", summary["loglik"])


if __name__ == "__main__":
    main()
