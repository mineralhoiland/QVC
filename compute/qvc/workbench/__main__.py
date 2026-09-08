"""``python -m qvc.workbench [--port 5006] [--fields-grid 40] [--gpe-grid 64]``."""

from __future__ import annotations

import argparse


def main() -> None:
    p = argparse.ArgumentParser(description="TEVC / QVC–TDVT Panel workbench")
    p.add_argument("--port", type=int, default=5006)
    p.add_argument("--show", action="store_true", help="open a browser tab")
    p.add_argument("--fields-grid", type=int, default=40, help="hopfion grid per axis (40 interactive, 64 acceptance)")
    p.add_argument("--gpe-grid", type=int, default=64, help="GPE grid (64 interactive, 128 lock run)")
    p.add_argument("--gpe-t-ps", type=float, default=0.5)
    p.add_argument("--gpe-frames", type=int, default=12)
    a = p.parse_args()
    from qvc.workbench.app import serve

    serve(port=a.port, show=a.show, fields_grid=a.fields_grid, gpe_grid=a.gpe_grid, gpe_t_ps=a.gpe_t_ps, gpe_frames=a.gpe_frames)


if __name__ == "__main__":
    main()
