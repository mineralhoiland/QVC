# Manifest

Compiled 8 September 2026 for the private GitHub repository `mineralhoiland/QVC`.

## Included

| Path | Source | License |
| --- | --- | --- |
| `paper/` | `QVCCursor/tex/QVC_arXiv_v4/` | CC BY 4.0 |
| `locks/` | `macros.tex` + `QVCCompute/qvc/locks/` | CC BY 4.0 |
| `compute/qvc/` | `QVCCompute/qvc/` | MIT |
| `compute/qvc_l1/` | `QVCCursor/qvc_l1/` | MIT |
| `compute/qvc_bkt/` | `QVCCursor/qvc_bkt/` | MIT |
| `compute/scripts/` | `QVCCompute/scripts/` | MIT |
| `compute/tests/` | both trees, namespaced | MIT |
| `compute/notebooks/` | both trees | MIT |
| `verification/` | `QVCCompute/output/` + `QVCCursor/output/l1_*` | CC BY 4.0 |
| `structures/` | `QVCCursor/TEVC_structures/` | CC BY 4.0 |
| `sim/` | `sim-suites/` (no `node_modules`) + `QVCCompute/web/tdvt3d/` | MIT |
| `docs/compute-notes/` | `QVCCompute/docs/` | CC BY 4.0 |
| `docs/companions/QVC_TDVT_Convergence_Hypotheses.tex` | `QVCCursor/tex/` | CC BY 4.0 |

VTK/PVD GPE replay dumps from `QVCCompute/output/crystal_suite/` are omitted
(size). CIF/JSON reports are kept under `verification/`.

## Excluded (on purpose)

| Material | Why |
| --- | --- |
| `garage/fusor/` | IEC safety notes, not the preprint |
| `releases/tabd-tdvt-warp/` | Separate public repo `TABD-TDVT-Warp` |
| `docs/preprint_revC/` | Unified Hopfion Spacetime, not QVC v4 |
| `QVCCompute/tex/` v3 Tier A | Superseded by `paper/` v4 |
| `sim-suites/node_modules/` | Reinstall with `npm install` |
| `QVCCompute/.venv/` | Local environment |
| Secrets, `.env`, credentials | Never ship |

Some files in `docs/compute-notes/` still say “provisional freeze / do not
quote meV.” That language is historical. The v4 lock is Option A* in
`locks/PARAMETERS.md`.
