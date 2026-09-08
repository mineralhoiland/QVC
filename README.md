# Quantum Vacuum Catalyzation (QVC)

Private research compilation of the **QVC / TEVC** arXiv v4 preprint: manuscript
sources, locked Option A* parameters, Python verification, L1 dictionary and
vortex-sector codes, crystal structures, and visualization suites.

This archive is the freeze to cite in the data-availability statement. Visual
simulations under `sim/` illustrate the geometry; they are **not** the numerical
proof. The lock is `locks/optionA_star.json` plus `compute/qvc/` and
`scripts/run_verification.py`.

**License:** [MIT](LICENSE-MIT) (code) · [CC BY 4.0](LICENSE-CC-BY-4.0) (paper,
figures, docs, locks). See [LICENSE](LICENSE). Snapshot version: [`VERSION`](VERSION).

Cite with [`CITATION.cff`](CITATION.cff).

## What this archive is

| Layer | Contents | Role |
| --- | --- | --- |
| **Paper** | `paper/` | arXiv v4 TeX, bibliography, figures, supplements |
| **Locks** | `locks/` | Frozen Option A* numbers matching `paper/macros.tex` |
| **Compute** | `compute/` | `qvc`, `qvc_l1`, `qvc_bkt`, tests, scripts, notebooks |
| **Verification** | `verification/` | Last verification report, TQGL numbers, L1 summaries |
| **Structures** | `structures/` | Graphene / hopfion CIF and cube files |
| **Sim** | `sim/` | Three.js TDVT viewers (visualization only) |

TABD analogue-warp, IEC fusor notes, and the Unified Hopfion Spacetime
Revision C manuscript are **not** bundled. TABD lives in its own repository:
[TABD-TDVT-Warp](https://github.com/mineralhoiland/TABD-TDVT-Warp).

## Locked Option A* (do not retune)

Fat torus \((r,R)=(8,50)\) nm, \(\hbar c_s=0.07\,\mathrm{meV\cdot\mu m}\):

| Quantity | Value |
| --- | --- |
| \(E_{\mathrm{vac}}\) | \(0.1287\) meV |
| \(F\) | \(0.58068\) |
| \(\ell\) | \(-2.7958\) |
| \(\alpha_c=\delta_c\) | \(0.1818\) |
| \(\lambda_{\mathrm{fold}}\) | \(0.8475\) |
| \(\lambda^\star\) | \(0.7628\) |
| \(\Delta^\star\) | \(0.2324\) meV |
| \(\lambda_{\min}\) | \(-g_0\) (multiplicity \(N-1\)) |
| \(Q_H\) | \(1/N\) |
| \(T^\star\) | \(118.7\) mK |
| Live gap (1 ps GPE) | \(0.6305\) meV |

Cavity moiré and Schwinger numbers are in `locks/PARAMETERS.md`. They are
parallel channels, not additive corrections to the static fold.

## Reproduce

### Paper

```bash
cd paper
latexmk -pdf QVC_arXiv_v4.tex
```

Requires a TeX distribution with `amsmath`, `tikz`, `pgfplots` (compat 1.18),
and `hyperref`. Committed PDFs are the September 2026 compile.

### Numerics

```bash
cd compute
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
python scripts/run_verification.py
python -m pytest tests -q
```

`qvc_l1` and `qvc_bkt` import the sibling `qvc` package. The bootstrap looks
in this `compute/` directory first.

### Visualizations

```bash
cd sim/standalone
python3 -m http.server 8765
```

Open `tdvt_hybrid_showcase.html`. These pages are illustrations, not the
verification suite.

## Repository map

```
locks/           Canonical JSON + PARAMETERS.md (paper numbers)
paper/           QVC_arXiv_v4.tex, macros.tex, qvc_refs.bib, figures/, tikz/
paper/           QVC_TQGL_supplement_v4.tex, QVC_Supplemental_Physics_v4.tex
compute/qvc/     Gap, fold, Casimir, Schwinger, TQGL, workbench
compute/qvc_l1/  L1 dictionary, Mexican hat, Hall, Lindblad, closure
compute/qvc_bkt/ Vortex / BKT C1–C6 suite
compute/tests/   compute_suite/ and l1_suite/
verification/    verification_report.md, tqgl_v4/, l1_* summaries
structures/      CIF / cube / generation scripts
sim/             GPU viewers (not the lock)
docs/            Closed-form identities and this manifest
```

## Provenance

Compiled 8 September 2026 from:

- `QVCCursor/tex/QVC_arXiv_v4/`
- `QVCCursor/qvc_l1/`, `qvc_bkt/`, `tests/`, `TEVC_structures/`, `sim-suites/`
- `QVCCompute/qvc/`, `tests/`, `scripts/`, `notebooks/`, `docs/`, `output/`

Working trees remain at `/Users/mineralhoiland/Code/QVCCursor` and
`/Users/mineralhoiland/Code/QVCCompute`. This snapshot is the publication
bundle.
