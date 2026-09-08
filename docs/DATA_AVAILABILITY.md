# Data and code availability

Private GitHub repository (this snapshot):

https://github.com/mineralhoiland/QVC

Until the repository is made public and a Zenodo DOI is minted, cite the
GitHub URL and `VERSION`. After a public release, add the version DOI to
`CITATION.cff` and to Appendix `app:data` in `paper/QVC_arXiv_v4.tex`.

## What to cite for each claim class

| Claim | Reproduce from |
| --- | --- |
| Vacuum \(E_{\mathrm{vac}}\), Bessel \(F\) | `compute/qvc/vacuum_freeze.py`, `locks/frozen_optionA_converged.json` |
| Fold \(\lambda^\star,\Delta^\star\) | `compute/qvc/gap_fold.py`, `gap_reduction.py`, `tqgl/locks.py` |
| Spec(\(G\)), \(Q_H=1/N\) | `compute/qvc/coupling_matrix.py`, `lens_chern_simons.py` |
| Cavity moiré / Schwinger | `compute/qvc/cavity_moire.py`, `schwinger_rate.py` |
| Live GPE / DCE | `compute/qvc/tqgl/`, `verification/tqgl_v4/` |
| L1 dictionary / Mexican hat / Hall | `compute/qvc_l1/` |
| Vortex C1–C6 | `compute/qvc_bkt/`, `paper/vortex_sector.tex` |
| Full PASS/FAIL suite | `compute/scripts/run_verification.py` → `verification/verification_report.md` |
| Figures | `paper/figures/*.py`, `paper/tikz/` |
| Bibliography | `paper/qvc_refs.bib` |
| Crystal files | `structures/cif/`, `verification/workbench_exports/` |

Visual HTML under `sim/` must not be cited as numerical verification.
