# Next steps — mathematical rigor, consistency, closed gaps (QVC ↔ TDVT)

## Status after figure regeneration

| Layer | Status |
|-------|--------|
| L1 QVC math language (preprint + supplemental) | Minimal Tier A corrections; bibliography restored |
| L1 figures (spectrum, fold, Q_H, cavity) | Regenerated from locked `qvc/` solvers + matplotlib |
| Vacuum energy meV | **Still provisional** — geometry / ħc_s / formula not frozen |
| L2 torsion → G_ij from action | Kernel stated; **not** variationally computed |
| L3 spacetime / TDVT Einstein limit | Ideation only (`docs/`) |

---

## P0 — close numerical & identity gaps

1. **Vacuum freeze (blocking)**  
   Choose one triple `(formula, (r,R), ħc_s)` in `qvc/vacuum_freeze.py`, promote lock `provisional` → `frozen`, re-run self-consistency table, update preprint \(E_{\mathrm{Cas}}\) / \(\Delta_{\mathrm{QVC}}\) only after that.

2. **Single Spec(G) everywhere**  
   Done in text + regenerated figures for democratic \(G\). Scrub any remaining prose that says \(-4g_0\) or \(-5g_0\) as the democratic eigenvalue.

3. **Reproduce fold numbers in preprint**  
   Commit \(\delta_c,\alpha_c,A_c,\Delta_c\) from `qvc/gap_fold.py` into a small appendix table; ensure PDF and text match to <1%.

4. **Explicit Hopf / CS check**  
   Numerical integral of \(F\wedge A\) on an \(N\)-layer ansatz → \(Q_H\approx 1/N\).

---

## P1 — microscopic closure (QVC)

5. **\(G_{ij}\) from overlaps**  
   Linearize TQGL + ZPF correlators on a Wigner–Seitz cell; output a real \(5\times 5\) Hermitian matrix; compare Spec to democratic lock (accept only if \(\lambda_{\min}<0\) emerges).

6. **Vortex RG (optional BKT upgrade)**  
   Follow `docs/VORTEX_RG_FRAMEWORK.md` checklist C1–C7. Until done, keep fold language; never claim essential singularity as derived.

7. **Lindblad TQGL**  
   Re-run 1 ps dynamics with \(\tau_{\mathrm{dec}}\sim 5\)–50 ps; report whether 92% collapse survives.

8. **Retardation / dynamic gap equation**  
   Promote \(\eta_{\mathrm{ret}}=3.1\) from static caveat to frequency-dependent \(V_{\mathrm{ret}}(q,\omega)\).

---

## P2 — TDVT / hopfion spacetime sharing (honest layering)

9. **Shared (keep)**  
   Hopf charge, CS/linking \(Q_H=1/N\), acoustic Unruh–Visser metric as **analogue**, five Goldstones (3+2) as CM counting.

10. **Not shared until derived**  
    Einstein equations from hopfion phonons; cosmological \(\Lambda\); “spacetime = lattice” ontology.

11. **Torsion program**  
    Vary \(S_{\mathrm{FN}}+S_{\mathrm{sf}}+S_{\mathrm{int}}\); extract \(T^\mu{}_{\nu\rho}\) and \(G_{ij}\) without inserting \(g_0(J-I)\) by hand (`docs/TORSION_TDVT_HOPFION.md`).

12. **Vortex–torsion coupling (conjecture)**  
    Ideation hypothesis H-VTOR1: vortex worldlines as torsion sources — write acceptance test before numerics.

---

## P3 — manuscript hygiene

13. Regenerate remaining structural PDFs if burned-in labels conflict (`hopfion_vortex_3d`, `laser_driven_qvc_setup`).

14. One **canonical equations** appendix (locked identities only) shared by preprint + supplemental + `docs/CLOSED_FORM_LOCKS.md`.

15. Separate TDVT spacetime note from MATBG Hall paper (do not cross-cite as proven).

---

## Suggested session order

| Session | Deliverable |
|---------|-------------|
| Next | Vacuum freeze decision + updated self-consistency numbers |
| +1 | \(G_{ij}\) quadrature notebook |
| +2 | Remaining structural figure scrub |
| +3 | Vortex RG toy flow (even 0+1D) or explicit deferral note |
| +4 | Canonical appendix + arXiv-ready bundle |

## Commands

```bash
cd /Users/mineralhoiland/Code/QVCCompute
MPLCONFIGDIR=/tmp/mpl .venv/bin/python scripts/run_verification.py
MPLCONFIGDIR=/tmp/mpl .venv/bin/python scripts/generate_corrected_figures.py
MPLCONFIGDIR=/tmp/mpl .venv/bin/python scripts/freeze_vacuum_demo.py
cd tex/tikz && pdflatex spectrum_democratic.tex && pdflatex gap_fold_pgfplots.tex
```
