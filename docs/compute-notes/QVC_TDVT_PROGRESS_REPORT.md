# QVC / TDVT progress report (this chat)

**Date:** 15 August 2026  
**Purpose:** Single index of what changed, what was computed, where the files are, and what to do next.  
**Repos:** `QVCCompute` (solvers, TeX, locks) · `QVCCursor` (this workspace, sim/docs, suite output)

Do not flatten status tags. **Proven algebra ≠ computed meV ≠ preprint folklore.**

---

## 1. One-page status

| Layer | Status now |
|-------|------------|
| Democratic Spec\((G)\): \(\lambda_{\max}=(N-1)g_0\), \(\lambda_{\min}=-g_0\) (mult. \(N-1\)) | **Proven** (algebra + numeric) |
| Fold map \(\delta=1+\alpha(\ln\delta+\ell)\), \(\alpha_c=\delta_c\) | **Proven** |
| \(Q_H=1/N\) from CS/linking; \(H_2(L(N,1))=0\) | **Proven** |
| Vacuum **form** \(E_{\mathrm{vac}}=(1/4\pi^2)(\hbar c_s/r)F\), \(F=\sum K_0(nx)\) | **Proven form** |
| Bridge \(\alpha=\lambda E_{\mathrm{vac}}/\Delta_0\), \(\ell(r)=\ln(\Delta_0 r/(2\hbar c_s))+\gamma_E\) | **Proven algebra** |
| Geometry + \(\hbar c_s\) | **Frozen working** (fat torus \(8,50\) nm, \(0.07\,\mathrm{meV\cdot\mu m}\)) |
| \(E_{\mathrm{vac}}\) meV | **Computed** \(0.1287\) — **replaces \(0.244\)** |
| Static \(\Delta^*\) at \(\lambda=1\) | **Fold-out** (theorem) |
| Participation \(\lambda\) | **Derived bound** \(\lambda_\alpha^{\max}=0.8475\); recommended \(\lambda^\star=0.90\lambda_{\mathrm{fold}}=0.7628\) |
| 35% / 52% / 92% as one stack | **Withdrawn** (three mechanisms, mixed baselines) |
| TQGL \(G_{ij}\) from WS \(K_0\) overlaps | **Program-computed** (not torsion variation) |
| TDVT Einstein / \(\Lambda\) / “spacetime = lattice” | **Ideation only** (L3) |
| HTML particle / atlas suites | **Stopped** (not physics source) |

**Freeze option chosen:** **A\*** — frozen-cavity fat torus + Tier-C form + derived bridge. Not Option B (moiré). Option C was the interim two-track; A\* is the working standard. meV is computed, not a lab measurement.

---

## 2. What changed in this chat (physics)

### 2.1 Vacuum energy

Old Casimir formula \((1/4\pi)(\hbar c_s/r)(R/r)F\) **does not** give \(0.244\,\mathrm{meV}\) on \((r,R)=(8,50)\,\mathrm{nm}\). Audit PASS (claim fails).

| Quantity | Preprint claim | This chat |
|----------|----------------|-----------|
| \(F(x)\) at \(x\approx 1.005\) | \(0.596\) | **\(0.58068\)** |
| \(E\) | \(0.244\,\mathrm{meV}\) | **\(E_{\mathrm{vac}}=0.1287\,\mathrm{meV}\)** (Tier-C) |
| Old Casimir on same geometry | (implied \(0.244\)) | \(2.527\,\mathrm{meV}\) |
| \(A\xi\) | \(0.296\) | \(AR=0.0919\) |
| Jacobi \(\theta_3=F\) | used | **Withdrawn** (\(\theta_3\neq F\)) |

Inverse diagnostics (not a freeze): to hit \(0.244\) under Tier-C on fat torus would need \(\hbar c_s\approx 0.133\) or \(r^*\approx 6.35\,\mathrm{nm}\).

### 2.2 Derived \(\alpha\)–\(E_{\mathrm{vac}}\) bridge

Frozen cavity \(a_H=r=\mathrm{const}\) (not running \(a_H=\hbar c_s/\Delta\)):

\[
\ell(r)=\ln\frac{\Delta_0 r}{2\hbar c_s}+\gamma_E,\qquad
A=\frac{E_{\mathrm{vac}}}{\hbar c_s},\qquad
\alpha=\lambda\frac{E_{\mathrm{vac}}}{\Delta_0}.
\]

- \(\ell(8\,\mathrm{nm})=-2.7958\) vs Corrections \(\ell=-2.791\) (**\(0.17\%\)**). Corrections \(\ell\) **is** this logarithm.
- Running healing makes \(\ell_{\mathrm{tree}}\approx-0.116\) and \(F\) exponentially small; it is **not** the Corrections theory.
- Do **not** put \(\sqrt{n}\) into the static loop by default (\(F\) already sums windings). The early suite’s \(\sqrt{n}\) bridge (\(\alpha/\alpha_c\approx 4.41\)) is a diagnostic, not the freeze.

### 2.3 Fold / self-consistency

| Quantity | Value |
|----------|--------|
| \(\delta_c=\alpha_c\) | \(0.1818\) (geometry \(\ell\)) / \(0.1820\) (Corrections \(\ell\)) |
| \(\Delta_c\) at \(\Delta_0=0.6\) | \(0.109\,\mathrm{meV}\) |
| \(\alpha(\lambda=1)\) | \(0.2145\) → \(\alpha/\alpha_c=1.180\) → **no \(\Delta^*\)** |
| \(\lambda_\alpha^{\max}=\lambda_{\mathrm{fold}}\) | **\(0.8475\)** |
| Max **static** reduction | **\(81.8\%\)** at \(\Delta_c\) (not \(92\%\)) |
| \(\lambda\) that recovers withdrawn \(0.389\,\mathrm{meV}\) | **\(0.508\)** |
| \(E_{\mathrm{vac}}\) that would give withdrawn \(35\%\) | \(0.0653\,\mathrm{meV}\) |
| \(\hbar c_s\) at the fold (\(\lambda=1\)) | \(0.0593\,\mathrm{meV\cdot\mu m}\) |

**Theorem:** \(\delta=0.032/0.6\approx 0.053<\delta_c\). No physical fold branch can produce the preprint \(92\%\) as a static gap.

Running-healing diagnostic: fat torus \(\Delta^*\approx 0.6\), \(\alpha/\alpha_c\sim 10^{-7}\).

### 2.4 Participation \(\lambda\) and the three channels

**Do not stack 35+52+92.** Parallel channels:

1. **Fold (static):** \(\alpha=\lambda_\alpha E_{\mathrm{vac}}/\Delta_0\). Branch iff \(\lambda_\alpha\le 0.8475\).
2. **Catalyzon:** \(\Delta_{\mathrm{cat}}=\Delta_{\mathrm{ex}}-g_0\). If the same vacuum feeds \(G_{ij}\), \(g_0=\lambda_G E_{\mathrm{vac}}\). Working \(g_0=0.21\,\mathrm{meV}>\ E_{\mathrm{vac}}\) — **withdrawn as a vacuum scale**. Scaled preprint \(|\lambda_{\min}|\sim 0.1\) → \(g_0\approx 0.0527\,\mathrm{meV}\), \(\lambda_{\mathrm{cat}}\approx 0.410\).
3. **TQGL 1 ps:** \(V_{\mathrm{ZPF}}\propto E_{\mathrm{Cas}}\). **No split-step solver in-repo.** Drive scale \(E_{\mathrm{vac}}/0.244\approx 0.527\). Proxy remaining \(\Delta\sim 0.128\,\mathrm{meV}\) (\(\sim 79\%\) of \(\Delta_0\)) is **not** a new GPE. Lindblad still omitted.

**Recommended \(\lambda^\star=0.90\lambda_{\mathrm{fold}}=0.7628\):**

| Object | Value |
|--------|--------|
| Fold \(\Delta^\star\) | \(0.232\,\mathrm{meV}\) (**\(61.3\%\)**) |
| Parallel catalyzon \(\Delta_{\mathrm{cat}}\) | \(0.502\,\mathrm{meV}\) (**\(16.4\%\)**) |
| Max static (at fold) | \(81.8\%\) |

Partition (unitarity \(\lambda_\alpha+\lambda_G\le 1\)) examples:

| Scenario | \(\lambda_\alpha\) | \(\lambda_G\) | \(\Delta^*\) | \(\Delta_{\mathrm{cat}}\) | 1 ps proxy |
|----------|-------------------|---------------|--------------|---------------------------|------------|
| \(\lambda_\alpha=1\) | \(1\) | \(0\) | fold-out | \(0.600\) | \(\sim 62\%\) |
| Fold-sat. leftover \(G\) | \(0.847\) | \(0.153\) | \(0.120\) | \(0.100\) | \(\sim 94\%\) |
| Equal split | \(0.50\) | \(0.50\) | \(0.393\) | \(0.328\) | \(\sim 79\%\) |
| Catalyzon only | \(0\) | \(1\) | \(0.600\) | \(0.471\) | \(\sim 70\%\) |

Equal split is a clean subcritical working point. Fold-saturated 1 ps \(>92\%\) can happen because the catalyzon floor is lower than the old \(0.20\,\mathrm{meV}\) stage-2 floor, even with a weaker drive. That is a **proxy**, not a lock.

### 2.5 \(G_{ij}\) (TQGL + Wigner–Seitz)

Five Gaussians on the moiré hexagon, kernel \(K_0(|r-r'|/R)\), \(G_{ii}=0\). Democratic \(g_0(J-I)\) **not inserted**.

| Test | Result |
|------|--------|
| Hermitian, \(\mathrm{Tr}\,G=0\), \(\lambda_{\min}<0\) | PASS |
| Rel. Frobenius to democratic at \(\langle G_{\mathrm{off}}\rangle\) | **\(7.0\%\)** (near-democratic) |
| Dark multiplicity | Split (pentagon in hexagon), not exact \(N-1\) |

**Not done:** variation of \(S_{\mathrm{FN}}+S_{\mathrm{sf}}+S_{\mathrm{int}}\) (H-TOR1 / H-TEGR1).

### 2.6 Easy math verification (suite)

| Check | Result |
|-------|--------|
| Spec\((G)\) numeric vs closed form | PASS: \(\lambda_{\max}=0.84\), \(\lambda_{\min}=-0.21\times 4\) (\(g_0=0.21\) example) |
| Berry \(c_1\) | \(0.999793\) → PASS |
| \(Q_H\) | \(0.2\) (\(N=5\)) |
| Old Casimir \(\to 0.244\) | FAIL_CLAIM (audit PASS) |
| Jacobi \(\theta_3\neq F\) | PASS withdrawal |
| `tests/test_locks.py` | **10 passed** (after gap-reduction tests) |

---

## 3. TDVT / hopfion spacetime (honest layering)

Unchanged in substance this chat; still **not** derived from QVC numerics.

| Layer | Content | Claim strength |
|-------|---------|----------------|
| **L1 QVC** | Gap, fold, catalyzon, Casimir/Tier-C, Hall \(Q_H\) | Partially locked as above |
| **L2 topological CM** | Hopf, CS/linking, acoustic Unruh–Visser **analogue**, 3+2 Goldstones | Identities + analogues |
| **L3 spacetime** | Einstein from hopfion phonons, cosmological \(\Lambda\), lattice = spacetime | **Ideation only** |

**Shared (keep):** \(Q_H=1/N\), linking/CS, analogue acoustic metric.  
**Not shared until derived:** Einstein–Cartan limit, \(\Lambda\), ontology.  
**Hypotheses register:** H-VRG1, H-TOR1, H-SHARE1, H-HAWK1, H-VTOR1 in `HYPOTHESES_PHASE1.md`.

Fusion / FRC notes remain speculative; no meV→MeV Coulomb-barrier extrapolation.

---

## 4. Computations and how to re-run

All from `QVCCompute` unless noted. Outputs also under `QVCCursor/output/qvc_freeze/`.

```bash
cd /Users/mineralhoiland/Code/QVCCompute
MPLBACKEND=Agg MPLCONFIGDIR=/tmp/mpl \
  QVC_OUTPUT_DIR=/Users/mineralhoiland/Code/QVCCursor/output/qvc_freeze \
  QVC_LOCKS_DIR=/Users/mineralhoiland/Code/QVCCursor/output/qvc_freeze/locks \
  .venv/bin/python scripts/run_vacuum_freeze_suite.py
MPLBACKEND=Agg .venv/bin/python scripts/run_bridge_and_gij.py
MPLBACKEND=Agg .venv/bin/python scripts/run_gap_reduction.py
.venv/bin/python -m pytest tests/test_locks.py -q
.venv/bin/python scripts/run_verification.py
```

Jupyter: `notebooks/vacuum_bridge_lab.ipynb`.

| Computation | Script / module | Principal results |
|-------------|-----------------|-------------------|
| Vacuum freeze suite | `scripts/run_vacuum_freeze_suite.py` | Fold table; two geometries; inverse freeze; Spec/\(c_1\)/\(Q_H\) |
| Fat-torus \(F\) | `qvc/fat_torus_corrections.py` | \(x\approx 1.005\), \(F=0.581\), n=1 rel. err. \(\sim 28\%\) |
| Bridge A\* | `qvc/bridge.py` | \(\ell(r)\), \(\alpha=E/\Delta_0\), fold-out at \(\lambda=1\) |
| \(G_{ij}\) TQGL–WS | `qvc/gij_tqgl_ws.py` | Hermitian, \(\lambda_{\min}<0\), rel. Frob. \(7\%\) |
| Gap reduction / \(\lambda\) | `qvc/gap_reduction.py` | \(\lambda_{\mathrm{fold}}=0.8475\), \(\lambda^\star=0.7628\), 92% not static |
| Fold solver | `qvc/gap_fold.py` | \(\delta_c,\alpha_c\), branches |
| Spec\((G)\) | `qvc/coupling_matrix.py` | Democratic lock |
| Hopf/CS numeric | `qvc/hopf_cs_check.py` | Berry \(c_1\approx 1\), \(Q_H=1/N\) |

---

## 5. File index

### New / updated in this chat (`QVCCompute`)

| Path | Role |
|------|------|
| `qvc/bridge.py` | Proven \(\alpha\)–\(E_{\mathrm{vac}}\)–\(\ell(r)\) map |
| `qvc/gap_reduction.py` | \(\lambda\) fold/catalyzon/TQGL proxy |
| `qvc/gij_tqgl_ws.py` | WS overlap \(G_{ij}\) |
| `qvc/gij_overlap_skeleton.py` | Earlier Spec acceptance skeleton |
| `qvc/self_consistency.py` | Bare vs \(\sqrt{n}\) bridges, fold-in bounds |
| `qvc/vacuum_suite.py` | Fat-torus + inverse freeze |
| `qvc/fat_torus_corrections.py` | \(F\) asymptotics, aspect scan |
| `scripts/run_vacuum_freeze_suite.py` | Suite orchestrator |
| `scripts/run_bridge_and_gij.py` | Bridge freeze + \(G_{ij}\) figures |
| `scripts/run_gap_reduction.py` | \(\lambda\) / 35/52/92 recompute |
| `tests/test_locks.py` | Lock + bridge + \(G_{ij}\) + \(\lambda\) tests |
| `notebooks/vacuum_bridge_lab.ipynb` | Interactive lab |
| `docs/BRIDGE_DERIVATION.md` | Bridge derivation |
| `docs/GAP_REDUCTION_LAMBDA.md` | \(\lambda\) and channel split |
| `docs/VACUUM_FREEZE_SPEC.md` | Decision log (A\* added) |
| `docs/VACUUM_FREEZE_SUITE.md` | Suite notes |
| `docs/CLOSED_FORM_LOCKS.md` | Locks + bridge subsection |
| `docs/P0_P1_DERIVATIONS.md` | Workstream (partly superseded by A\*) |
| `tex/appendix_fold_critical.tex` | Fold appendix stub |
| `docs/QVC_TDVT_PROGRESS_REPORT.md` | **This file** |

### Pre-existing theory docs (still current)

| Path | Role |
|------|------|
| `docs/NEXT_STEPS_RIGOR.md` | P0–P3 roadmap (update freeze item: A\* chosen) |
| `docs/TORSION_TDVT_HOPFION.md` | L1/L2/L3, torsion kernel |
| `docs/VORTEX_RG_FRAMEWORK.md` | Fold → BKT checklist |
| `docs/HYPOTHESES_PHASE1.md` | H-VRG1, H-TOR1, … |
| `tex/QVC_VortexRG_Torsion_Ideation.tex` | Ideation TeX |
| `tex/QVC_arXiv_preprint_v3_tierA.tex` | Preprint (still contains \(0.244\), \(0.389\), \(92\%\)) |
| `qvc/locks/` or output locks | `frozen_working_fat_torus.json`, `provisional_corrections_self_energy.json` |

### Outputs (`QVCCursor/output/qvc_freeze/`)

| Path | Role |
|------|------|
| `vacuum_freeze_suite_report.md/.json` | Full freeze suite |
| `bridge_gij_report.md/.json` | Bridge + \(G_{ij}\) |
| `gap_reduction_report.md/.json` | \(\lambda\) and channel table |
| `locks/frozen_working_fat_torus.json` | Working freeze |
| `locks/provisional_corrections_self_energy.json` | Old provisional (moiré) |
| `figures/` | Fold, \(F(x)\), \(G_{ij}\) heatmap/spectrum (if generated) |

### Canvases (open beside chat)

- [QVC bridge freeze](/Users/mineralhoiland/.cursor/projects/Users-mineralhoiland-Code-QVCCursor/canvases/qvc-bridge-freeze.canvas.tsx)
- [Gap reduction λ](/Users/mineralhoiland/.cursor/projects/Users-mineralhoiland-Code-QVCCursor/canvases/qvc-gap-lambda.canvas.tsx)

### Stopped this chat (do not treat as physics)

- `QVCCursor/sim-suites/` particle lab / `tdvt_atlas.html` — user rejected; not a source of numbers.

Copy of this report also: `QVCCursor/docs/QVC_TDVT_PROGRESS_REPORT.md`.

---

## 6. What is *not* complete

- Preprint TeX still quotes \(E_{\mathrm{Cas}}=0.244\), \(\Delta_{\mathrm{QVC}}=0.389\), \(A\xi=0.296\), \(F=0.596\), stacked 35/52/92%.
- No Lindblad TQGL; no in-repo 1 ps split-step.
- \(G_{ij}\) not from torsion action (H-TOR1 open).
- Vortex RG not done (BKT still contingent).
- Dynamic \(V_{\mathrm{ret}}(q,\omega)\) not done (\(\eta_{\mathrm{ret}}\approx 3.1\) still a static caveat).
- Absolute meV not a measured MATBG energy.
- TDVT L3 not derived.

---

## 7. Next big steps (ordered)

Same ladder as before, updated for A\* and \(\lambda^\star\).

### Do next (highest leverage)

1. **Wire A\* + \(\lambda^\star\) into the preprint**  
   Replace \(0.244\), \(F=0.596\), \(A\xi=0.296\), \(\Delta_{\mathrm{QVC}}=0.389\) with computed values and **provisional/computed** tags. State fold-out at \(\lambda=1\) and \(\lambda^\star=0.7628\) as the working participation. Commit `appendix_fold_critical.tex`. This is manuscript P0 and unblocks honest quoting.

2. **Scale \(G_{ij}\) to meV**  
   The WS kernel is in arbitrary units (\(g_{0,\mathrm{eff}}\sim 10^{-4}\)). Fix the prefactor so \(g_0=\lambda_G E_{\mathrm{vac}}\) and re-check Spec vs democratic. That is the next **theoretical** step after the bridge.

3. **Lindblad / actual TQGL 1 ps**  
   The 92% proxy is the weakest number. Either implement a minimal 1D/2D dissipative run (\(\tau_{\mathrm{dec}}\sim 5\)–\(50\,\mathrm{ps}\)) or **formally defer** Stage-3 and never quote 92% as derived.

### Then (P1)

4. Vortex RG toy (0+1D fugacity) **or** explicit BKT deferral paragraph in the preprint.  
5. Dynamic \(V_{\mathrm{ret}}(q,\omega)\) for \(\eta_{\mathrm{ret}}\).  
6. Torsion variation H-TOR1 (same Spec target, no hand-inserted \(J-I\)).

### Later (P2/P3 TDVT + hygiene)

7. Keep L3 spacetime claims out of the MATBG paper.  
8. Scrub leftover \(-4g_0/-5g_0\) if any remain in figures.  
9. Canonical equations appendix shared by preprint + supplemental + `CLOSED_FORM_LOCKS.md`.

---

## 8. Suggested next session (concrete)

If only one session: **(1) TeX/table patch** — fold appendix + withdrawn 0.244/0.389/92% language using this report’s numbers.  
If two: **(1)** plus **(2)** dimensional \(G_{ij}\) so catalyzon \(g_0\) matches \(\lambda_G E_{\mathrm{vac}}\).

Do not promote a lock to `status: frozen` as a laboratory energy. A\* is **frozen working**: form + bridge + geometry. Participation \(\lambda^\star\) is a derived operating choice, not a measured overlap.
