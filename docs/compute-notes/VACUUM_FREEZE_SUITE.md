# Vacuum freeze suite + P0/P1 brainstorm

## What was built

| Module / script | Role |
|-----------------|------|
| `qvc/vacuum_suite.py` | Fat-torus diagnostics, PFA ratio, F convergence, inverse freeze |
| `qvc/fat_torus_corrections.py` | Fat-torus F asymptotics + aspect-ratio scan |
| `qvc/self_consistency.py` | Rewritten Δ↔A_ZPF↔fold loop; bare/√n bridge scan |
| `qvc/gap_fold.py` | Branch solver + fold critical table (ℓ≈−2.791) |
| `qvc/hopf_cs_check.py` | Berry c₁ + Spec(G) + Q_H easy checks |
| `qvc/gij_overlap_skeleton.py` | P1 G_ij overlap skeleton vs democratic Spec |
| `scripts/run_vacuum_freeze_suite.py` | One-shot suite → `output/vacuum_freeze_suite_report.*` |
| `docs/P0_P1_DERIVATIONS.md` | Post-suite P0/P1 derivation workstream |

```bash
cd /Users/mineralhoiland/Code/QVCCompute
.venv/bin/python scripts/run_vacuum_freeze_suite.py
.venv/bin/python scripts/run_verification.py
```

## Locked vs provisional (do not conflate)

| Item | Status |
|------|--------|
| \(F=\sum K_0(nx)\), Tier-C \(E_{\mathrm{vac}}\) form | Locked |
| Fold map + \(\delta_c,\alpha_c\) at ℓ≈−2.791 | Locked |
| Spec\((G)\), \(Q_H=1/N\) | Locked |
| meV \(E_{\mathrm{vac}}\), \(\Delta^*\), \(T^*\) | **Provisional** |
| Inverse freeze (ℏc_s or r to hit 0.244) | Diagnostic only |

## Fat-torus findings (suite)

On preprint \((r,R)=(8,50)\,\mathrm{nm}\), \(x=2\pi r/R\approx 1\) (fat torus):

- Old Casimir formula **does not** reproduce 0.244 meV (audit PASS).
- Tier-C \(E_{\mathrm{vac}}\approx 0.129\,\mathrm{meV}\) at package \(\hbar c_s=0.07\) (provisional).
- Inverse Tier-C freeze on fat torus: either raise \(\hbar c_s\) or shrink \(r\) to hit 0.244 (suite prints both).
- **Fixed-geometry bridge** \(\alpha=\sqrt{n}\,E_{\mathrm{vac}}/\Delta_0\) is the primary self-consistency feed for freeze candidates.
- **Running-healing** \(A_{\mathrm{ZPF}}(\Delta)\) is exponentially small when \(a_H\gg R\) — diagnostic only, not the freeze loop.
- **Fold-in bound:** \(E_{\mathrm{vac}}\le\alpha_c\Delta_0\approx 0.109\,\mathrm{meV}\) (bare) or \(\approx 0.029\,\mathrm{meV}\) with \(\sqrt{n}\) at \(n=14\). Package fat Tier-C **exceeds both** → folded out until bridge/\(\hbar c_s\)/geometry change.
- Corrections moiré \(E_{\mathrm{vac}}\) is far past \(\alpha_c\) under either bridge.

## Self-consistency rewrite

Old prose mixed: static 35% vacuum loop, Spec(G) 52% catalyzon, TQGL 92%.

New table:

1. Geometry candidate + Tier-C \(E_{\mathrm{vac}}\)
2. Fixed-geometry fold loop \(\Delta^*\) with ℓ Corrections (bare **and** √n bridges in suite)
3. Running-healing diagnostic row
4. Separate catalyzon \(\Delta_{\mathrm{cat}}=\Delta_{\mathrm{ex}}-g_0\)

## Freeze decision brainstorm (solutions)

### Option A — Freeze fat torus + Tier-C, retune \(\hbar c_s\)
- Keep \((r,R)=(8,50)\,\mathrm{nm}\) narrative.
- Set \(\hbar c_s\) from inverse freeze so \(E_{\mathrm{vac}}=0.244\) **or** abandon 0.244 and quote Tier-C meV honestly.
- **Pros:** matches paper geometry language. **Cons:** \(\hbar c_s\) may disagree with phonon calibration.

### Option B — Freeze Corrections moiré + Tier-C, drop 0.244
- Use \((r,R)=(0.335,13.4)\,\mathrm{nm}\).
- Recompute self-consistency; update all **(prov.)** numbers.
- **Pros:** consistent with Corrections fold programme. **Cons:** breaks preprint 0.244 folklore.

### Option C — Two-track (recommended interim)
- **Track 1 (locked form):** Tier-C + F always.
- **Track 2 (geometry):** report both candidates in tables until platform \(\hbar c_s\) is measured.
- Promote to frozen only when (formula, geometry, \(\hbar c_s\), \(n_{\mathrm{modes}}\)) are fixed and verification reproduces to &lt;1%.

### Option D — Abandon absolute meV until SI calibration
- Keep dimensionless fold / Spec / \(Q_H\) as paper backbone.
- Quote only ratios \(\Delta/\Delta_{\mathrm{ex}}\), \(\alpha/\alpha_c\), \(F(x)\).

**Recommendation:** Option C now; choose A or B after one \(\hbar c_s\) lab/phonon constraint. Do not promote frozen this session.

## P0 next (after this suite)

1. Pick Option A/B/C explicitly in `VACUUM_FREEZE_SPEC.md` decision log.  
2. Commit fold critical table into Tier A TeX appendix from `fold_critical_table()`.  
3. Spec(G) prose scrub for any remaining −4g₀/−5g₀.  
4. Wire suite numbers into `CLOSED_FORM_LOCKS.md` “provisional numerics” subsection.

## P1 next (derivations)

1. **\(G_{ij}\) overlaps** — discretize TQGL+ZPF on WS cell; compare Spec to democratic lock.  
2. **H-TEGR1** — torsion-overlap kernel → same Spec target.  
3. **Vortex RG toy** — 0+1D fugacity flow *or* formal deferral.  
4. **Lindblad TQGL** — fate of 92% with \(\tau_{\mathrm{dec}}\sim5\)–50 ps.  
5. **Dynamic \(V_{\mathrm{ret}}(q,\omega)\)** — promote \(\eta_{\mathrm{ret}}\approx 3.1\).

## Acceptance tests before freeze promotion

- [ ] `assert_preprint_0p244_not_reproduced` still PASS for old formula  
- [ ] Jacobi θ₃ ≠ F PASS  
- [ ] Fold \(\delta_c\) matches appendix to &lt;1%  
- [ ] Spec(G) numeric = closed form  
- [ ] Chosen geometry + ℏc_s + n_modes written in lock JSON with `status: frozen`  
- [ ] Self-consistency table regenerated and TeX updated in the same commit  
