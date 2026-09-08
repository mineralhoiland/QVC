# P0 / P1 derivation workstream (post vacuum-freeze suite)

Status: vacuum freeze suite is **complete and provisional**. This note starts
the next derivation layer. No meV vacuum number is frozen.

## Suite verdict (inputs to P0)

| Lock | Value / status |
|------|----------------|
| Fold ℓ | −2.791 (locked) |
| δ_c = α_c | ≈ 0.182 (locked) |
| Δ_c at Δ₀=0.6 | ≈ 0.109 meV (locked map) |
| Fat torus x | ≈ 1.005, F ≈ 0.581 |
| E_TierC (package ℏc_s) | ≈ 0.129 meV **provisional** |
| Old Casimir → 0.244 | **FAIL** (audit PASS) |
| Fold-in bare | E_vac ≤ α_c Δ₀ ≈ **0.109** meV |
| Fold-in √n (n=14) | E_vac ≤ ≈ **0.029** meV |
| Package fat Tier-C | **past fold** on both bare and √n bridges |
| Running healing | Δ*≈Δ₀ (α/α_c ∼ 10⁻⁷) — not the freeze loop |
| Spec(G), Q_H, Berry c₁ | verification PASS |

**Implication:** freezing preprint geometry + package ℏc_s + Tier-C + √n
enhancement is **inconsistent with a gapped fold branch**. Freeze Option C
(two-track) remains mandatory until ℏc_s / n / bridge are chosen.

## P0 — begin now

### P0.1 Vacuum freeze decision log
- Keep Option C interim in `VACUUM_FREEZE_SPEC.md`.
- Do **not** promote `status: frozen` this session.
- Candidate resolutions (brainstorm):
  1. Drop √n for static bridge; still need E_vac ≲ 0.109 meV (shrink ℏc_s or grow r).
  2. Abandon absolute meV (Option D): quote only α/α_c, F(x), δ.
  3. Reinterpret 0.244 as a **target** only via inverse freeze (ℏc_s≈0.133 or r*≈6.35 nm).

### P0.2 Fold critical TeX appendix
Source of truth: `qvc.gap_fold.fold_critical_table()`.

Committed numbers for appendix (regenerate if ℓ changes):

```
ℓ = -2.791
ℓ_tree = ln(1/2)+γ_E ≈ -0.115932
δ_c = α_c ≈ 0.181990
Δ₀ = 0.6 meV  →  Δ_c ≈ 0.109194 meV
scaling: δ − δ_c ∼ ±√(α_c − α)
BKT: contingent_on_vortex_RG
```

TeX stub: `tex/appendix_fold_critical.tex` (draft).

### P0.3 Spec(G) scrub
Democratic lock only: λ_max=(N−1)g₀, λ_min=−g₀ (mult. N−1).
Never −4g₀ / −5g₀ as democratic eigenvalues.

### P0.4 Provisional numerics subsection
Wire suite report into `CLOSED_FORM_LOCKS.md` under “provisional — not frozen”.

## P1 — begin skeleton

### P1.1 G_ij overlaps
Module: `qvc/gij_overlap_skeleton.py`

Acceptance (before claiming derivation):
- [ ] Hermitian 5×5 from overlaps (no hand-inserted g₀(J−I))
- [ ] λ_min < 0
- [ ] Distance to democratic Spec documented (Frobenius / eigenvalue errs)

### P1.2 H-TEGR1
Torsion-overlap kernel → same Spec target (`docs/TORSION_TDVT_HOPFION.md`).

### P1.3 Vortex RG
Either 0+1D fugacity toy **or** explicit deferral note in preprint.

### P1.4 Lindblad TQGL
Fate of 92% with τ_dec ∼ 5–50 ps.

### P1.5 Dynamic V_ret
Promote η_ret ≈ 3.1 from static caveat.

## Easy verification programs (keep green)

```bash
cd /Users/mineralhoiland/Code/QVCCompute
QVC_OUTPUT_DIR=.../QVCCursor/output/qvc_freeze \
  .venv/bin/python scripts/run_vacuum_freeze_suite.py
.venv/bin/python -m qvc.gij_overlap_skeleton
.venv/bin/python scripts/run_verification.py
```

## Recommended next coding session

1. User picks freeze Option A/B/C/D explicitly.
2. Fill TeX fold appendix + provisional tags on 0.244 in preprint.
3. Replace `toy_ws_overlap_G` with real TQGL mode overlaps.
