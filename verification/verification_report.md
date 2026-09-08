# QVC Tier B/C Verification Report

Summary: `{'PASS': 23, 'FAIL': 0, 'WARN': 1, 'INFO': 5}`

| Status | Check | Detail | Value | Expected |
|--------|-------|--------|-------|----------|
| PASS | democratic Spec(G) exact match | eigs=[-1.0000000000000004, -1.0000000000000002, -1.0, -0.9999999999999993, 3.9999999999999996] | -1.0000000000000004 | -1.0 |
| PASS | preprint claim lambda_min=-4g0 is FALSE | Democratic matrix does not have lambda_min=-4g0 | -1.0000000000000004 | -4.0 |
| PASS | Tr(G)=0 | Tr=0.0 | None | None |
| PASS | catalyzon dark basis is sum-zero | row sums=[0.0, 0.0, 1.1102230246251565e-16, 0.0] | None | None |
| PASS | catalyzon dark basis orthonormal (Helmert) | B B^T ≈ I | None | None |
| PASS | catalyzon gap shift Δ_cat = Δ_ex - g0 | Identity from λ_min | -0.4 | -0.4 |
| PASS | supplemental lambda_min=-5g0 != democratic | Octahedral adjacency claim needs a different matrix; democratic is -g0 | -1.0000000000000004 | -5.0 |
| PASS | CS[A_1] = 1/N | Chern-Simons of fundamental flat connection | 0.2 | 0.2 |
| PASS | linking self-number = 1/N | Rational linking form on H1=Z_N | 0.2 | 0.2 |
| PASS | Q_H identification = 1/N | Physical Hall fractionalization target | 0.2 | 0.2 |
| INFO | anyonic theta = pi/N | Programmatic prediction from same linking form | 0.6283185307179586 | 0.6283185307179586 |
| PASS | H2(L(N,1))=0 — no fundamental 2-cycle | Preprint 2-cycle fractionalization proof path is INVALID; use CS/linking on H1 instead. | 0 | 0 |
| PASS | fold δ_c recovers Corrections ~0.182 | ell=-2.790757 | 0.182 | 0.182 |
| PASS | fold α_c = δ_c | Saddle-node condition | 0.182 | 0.182 |
| PASS | fold has FINITE Delta_c (not BKT vanishing) | Without vortex RG, transition is saddle-node | 0.10919999999999999 | 0.109 |
| PASS | two real branches exist below alpha_c | max sign-change pairs=2 | 2 | 2 |
| PASS | tree-level log coefficient ln(1/2)+γ_E | Negative ⇒ vacuum amplitude suppresses gap | -0.11593151565841242 | -0.11593 |
| PASS | preprint formula does NOT yield 0.244 meV on stated (r,R) | FAIL_CLAIM: preprint 0.244 meV is NOT reproduced by old Casimir formula on stated preprint geometry | 2.5270362279386474 | 0.244 |
| INFO | Bessel F(x) at preprint aspect ratio | x=1.0053 | 0.5806751779874608 | 0.596 |
| INFO | corrected prefactor on preprint geometry (no R/r) | Shows sensitivity to formula choice | 0.12870089825559464 | None |
| INFO | corrected self-energy on Corrections geometry (package units) | F(0.157)=8.097; provisional — not a paper meV claim | 42.85874528602627 | None |
| WARN | Corrections SI self-energy order ~0.25 meV | F(0.157)=8.097; sensitive to K0 accuracy / c_s | 8.467411287782342 | 0.256 |
| PASS | Jacobi θ₃ alternating does NOT equal F≈0.596 | theta=-2.415966e-16 at alpha_c=0.0442 | -2.4159660232755273e-16 | 0.596 |
| PASS | provisional vacuum lockfile written | qvc/locks/provisional_corrections_self_energy.json | provisional_corrections_self_energy.json | provisional_corrections_self_energy.json |
| PASS | single-mode A_ZPF at bare gap is negligible | Matches Corrections iteration-0 exponential suppression | 6.411531954925066e-26 | None |
| PASS | √14 enhancement ratio | Quadrature multiplicity factor | 3.7416573867739418 | 3.7416573867739413 |
| PASS | η_ret = 2 when ξ = a_H = ħc_s/Δ₀ | Causal-propagator algebra identity | 2.0 | 2.0 |
| PASS | η_ret = 2 Δ₀ ξ / (ħ c_s) for ξ=λ_M | General retardation factor | 0.22285714285714284 | None |
| INFO | acoustic R_WS matching helper (not duality proof) | MATCHING_NOT_PROOF: invert T_H=T* for R_WS at T*=1 K | 0.25856829897854017 | None |

## Canonical lock notes

- **lambda_min**: Use -g0 (mult N-1), not -4g0 or -5g0
- **Q_H**: Use CS/linking; never H2-cycle language
- **gap_transition**: Fold unless vortex RG is derived
- **vacuum**: Freeze one (formula, geometry) before quoting meV; see locks/
- **retardation**: η_ret=2Δ0 ξ/(ħ c_s) — algebra lock
- **mode_enhancement**: A_eff=√n A_single — quadrature lock

## Provenance

- Closed-form locks: `docs/CLOSED_FORM_LOCKS.md`
- Vacuum freeze decisions: `docs/VACUUM_FREEZE_SPEC.md`
- Ported/improved from `QVCCursor/rigor/`
