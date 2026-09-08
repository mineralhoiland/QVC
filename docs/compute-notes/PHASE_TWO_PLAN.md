# QVC Phase Two Computation Plan

## Ground Truth (Option A* Lock — Authoritative)

| Quantity | Locked Value | Status |
|----------|--------------|--------|
| (r, R) | (8, 50) nm | Frozen |
| ℏc_s | 0.07 meV·μm | Provisional (platform calibration needed) |
| Δ₀ | 0.600 meV | Locked |
| F(x) = Σ K₀(nx), x=2πr/R | **0.58068** | Locked |
| E_vac = (1/4π²)(ℏc_s/r)F | **0.1287 meV** | Locked under Option A* |
| ℓ (fold log param) | −2.7958 | Locked |
| λ_fold | 0.8475 | Locked |
| λ* | 0.7628 (= 0.90 × λ_fold) | Locked |
| Δ* | 0.2324 meV | Locked |
| α_c = δ_c | 0.1818 | Locked |
| Δ_c (at fold) | 0.109 meV | Locked |
| Max static reduction | 81.8% (at fold) | Locked |
| g₀ = λ* E_vac | 0.0982 meV | Locked |
| λ_max = (N−1)g₀ | 0.3928 meV | Locked |
| λ_min = −g₀ | −0.0982 meV (mult. 4) | Locked |
| Q_H = 1/N | 0.2 | Locked (CS/linking) |
| θ = π/(2N) | π/10 | Locked |
| δσ_xy | e²/(10h) ≈ 3.874×10⁻⁶ S | Locked |
| A*/A_c^fold | 0.9117 (resonant) | Locked |
| T* | 118.7 mK (diagnostic) | Locked |

### DEPRECATED (do NOT use)
- E_Cas = 0.244 meV
- 92% gap reduction / Δ→0.032 meV
- T* = 285–290 mK
- A_ZPF·ξ = 0.296 as coupling constant
- λ_min = −4g₀ or −5g₀

---

## Workstream 1: Vacuum Lock Convergence Analysis

### Physics Goal
Demonstrate that F(x=1.0053) = 0.58068 is converged to <0.01% under variation of:
(a) truncation order N_max in Σ K₀(nx)
(b) alternative boundary conditions (Dirichlet vs periodic on torus)
(c) continuum limit vs discrete lattice
(d) sensitivity to ℏc_s ∈ [0.05, 0.10] meV·μm

### Mathematical Specification
```
F(x) = Σ_{n=1}^{N_max} K₀(nx)

Convergence: |F(N_max) − F(N_max−1)| / |F(N_max)| < ε

Sensitivity: ∂E_vac/∂(ℏc_s) |_{ℏc_s=0.07} = F/(4π² r) = 1.838 meV/(meV·μm)

Alternative spectral sum (Dirichlet on cylinder):
F_D(x) = Σ_{n=1}^∞ K₀(j_{0,n} x / π)   where j_{0,n} = Bessel zeros
```

### Acceptance Criteria
1. F converges to 6 significant figures by N_max = 500
2. Dirichlet vs periodic gives ΔF/F < 5% (boundary-independent at this aspect ratio)
3. Error budget: dominant uncertainty is ℏc_s, not truncation
4. Lockfile upgraded from "provisional" to "frozen_converged"

### Deliverables
- `notebooks/01_vacuum_convergence.ipynb` — full analysis + figures
- `qvc/vacuum_convergence.py` — production convergence checker
- Figures: convergence plot, sensitivity band, BC comparison
- Updated lockfile: `qvc/locks/frozen_optionA_converged.json`

---

## Workstream 2: G_ij Wigner-Seitz Cell Quadrature

### Physics Goal
Verify that the 5×5 coupling matrix emerging from TQGL mode overlaps on the
hexagonal WS cell approaches the democratic form G = g₀(J−I) when the kernel
range ξ_K ≫ inter-site spacing d, and characterize deviations at physical ξ_K.

### Mathematical Specification
```
G_ij = ∫∫_{WS} Φ_i(r) K₀(|r−r'|/ξ_K) Φ_j(r') d²r d²r'   (i≠j)
G_ii = 0

Φ_i: Gaussians of width σ centered on regular N-gon of radius c·R_hex
K₀: modified Bessel function of second kind

Controls:
  σ/λ_M ∈ [0.10, 0.30]  (mode width)
  c = r_centre/R_hex ∈ [0.30, 0.60]  (how close to edge)
  ξ_K/λ_M ∈ [1, 10]  (kernel range relative to moiré)
  n_grid ∈ [32, 128]  (convergence in spatial quadrature)
```

### Acceptance Criteria
1. Hermitian: ||G − G^T|| < 10⁻¹² (machine precision)
2. Trace-zero: |Tr G| < 10⁻¹² (no self-coupling)
3. λ_min < 0 (catalyzon exists) WITHOUT inserting democratic form
4. At ξ_K/λ_M = 4 (physical): ||G − g₀(J−I)||_F / ||g₀(J−I)||_F < 15%
5. Democratic limit: at ξ_K/λ_M → ∞, ratio → 0
6. Grid convergence: results stable n_grid = 64 vs 128 to 1%

### Deliverables
- `notebooks/02_gij_wigner_seitz.ipynb` — parameter sweep + convergence
- Publication figure: Spec(G) vs ξ_K/λ_M with democratic limit shown
- Heatmap of G_ij matrix at physical parameters
- `output/gij_sweep.csv` — tabulated eigenvalues

---

## Workstream 3: Vortex RG Numerical Flow

### Physics Goal
Solve the Kosterlitz-Thouless RG flow equations adapted to moiré excitonic
order. Determine whether the fold saddle-node transitions to BKT-class
behavior for realistic stiffness parameters.

### Mathematical Specification
```
Textbook KT flow (Kosterlitz normalisation):
  dK⁻¹/dℓ = 4π³ y²
  dy/dℓ = (2 − πK) y

Diagnostic stiffness on fold branch:
  K(α) = K₀ [δ_+(α)]²,   K₀ = (2/π)(T_BKT^bare / T_eff)

With moiré anyon correction (conjecture H-VRG2):
  K_c^anyon = K_c^boson × (1 − θ/π)²,   θ = π/N

Critical questions (C1–C7):
  C1. Does the flow reach K < K_c before α reaches α_c?
  C2. Does y diverge (vortex proliferation)?
  C3. What is the length scale ℓ* at which BKT kills the fold?
  C4. Core energy E_core ~ π κ_s ln(L/a_H) with a_H = 8 nm
  C5. Fugacity y_0 = exp(−E_core / k_B T_eff)
  C6. Does phonon damping suppress vortex nucleation?
  C7. Experimental signature: Δ(A) closing shape near A_c
```

### Acceptance Criteria
1. RG solver converges (adaptive RK45, rel tol 10⁻⁸)
2. Phase diagram: (K₀, y₀) → fold vs BKT outcome
3. Clear identification of crossover regime
4. Quantitative prediction: T_BKT(α) curve on the fold branch
5. Honest conclusion on whether BKT is activated at T* = 118.7 mK

### Deliverables
- `notebooks/03_vortex_rg_flow.ipynb` — flow diagrams + phase portrait
- Publication figure: RG flow lines in (K⁻¹, y) plane
- Phase diagram: K₀ vs T_eff showing fold/BKT boundary
- `output/vortex_rg_results.json` — computed crossover parameters

---

## Workstream 4: Hopf Integral on N-Layer Ansatz

### Physics Goal
Compute the Hopf invariant Q_H numerically from an explicit 5-layer texture
field on S³ (or equivalently on R³ with boundary conditions). Verify Q_H = 1/N
for the proposed configuration.

### Mathematical Specification
```
Hopf invariant (Whitehead integral):
  Q_H = (1/4π²) ∫_{R³} A ∧ F

where F = dA is the pulled-back curvature of the Hopf map n: S³ → S²,
and A satisfies dA = F (the Chern-Simons gauge field on S³).

For N-layer configuration with linking number 1/N:
  n(r) = (sin f(r) cos(φ/N), sin f(r) sin(φ/N), cos f(r))
  where f(r) = π(1 − r²/R²)^k for r ≤ R, and f→0 as r→∞.

Numerical evaluation:
  Q_H = (1/4π²) ∫ ε^{ijk} A_i ∂_j A_k d³x
  with A_i computed from n via: A_i = ε^{abc} n_a ∂_i n_b n_c / (1+n_z)
```

### Acceptance Criteria
1. On standard Hopf texture (N=1): Q_H = 1.000 ± 0.01
2. On N=5 ansatz: Q_H = 0.200 ± 0.01
3. Grid convergence: stable at n_grid³ = 64³ vs 128³ to 1%
4. Topological robustness: smooth deformation preserves Q_H

### Deliverables
- `notebooks/04_hopf_integral.ipynb` — computation + 3D visualization
- `qvc/hopf_integral.py` — reusable computation module
- Figure: Q_H vs N for N=1..10, showing 1/N law
- Figure: 3D texture visualization (cross-sections of n-field)

---

## Workstream 5: Lindblad TQGL Dynamics

### Physics Goal
Determine whether the gap dynamics under open-system (Lindblad) evolution
preserves the fold-predicted Δ* = 0.2324 meV, or whether decoherence
destroys the self-consistent fixed point on experimentally relevant timescales.

### Mathematical Specification
```
Lindblad master equation for gap-mode density matrix:
  dρ/dt = −(i/ℏ)[H_eff, ρ] + Σ_k γ_k (L_k ρ L_k† − ½{L_k† L_k, ρ})

Effective Hamiltonian in number basis (n̂ = Δ†Δ):
  H_eff = ω_0 n̂ + U n̂(n̂−1) + V_pump(t)
  ω_0 = 2Δ*/ℏ,  U = g₀

Lindblad operators:
  L₁ = √γ_loss · â         (cavity loss; γ_loss ~ 1/τ_cav ~ 1/(30 ps))
  L₂ = √γ_deph · n̂        (pure dephasing; γ_deph ~ 1/T₂* ~ 1/(50 ps))
  L₃ = √γ_phonon · â†â     (phonon bath at T=20 mK)

Compare:
  τ_coherent ~ 1 ps (GPE gap formation)
  τ_decoherence ~ 30–50 ps (cavity + dephasing)
  τ_DCE ~ 80 ps (self-limiting)

Key question: Is τ_coherent ≪ τ_decoherence so that gap forms before
dissipation destroys it?
```

### Acceptance Criteria
1. Steady-state ⟨Δ⟩ from Lindblad = Δ* from fold within 20%
2. Clear separation of formation vs decoherence timescales
3. Decoherence rate sensitivity: how much γ can increase before fixed point dies
4. Comparison: coherent (von Neumann) vs Lindblad trajectory for same IC
5. Temperature dependence: stable up to T* = 118.7 mK

### Deliverables
- `notebooks/05_lindblad_tqgl.ipynb` — dynamics + steady state
- `qvc/tqgl/lindblad.py` — Lindblad evolution module
- Figure: Δ(t) coherent vs Lindblad
- Figure: Steady-state Δ vs γ_loss (robustness curve)
- Figure: Wigner function at steady state

---

## Workstream 6: TikZ/PGFPlots Figures (Locked Option A*)

### Physics Goal
All publication figures reflect Option A* numbers. No deprecated values appear.

### Figures to Generate/Update

| Figure | Content | Key Data |
|--------|---------|----------|
| fig_fold_bifurcation | Gap δ vs α with fold at α_c=0.1818 | ℓ=−2.7958, two branches |
| fig_spectrum | Spec(G) showing λ_min=−g₀, λ_max=4g₀ | N=5, g₀=0.0982 meV |
| fig_vacuum_F | F(x) convergence + value F=0.58068 at x=1.005 | Bessel sum |
| fig_dce_selflimit | A(t) approaching A*/A_c=0.9117 | DCE ODE timeseries |
| fig_hall_staircase | σ_xy steps at 1/N = 0.2 | δσ_xy = e²/(10h) |
| fig_lindblad_gap | Δ(t) coherent vs Lindblad | τ comparison |
| fig_gij_convergence | Spec(G) vs ξ_K/λ_M approaching democratic | Quadrature results |
| fig_vortex_flow | RG flow in (K⁻¹, y) plane | KT flow solver |
| fig_hopf_Q | Q_H vs N showing 1/N law | Numerical Hopf integral |

### Format Requirements
- PGFPlots reading .csv data files from `output/`
- Colorblind-safe (Okabe-Ito palette)
- Single-column width: 86 mm; double: 178 mm
- Font: 8 pt minimum at final size
- Vector output (PDF via pdflatex)

### Deliverables
- `tex/figures/data/` — CSV files for each plot
- `tex/figures/*.tex` — standalone TikZ source
- `notebooks/06_generate_figures.ipynb` — data generation + matplotlib preview
- Compiled PDF figure set

---

## Implementation Order

1. **Vacuum convergence** (Workstream 1) — fastest to complete; upgrades lock status
2. **G_ij quadrature** (Workstream 2) — infrastructure exists; needs convergence sweep
3. **TikZ figures** (Workstream 6) — can proceed in parallel with data from 1+2
4. **Vortex RG** (Workstream 3) — self-contained diagnostic
5. **Hopf integral** (Workstream 4) — requires careful 3D grid implementation
6. **Lindblad TQGL** (Workstream 5) — builds on existing coupled.py

## Technical Infrastructure

```
notebooks/
  01_vacuum_convergence.ipynb
  02_gij_wigner_seitz.ipynb
  03_vortex_rg_flow.ipynb
  04_hopf_integral.ipynb
  05_lindblad_tqgl.ipynb
  06_generate_figures.ipynb
```

All notebooks import from `qvc/` package. Figures saved to `output/figures/`.
Data saved to `output/data/`. JSON results to `output/results/`.
