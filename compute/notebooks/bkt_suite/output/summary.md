# QVC / TEVC vortex RG — C1–C6 summary

BKT is treated as a physically interesting vortex sector. C1 is partly ansatz;
this suite is **not** a BKT theorem. Hopfion ≠ 2D phase vortex ≠ CS fluxon.
Chern–Simons is truncated unless the anyon $K_c$ option is on. Option A* freezes $\alpha$.

## Locked numbers (Option A*, fat torus)

- $(r,R)=(8,50)$ nm, $\hbar c_s=0.07$ meV·µm, $F=0.58068$, $E_{\mathrm{vac}}=0.1287$ meV
- $\alpha_c=\delta_c=0.1818$, $\Delta_c=0.1091$ meV, $\Delta_0=0.600$ meV
- $\lambda^\star=0.7628$, $\Delta^\star=0.2324$ meV, $T^*=118.7$ mK, $E_{\mathrm{pair}}=0.1963$ meV
- Democratic $\lambda_{\min}=-g_0=-0.0982$ meV, $Q_H=1/N=0.2$

## Competing-scale test (H-VRG1)

- Primary closure (healing-disk $n_0$ + frozen $\kappa_{\mathrm{geom}}$) at $T^*$: **KT before fold = False**
- Conventional-only (no geometric floor) at $T^*$: **KT before fold = True** (IR plasma at $\alpha/\alpha_c\approx 0.945$, finite-$y$ separatrix; bare $K$ stays above $2/\pi$)
- $K(\alpha=0)=46.694$, $K(\alpha_c)=10.573$, $K_c=2/\pi=0.6366$
- $\kappa_{\mathrm{UV}}=0.4775$ meV, $\kappa_{\mathrm{fold}}=0.1081$ meV

H-VRG1 is a competing-scale test on a computed-or-ansatz κ(α), not a BKT theorem. Option A* keeps α frozen. CS is truncated in the default scan.

## C2 core

- $E_{\mathrm{core}}/\kappa=0.3427$ (cutoff: healing_length_ansatz)
- $a_c=0.008000$ µm, $\xi=0.008000$ µm

## C5 fold-branch comparison (not a BKT disproof)

- free slope of $\ln(\delta-\delta_c)$ vs $\ln t_A$: 0.5421 (fold theorem: $1/2$)
- $\Delta\ell\ell$ (√ minus essential): 431.329

## C6 H-VRG2

- sharp_jump_possible

## Computed vs assumed

### C1
- **computed:** physical-branch δ_+(α) from locked gap map; ⟨ρ²⟩=δ² via live-gap identification; uniform phase-twist consistency of 2 n0 κ_GL ⟨ρ²⟩; hopfion ring-weighted ⟨|ψ|²⟩ (amplitude sector only)
- **derived:** two-sector identity κ=2 κ_GL ⟨ρ²⟩+κ_geom (algebra); κ_GL=Δ₀ a_H² from GPE kinetic match; K=κ/T_eff in Kosterlitz convention
- **ansatz:** n0=1/(π a_H²) healing-disk density (primary); Peotta–Törmä κ_geom=(C_eff/2π)Δ₀ frozen geometric floor; alternate moiré-cell n0 and C_eff=Q_H variants

### C2
- **computed:** radial GP winding-1 profile; E_v on a disk; half-density core radius; E_core = E_v − πκ ln(R/a_c) given the cutoff
- **derived:** GP energy functional matching the truncated EOM
- **ansatz:** a_c=ξ healing-length cutoff used in the IR logarithm unless noted

### C3
- **computed:** Kosterlitz ODE at frozen α; physical-branch competing-scale scan; T_eff grid winner table
- **derived:** H-VRG1 criterion (KT surface before fold); y=exp(−E_core/T_eff)
- **truncated:** Chern–Simons (default)

### C4
- **computed:** ξ_KT from flow ℓ* when unbinding occurs
- **conjecture:** Δ∼ħ c_s/ξ_KT as a gap identification

### C5
- **computed:** OLS slopes and Gaussian log-likelihood on ≳1.5 decades of t_A; locked fold branch vs √t_A and vs essential singularity

### C6
- **computed:** AHNS ratio (ħ/τ_v)/(πκ) along the branch given the estimates
- **labeled_estimate:** ħ/τ_v ∼ g_ph=0.1Δ₀ from existing GPE phonon potential; ħ/τ_v ∼ 2π λ_ep kT with λ_ep=1 Gao–Khalaf-type OOM

Withdrawn Casimir / gap-schedule folklore is not used as physics.

