# Track A L1 parameter dictionary

Fat-torus Option A* working point. Three lengths are distinct. Catalyzon is a parallel channel (not stacked with fold or live-gap). Hall algebra uses CS level $k=2N$ with $Q_H=1/N$.

## Parameter card

- $(r,R)=(8,50)$ nm, $\hbar c_s=0.07$ meV·µm, $\Delta_0=0.600$ meV  **[locked]**
- $F=0.58068$, $E_{\mathrm{vac}}=0.1287$ meV, $\alpha_c=\delta_c=0.1818$, $\Delta_c=0.1091$ meV  **[locked / derived]**
- $\lambda_{\mathrm{fold}}=0.8475$, $\lambda^\star=0.7628$, $\Delta^\star=0.2324$ meV  **[derived / locked]**
- $T^*=118.7$ mK, $E_{\mathrm{pair}}=0.1963$ meV, $N=5$, $Q_H=0.2$, $k=10$  **[computed / locked]**

### Lengths

- $a_H=r=0.008$ µm $=8$ nm, $\eta_{\mathrm{ret}}(a_H)=0.137$  **[locked]**
- $\xi_{\mathrm{ph}}=\hbar c_s/\Delta_0=0.1167$ µm, $\eta_{\mathrm{ret}}(\xi_{\mathrm{ph}})=2$ identically  **[derived]**
- $\xi_{181}=0.181$ µm (alternate; preprint coherence $\hbar c_s/\Delta_{\mathrm{QVC}}$), $\eta_{\mathrm{ret}}=3.10$  **[alternate, not the lock]**

### GL coefficients (matched, not Hubbard)

- $\alpha=\lambda^\star E_{\mathrm{vac}}/\Delta_0=0.1636$  **[derived]**
- $\alpha_{\mathrm{GL}}=-\Delta_0=-0.600$ meV  **[matched]**
- $\kappa_{\mathrm{GL}}=\Delta_0 a_H^2=3.8400e-05$ meV·µm$^2$  **[derived]**
- $\lambda=\lambda^\star=0.7628$  **[locked]**
- $n_{\mathrm{ref}}=0.86435$ on hopfion_initial  **[computed]**
- $\beta=\Delta_0/n_{\mathrm{ref}}=0.6942$ meV  **[matched]**
- $\kappa_S\approx\Delta_0 a_H^4=2.4576e-09$ meV·µm$^4$  **[ansatz]**

### Catalyzon (parallel channel)

- $g_0^\star=\lambda^\star E_{\mathrm{vac}}=0.0982$ meV, $\Delta_{\mathrm{cat}}=0.5018$ meV (16.4%)  **[derived]**
- at $\lambda_{\mathrm{fold}}$: $g_0=\Delta_c=0.1091$ meV (18.2%=$\alpha_c$)  **[derived]**
- Spec$(G)$: $\lambda_{\min}=-g_0^\star$, multiplicity $N-1=4$  **[locked algebra]**

### Mexican hat

- GPE-sign well-depth change: 35.4% (attractive $V_{\mathrm{ZPF}}$, well deepens)  **[computed]**
- Catalyzon-dressed barrier reduction: 30.0% $=1-(\Delta_{\mathrm{cat}}/\Delta_0)^2$ (replaces cartoon 42%)  **[computed]**

### Hall and matching

- $\delta\sigma_{xy}=e^2/(k h)=3.8740e-06$ S, parent $\sigma_{xy}=(k/2)e^2/h=1.9370e-04$ S  **[derived]**
- $R_{\mathrm{WS}}=\hbar c_s/(\pi k_B T^*)=2.18$ µm  **[matched]**
- Berry $R_{\mathrm{mean}}=0.908$ is a ratio, not a radius  **[computed]**

### Coherent GPE baseline (Lindblad is another agent)

- $\Delta_{\mathrm{live}}(1\,\mathrm{ps})=0.6305$ meV, $g^{(1)}=0.8832$  **[computed]**

## Status ledger

### locked
- fat torus (r,R)=(8,50) nm, ħ c_s=0.07 meV·µm, Δ₀=0.600 meV
- a_H = r (Option A* GL/GPE healing freeze)
- η_ret = 2 Δ₀ ξ / (ħ c_s) as algebra (ξ remains a choice)
- λ★ = 0.90 λ_fold, N=5, Q_H=1/N, k=2N
- F, E_vac, α_c=δ_c, Δ_c from the frozen-cavity bridge

### derived
- ξ_ph = ħ c_s / Δ₀ and η_ret(ξ_ph)=2
- η_ret(a_H)=2 Δ₀ a_H/(ħ c_s)≈0.137
- κ_GL = Δ₀ a_H²
- α = λ E_vac / Δ₀ (gap-map coupling)
- g0 = λ E_vac, Δ_cat = Δ₀ − g0, Spec(G): λ_min=−g0
- λ_fold = α_c Δ₀ / E_vac
- Hall: k=2N, δσ_xy = e²/(k h) = Q_H e²/(2h)
- pairing k Q_H = 2

### matched
- α_GL = −Δ₀ (Landau mass |α_GL|=Δ₀)
- β = Δ₀ / n_ref so the Mexican-hat vacuum sits on hopfion n_ref
- R_WS = ħ c_s/(π k_B T*) at T* (analogue Hawking matching, not duality)

### ansatz
- κ_S ≈ κ_GL a_H² = Δ₀ a_H⁴ (Derrick balance on the locked tube)

### computed
- n_ref from hopfion_initial (ring-weighted |ψ|²)
- Mexican-hat well depths and catalyzon barrier reduction
- T*, E_pair from Schwinger twin-channel at g0★
- live gap 0.6305 meV / g(1)=0.8832 coherent GPE baseline (not rerun here)
- Berry R_mean≈0.908 (ratio, not a radius)

### alternate_not_lock
- ξ=0.181 µm → η_ret≈3.1 (preprint coherence; ħ c_s/Δ_QVC withdrawn)
- R_WS(285 mK)≈0.91 µm (same matching formula, withdrawn T*)

### retired
- g0~0.21 as a derived coupling
- κ_tor≈0.25 as an insertion
- cartoon α0=−0.5, α_QVC=−0.35, 42% barrier drop
- holographic γ_QM≈3.73 meV·nm² as κ_GL
- folklore 0.244 / 0.389 / 0.032 / 92%

### out_of_scope
- Hubbard microscopic derivation of GL coefficients
- Chern–Simons in the GPE EOM
- H-TOR1 torsion vertex
- holographic dictionary

## Notes

- Democratic $\lambda_{\min}=-g_0$ check: True.
- $n_{\mathrm{ref}}$ formula: n_ref = ⟨|ψ|²⟩_{ring} on hopfion_initial.
- $\kappa_S$ ansatz: Faddeev–Skyrme 3D Derrick on locked tube a_H; 2πR cancels.
- Do not stack catalyzon, fold, and live-gap percentages.
- Do not box $\eta_{\mathrm{ret}}=3.1$ as the unique lock.
- Retired insertions: $g_0\sim0.21$ meV, $\kappa_{\mathrm{tor}}\approx0.25$.

