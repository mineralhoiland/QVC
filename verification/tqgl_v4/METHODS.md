# TQGL v4 methods (equations actually solved)

Status: truncated 2D GPE + DCE ODE. Option A* fat torus (r,R)=(8,50) nm.

## Truncation (not in the EOM)

Chern–Simons, Berry connection, and Skyrme/Hopf terms are **not** discretised.
Democratic Spec(G) is a static lock (λ_min=−g₀, multiplicity N−1), not a GPE term.
The 1 ps window does **not** include Lindblad dissipation.

## GPE (split-step Fourier)

$$
i\hbar\partial_t\psi = \Bigl[-\tfrac{\hbar^2}{2m^*}\nabla^2 + \alpha|\psi|^2 + \beta|\psi|^4 + V_{\rm ZPF}(\mathbf r,t) + V_{\rm ph}(\mathbf r,t)\Bigr]\psi
$$

- Box: 200 nm, grid 128×128, dt=2.0 fs, T=1.0 ps.
- ℏ²/(2m*)=Δ₀ a_H² with a_H=r=8 nm; α=−Δ₀, β=Δ₀/n₀, n₀=1.
- Hopfion IC: 10% core depletion, phase Q_H θ with Q_H=1/N, noise 3%.
- $V_{\rm ZPF}=-E_{\rm vac}\,(A(t)/A(0))\,\mathrm{Gaussian}_{\rm ring}$. Drive is locked $E_{\rm vac}=0.1287$ meV, **not** 0.244. Scale vs withdrawn Casimir: 0.527.
- $V_{\rm ph}=0.1\Delta_0\sin(\omega_0 t)\,\mathrm{Gaussian}_{\rm ring}$, $\omega_0=2\Delta_0/\hbar$.

## Live gap (uses ψ; no time schedule)

$$
\Delta_{\rm live}=\Delta_0\sqrt{n_w/n_{\rm ref}},\qquad n_w=\langle|\psi|^2\rangle_{\rm ring}.
$$

Compared to the static physical fold branch Δ(A) at the running DCE amplitude, and to Δ★(λ★)=0.2324 meV. If Δ_live<Δ_c=0.1091 meV the run is flagged unphysical vs the static fold. There is **no** 0.600→0.389→0.032 prescription. Do not quote 92%.

## DCE ODE (RK4)

$$
\dot A=\bigl[\Lambda(\Delta)\,\mathcal L(\delta\omega/\Gamma)-\Lambda_{\rm loss}\bigr]A,\quad
\Lambda(\Delta)=(m\Delta/2\hbar)\,Q/(2\pi),\quad \mathcal L(x)=1/(1+x^2).
$$

- Resonant: ω_pump/2 = 2Δ★/ℏ. Detuning-feedback: ω_pump/2 = 2Δ★/ℏ + 2Γ.
- Λ_loss = Γ_cav(Δ★) (cavity linewidth), **not** fitted to A*/A_c=0.52.
- A_c = α_c Δ₀/(ℏ c_s) = 1.5582 µm⁻¹ (fold, not BKT 1.2233).
- Coupled 1 ps: Δ = Δ_live(ψ). Long window (~80 ps): Δ slaved to the physical fold branch Δ(A). Two-timescale justification: Γ_cav⁻¹ ~ tens of ps ≫ 1 ps GPE window.

Fold-slaved resonant A*/A_c (computed, not imposed): 0.9117.

## Schwinger / Berry (diagnostics only)

E_pair = 2|λ_min| with |λ_min|=g₀=λ★ E_vac (democratic). Withdrawn E_pair=0.172 is not used unless that number is re-derived. Berry GOE/GUE Monte Carlo is a random-matrix probe, **separate** from democratic G.
