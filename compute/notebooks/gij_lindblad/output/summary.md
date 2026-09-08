# L1 dictionary — B1 $G_{ij}$ and B2 Lindblad GPE

Computed kernel model + truncated GPE with phenomenological dissipation.
Not torsion, not a linearized TQGL theorem, not Hubbard / CS-in-EOM / H-TOR1 / holography.
Hopfion ≠ 2D phase vortex ≠ CS fluxon. CS/Berry/Skyrme truncated.

## Locks (Option A*)

- $(r,R)=(8,50)$ nm, $\hbar c_s=0.07$ meV·µm, $\Delta_0=0.600$ meV
- $E_{\mathrm{vac}}=0.1287$ meV, $\lambda^\star=0.7628$
- **BRIDGE algebra:** $g_0^\star=\lambda^\star E_{\mathrm{vac}}=0.0982$ meV
- Democratic Spec (algebra): $\lambda_{\max}=4g_0$, $\lambda_{\min}=-g_0$ (mult 4). **Not inserted into $G_{ij}$.**

## B1  $G_{ij}$  (computed kernel, K0 units ≠ meV)

### WS Gaussian (finer grid)

- grid 80, $g_{0,\mathrm{eff}}=0.000136454$ (kernel)
- Spec $= [-0.00016 -0.00016 -0.00012 -0.00012  0.00055]
- $\lambda_{\min}=-0.00015685$, $\lambda_{\max}=0.00054582$
- rel Frobenius to democratic at $g_{0,\mathrm{eff}}$: **0.0667**
- near-democratic (rel Fro < 0.15): True

### Ring Fourier / Wannier

- grid 96, $g_{0,\mathrm{eff}}=0.000547169$ (kernel)
- Spec $= [-0.00074 -0.00072 -0.00043 -0.00029  0.00219]
- $\lambda_{\min}=-0.0007385$, $\lambda_{\max}=0.0021888$
- rel Frobenius to democratic at $g_{0,\mathrm{eff}}$: **0.1564**
- near-democratic: False

$g_{0,\mathrm{eff}}$ is **not** forced equal to $g_0^\star$. Different units, different objects.

## B2  Lindblad scan

Coherent 1 ps target (published truncated GPE+DCE): $\Delta_{\mathrm{live}}\approx 0.6305$ meV, $g^{(1)}\approx 0.8832$.

### Table, 1 ps (coupled DCE)

| label | channel | $\tau_{\mathrm{dec}}$ (ps) | $\Delta_{\mathrm{live}}$ (meV) | $g^{(1)}$ |
|---|---|---:|---:|---:|
| coherent | coherent | — | 0.6305 | 0.8832 |
| tau_5 | number+dephasing | 5 | 0.5720 | 0.7102 |
| tau_20 | number+dephasing | 20 | 0.6152 | 0.8351 |
| tau_50 | number+dephasing | 50 | 0.6248 | 0.8639 |
| number_tau_5 | number | 5 | 0.5697 | 0.8806 |
| dephasing_tau_5 | dephasing | 5 | 0.6309 | 0.7031 |

### Table, 20 ps window (V_ZPF frozen, no DCE)

| label | channel | $\tau_{\mathrm{dec}}$ (ps) | $\Delta_{\mathrm{live}}$ (meV) | $g^{(1)}$ |
|---|---|---:|---:|---:|
| coherent | coherent | — | 0.6623 | 0.1075 |
| tau_5 | number+dephasing | 5 | 0.0872 | 0.0102 |
| tau_20 | number+dephasing | 20 | 0.3920 | 0.0812 |
| tau_50 | number+dephasing | 50 | 0.4761 | 0.0032 |

## Computed vs ansatz

- **Computed:** K0 overlaps; Spec(G); Δ_live(ψ); g^{(1)}(ψ); Lindblad trajectories.
- **Locked algebra:** democratic (4g0, −g0×4); g0★=λ★ E_vac.
- **Phenomenological:** γ=1/τ_dec; local phase kicks.
- **Not:** J−I insertion, torsion theorem, 92%, prescribed Δ(t).

