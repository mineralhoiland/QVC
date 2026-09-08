# Locked Option A* parameters

These numbers match `paper/macros.tex` and `compute/qvc/tqgl/locks.py`.
Do not retune for the v4 preprint. Machine-readable copy:
[`optionA_star.json`](optionA_star.json). Vacuum-only freeze:
[`frozen_optionA_converged.json`](frozen_optionA_converged.json).

## Geometry and vacuum

| Symbol | Value | Notes |
| --- | --- | --- |
| \(r, R\) | \(8,\ 50\) nm | Fat torus |
| \(\hbar c_s\) | \(0.07\) meV·µm | Dominant \(E_{\mathrm{vac}}\) uncertainty (~10%) |
| \(E_{\mathrm{vac}}\) | \(0.1287\) meV | Corrected Bessel self-energy |
| \(F\) | \(0.58068\) | \(\sum_n K_0(nx)\), \(x=R/r\) |
| \(\ell\) | \(-2.7958\) | Frozen log coefficient |

Withdrawn folklore (not physics): \(E_{\mathrm{Cas}}=0.244\) meV,
\(\Delta_{\mathrm{QVC}}=0.389\) meV, \(\Delta_{1\mathrm{ps}}=0.032\) meV, 92%.

## Fold and star point

| Symbol | Value |
| --- | --- |
| \(\alpha_c=\delta_c\) | \(0.1818\) |
| \(\Delta_c\) | \(0.109\) meV |
| \(\lambda_{\mathrm{fold}}\) | \(0.8475\) |
| \(\lambda^\star\) | \(0.7628\) |
| \(\Delta^\star\) | \(0.2324\) meV |
| \(\lambda_{\min}\) | \(-g_0\) (mult. \(N-1\)) |
| \(Q_H\) | \(1/N\) |

## Bridge, Hall, temperature

| Symbol | Value |
| --- | --- |
| \(g_0^\star=\lambda^\star E_{\mathrm{vac}}\) | \(0.0982\) meV |
| \(E_{\mathrm{pair}}\) | \(0.1963\) meV |
| \(T^\star\) | \(118.7\) mK |

The overlap mean \(g_{0,\mathrm{eff}}\) is a computed kernel number. Do not
force \(g_{0,\mathrm{eff}}=g_0^\star\).

## Parallel channels (not additive to the static fold)

| Channel | Locked numbers |
| --- | --- |
| Catalyzon Mexican hat | ~30% barrier reduction |
| Live GPE (1 ps) | \(\Delta=0.6305\) meV, \(g_1=0.8832\) |
| Cavity moiré \(N=5\) | \(\mathcal{F}=1.119\), \(A_{\mathrm{geom}}=1.839\,\mu\mathrm{m}^{-1}\), \(A_{\mathrm{eff}}=2.058\), \(\lambda_{\mathrm{fold}}^{(F)}=0.7573\), \(\lambda^\star_F=0.6816\), \(\Delta^\star_F=0.2914\) meV |
| Schwinger | \(\alpha^\star=0.1636\), exponent \(19.20\), \(\Gamma=530\) Hz. No \(\mathcal{F}_{\mathrm{moiré}}\) in the exponent |

Frozen Lindblad \(\tau=5\) ps at \(0.0872\) meV \(<\Delta_c\) is unphysical on
the static branch. H-FN1 is not realized (one Goldstone). H-TOR1 is not
evaluated.

## Phonon-loop IR cutoff

\(q_{\mathrm{IR}}=3.65\times10^3\,\mu\mathrm{m}^{-1} > \Lambda=125\,\mu\mathrm{m}^{-1}\).
The logarithm is a regulated phenomenological map, not a one-loop theorem.

## How to regenerate

```bash
cd compute
python scripts/run_verification.py
python -c "from qvc.tqgl.locks import build_frozen_tqgl; print(build_frozen_tqgl().as_dict())"
python -c "from qvc_l1.locks import build_l1_locks; print(build_l1_locks().as_dict())"
```
