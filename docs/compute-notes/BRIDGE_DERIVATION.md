# Derived vacuum bridge and freeze decision

MCP scientific-agent servers were unavailable in this session
(`user-claude-skills` error; Hugging Face MCP needs auth). Derivation
uses the locked QVC gap equation, Bessel form, and fold theorem, then
numeric solvers in `qvc/`.

## 1. Problem

The preprint identified \(E_{\mathrm{Cas}}=0.244\,\mathrm{meV}\) and a
self-consistent \(\Delta_{\mathrm{QVC}}=0.389\,\mathrm{meV}\) (35%).
Audit: the old Casimir formula does **not** produce 0.244 on the stated
geometry, and \(\alpha=E/\Delta_0\) was never derived as a map into the
fold equation.

## 2. Proven algebra (gap equation)

Locked gap equation:

\[
\Delta = \Delta_0 + \hbar c_s A\left[\ln\frac{\Delta a_H}{2\hbar c_s}+\gamma_E\right],
\qquad \Delta_0=2J.
\]

**Frozen cavity** (preprint torus): \(a_H=r=\mathrm{const}\). Set
\(\delta=\Delta/\Delta_0\):

\[
\ell(r)=\ln\frac{\Delta_0 r}{2\hbar c_s}+\gamma_E,
\qquad
\alpha=\frac{\hbar c_s A}{\Delta_0},
\qquad
\delta=1+\alpha(\ln\delta+\ell).
\]

This is the Corrections fold map. Numerically
\(\ell(8\,\mathrm{nm})=-2.7958\) vs \(\ell_{\mathrm{Corr}}=-2.791\)
(relative error \(0.17\%\)). So \(\ell_{\mathrm{Corr}}\) **is** the
frozen-cavity logarithm, not a separate IR fit.

**Running healing** \(a_H=\hbar c_s/\Delta\) forces \(\ell_{\mathrm{tree}}=\ln\tfrac12+\gamma_E\approx-0.116\)
and exponentially small \(F\) when \(a_H\gg R\). That is a different
theory; it is not what \(\ell\approx-2.79\) describes.

## 3. Proven map \(A\leftrightarrow E_{\mathrm{vac}}\)

Locked Tier-C form:

\[
E_{\mathrm{vac}}=\frac{1}{4\pi^2}\frac{\hbar c_s}{r}F\!\left(\frac{2\pi r}{R}\right),
\qquad
A_{\mathrm{geom}}=\frac{1}{4\pi^2}\frac{1}{r}F
=\frac{E_{\mathrm{vac}}}{\hbar c_s}.
\]

Unit participation \(\lambda=1\) (condensate sees the full geometric
amplitude):

\[
\boxed{\alpha=\lambda\,\frac{E_{\mathrm{vac}}}{\Delta_0}.}
\]

\(\sqrt{n}\) is **not** inserted: \(F=\sum K_0(nx)\) already sums
windings. Extra quadrature would be a distinct channel.

## 4. Freeze option chosen: **A\***

| Item | Choice |
|------|--------|
| Formula | Tier-C \(E_{\mathrm{vac}}=(1/4\pi^2)(\hbar c_s/r)F\) |
| Geometry | Fat torus \((r,R)=(8,50)\,\mathrm{nm}\) |
| \(\hbar c_s\) | Package \(0.07\,\mathrm{meV\cdot\mu m}\) |
| Bridge | \(\alpha=E_{\mathrm{vac}}/\Delta_0\), \(\ell=\ell(r)\) |
| Status | **FROZEN_WORKING** (form+bridge+geometry). meV is computed, not folklore. |

**Not chosen:** moiré \((0.335,13.4)\,\mathrm{nm}\) (far past fold);
retuning \(\hbar c_s\) to force 0.244; \(\sqrt{n}\) in the static loop.

## 5. Computed replacements for withdrawn numbers

| Withdrawn | Computed (A\*) |
|-----------|----------------|
| \(F=0.596\) | \(F=0.58068\) |
| \(E=0.244\,\mathrm{meV}\) | \(E_{\mathrm{vac}}=0.1287\,\mathrm{meV}\) |
| \(A\xi=0.296\) | \(AR=0.0919\) |
| \(\Delta^*=0.389\,\mathrm{meV}\) | **no real branch** (\(\alpha/\alpha_c\approx 1.180\)) |

Fold theorem: \(\alpha_c=\delta_c\approx 0.1818\), \(E_{\mathrm{vac}}^{\max}=\alpha_c\Delta_0\approx 0.109\,\mathrm{meV}\).
Package fat-torus energy exceeds this by 18%. That is a **proof of
absence** of a static mean-field \(\Delta^*\) at \(\lambda=1\), not a
solver failure.

The withdrawn 35% reduction would require \(\alpha\approx 0.109\) and
\(E_{\mathrm{vac}}\approx 0.065\,\mathrm{meV}\) — a different energy,
not 0.244 and not 0.129.

Companion diagnostics (not the freeze):
- \(\hbar c_s\) at fold \(\approx 0.0593\,\mathrm{meV\cdot\mu m}\)
- \(\lambda\approx 0.76\) would sit at \(0.90\,\alpha_c\)

## 6. What is proven vs computed vs open

| Claim | Grade |
|-------|-------|
| \(\delta=1+\alpha(\ln\delta+\ell)\), \(\alpha_c=\delta_c\) | Proven |
| \(\ell(r)\) formula and numerical match to Corrections | Proven |
| \(\alpha=E_{\mathrm{vac}}/\Delta_0\) at \(\lambda=1\) | Proven identification |
| Fat-torus \(F,E_{\mathrm{vac}}\) | Proven given frozen triple |
| No \(\Delta^*\) at \(\lambda=1\) | Proven by fold theorem |
| Absolute meV as a measured MATBG energy | **Not** claimed |
| \(G_{ij}\) from torsion action variation | Open (H-TEGR1) |
| \(G_{ij}\) from TQGL–WS \(K_0\) overlaps | Program-computed; \(\lambda_{\min}<0\) without inserting \(J-I\) |

## 7. Self-consistency statement

The self-consistent gap equation **is** the fold map with derived
\((\alpha,\ell)\). Its solution at the frozen triple is: the physical
branch does not exist. Remaining gapped physics at these parameters,
if any, is **not** this static vacuum loop (catalyzon Spec\((G)\) is a
separate channel).
