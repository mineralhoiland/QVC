# Closed-form locks (Tier C)

These identities are implemented in `qvc/` and checked by
`scripts/run_verification.py`. Distinguish **proven identities** from
**provisional numerics**.

## 1. Democratic coupling matrix

For \(G = g_0(J_N - I_N)\):

\[
\lambda_{\max} = (N-1)g_0 \quad (\mathrm{mult.\ }1),\qquad
\lambda_{\min} = -g_0 \quad (\mathrm{mult.\ }N-1).
\]

Helmert (sum-zero) basis spans the dark / catalyzon subspace.

**False claims (audit):** \(\lambda_{\min}=-4g_0\) (preprint, \(N=5\)) and
\(\lambda_{\min}=-5g_0\) (supplemental octahedral) are **not** spectra of
the democratic matrix.

**Catalyzon gap shift (proven on Spec(G)):**

\[
\Delta_{\mathrm{cat}} = \Delta_{\mathrm{ex}} - g_0.
\]

Module: `qvc/coupling_matrix.py`.

## 2. Gap equation fold

\[
\delta = 1 + \alpha(\ln\delta + \ell).
\]

At the saddle-node:

\[
\alpha_c = \delta_c,\qquad
\delta_c(1 - \ln\delta_c - \ell) = 1.
\]

Near the fold:

\[
\Delta - \Delta_c \sim \pm\sqrt{A_c - A}.
\]

This is a **fold**, not BKT vanishing, unless a vortex RG is derived.

**Tree-level log coefficient** (log argument \(1/2\)):

\[
\ell_{\mathrm{tree}} = \ln\frac12 + \gamma_E \approx -0.11593.
\]

Module: `qvc/gap_fold.py`.

## 3. Chern–Simons / linking on \(L(N,1)\)

\[
\mathrm{CS}[A_k] = \frac{k^2}{N},\qquad
\mathrm{lk}([\gamma],[\gamma]) = \frac1N,\qquad
Q_H = \frac1N,\qquad
\theta = \frac\pi N.
\]

**Topology lock:** \(H_2(L(N,1))=0\). There is no fundamental 2-cycle;
any preprint proof that fractionalizes via an \(H_2\) class is **invalid**.
Use CS / linking on \(H_1\cong\mathbb{Z}_N\).

Module: `qvc/lens_chern_simons.py`.

## 4. Vacuum spectral form

Locked **form** (2D phase-space measure):

\[
E_{\mathrm{vac}} = \frac{1}{4\pi^2}\frac{\hbar c_s}{r}\, F\!\left(\frac{2\pi r}{R}\right),\qquad
F(x)=\sum_{n=1}^\infty K_0(nx).
\]

Mode quadrature:

\[
A_{\mathrm{eff}} = \sqrt{n}\, A_{\mathrm{single}}.
\]

**Do not insert \(\sqrt{n}\) into the static gap loop by default:** \(F=\sum K_0(nx)\) already sums toroidal windings.

**Derived bridge (proven algebra, frozen-cavity \(a_H=r\)):**

\[
\ell(r)=\ln\frac{\Delta_0 r}{2\hbar c_s}+\gamma_E,\qquad
A=\frac{E_{\mathrm{vac}}}{\hbar c_s},\qquad
\alpha=\lambda\frac{E_{\mathrm{vac}}}{\Delta_0}.
\]

At \(\lambda=1\) on the fat torus, \(\ell(r)=-2.7958\) matches Corrections \(\ell=-2.791\) to \(0.17\%\). Computed \(E_{\mathrm{vac}}=0.1287\,\mathrm{meV}\) **replaces** \(0.244\,\mathrm{meV}\). Unit-participation \(\alpha/\alpha_c\approx 1.18\) is past the fold (no \(\Delta^*\)). See `docs/BRIDGE_DERIVATION.md`.

**Geometry \((r,R)\) is frozen working (fat torus).** \(\hbar c_s\) remains a platform scale. The old Casimir formula
\((1/4\pi)(\hbar c_s/r)(R/r)F\) does **not** reproduce the claimed
\(0.244\,\mathrm{meV}\) on the preprint geometry.

Modules: `qvc/vacuum_freeze.py`, `qvc/mode_enhancement.py`.

## 5. Retardation

\[
\eta_{\mathrm{ret}} = \frac{2\Delta_0\,\xi}{\hbar c_s}.
\]

With \(\xi=a_H=\hbar c_s/\Delta_0\), one obtains \(\eta_{\mathrm{ret}}=2\).

Module: `qvc/retardation.py`.

## 6. Acoustic \(R_{\mathrm{WS}}\) matching (not a proof)

\[
T_H = \frac{\hbar c_s}{\pi k_B R_{\mathrm{WS}}}
\quad\Rightarrow\quad
R_{\mathrm{WS}} = \frac{\hbar c_s}{\pi k_B T^*}.
\]

**Label:** matching condition from analogue gravity — **not** a QED/gravity
duality derivation.

Module: `qvc/acoustic_matching.py`.

## Lock table

See `CLOSED_FORM_LOCKS` in `qvc/canonical_formalism.py` for machine-readable
status tags (`PROVEN`, `PROGRAM_PREDICTION`, `PROVEN_FORM_GEOMETRY_OPEN`,
`MATCHING_NOT_PROOF`).
