# Gap reduction recomputation and participation λ

## Original three numbers (do not stack)

| Label | Claim | Actual arithmetic | Baseline mix |
|-------|--------|-------------------|--------------|
| 35% | 0.600 → 0.389 | 35.2% | static loop at \(E_{\mathrm{Cas}}=0.244\) |
| 52% | 0.600 → 0.29 | 51.7% | \((\Delta_{\mathrm{QVC}}-0.1)/\Delta_{\mathrm{ex}}\) |
| 92% | \(\Delta(1\,\mathrm{ps})=0.032\) | **91.8% of 0.389**, **94.7% of 0.600** | coherent TQGL; \(V_{\mathrm{ZPF}}=-0.244\) |

The 92% label is \((0.389-0.032)/0.389\), not a 92% drop from \(\Delta_0\).

## Fold channel → \(\lambda_{\mathrm{fold}}\)

\[
\alpha=\lambda\frac{E_{\mathrm{vac}}}{\Delta_0},\qquad
\lambda_{\mathrm{fold}}=\frac{\alpha_c\Delta_0}{E_{\mathrm{vac}}}.
\]

At unit participation the frozen fat torus is **past the fold** (\(\lambda_{\mathrm{fold}}<1\)).
The physical branch exists only for \(\lambda\le\lambda_{\mathrm{fold}}\).
Its deepest gap is \(\Delta_c=\delta_c\Delta_0\), a **maximum static reduction** \(1-\delta_c\approx 81.8\%\).

**Theorem:** 92% reduction (\(\delta=0.032/0.6\approx 0.053<\delta_c\)) is **not** a static physical solution for any \(\lambda\). It is either the unstable branch or a dynamical overshoot of a simulation that used withdrawn 0.244 meV.

## Catalyzon channel → \(\lambda_{\mathrm{cat}}\)

Keep Spec\((G)\) separate: \(\Delta_{\mathrm{cat}}=\Delta_{\mathrm{ex}}-g_0\).
If the **same** participation feeds off-diagonal hybridisation,

\[
g_0=\lambda E_{\mathrm{vac}},\qquad \lambda_{\mathrm{cat}}=\frac{g_0}{E_{\mathrm{vac}}}.
\]

At \(\lambda_{\mathrm{fold}}\) one has \(g_0=\alpha_c\Delta_0=\Delta_c\), so catalyzon reduction is \(\alpha_c\approx 18\%\), **not** 52%.
The 52% mixed \(\Delta_{\mathrm{QVC}}\) with a 0.1 meV guess and \(\Delta_{\mathrm{ex}}\).

Scaling the 0.1 meV guess with \(E_{\mathrm{vac}}/0.244\) gives a smaller \(g_0\) and \(\lambda_{\mathrm{cat}}=0.1/0.244\approx 0.41\).

## Recommended \(\lambda^\star\)

\[
\lambda^\star=0.90\,\lambda_{\mathrm{fold}}=0.90\,\frac{\alpha_c\Delta_0}{E_{\mathrm{vac}}}.
\]

Numerics (fat torus, Tier-C \(E_{\mathrm{vac}}=0.1287\,\mathrm{meV}\)):

| Object | Value |
|--------|-------|
| \(\lambda_{\mathrm{fold}}\) | 0.8475 |
| \(\lambda^\star=0.90\lambda_{\mathrm{fold}}\) | **0.7628** |
| Fold \(\Delta^\star(\lambda^\star)\) | 0.232 meV (**61.3%**) |
| Max static reduction | **81.8%** at \(\Delta_c=0.109\) meV |
| \(\lambda_{35}\) (withdrawn 0.389) | 0.508 |
| Catalyzon at \(\lambda^\star\) (parallel) | \(\Delta_{\mathrm{cat}}=0.502\) meV (16.4%) |
| \(\lambda_{\mathrm{cat}}\) from scaled \(g_0=0.1\times E_{\mathrm{vac}}/0.244\) | 0.410 |

Partition (unitarity \(\lambda_\alpha+\lambda_G\le 1\)): \(\lambda_\alpha\) feeds the fold, \(\lambda_G\) feeds \(g_0=\lambda_G E_{\mathrm{vac}}\). Equal split \((0.5,0.5)\) is subcritical. Unit \(\lambda_\alpha=1\) is fold-out.

Evaluate fold \(\Delta^\star(\lambda^\star)\) and catalyzon \(\Delta_{\mathrm{cat}}\) **in parallel**. Do not add 35+52+92.

## TQGL 92%

No split-step solver is in the repository. Drive scale \(E_{\mathrm{vac}}/0.244\approx 0.527\). Phenomenological exponential remaining \(\Delta\sim 0.128\,\mathrm{meV}\) (78.7%) is **not** a TDGL rerun. Do not quote a new 92%.