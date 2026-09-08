# Vacuum freeze specification

## Purpose

Prevent quoting a vacuum energy in meV until **one** spectral problem is
frozen: formula + geometry + scale + \(F\) definition.

## Locked form (Tier C)

\[
E_{\mathrm{vac}}
  = \frac{1}{4\pi^2}\frac{\hbar c_s}{r}\,
    F\!\left(\frac{2\pi r}{R}\right),\qquad
F(x)=\sum_{n=1}^\infty K_0(nx).
\]

Old (audit-failed) Casimir formula:

\[
E_{\mathrm{old}}
  = \frac{1}{4\pi}\frac{\hbar c_s}{r}\frac{R}{r}\, F.
\]

**Audit result (enforced in code):** on preprint geometry
\((r,R)=(8\,\mathrm{nm},50\,\mathrm{nm})\) with \(\hbar c_s=0.07\,\mathrm{meV\cdot\mu m}\),
\(E_{\mathrm{old}}\) is **not** \(0.244\,\mathrm{meV}\).

## Candidate geometries

| Label | \(r\) | \(R\) | Source |
|-------|-------|-------|--------|
| `preprint_fat_torus` | 8 nm | 50 nm | Preprint healing / coherence |
| `corrections_moire` | 0.335 nm | 13.4 nm | Corrections local-curvature / moiré |

Do not mix silently.

## API

```python
from qvc.vacuum_freeze import (
    FrozenVacuumSpec,
    freeze,
    provisional_corrections_self_energy_spec,
    assert_preprint_0p244_not_reproduced,
)

spec = provisional_corrections_self_energy_spec()  # marked provisional
path = freeze(spec)  # writes qvc/locks/*.json
```

`FrozenVacuumSpec` fields: `formula_id`, `r`, `R`, `hbar_cs`,
`F_definition`, `notes`, plus `status` / `geometry_label`.

Default lockfile:
`qvc/locks/provisional_corrections_self_energy.json`
— **clearly marked provisional**.

## What must still be decided before quoting meV in papers

1. **Geometry:** preprint fat-torus vs Corrections moiré (or a third,
   derived scale). State the operator, BC, and measure.
2. **Formula:** corrected self-energy (recommended lock) vs any remaining
   Casimir-like variant — document the 2D measure argument.
3. **Scale \(\hbar c_s\):** confirm \(0.07\,\mathrm{meV\cdot\mu m}\) against
   the platform (phonon / optical) calibration; SI vs package units
   currently disagree at O(1)–O(10) unless \(c_s\) is retuned.
4. **Mode content of \(F\):** truncation `nmax`, continuum vs discrete
   WS modes, and whether Jacobi-\(\theta\) language is retired (it fails
   the \(F\approx 0.596\) claim).
5. **Mode multiplicity \(n\)** entering \(A_{\mathrm{eff}}=\sqrt{n}A_{\mathrm{single}}\)
   if vacuum amplitude feeds the gap equation.
6. **Status promotion:** change lock `status` from `provisional` to
   `frozen` only after (1)–(5) are fixed and verification reproduces
   the target energy to an agreed tolerance (e.g. &lt;1%).

Until then: **do not quote a vacuum meV as a final paper result.**

## Demo (no final claim)

```bash
.venv/bin/python scripts/freeze_vacuum_demo.py
.venv/bin/python scripts/run_vacuum_freeze_suite.py
```

Full suite docs: [`VACUUM_FREEZE_SUITE.md`](VACUUM_FREEZE_SUITE.md) (fat-torus diagnostics, rewritten self-consistency table, inverse freeze, P0/P1 brainstorm).

## Decision log (interim)

| Date | Decision |
|------|----------|
| 2026-08 | **Option C (two-track)** until platform \(\hbar c_s\) is fixed. Locked form always; both geometries in tables. |
| 2026-08 | Package Tier-C on fat torus (\(E\approx0.129\,\mathrm{meV}\)) exceeds fold-in bare bound \(\alpha_c\Delta_0\approx0.109\,\mathrm{meV}\) → no physical \(\Delta^*\) under \(\alpha=E/\Delta_0\) or \(\sqrt{n}\) bridge. |
| 2026-08 | Do **not** promote `status: frozen` until geometry + \(\hbar c_s\) + \(n\) + bridge are chosen together. |
| 2026-08-15 | **Option A\*** chosen: frozen-cavity fat torus + derived bridge \(\alpha=E_{\mathrm{vac}}/\Delta_0\), \(\ell=\ell(r)\). Form+bridge **proven**; geometry **frozen working**; \(E_{\mathrm{vac}}=0.1287\,\mathrm{meV}\) **computed** (0.244 withdrawn); \(\Delta^*\) at \(\lambda=1\) **fold-out (theorem)**. See [`BRIDGE_DERIVATION.md`](BRIDGE_DERIVATION.md). |

See also [`P0_P1_DERIVATIONS.md`](P0_P1_DERIVATIONS.md) and [`BRIDGE_DERIVATION.md`](BRIDGE_DERIVATION.md).
