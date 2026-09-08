# TDVT 3D

Browser visualizer for the Topological Dynamic Vacuum (Hopfion / Whitehead / analogue geometry). Live WebGL scenes sit on top of the same Option A\* laboratory locks used by the Python workbench. The HUD labels every number as a **lock**, a **live compute**, a **diagnostic**, or **artistic** kinematics.

## How to run

From the `web/` directory (ES modules need HTTP, not `file://`):

```bash
cd /Users/mineralhoiland/Code/QVCCompute/web && python3 -m http.server 8765
```

Then open:

- Play mode: http://localhost:8765/tdvt3d/
- Self-test (`?selftest=1`): http://localhost:8765/tdvt3d/?selftest=1

Whitehead 48³ in the self-test can take several seconds. Leave the tab visible; `document.hidden` skips the render loop.

## Tabs

| Tab | Layers | What you see |
|-----|--------|----------------|
| Hopf fibration | L2 | Hopf bundle π: S³ → S². Colour = base point. Two fibres: live Gauss linking → ±1. |
| Whitehead core | L2 | Stereographic Hopf ansatz. Whitehead Q = (1/16π²)∫A·B, distinct from Q_H = 1/N. |
| BCC lattice | L2 | Schematic Im-3m hopfion crystal, a = ξ = 50 nm. N slider updates lens-space identities only. |
| Acoustic horizon | L2 | Unruh–Visser g_tt = −(c_s² − v²). Gold ring = sonic surface when an exterior horizon exists. |
| Teleparallel torsion | L2, L3 | Torsion glyphs ∝ ρ_H. H-TOR1 is a parked experiment and never feeds G_ij. |
| Solitons & spectra | L1, L2 | Fold / Δ\* well, VK bound, Spec(G) stems. SYK is an analogy, not an identification. |
| S³ flow | L2, L3 | Hopf U(1) fibre flow cartoon. H-HS1 (torsion ∝ mean curvature) is a conjecture. |

## Computed vs artistic

- **Computed / live:** Gauss linking, Whitehead Q and FN energy on the displayed grid, CS-flux → Q_H, acoustic ρ_h / κ / T_H, Spec(G) eigenvalues, fold α_c = δ_c.
- **Locks (never re-fitted):** Option A\* table below. HUD rows with kind `lock` are these numbers.
- **Diagnostic:** analogue Hawking T_H vs locked T\* (H-HAWK1). T_H is Unruh–Visser kinematics, not T\*.
- **Artistic / schematic:** fibre-phase breathing, lattice site glyphs, Goldstone arrows, S³ flow, H-TOR1 overlay, H-HS1 cartoon. These do not replace a dynamical solver.

Play-mode verification dots stay empty. `?selftest=1` fills them from the live JS suite (`js/math/verify.js`) against `data/whitehead_reference.json` and `data/analogue_reference.json`.

## Layers

| Layer | Meaning |
|-------|---------|
| **L1** | Laboratory lock (QVC / TEVC). Frozen Option A\* numbers. |
| **L2** | Analogue geometry (TDVT). Diagnosed from the condensate; **not** Einstein gravity. |
| **L3** | Spacetime ontology. Schematic / parked; not locked. |

## Option A\* locks

Mirrored in `js/math/locks.js` from `qvc/workbench/locks.py` (`data/locks.json`). Displaying a withdrawn number throws.

| Quantity | Value | Do not use |
|----------|-------|------------|
| Geometry (r, R) | 8 nm, 50 nm fat torus | — |
| ħ c_s | 0.07 meV·μm | — |
| E_vac | **0.1287 meV** | 0.244 meV (withdrawn Casimir) |
| Δ₀ | 0.6 meV | — |
| α_c = δ_c | 0.1818 | 92% static reduction |
| Δ\* | 0.2324 meV | 0.389 meV, 0.032 meV |
| T\* | **118.7 mK** | 285–290 mK |
| Q_H = 1/N | 0.2 (N = 5) | — |
| θ_anyon | π/N | Paper C π/(2N) is a different CS level |
| Hall step | e²/(2Nh) = 0.1 | — |
| Spec(G) | λ_max = (N−1)g₀\*, λ_min = −g₀\* ×(N−1) | −4g₀ / −5g₀ |

T_H at a chosen drain is an analogue diagnostic, not T\*. H-EIN1: TA vs TT overlap is not graviton-like.
