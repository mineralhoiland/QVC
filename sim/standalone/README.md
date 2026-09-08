# Standalone simulation suites

## TDVT Hybrid Showcase (`tdvt_hybrid_showcase.html`)

**Open the file in a browser — no localhost.** Classic UMD Three.js in `vendor-umd/` (same `file://` path as the processor chip). Default **524,288** additive GPU points (slider to 1,048,576). Planetarium cabinet UI (Cormorant / Figtree / ivory–gold), not the neon HUD and not the clay chip.

Eight particle families share one vertex shader. Gallery I is the full orrery; II–IX isolate a mechanic. **Overlay chips** add/remove families; **Stacked** vs **Cabinet** control spatial layout. **Listening** is a 16-step Φ₀ sequencer tuned to locked parameters:

- **f₀** scales with **g₀** (democratic coupling); **4f₀** bright partial is the Spec(\(G\)) antenna arpeggio
- **Register** arpeggio: just-intonation cluster \([4,1,5/4,3/2,9/8,\ldots]\) f₀; island toms on pentamer bits
- **SSH** arpeggio: eight scale degrees walk with soliton; snare on soliton site; hat every 16th
- **Hopf** dyad: θ beat at **π/(2N)**; **Horizon** dual arpeggio with Doppler detune from **v_s/κ**
- **Cavity** arpeggio gated by **A\***; **E_vac/g₀** ratio as coincidence bell (not a detection claim)
- **Helix** staccato **Φ₀** plucks + kick on quarters; **Torsion** saw stabs ∝ **α_T Q**
- **Kick** = clock; **clap** on active coupling identities; **Arpeggio** slider scales note/perc level

Locked couplings when overlays combine:

- Helix ∩ SSH — soliton snaps to 16ths
- SSH ∩ Register — soliton writes island \(i \bmod 5\)
- Horizon ∩ Cavity — Doppler stretches cavity partials
- Cavity ∩ Register — DCE \(A^\star\) gates the bright \(4f_0\) bit
- Hopf ∩ Lattice — linking delay beats BCC tubes
- Register ∩ Torsion — Helmert word on the torsion stream

L3 ontology stays off.

| Gallery | Particles | Physics |
|---------|-----------|---------|
| I Orrery | all eight families | spatial cabinet of the processor branch + analogue spacetime |
| II Hopf S³ | stereographic \(\pi(z_1,z_2)\) | linking / \(Q_H=1/N\) |
| III Hopfion lattice | BCC tubes + filaments | 5 Goldstones; FeGe is a later module |
| IV Acoustic horizon | Unruh–Visser \(v_s=\kappa/(2\pi r)\) | \(T^\star=118.7\) mK Schwinger, **not** Hawking duality |
| V Josephson cavity | TE standing wave | DCE self-limit \(A^\star/A_c^{\mathrm{fold}}=0.9117\); \(E_{\mathrm{vac}}/h\approx 31.1\) GHz is a **coincidence** with 32 GHz bolometers |
| VI SSH solitons | dimerized chain | adiabatic shuttle; dump is the falsifier; \(\theta=\pi/(2N)\) abelian |
| VII Dark register | pentamer \(N=5\) | \(\lambda_{\max}=4g_0\), four-fold Helmert dark; CS \(N\) moves \(\theta\)/Hall only |
| VIII Flux helix | WC helix, \(R=50\) nm | \(\Phi_0\) beads; TA phonons \(\neq\) gravitons |
| IX Torsion analogue | \(S^r\sim\alpha_T Q/r^2\) | H-TEGR1 open; Einstein not claimed |

Double-click `tdvt_hybrid_showcase.html`. Click **Enable listening** once. Keys: `1–9` galleries, `H` hide plaque, `M` mute, space pause.

Locked: \((r,R)=(8,50)\) nm, \(F=0.58068\), \(E_{\mathrm{vac}}=0.1287\) meV, \(g_0=0.0982\) meV, \(\lambda^\star=0.7628\), \(\lambda_{\mathrm{fold}}=0.8475\), \(A^\star/A_c=0.9117\).

---

## TDVT Astrophysical Cosmos (`tdvt_astro_cosmos.html`)

Three-tab **cinematic** GPU environment (not a lab overlay). Full-viewport Painlevé-style Schwarzschild geodesics, volumetric Keplerian disk, dipole magnetosphere, and SGR pair fireball. Instruments are a side console; hide it for a clean view.

| Tab | Object | Physics |
|-----|--------|---------|
| A Magnetar / Superfluid | 1E 1547, 4U 0142+61 | Hopfion/vortex core, Adler birefringence, E_ind/E_c |
| B Accretion / Horizon | Sgr A*, M87*, Schwarzschild | Null geodesics, Kerr ISCO, photon ring, T_AH overlay |
| C Stellar / SGR Flare | SGR 1806−20 | Reconnection fireball, pair opacity, ejecta ~c/3 |

```bash
cd sim-suites/standalone
python3 -m http.server 8765
# http://localhost:8765/tdvt_astro_cosmos.html
```

Presets: 1E 1547.0−5408, 4U 0142+61, SGR 1806−20, Sgr A*, M87*, Schwarzschild, Sonic overlay, Millisecond gap.

---

## TDVT Spacetime Evolution (`tdvt_spacetime.html`)

Paper-faithful **time-evolving** GPU suite for Unified Hopfion Spacetime, Revision B. Four locked simulations only. S³ Clifford rotation visualizes the Hopf total space (not a tesseract ontology). Fractal $d_H$ is omitted. Holography is a 2+1D boundary overlay on the lattice tab.

| Tab | Paper | GPU field |
|-----|--------|-----------|
| Hopf S³ | §2 | Exact π(z₁,z₂), 4π fiber, linking |
| Acoustic BH | §4 | v_s=κ/(2πr), Unruh–Visser horizon |
| TEGR Torsion | §5 | S^r∼α_T Q/r², TA roton |
| Lattice / Boundary | §3 | BCC + 5 Goldstones + Q_H=1/N |

```bash
cd sim-suites/standalone
python3 -m http.server 8765
# http://localhost:8765/tdvt_spacetime.html
```

Default 524,288 points (slider to 1,048,576). Three.js additive GPU vertices + Babylon.js inset. Presets: Fiber lock, Linked Q=2, Sonic core, Torsion steer, BCC crystal, Hypersphere.

---

## TDVT Particle Lab (`tdvt_particle_lab.html`) — GPGPU playground

GPU computational + visual suite for the hopfion-crystal TDVT programme (Revision B).

- **Three.js** additive GPU points (default 524,288; up to 1,048,576)
- **GPUComputationRenderer** 1024² ping-pong position/velocity (field forces per world)
- **Babylon.js** inset of hull / horizon ring / torsion axis
- **Plotly** live analysis + sparkline time series
- Orbit camera, presets, Kick Field, Reseed GPGPU
- Claim chips: lock / prov / conj / analogue (`T_H` matching ≠ duality; Einstein not claimed)

| Tab | Physics |
|-----|---------|
| Hopfion Crystal | BCC tubes, VK binding, 5 Goldstones |
| Acoustic Horizon | Unruh–Visser, critical slowdown, T_H matching |
| Hypersphere / Vortices | S³→S² flow, linking, H-VTOR1 |
| TEGR Torsion | Weitzenböck T∝∇q_H, H-TEGR1 open |
| Quantum Geometry | Spec(G), fold, Hall step, D_s language |
| Lattice / Audit | CS Q_H=1/N, Bessel F, freeze badge |

```bash
cd sim-suites/standalone
python3 -m http.server 8765
# http://localhost:8765/tdvt_particle_lab.html
```

CDN required (importmap): three@0.175, Babylon, Plotly. Serve over http (not `file://`) so ES modules load.

---

## TABD · QVC Drive (`tabd_qvc_drive.html`)


Sci-fi GPU particle suite for the **TDVT Acoustic Bubble Drive** and **QVC** catalysis papers. Three.js renders up to 1,048,576 GPU points; Babylon.js drives the hull inset. Tabs: TABD bubble, hopfions, torsion, acoustic horizon, QVC catalysis, TDVT lattice. Locked numerics: ν₀ = 0.188 THz, Γ_Sch = 4.73×10⁸ Hz, fold δ_c = 0.182, Hall step e²/(2Nh). L3 GR-FTL is labeled open.

```bash
cd sim-suites/standalone
python3 -m http.server 8765
# http://localhost:8765/tabd_qvc_drive.html
```

Orbit drag · scroll zoom · sliders · presets · Nucleate / Reverse ŝ · live plots.

---

## TDVT Processor Chip (`tdvt_processor_chip.html`)

Open the file in a browser — **no localhost**. Classic UMD Three.js lives in `vendor-umd/` (ES modules are what failed on `file://`). Clay/copper instrument UI, not the neon HUD. Web Audio maps Spec(G), the SSH grain, DCE limiter, Q_H subharmonic, and hopfion beating into a live FFT.

| View | Analog |
|------|--------|
| Register | Bright partial at 4f₀, four darks at f₀ (stereo) |
| Shuttle | Bandpass grain, panned with the soliton |
| Gate | Limited subharmonic ∝ min(A, A*) |
| Memory | Linked two-tone beat if FeGe is on |
| Clock | Φ₀ clicks; fold as low-pass cutoff ∝ Δ |

Double-click `tdvt_processor_chip.html`. Click **Enable listening** once (browser policy). Keep `vendor-umd/three.min.js` and `OrbitControls.js` beside the HTML.

---

# TDVT Atlas (standalone)

Paper-structured interactive verification suite for **Unified Hopfion Spacetime, Revision B**.

**Not** a clone of the earlier Claude HTML labs. Those were used only as inspiration (GPU/interaction energy). This atlas follows the manuscript section order and locked math.

**File:** [`tdvt_atlas.html`](./tdvt_atlas.html)

## Open

```bash
cd sim-suites/standalone
python3 -m http.server 8765
# http://localhost:8765/tdvt_atlas.html
```

Or open the HTML file directly in a browser.

## Chapters (paper §§)

| Tab | Content |
|-----|---------|
| §2 Hopf | Exact π(z₁,z₂), linked fibers, Berry A/F |
| §2 VK/FN | Toroidal hopfion + E≥C\|Q\|^{3/4} |
| §3 Lattice | BCC + five Goldstone channels |
| §4 Metric | Unruh–Visser g_eff, horizon, g_tt(r) |
| §5 Torsion | Hopf-sourced axial torsion quiver |
| §6 Vacuum | F=ΣK₀(nx), Tier-C E_vac, freeze badge |
| §7 Fold | Saddle-node vs dashed BKT (rejected) |
| §8 Spec(G) | Democratic G heatmap + λ spectrum |

## Epistemic chips

- **lock** — Revision B / QVC lock  
- **prov** — provisional meV until vacuum freeze  
- **conj** — H-TEGR1 open  
- **contingent** — BKT only after vortex RG  
