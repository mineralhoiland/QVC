# QVC / TEVC Experimental Program

Lab-testable experiments, engineering developments, and a research timeline for Quantum Vacuum Catalyzation (QVC) on Topological Excitonic Vacuum Crystal (TEVC) platforms.

**Governing rule:** QVC claims *rate / barrier reshaping*, not net energy extraction from the vacuum. External drives (bias, laser, gate, chemical potential) supply energy.

**Confidence tiers**
| Tier | Meaning |
|------|---------|
| **A** | Direct condensed-matter test of a locked or nearly locked QVC identity |
| **B** | Testable but depends on provisional vacuum freeze / phenomenological catalyzon |
| **C** | Constructible engineering demonstrator; theory incomplete |
| **D** | Diagnostic analogue (IEC RF/optical spectra vs geometry); not fusion catalysis |

Numbers marked **(prov.)** inherit the provisional Casimir / \(\hbar c_s\) freeze and must be re-anchored after `qvc/vacuum_freeze.py` promotion.

---

## 0. Reference parameters (MATBG operating point)

| Symbol | Value | Role |
|--------|-------|------|
| \(\theta\) | \(\approx 1.1^\circ\) | Magic-angle twist |
| \(N\) | 5 (vary 3–9) | Layer count |
| \(\lambda_M\) | \(\sim 13\) nm | Moiré period / intrinsic cavity scale |
| \(d_\perp\) | 0.34 nm | Interlayer spacing |
| \(\Delta_{\mathrm{ex}}=2J\) | 0.600 meV | Bare exchange gap |
| \(\Delta_{\mathrm{QVC}}\) | 0.389 meV **(prov.)** | Self-consistent gap (~35% static loop) |
| \(E_{\mathrm{Cas}}\) | 0.244 meV **(prov.)** | Toroidal Casimir energy |
| \(A_{\mathrm{ZPF}}\times\xi\) | 0.296 **(prov.)** | Dimensionless ZPF |
| \(A_{\mathrm{eff}}/A_c\) | \(\approx 0.45\) **(prov.)** | Operating point on fold map |
| \(a_H\), \(\xi\), \(R\) | 8, 50, 50 nm | Healing / coherence / hopfion major radius |
| \(q^*\), \(\lambda^*\) | 0.056 nm\(^{-1}\), 113 nm | Phonon resonance |
| \(Q_H\) | \(1/N=1/5\) | Fractional Hopf (CS/linking) |
| \(\delta\sigma_{xy}\) | \(e^2/(2Nh)\approx 7.75\times10^{-6}\) S | Hall step |
| \(\omega_0=2\Delta/\hbar\) | \(\sim 10^{12}\) rad/s; IR \(\sim\) few THz–mid-IR depending on \(\Delta\) | Drive / resonance |
| \(T_{\mathrm{BEC}}\) | \(\sim 1\) K | Condensate scale |
| \(T^*\) | 285 mK **(prov.)** | Thermal–Schwinger crossover |

**Mechanism bookkeeping (do not conflate):**
1. Static \(A_{\mathrm{ZPF}}\) loop ≈ 35% gap reduction  
2. Catalyzon hybridization ≈ 52% (phenomenological; \(\lambda_{\min}=-g_0\))  
3. Coherent TQGL dynamics ≈ 92% at 1 ps (no Lindblad; upper bound)

---

## 1. Core lab experiments (Tier A–B)

### E1 — Hall staircase vs hopfion filling (Tier A)
**Concept:** Each hopfion nucleation changes Hall conductance by \(\delta\sigma_{xy}=e^2/(2Nh)=Q_H e^2/(2h)\).

| | |
|--|--|
| **Platform** | Dual-gated MATBG or TMD bilayer TEVC, \(N=3,5,7\); dilution fridge \(T\ll T_{\mathrm{BEC}}\) |
| **Controls** | Gate density \(n\), displacement field \(D\), \(B_\perp\), pressure / twist (proxy for \(A_{\mathrm{ZPF}}\)) |
| **Measure** | \(\sigma_{xy}(n,B)\), noise, simultaneous STM/SNOM if available |
| **Prediction** | Steps of height \(\propto 1/(2N)\); \(N\)-dependence of step height |
| **Falsify** | No \(1/N\) scaling of step height; steps only at integer Chern without fractional substructure when resolution allows |
| **Equipment** | Lock-in Hall bar, vector magnet, low-\(T\) cryostat |

### E2 — Layer-number gap scan (Tier A/B)
**Concept:** Moiré self-cavity \(\mathcal{F}_{\mathrm{moir\'e}}(N)\) boosts \(A_{\mathrm{ZPF}}^{\mathrm{cav}}\); gap ratio vs \(N\) is non-monotonic toward fold.

| | |
|--|--|
| **Platform** | Identical twist series \(N=3,5,7,9\) |
| **Measure** | Excitonic / pairing gap via tunneling spectroscopy or ARPES; \(T_c\) if superconducting |
| **Prediction** | \(\Delta(N)/\Delta_{\mathrm{ex}}\) decreases as \(\mathcal{F}(N)\) rises; strongest catalyzation at thicker stacks before other decoherence wins |
| **Falsify** | Gap ratio independent of \(N\) at fixed \(\theta\), or opposite trend |

### E3 — Fold criticality (not BKT until vortex RG) (Tier A)
**Concept:** Self-consistent gap map ends in a saddle-node fold; \(\delta-\delta_c\sim\pm\sqrt{\alpha_c-\alpha}\).

| | |
|--|--|
| **Tune** | Pressure, dielectric environment, cavity enhancement, or chemical potential as proxy for \(A_{\mathrm{ZPF}}/A_c\) |
| **Measure** | Gap vs tuning parameter; hysteresis / bistability near criticality |
| **Prediction** | Square-root collapse of physical branch; possible bistable window |
| **Falsify** | Essential singularity \(\exp(-c/\sqrt{t})\) *without* independent evidence of vortex unbinding (would require vortex RG upgrade) |

### E4 — Laser-driven triple-coincidence (Tier B)
**Concept:** Circular drive at \(\omega_L=\omega_0=2\Delta/\hbar\) triggers L1+L2+L3 *only* on resonance.

| Channel | Observable | Prediction |
|---------|------------|------------|
| **L1** | Stroboscopic / Faraday-wave imaging | Sidebands at \(\omega_L/2\); hopfion count \(\propto A^2\) above Mathieu threshold |
| **L2** | Hall transport during drive | Staircase steps lock to nucleation events |
| **L3** | IR photoluminescence | Superlinear pair production past instability tongue |

| | |
|--|--|
| **Setup** | Cryostat optical access; \(\sigma^\pm\) laser / THz–mid-IR source tunable through \(\omega_0\); simultaneous Hall + PL + camera |
| **Key knobs** | \(\omega_L\), amplitude \(A\), helicity, detuning \(\delta\omega\) |
| **Falsify** | Any channel present far off \(\omega_0\); L2 steps without L1 nucleation; no helicity dependence if theory requires circular drive |

**Engineering note:** Preprint figure quotes \(\sim 290\) THz for a particular \(\Delta\); retune \(\omega_L\) to measured \(\Delta\) on *that* device. Do not hard-code optical frequency.

### E5 — Thermal–Schwinger crossover \(T^*\) (Tier B)
**Concept:** Twin-channel pair creation: Arrhenius above \(T^*\), \(T\)-independent Schwinger plateau below.

| | |
|--|--|
| **Measure** | Pair-creation luminescence vs \(1/T\) from \(\sim 50\) mK to \(\sim 1\) K |
| **Prediction** | Kink at \(T^*\approx 285\) mK **(prov.)**; plateau \(\Gamma\sim 4\times10^8\) Hz per cell × \(N_{\mathrm{cells}}\) **(prov.)** |
| **Falsify** | Pure Arrhenius to base \(T\); plateau at wrong scale after vacuum freeze update |

### E6 — Catalyzon spectral function (Tier B)
**Concept:** Dark multiplet at \(\lambda_{\min}=-g_0\) (mult. \(N-1\)); bright mode at \((N-1)g_0\); resonance at \(q^*\).

| | |
|--|--|
| **Measure** | Inelastic light scattering / Raman / tunneling DOS; momentum-resolved EELS if available |
| **Prediction** | Soft mode / spectral weight at \(\omega_0\); dark-mode multiplicity signature under mode-selective probes |
| **Falsify** | Spectrum matching sparse adjacency (\(\lambda_{\min}=-5g_0\)) rather than democratic \(-g_0\); no resonance near \(q^*\approx 0.056\) nm\(^{-1}\) |

### E7 — Catalyzon-polariton in external cavity (Tier B)
**Concept:** Six-mode hybridization with split-ring / Fabry–Pérot cavity (Enkner-style vacuumronics).

| | |
|--|--|
| **Setup** | MATBG in tunable resonator; vary sample–resonator distance / \(\omega_c\) |
| **Measure** | Gap vs \(g_3/\omega_c\); avoided crossing in microwave/THz response |
| **Prediction** | Extra 10–30% gap reduction in ultrastrong regime; nonlinear (not pure Rabi) gap vs coupling |
| **Falsify** | Only linear vacuum-Rabi splitting with no gap softening |

### E8 — Devil’s staircase fractal Hall response (Tier B)
**Concept:** Rational hopfion fillings \(\nu_H=p/q\) give plateaux \(\sigma_{xy}\propto p/(2qN)\); Hausdorff \(d_H=3/4\) for dipolar \(\alpha=3\).

| | |
|--|--|
| **Measure** | Fine \(B\) or \(\mu_{\mathrm{eff}}\) sweeps (\(\delta B/B\sim10^{-4}\)) |
| **Prediction** | Sub-steps between main \(N\)-steps; plateau widths narrow with \(A_{\mathrm{ZPF}}\) |
| **Falsify** | Only integer steps at all resolutions; \(d_H\) inconsistent with dipolar model |

### E9 — \(C\to C'\) transition field reduction (Tier B)
**Concept:** Vacuum-catalyzed topological transition: \(B_c^{\mathrm{QVC}}\approx 0.72\,B_c^{(0)}\) at operating point **(prov.)**.

| | |
|--|--|
| **Measure** | Critical field for \(C=1\to 2\) (or valley Chern switch) vs cavity enhancement / \(N\) |
| **Prediction** | \(\sim 28\%\) reduction when vacuum coupling increased toward \(A_c\) |
| **Falsify** | \(B_c\) independent of cavity/\(N\) proxies for \(A_{\mathrm{ZPF}}\) |

### E10 — Coherence length / \(\xi_{\mathrm{QVC}}\) (Tier A/B)
**Concept:** Catalyzon-mediated pairing uses \(v_{\mathrm{cat}}\sim c_s\), not \(v_F\); \(\xi\sim 40\)–80 nm.

| | |
|--|--|
| **Measure** | Vortex lattice spacing, STM coherence maps, critical-current vs size |
| **Prediction** | \(\xi\) two orders below naive BCS \(v_F\) estimate |
| **Falsify** | \(\xi\sim\hbar v_F/k_BT_c\) with graphene \(v_F\) |

---

## 2. TEVC engineering developments (Tier A–C)

### T1 — Twist-angle and encapsulation stack
- hBN-encapsulated MATBG / twisted TMD; AFM-assisted tear-and-stack  
- Target \(\theta=1.05^\circ\)–\(1.15^\circ\); dual graphite gates  
- Metrology: STEM twist map, Raman, Landau fan  

### T2 — Multilayer TEVC library
- Alternating \(\pm\theta\) stacks for \(N=3,5,7,9\)  
- Controlled \(d_\perp\) via pressure cells or intercalants  
- Deliverable: gap vs \(N\) dataset for E2  

### T3 — On-chip Hall + optical hybrid package
- Transparent top gate or photonic crystal window  
- Side-contact Hall bar + fiber / free-space ports  
- Enables E4 triple-coincidence without sample swaps  

### T4 — Integrated split-ring / superconducting resonator
- Enkner-like vacuumronics on the same chip as TEVC  
- Tunable \(g_c\) for E7  

### T5 — Strain / pressure \(A_{\mathrm{ZPF}}\) knob
- Diamond anvil or MEMS pressure cell at mK  
- Primary fold-map tuner for E3  

---

## 3. Catalyzon generation protocols (Tier B–C)

### C1 — Passive (equilibrium) catalyzon
Cool below \(T_{\mathrm{BEC}}\); verify gap \(\Delta_{\mathrm{QVC}}\) and spectral weight at \(\omega_0\) (E6). No drive.

### C2 — Parametric (Mathieu) nucleation
Modulate \(\Delta(t)=\Delta_0+\delta\Delta\cos(\omega_{\mathrm{pump}}t)\) or optical \(A(t)\) at \(\omega_L\approx\omega_0\).  
Threshold \(A_{\mathrm{thr}}\); subharmonic \(\omega_L/2\); hopfion count vs \(A^2\).

### C3 — Cavity-assisted catalyzon-polariton
Park \(\omega_c\) on catalyzon; increase \(g_3\) into ultrastrong; track lowest eigenvalue (E7).

### C4 — Spatially seeded hopfion lattice
Use structured light / AFM tip / local gates to seed rings of radius \(R\sim 50\) nm; measure Hall staircase as lattice filling increases (E1/E8).

### C5 — Lindblad-aware pulsed generation
Repeat C2 with measured \(\tau_{\mathrm{dec}}\); accept Stage-3 “92%” only if coherence survives decoherence — research deliverable, not a guaranteed product metric.

---

## 4. Laser-driven and photonic setups (constructible)

| Setup | Purpose | Build notes |
|-------|---------|-------------|
| **LD-1 Free-space cryo optics** | E4 L1+L3 | Dilution/optical cryostat, \(\sigma^\pm\) control, mid-IR/THz OPO or QCL tuned to measured \(\Delta\) |
| **LD-2 Fiber-coupled PL** | E5 \(T^*\) | High-NA collection; calibrated count rate per moiré cell |
| **LD-3 On-chip photonic crystal** | Local \(A_{\mathrm{ZPF}}\) boost | Defect mode at \(\lambda^*\) or \(\omega_0\); couples to E3 fold |
| **LD-4 Stroboscopic Faraday imaging** | L1 | Pulsed probe at \(\omega_L/2\); ns–ps resolution optional |
| **LD-5 Dual-tone pump–probe** | Separate catalyzon vs Schwinger | Detune one tone off \(\omega_0\); expect L3 residual if Schwinger-dominated |

---

## 5. Fusor / plasma analogues (Tier D — diagnostic, not catalysis)

These are **not** claims that QVC produces fusion energy or extracts vacuum energy. The preprint’s fusion remark remains L3 / conjectural. An IEC chamber already under construction is allowed as a **spectral diagnostic** (geometry-dependent RF), not as a MATBG confirmation path.

Shop-floor protocols (Tormach PCNC 440 envelope, electrode CAD, F-S1–F-S3): `QVCCursor/tex/QVC_Garage_Lab_Manual.tex`.

### F1 — Barrier-catalysis analogue bench (cold, CM)
Map Mexican-hat barrier reduction (E3/E6) onto a *classical* tunable barrier (Josephson junction array / SQUID).  
**Goal:** Measure rate enhancement vs barrier height with a “catalytic” mode hybrid — pedagogical, not nuclear.

### F2 — IEC chamber as RF / optical diagnostic (F-S series)
Use the existing IEC vessel as a driven cavity. Swap inner electrodes (spherical-approximation rings vs fat torus, both millable on a PCNC 440) at matched voltage, current, and fill pressure. Headline observable is the **RF spectrum** (VNA / pickup loop), not neutrons and not Balmer centroids as \(\Delta_{\mathrm{QVC}}\). Optical lines are covariates for plasma density/Stark. A hovering SRR at a viewport is an Enkner-geometry analogue (F-S3).  
**Goal:** Falsify whether fat-torus mode-sum shifts are visible in this chamber. A null does not falsify MATBG QVC. **No** fusion-rate catalysis claim.

### F3 — Plasma–cavity RF coupling study
Couple the G4 6061 cavity / viewport SRR to the IEC RF envelope; compare empty-cavity \(\Delta f(d)\) to plasma-on pull.  
**Goal:** Cross-domain scaling for “drive + boundary geometry” — still Tier D.

### F4 — Do-not-build / do-not-claim list
- Devices claiming “vacuum energy → net power”
- Fusors marketed as QVC-powered reactors
- Neutron yield, deuterium recipes, or grid retuning as QVC signatures
- Identifying a visible-line shift with \(\Delta_{\mathrm{QVC}}\approx 0.4\,\mathrm{meV}\) (THz, not optical)
- Any write-up that treats F-S as a substitute for E1–E5 condensed-matter falsification  

---

## 6. Measurement matrix (quick reference)

| Parameter | Method | Experiments |
|-----------|--------|-------------|
| \(\sigma_{xy}\), \(\delta\sigma_{xy}\) | Lock-in Hall | E1, E4-L2, E8 |
| \(\Delta\), \(\Delta_{\mathrm{cat}}\) | Tunneling / ARPES / optical | E2, E3, E6, E7 |
| \(T_c\), \(\xi\) | Transport / STM | E10 |
| \(A_{\mathrm{ZPF}}/A_c\) proxy | Pressure, \(N\), cavity \(Q\) | E2, E3, E9 |
| Hopfion count / Faraday | Imaging | E4-L1, C2, C4 |
| Pair PL rate \(\Gamma(T)\) | Calibrated photon counting | E4-L3, E5 |
| Mode spectrum | Raman / EELS / microwave | E6, E7 |
| \(B_c\) for Chern switch | Magnetotransport | E9 |

---

## 7. Development & research timeline

### Phase 0 — Theory freeze & facility prep (Months 0–6)
- [ ] Vacuum freeze decision; update all **(prov.)** numbers  
- [ ] Commit fold critical table from `qvc/gap_fold.py`  
- [ ] \(G_{ij}\) from overlaps (accept/reject democratic Spec)  
- [ ] Secure optical dilution fridge + MATBG fab partner  
- [ ] Write E1–E5 SOPs and falsification checklists  

### Phase 1 — TEVC baseline library (Months 4–18)
- [ ] T1–T2: fabricate \(N=3,5,7\) devices  
- [ ] E1 Hall staircase; E2 gap vs \(N\); E10 \(\xi\)  
- [ ] Publish baseline (even null results)  

### Phase 2 — Dynamics & twin channel (Months 12–30)
- [ ] T3 hybrid package; LD-1/LD-2  
- [ ] E4 triple-coincidence; E5 \(T^*\) kink  
- [ ] E3 fold criticality with pressure/cavity knob  
- [ ] Lindblad TQGL re-sim with measured \(\tau_{\mathrm{dec}}\)  

### Phase 3 — Catalyzon & cavity engineering (Months 24–42)
- [ ] T4 resonator integration; E7 polariton  
- [ ] C2–C4 generation protocols; mode spectroscopy E6  
- [ ] E8 devil’s staircase high-resolution scan  
- [ ] E9 \(B_c\) reduction test  

### Phase 4 — Systems & speculative analogues (parallel; does not gate MATBG)
- [ ] On-chip photonic \(A_{\mathrm{ZPF}}\) control (LD-3)
- [ ] F1 classical barrier-catalysis demonstrator
- [ ] F-S1–F-S3 IEC spectral diagnostics (RF headline; optical covariate) — shop work, not a Phase-2 gate
- [ ] Separate TDVT spacetime program (do not gate MATBG papers on it)  

### Parallel theory track (continuous)
- Vortex RG → optional BKT language upgrade  
- Torsion → \(G_{ij}\) derivation  
- Retarded \(V(q,\omega)\) replacing \(\eta_{\mathrm{ret}}=3.1\) static caveat  

---

## 8. Success criteria (program-level)

| Gate | Requirement to proceed |
|------|------------------------|
| **G0→G1** | Vacuum freeze + ≥1 working TEVC Hall bar at \(T<1\) K |
| **G1→G2** | E1 or E2 positive *or* well-documented null with calibrated sensitivity |
| **G2→G3** | E4: ≥2 of 3 coincidence channels on resonance, off-resonance null |
| **G3→G4** | E5 kink *or* E7 nonlinear gap softening reproduced in second lab |
| **Fusor spectral diagnostics (F-S)** | Allowed in parallel as Tier D analogue; RF geometry scan only |
| **Fusor as QVC fusion / yield** | Never a program gate; L3 remains quarantined |

---

## 9. Related docs
- `docs/NEXT_STEPS_RIGOR.md` — math closure  
- `docs/VACUUM_FREEZE_SPEC.md` — provisional energy  
- `docs/CLOSED_FORM_LOCKS.md` — Spec(G), fold, \(Q_H\)  
- `docs/VORTEX_RG_FRAMEWORK.md` — BKT contingency  
- Preprint §Experimental Signatures, §Cavity QVC, §Schwinger Twin Channel  

```bash
cd /Users/mineralhoiland/Code/QVCCompute
MPLCONFIGDIR=/tmp/mpl .venv/bin/python scripts/run_verification.py
MPLCONFIGDIR=/tmp/mpl .venv/bin/python scripts/freeze_vacuum_demo.py
```
