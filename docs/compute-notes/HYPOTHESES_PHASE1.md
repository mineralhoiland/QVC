# Phase-1 Hypotheses Register (Vortex RG · Torsion · TDVT–QVC Sharing)

**Purpose:** Structured, falsifiable conjectures for the Phase-1 ideation package.  
**Tone:** Ambitious but honest — every entry is **Status: conjecture** unless noted.  
**Companions:** `VORTEX_RG_FRAMEWORK.md`, `TORSION_TDVT_HOPFION.md`  
**Related locks:** `CLOSED_FORM_LOCKS.md`, `VACUUM_FREEZE_SPEC.md`

---

## How to read this file

Each hypothesis uses:

- **Hypothesis statement** — precise claim  
- **Status** — conjecture (unless upgraded later)  
- **Depends on** — prior locks / data / computations  
- **Acceptance / falsification test** — what would upgrade or kill it  
- **Priority** — P0 (blocks preprint language), P1 (next paper), P2 (exploratory)

---

## H-VRG1 — Vortex RG converts fold → BKT

**Hypothesis statement:**  
Including a vortex fugacity sector in the RG of the moiré excitonic phase \(\theta=\arg\Psi\) converts the locked saddle-node fold of \(\delta=1+\alpha(\ln\delta+\ell)\) into a BKT essential singularity \(\Delta\sim\exp(-C/\sqrt{A_c-A})\) with continuous gap vanishing.

**Status:** conjecture

**Depends on:** Locked fold (V1)–(V3); microscopic \(\kappa(\Delta)\); KT flow (V7)–(V8)

**Acceptance / falsification test:**  
- **Accept:** C1–C7 in `VORTEX_RG_FRAMEWORK.md` pass; evidence ratio favors essential singularity over \(\sqrt{\,}\) for \(\Delta(A)\) and \(\eta\to 1/4\) in phase correlations.  
- **Falsify:** Physical drive path hits the amplitude fold at finite \(\Delta_c\) before \(\kappa=\kappa_{\mathrm{BKT}}\); or phonon damping destroys the KT separatrix with no residual (V9) scaling.

**Priority:** P0

---

## H-TOR1 — Torsion overlaps produce democratic Spec\((G)\)

**Hypothesis statement:**  
Evaluating \(G_{ij}=\int_{\mathrm{WS}} T^\mu{}_{\nu\rho}\Phi_i^{*\nu}\Phi_j^\rho\,d^3x\) with \(T\propto\varepsilon\,\partial q_H\) on an \(O_h\)-symmetric hopfion/WS cell yields (up to \(O_h\)-allowed perturbations) the democratic spectrum \(\lambda_{\max}=(N-1)g_0\), \(\lambda_{\min}=-g_0\) (mult.\ \(N-1\)) **without** inserting \(G=g_0(J-I)\) by hand.

**Status:** conjecture

**Depends on:** Mode profiles \(\Phi_i\); \(\kappa_{\mathrm{tor}}\); freeze of WS geometry (`VACUUM_FREEZE_SPEC.md`)

**Acceptance / falsification test:**  
- **Accept:** Numerical quadrature produces Hermitian traceless \(G\) with \(G_{ii}=0\), \(\lambda_{\min}<0\), and eigenvector structure matching the Helmert dark basis within a stated tolerance (e.g. relative eigenvalue error \(<10\%\) after scale match of \(g_0\)).  
- **Falsify:** Generic \(\nabla q_H\) produces \(\lambda_{\min}\ge 0\), large diagonal \(G_{ii}\), or spectrum incompatible with multiplicity \(N-1\) dark subspace.

**Priority:** P0

---

## H-SHARE1 — \(Q_H=1/N\) is the correct shared topological invariant

**Hypothesis statement:**  
The physically shared topological invariant between QVC fractionalization and hopfion/\(L(N,1)\) spacetime ideation is \(Q_H=1/N\) descending from \(\mathrm{CS}[A_1]=1/N\) and \(\mathrm{lk}=1/N\), not a putative fundamental 2-cycle of \(L(N,1)\).

**Status:** conjecture *(CS/lk identities themselves are locked; the “correct shared invariant” interpretation is the conjecture)*

**Depends on:** Locked CS/linking identities; N-layer \(\mathbb{Z}_N\) identification

**Acceptance / falsification test:**  
- **Accept:** Anyonic/Hall sector of TEVC shows fractionalization consistent with \(1/N\) (or abelian \(\theta=\pi/N\)) and no independent H\(_2\)-cycle observable is required.  
- **Falsify:** Experiment or exact solvable limit requires a different primitive invariant (e.g. integer Hopf only, or filling-dependent fraction unrelated to \(N\)-layers).

**Priority:** P0

---

## H-HAWK1 — \(T_H=T^*\) is matching, not duality

**Hypothesis statement:**  
Equality of acoustic Hawking temperature \(T_H=\hbar c_s/(\pi k_B R_{\mathrm{WS}})\) with Schwinger crossover \(T^*\) is a **geometric matching condition** that fixes \(R_{\mathrm{WS}}\), not a derived bulk/boundary duality from a single spectrum, unless and until both are obtained from one eigenvalue problem.

**Status:** conjecture *(as duality)* — the matching interpretation is the **disciplined default**

**Depends on:** Analogue gravity formula; independent determination of \(T^*\) and \(R_{\mathrm{WS}}\)

**Acceptance / falsification test:**  
- **Accept as duality:** Derive \(T^*\) and \(T_H\) from one spectral object (same modes, same measure) with no free \(R_{\mathrm{WS}}\) insert.  
- **Falsify duality (retain matching):** Measured \(R_{\mathrm{WS}}\) (or coherence geometry) disagrees with \(R=\hbar c_s/(\pi k_B T^*)\) by \(\gg\) systematic error, or \(T^*\) mechanism proven unrelated to acoustic horizons.  
- **Language rule:** Until acceptance, write “matched” / “identified by matching,” never “dual.”

**Priority:** P0

---

## H-VRG2 — Phonon damping selects overdamped KT crossover

**Hypothesis statement:**  
Gao–Khalaf-type electron–phonon coupling in MATBG/TEVC places vortices in the overdamped regime (\(\hbar/\tau_v \gtrsim \pi\kappa\)), converting textbook BKT into a rounded essential-singularity-like crossover rather than a sharp universal jump.

**Status:** conjecture

**Depends on:** H-VRG1 scaffolding; phonon spectral function; vortex mobility estimate

**Acceptance / falsification test:**  
- **Accept:** Computed \(\eta_v\) implies overdamped dynamics; stiffness jump washed out while \(\Delta(A)\) still closer to (V9) than to MF \(\sqrt{\,}\).  
- **Falsify:** Underdamped vortices with clear universal jump in kinetic inductance at the gap collapse.

**Priority:** P1

---

## H-CAT1 — Dark-mode catalyzon is dynamically protected relative fluctuation

**Hypothesis statement:**  
The Helmert dark subspace of democratic \(G\) (sum-zero modes with shift \(\Delta_{\mathrm{cat}}=\Delta_{\mathrm{ex}}-g_0\)) remains long-lived under uniform drives and dominates vacuum-catalyzed gap reduction relative to the bright eigenvalue \(\lambda_{\max}\).

**Status:** conjecture

**Depends on:** Locked Spec\((G)\); TQGL/Lindblad or quench numerics

**Acceptance / falsification test:**  
- **Accept:** Simulated or measured response shows dark-sector weight dominating gap renormalization; bright mode radiates/decays faster; coherence time hierarchy \(\tau_{\mathrm{dark}}>\tau_{\mathrm{bright}}\).  
- **Falsify:** Dark modes decohere as fast as bright; gap shift tracks \(\lambda_{\max}\) rather than \(\lambda_{\min}\).

**Priority:** P1

---

## H-VTOR1 — Vortex worldlines source measurable torsion-like mode mixing

**Hypothesis statement:**  
Increasing vortex density (by drive or temperature) increases off-diagonal mode mixing in the effective \(G_{ij}\) (or measured hybridization splittings) in a manner consistent with defect-sourced torsion/spin density in the Einstein–Cartan analogue — testable in CM without GR ontology.

**Status:** conjecture

**Depends on:** Vortex densitometry; mode spectroscopy; H-TOR1 kernel

**Acceptance / falsification test:**  
- **Accept:** Controlled vortex injection correlates with growth of off-diagonal \(G_{ij}\) proxies at fixed average density.  
- **Falsify:** No correlation after controlling for amplitude suppression and heating.

**Priority:** P2

---

## H-FN1 — Five Goldstones are the correct IR mode content for \(N=5\) QVC

**Hypothesis statement:**  
The IR theory relevant to catalyzon hybridization has five soft modes (3 acoustic + 2 fiber/orientational), matching \(N=5\) democratic channels, as the symmetry limit of an \(O_h\) hopfion crystal projected to the moiré boundary.

**Status:** conjecture

**Depends on:** Broken-generator count; pinning gaps; boundary projection dictionary

**Acceptance / falsification test:**  
- **Accept:** Spectral function shows five soft poles (or five dominant hybridization channels) before pinning gaps; democratic \(5\times 5\) \(G\) fits data better than \(N\neq 5\).  
- **Falsify:** Only 2–3 soft modes appear; or best-fit channel count \(\neq 5\) at operating point.

**Priority:** P1

---

## H-GEO1 — Vacuum geometry freeze selects Corrections scale, not preprint fat-torus

**Hypothesis statement:**  
Once a single vacuum spectral problem is frozen (`VACUUM_FREEZE_SPEC.md`), the self-consistent \((r,R)\) aligning Casimir/ZPF amplitude with the physical fold branch is the Corrections local-curvature/moiré scale, not the preprint fat-torus \((r,R)\approx(8,50)\,\mathrm{nm}\).

**Status:** conjecture

**Depends on:** Frozen operator, BC, measure; locked vacuum formula \(E_{\mathrm{vac}}\propto(\hbar c_s/r)F\)

**Acceptance / falsification test:**  
- **Accept:** Unique geometry reproduces target energy to \(<1\%\) and consistent \(A_{\mathrm{ZPF}}\) on the physical \(\delta\) branch.  
- **Falsify:** Fat-torus wins the same spectral problem; or no geometry meets \(<1\%\) without retuning forbidden parameters.

**Priority:** P1

---

## Priority board (Phase-1)

| ID | One-line | Priority |
|----|----------|----------|
| H-VRG1 | Vortex RG fold→BKT | P0 |
| H-TOR1 | Torsion→democratic Spec\((G)\) | P0 |
| H-SHARE1 | Shared invariant \(Q_H=1/N\) | P0 |
| H-HAWK1 | \(T_H=T^*\) matching ≠ duality | P0 |
| H-VRG2 | Phonon-overdamped KT | P1 |
| H-CAT1 | Dark catalyzon protection | P1 |
| H-FN1 | Five Goldstone channels | P1 |
| H-GEO1 | Geometry freeze → Corrections scale | P1 |
| H-VTOR1 | Vortex→torsion mixing | P2 |

---

## Upgrade / demote rules

1. **No silent upgrades.** Moving conjecture → locked requires a named computation or experiment logged in `output/` or `CLOSED_FORM_LOCKS.md`.  
2. **Preprint language:** P0 conjectures may appear only with explicit “conjecture/program” wording.  
3. **L3 quarantine:** Spacetime ontology, \(\Lambda\), circular time never enter L1 abstracts.  
4. **Anti-overclaim:** Failure of H-VRG1 does **not** invalidate the locked fold; failure of H-TOR1 does **not** invalidate democratic phenomenology of \(G\).
