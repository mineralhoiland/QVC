# Vortex RG Framework: Fold → BKT Essential Singularity

**Status:** Framework sketch / rigorous speculation  
**Not claimed:** A completed derivation that QVC exhibits a Berezinskii–Kosterlitz–Thouless (BKT) transition  
**Locked baseline:** The QVC gap map is a saddle-node fold with finite \(\delta_c\) unless vortex renormalization-group (RG) physics is included  
**Companions:** `TORSION_TDVT_HOPFION.md`, `HYPOTHESES_PHASE1.md`  
**Canon:** `QVCCursor/rigor/canonical_formalism.py` (v5-lock); democratic Spec\((G)\) and fold locked in `CLOSED_FORM_LOCKS.md`

---

## 0. Epistemic labels

| Label | Meaning |
|-------|---------|
| **Locked** | Survives audit; usable without “unless” caveat |
| **Theorem-class (internal)** | Follows from stated assumptions by algebra |
| **Conjecture** | Plausible; needs computation/experiment before claim |
| **Analogy** | Structural resemblance; not identity |
| **Open** | Must be computed/measured before language upgrades |

This note is **ideation with equations**. Today’s locked theory is a **finite-gap fold**, not BKT.

---

## 1. Locked starting point: saddle-node fold

### 1.1 Gap map

Dimensionless gap \(\delta=\Delta/\Delta_0\) and vacuum amplitude \(\alpha\propto A_{\mathrm{ZPF}}\), with IR parameter \(\ell\):

\[
\delta = 1 + \alpha\bigl(\ln\delta + \ell\bigr).
\tag{V1}
\]

**Status:** Locked.

### 1.2 Fold conditions

Vertical tangent in the \((\alpha,\delta)\) plane:

\[
\alpha_c = \delta_c,\qquad
\delta_c\bigl(1 - \ln\delta_c - \ell\bigr) = 1.
\tag{V2}
\]

For Corrections \(\ell\approx -2.791\): \(\delta_c\approx 0.182\), \(\Delta_c\sim 0.11\,\mathrm{meV}\) (finite). **Status:** Locked.

### 1.3 Square-root scaling

Near the fold, saddle-node normal form:

\[
\Delta - \Delta_c \;\sim\; \pm\sqrt{A_c - A}.
\tag{V3}
\]

**Status:** Theorem-class for the amplitude truncation that produced (V1). Finite \(\Delta_c\) is **incompatible** with continuous BKT vanishing. Preprint language: “saddle-node fold unless vortex RG is included.”

### 1.4 Target BKT form (not locked)

\[
\Delta(A)\;\sim\; \Delta_*\exp\!\left(-\frac{C}{\sqrt{A_c-A}}\right)
\quad\text{as }A\to A_c^-.
\tag{V4}
\]

Earlier fits of the mean-field solver to (V4) are **artifacts of wrong functional form** unless a vortex sector is added.

---

## 2. Microscopic vortex in TEVC

### 2.1 Definition (proposed)

Let \(\Psi=|\Psi|e^{i\theta}\) be the excitonic/pairing condensate on the moiré scale. A vortex is a point defect with

\[
\oint\nabla\theta\cdot d\mathbf{l}=2\pi n,\quad n\in\mathbb{Z},
\]

core size \(\sim a_H=\hbar c_s/\Delta\) (or the excitonic healing length if distinct).

**Status:** Standard 2D superfluid identification; **application to vacuum-catalyzed QVC gap** is conjecture.

### 2.2 Energy and fugacity

With phase stiffness \(\kappa\),

\[
E_v \simeq \pi\kappa\,\ln\!\left(\frac{L}{a_H}\right)+E_{\mathrm{core}},
\tag{V5}
\]

\[
y=\exp\!\bigl(-E_{\mathrm{core}}/T_{\mathrm{eff}}\bigr).
\tag{V6}
\]

\(T_{\mathrm{eff}}\) may be thermal \(k_BT\) or an effective phonon/vacuum noise temperature. At true \(T=0\), quantum KT / dual gauge physics is required (**open**).

**Literature directions:** Kosterlitz–Thouless; Halperin–Nelson–Young; Minnhagen; MATBG stiffness/geometry (Tian–Gao–…, *Nature* 2023).

---

## 3. KT-like flow for moiré excitonic order

### 3.1 Classical KT equations

\[
\begin{aligned}
\frac{dK^{-1}}{d\ell_{\mathrm{RG}}} &= 4\pi^3 y^2+O(y^4),\\
\frac{dy}{d\ell_{\mathrm{RG}}} &= (2-\pi K)y+O(y^3),
\end{aligned}
\tag{V7}
\]

with \(K=\pi\kappa/(k_BT)\) and \(\ell_{\mathrm{RG}}=\ln(b/a_H)\).

### 3.2 QVC adaptation (conjecture)

1. Bare \(\kappa_0\to\kappa(\Delta)\) from quantum geometry / excitonic condensate.  
2. Couple \(\Delta\) (or \(\alpha\)) to running \(K\) so driving \(A\to A_c\) reduces \(\kappa\) and increases \(y\).  
3. IR cutoff = system / domain / screening length \(L_*\).

Schematic coupled flow (**not derived**):

\[
\begin{aligned}
\partial_\ell K^{-1} &= 4\pi^3 y^2,\\
\partial_\ell y &= (2-\pi K)y,\\
\partial_\ell \alpha &= \beta_\alpha(\alpha,K,y;\ell_{\mathrm{phonon}}).
\end{aligned}
\tag{V8}
\]

### 3.3 How RG can replace \(\sqrt{\,}\) by essential singularity

**Mechanism sketch (conjecture):** The fold describes the amplitude sector with phase fluctuations Gaussian-integrated. With vortices retained, the transition is controlled by the KT separatrix \(K\to 2/\pi\), \(y\to 0\). Identifying reduced drive \(t_A=(A_c-A)/A_c\) with the parameter that drives \(\kappa(A)\) through \(\kappa_{\mathrm{BKT}}\) yields

\[
\Delta\sim\frac{\hbar c_s}{\xi_{\mathrm{KT}}}\sim\exp\!\Bigl(-\frac{C}{\sqrt{A_c-A}}\Bigr).
\tag{V9}
\]

**Falsifier:** If measured \(\Delta(A)\) near collapse fits \(\sqrt{A_c-A}\) better than (V9) over \(\gtrsim 1.5\) decades in \(t_A\), vortex-RG promotion fails for that protocol.

---

## 4. Schematic: fold vs BKT

```
Amplitude-only (LOCKED)
  α ↑
  |     fold: branches merge at finite δ_c
  |        *  ← saddle-node (√ scaling)
  |       / \
  +-------------→ δ

With vortex RG (CONJECTURE)
  K^{-1}
  |     flows → y→∞ (plasma) or y→0 (bound)
  |        ╲   ╱  separatrix
  |         ╲ ╱
  |          *  ← essential singularity; Δ→0
  +-------------→ y
```

---

## 5. Gao–Khalaf phonons and vortex damping

Phonons in MATBG/moiré systems renormalize \(\kappa\) and \(E_{\mathrm{core}}\) and provide Bardeen–Stephen-like drag on vortices.

**Conjectural consequences:**

1. Overdamped vortices can destroy the sharp KT jump → rounded crossover.  
2. Vacuum/ZPF drive may raise \(T_{\mathrm{eff}}\) via phonon/Floquet heating → \(y\uparrow\).  
3. At mK scales, phonon-assisted quantum phase slips may dominate thermal KT.

**Open:** Extract vortex mobility \(\eta_v\); compare \(\hbar/\tau_v\) to \(\pi\kappa\).

**Literature directions:** Gao–Khalaf phonon/geometry line in TBG; AHNS dissipative KT; exciton-polariton BKT with pump heating.

---

## 6. Must-compute before claiming BKT in a preprint

| # | Task | Acceptance |
|---|------|------------|
| C1 | Microscopic \(\kappa(\Delta,n,T)\) for TEVC | OOM match to measured stiffness if available |
| C2 | \(E_{\mathrm{core}}\), \(a_H\) from BdG/GP on moiré cell | \(E_{\mathrm{core}}/\kappa=O(1)\)–\(O(10)\) |
| C3 | Integrate (V7)/(V8) with \(\kappa(A)\) | Separatrix intersects physical drive path |
| C4 | Continuous \(\Delta\to 0\) or \(\kappa\to\kappa_{\mathrm{BKT}}\) | No finite-\(\Delta_c\) remnant |
| C5 | Model comparison: \(\ln\Delta\) vs \(1/\sqrt{t_A}\) vs \(\frac12\ln t_A\) | Essential singularity preferred (AIC/evidence) |
| C6 | Include phonon damping | Separatrix survives or controlled crossover still yields (V9) |
| C7 | Finite-size \(g(r)\sim r^{-\eta}\), vortex density | \(\eta\to 1/4\) at putative \(A_c\) |

Until C1–C7: retain “saddle-node fold; BKT is a vortex-RG target, not locked.”

---

## 7. Experimental signatures: fold vs BKT

| Observable | Fold (locked) | BKT (conjecture) |
|------------|---------------|------------------|
| Gap vs drive \(A\) | Ends at finite \(\Delta_c\); \(\sqrt{\,}\) | Continuous \(\Delta\to 0\); \(\exp(-C/\sqrt{t})\) |
| Stiffness \(\kappa\) | Need not hit universal jump | Universal jump (thermal) or quantum analogue |
| Vortex density | Not required | Proliferates above \(A_c\) |
| Latent heat | Possible MF anomaly | None; essential singularity |
| Phase correlator | Amplitude death first | Algebraic → exponential |
| Hall/Chern | May jump with fold | Speculative staircase with vortex plasma |

**Protocol (ideation):** Gate/cavity-tune \(A_{\mathrm{ZPF}}\) while measuring gap (STS/ARPES) and stiffness (kinetic inductance).

---

## 8. Layer discipline

- **L1 QVC:** fold locked; vortex RG conjectural.  
- **L2 topological CM:** Hopf/CS fractionalization — partially locked as identities.  
- **L3 spacetime:** hopfion geometry — conjecture / analogue gravity.

Vortex worldlines as torsion sources: see `TORSION_TDVT_HOPFION.md`. Do not use L3 language to “prove” BKT.

---

## 9. Honest preprint sentence

> The present amplitude theory yields a saddle-node fold at finite gap; whether vacuum driving instead precipitates a BKT vortex-unbinding transition remains an open renormalization-group problem whose acceptance tests are listed in `HYPOTHESES_PHASE1.md`.
