# QVCCompute: Multi-Tier Computational Framework for TEVC Physics

**Version:** 0.2.0 (Architecture Specification)  
**Author:** Mineral Hoiland  
**Date:** August 2026  
**Status:** Design document — implementation in progress

---

## 1. Design Philosophy

The TEVC (Topological Excitonic Vacuum Crystal) is an intrinsically multi-scale object.
Three energy scales interact non-trivially (Section 1.3 of the v4 paper):

| Scale | Physics | Length | Energy | Computation |
|-------|---------|--------|--------|-------------|
| Crystallographic | Atomic lattice → moiré | 0.14–13 nm | ~eV (bonds) → ~10 meV (flat band) | DFT, tight-binding |
| Mesoscale | Condensate, hopfion, cavity | 8–50 nm | 0.1–0.6 meV (gap) | TQGL PDE, eigenvalue problems |
| Macroscopic | Transport, phase transitions | 50 nm–µm | µeV–mK (BKT, Schwinger) | RG flow, Monte Carlo |

The computational framework mirrors this separation. Each tier is **independently testable**, feeds
parameters upward, and validates against the tier below. The critical principle:

> **Ab initio sets parameters. Effective theory does the physics. Visualization communicates both.**

---

## 2. Tier Architecture

```
┌─────────────────────────────────────────────────────────────────────┐
│                    TIER 0: AB INITIO GROUND TRUTH                     │
│         Quantum ESPRESSO  ·  LAMMPS  ·  Wannier90  ·  phonopy        │
│                                                                       │
│  Purpose: Establish microscopic parameters for monolayer/bilayer      │
│  graphene and small-angle commensurate approximants. Run once (or     │
│  when new materials are considered). Results stored as JSON/HDF5.     │
└───────────────────────────────────┬─────────────────────────────────┘
                                    │ extracted parameters
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                 TIER 1: STRUCTURE & MODEL HAMILTONIANS                 │
│            pymatgen  ·  ASE  ·  PythTB  ·  kwant  ·  numpy           │
│                                                                       │
│  Purpose: Build moiré supercells, relax geometry, construct           │
│  tight-binding and continuum model Hamiltonians, compute bands,       │
│  Berry curvature, folded phonon dispersions.                          │
└───────────────────────────────────┬─────────────────────────────────┘
                                    │ matched coefficients (α_GL, κ_GL, ...)
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                  TIER 2: EFFECTIVE THEORY (QVC Core)                   │
│            qvc/ package  ·  scipy  ·  sympy  ·  numba                 │
│                                                                       │
│  Purpose: Gap equation, TQGL dynamics, DCE self-limiting,             │
│  vortex-RG sector, catalyzon hybridization. The existing              │
│  QVCCompute modules live here.                                        │
└───────────────────────────────────┬─────────────────────────────────┘
                                    │ observables, phase boundaries
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                   TIER 3: TOPOLOGY & TRANSPORT                         │
│              numpy  ·  scipy  ·  numba  ·  Monte Carlo                │
│                                                                       │
│  Purpose: Chern number, Hopf charge, Hall staircase, Schwinger        │
│  rate, BKT transition diagnostics, g(r) correlator.                   │
└───────────────────────────────────┬─────────────────────────────────┘
                                    │ figures, animations, interactives
                                    ▼
┌─────────────────────────────────────────────────────────────────────┐
│                      TIER 4: VISUALIZATION                            │
│     matplotlib  ·  Three.js  ·  Plotly  ·  pymatgen.vis  ·  VESTA    │
│                                                                       │
│  Purpose: Publication figures, interactive parameter explorers,        │
│  3D structure renderings, time-evolution animations.                   │
└─────────────────────────────────────────────────────────────────────┘
```

---

## 3. TIER 0: Ab Initio Ground Truth

### 3.1 Quantum ESPRESSO (QE)

**Role:** DFT electronic structure and DFPT phonon calculations for graphene-based systems.

**What we compute:**
- Monolayer graphene band structure (validates tight-binding parameterization)
- Bilayer graphene at AA, AB, SP stacking (interlayer coupling V_pp_σ, V_pp_π)
- Phonon dispersions from DFPT (LA, TA, ZA branches + optical modes)
- Electron-phonon coupling matrix elements g_ν(k,q)
- Small-angle commensurate approximants (e.g., θ = 1.08° with ~11,000 atoms)

**Key calculations:**

```
QE_calculations/
├── monolayer/
│   ├── scf/              # Self-consistent field (ground state)
│   ├── bands/            # Band structure along high-symmetry path
│   ├── phonon/           # DFPT phonon dispersions
│   └── epw/              # Electron-phonon (Wannier interpolation)
├── bilayer_AA/
│   ├── scf/
│   └── bands/            # Interlayer hybridization at AA stacking
├── bilayer_AB/
│   ├── scf/
│   └── bands/            # Interlayer hybridization at AB stacking
├── bilayer_SP/
│   └── scf/              # Saddle-point stacking (domain walls)
└── commensurate_1p08/
    ├── scf/              # Full moiré cell DFT (expensive!)
    └── bands/            # Flat-band verification
```

**Input generation:** pymatgen's `PWInput` class generates QE input files from structures.

**Output parsing:** pymatgen parses `pw.x` and `ph.x` output; for EPW we use custom parsers.

**Key parameters extracted:**

| Parameter | Source calculation | Value (graphene) | Use in QVC |
|-----------|-------------------|------------------|------------|
| v_F | monolayer bands | ~10⁶ m/s | BM model velocity |
| t_⊥ | bilayer bands | ~0.11 eV | Interlayer hopping |
| c_LA | monolayer phonon | ~21.4 km/s | Phonon dispersion |
| c_TA | monolayer phonon | ~13.6 km/s | Phonon dispersion |
| α_ZA | monolayer phonon | ~6.2×10⁻⁷ m²/s | Flexural mode |
| g_ep | EPW | ~0.1 eV (deformation) | Electron-phonon coupling |

**Pseudopotentials:** PAW (PBE) from PSlibrary; van der Waals: DFT-D3 or vdW-DF2.

**Computational cost:**
- Monolayer/bilayer: minutes on a workstation (2-4 atoms/cell)
- Commensurate approximant: ~10⁴ atoms → needs HPC (100–1000 CPU-hours)
- Alternative: use the BM model parameters from literature as Tier 0 proxy

### 3.2 LAMMPS

**Role:** Classical molecular dynamics for structural relaxation of the moiré superlattice.

**Why this matters:** At θ ≈ 1.1°, atoms in MATBG don't remain on ideal twisted-lattice
positions. They reconstruct: large triangular AB/BA domains form with sharp domain walls
(~1–2 nm width). This reconstruction affects:
- Local stacking registry → modulates interlayer coupling
- In-plane strain field → modifies phonon dispersion
- Out-of-plane corrugation → changes effective interlayer distance
- Domain wall network → defines the effective Casimir cavity geometry

**What we compute:**
- Fully relaxed atomic positions for N-layer moiré at given θ
- Corrugation profile h(x,y) for each layer
- In-plane displacement field u(x,y) (strain tensor)
- Stacking registry map (AA fraction, AB/BA domains, SP walls)
- Interlayer distance variation d_⊥(x,y)

**Force field:** DRIP (Dihedral-angle-corrected Registry-dependent Interlayer Potential)
for interlayer interactions + REBO for intralayer C-C bonds. Alternative: Kolmogorov-Crespi
potential for interlayer + AIREBO for intralayer.

```python
# LAMMPS input template structure
"""
units           metal
atom_style      atomic
boundary        p p p

# Read moiré supercell (generated by pymatgen)
read_data       moire_5layer_1p1deg.data

# Force fields
pair_style      hybrid/overlay rebo kolmogorov/crespi/z 16.0
pair_coeff      * * rebo CH.rebo C C C C C
pair_coeff      * * kolmogorov/crespi/z CC.KC C C C C C

# Relaxation
fix             1 all box/relax x 0.0 y 0.0
minimize        1.0e-8 1.0e-10 10000 100000

# Dump relaxed structure
dump            1 all custom 1 relaxed.dump id type x y z
"""
```

**Output → Tier 1 bridge:**
- Relaxed positions feed into tight-binding model (modulated hopping)
- Corrugation profile feeds into effective cavity geometry
- Domain wall width validates Casimir cavity scale assumptions

**Computational cost:** ~10⁴–10⁵ atoms per moiré cell for 5 layers.
Minutes to hours on a workstation. Very feasible without HPC.

### 3.3 Wannier90

**Role:** Extract maximally-localized Wannier functions from DFT, producing a
tight-binding Hamiltonian that reproduces the DFT bands exactly in the low-energy window.

**Workflow:**
1. QE `pw.x` → ground state
2. QE `pw2wannier90.x` → overlaps and projections
3. Wannier90 → spreads, centers, hopping matrices t_ij(R)
4. Result: real-space Hamiltonian H(R) that can be Fourier-interpolated to any k

**Use in QVC:** The Wannier TB model serves as the definitive parameterization
for the Bistritzer-MacDonald-like Hamiltonian in Tier 1. It gives you:
- Accurate interlayer tunneling as a function of stacking registry
- Correct Fermi velocity renormalization
- Berry connection matrix elements for topological invariants

### 3.4 phonopy / phono3py

**Role:** Phonon calculations via the finite-displacement method (supercell approach).

**Advantage over DFPT:** Can use any force calculator (DFT, classical potential, ML potential).
For the moiré phonon problem, we can use LAMMPS forces (classical potential) to get
the dynamical matrix of the full moiré cell — something DFPT cannot do at 10⁴ atoms.

**Workflow:**
1. Generate displaced supercells (phonopy)
2. Compute forces (LAMMPS or QE)
3. Construct dynamical matrix
4. Diagonalize → phonon bands, DOS, group velocities

**Key output for QVC:**
- Folded phonon dispersion in moiré BZ
- Flat phonon bands near zone boundary (from hybridization)
- Phonon density of states at resonance frequency ω₀ = 2Δ/ℏ
- Acoustic sound velocities (validates c_s = 0.07 meV·µm / ℏ)

---

## 4. TIER 1: Structure & Model Hamiltonians

### 4.1 pymatgen (Structure Management Hub)

**Role:** The Swiss army knife for crystal structure manipulation. Serves as the central
structure-handling layer that generates inputs for all other codes.

**Key operations for TEVC:**

```python
"""
TEVC structure generation with pymatgen
"""
from pymatgen.core import Structure, Lattice
from pymatgen.transformations.standard_transformations import (
    SupercellTransformation,
    RotationTransformation,
)
import numpy as np


def build_graphene_monolayer(a_cc: float = 1.42) -> Structure:
    """Monolayer graphene unit cell.
    
    Parameters
    ----------
    a_cc : float
        C-C bond length in Angstroms (1.42 Å for graphene).
        
    Returns
    -------
    Structure
        2-atom hexagonal unit cell with 20 Å vacuum.
    """
    a = a_cc * np.sqrt(3)  # lattice constant = 2.46 Å
    c = 20.0  # vacuum spacing
    
    lattice = Lattice.hexagonal(a, c)
    coords = [[1/3, 2/3, 0.5], [2/3, 1/3, 0.5]]
    
    return Structure(lattice, ["C", "C"], coords)


def build_twisted_bilayer(theta_deg: float, a_cc: float = 1.42) -> Structure:
    """Commensurate twisted bilayer graphene.
    
    For small angles, the commensurate approximant has moiré period
    λ_M = a / (2 sin(θ/2)). We find the smallest commensurate cell
    by matching (m, n) indices.
    
    Parameters
    ----------
    theta_deg : float
        Twist angle in degrees.
        
    Returns
    -------
    Structure
        Commensurate moiré supercell (may be large for small θ).
    """
    # Implementation uses commensurate angle algorithm
    # For θ = 1.1°: approximately (m,n) = (31,32), N_atoms ≈ 11,164
    ...


def build_tevc_stack(
    theta_deg: float = 1.1,
    n_layers: int = 5,
    d_perp: float = 3.35,
    alternating: bool = True,
) -> Structure:
    """Full TEVC multilayer structure.
    
    Builds N-layer stack with alternating twist angles ±θ.
    
    Parameters
    ----------
    theta_deg : float
        Twist angle magnitude (each layer alternates ±θ/2 from vertical).
    n_layers : int
        Number of graphene layers (default 5 for MATBG TEVC).
    d_perp : float
        Interlayer spacing in Angstroms (3.35 Å nominal, varies with relaxation).
    alternating : bool
        If True, twist angles alternate: +θ, -θ, +θ, -θ, +θ.
        If False, all layers twisted by θ relative to bottom.
        
    Returns
    -------
    Structure
        Full multilayer moiré supercell.
    """
    ...


def generate_lammps_input(struct: Structure, output_dir: str) -> None:
    """Write LAMMPS data file from pymatgen Structure."""
    ...


def generate_qe_input(struct: Structure, output_dir: str, calc_type: str = "scf") -> None:
    """Write Quantum ESPRESSO input from pymatgen Structure.
    
    Uses pymatgen.io.pwscf.PWInput with appropriate settings for
    graphene/MATBG calculations.
    """
    from pymatgen.io.pwscf import PWInput
    ...
```

**Additional pymatgen roles:**
- Symmetry analysis of moiré cell (SpacegroupAnalyzer)
- K-path generation for band structure plots (HighSymmKpath)
- Brillouin zone construction and high-symmetry point identification
- Format conversion: POSCAR ↔ CIF ↔ XYZ ↔ LAMMPS data
- Interface to Materials Project (retrieve graphene reference data)

### 4.2 Bistritzer-MacDonald Continuum Model

**Role:** The workhorse Hamiltonian for moiré flat bands. Operates in k-space on the
continuum (no need for the full 11,000-atom supercell).

**Hamiltonian structure (per valley, per spin):**

```
H_BM(k) = | H_1(k)    T(r)   |
           | T†(r)     H_2(k) |
```

where H_i(k) = ℏv_F (k - K_i) · σ is the Dirac Hamiltonian for layer i,
and T(r) is the interlayer tunneling matrix modulated by the moiré pattern:

T(r) = Σ_j T_j exp(i q_j · r)

with T_j encoding AA vs AB tunneling (w_AA ≈ 80 meV, w_AB ≈ 110 meV).

**Implementation:**

```python
"""
Bistritzer-MacDonald model for TEVC flat bands.
Extends to N-layer system with alternating twist.
"""
import numpy as np
from scipy.linalg import eigh


class BistritzerMacDonald:
    """N-layer continuum model for twisted graphene.
    
    Parameters extracted from Tier 0 (QE + Wannier90) or literature:
    - v_F: Fermi velocity (≈ 10⁶ m/s for graphene)
    - w_AA: interlayer tunneling at AA sites (≈ 80 meV)  
    - w_AB: interlayer tunneling at AB sites (≈ 110 meV)
    - theta: twist angle
    """
    
    def __init__(self, theta_deg, n_layers=5, v_F=1e6, w_AA=0.080, w_AB=0.110):
        self.theta = np.radians(theta_deg)
        self.n_layers = n_layers
        self.v_F = v_F
        self.w_AA = w_AA  # eV
        self.w_AB = w_AB  # eV
        
        # Moiré reciprocal vectors
        self.K_M = 4 * np.pi / (3 * 0.246e-9) * 2 * np.sin(self.theta / 2)
        
    def hamiltonian(self, k, n_shells=4):
        """Build Hamiltonian matrix at momentum k.
        
        Parameters
        ----------
        k : array (2,)
            Momentum in moiré BZ (units of K_M).
        n_shells : int
            Number of reciprocal lattice shells to include.
            
        Returns
        -------
        H : ndarray
            Hamiltonian matrix (size depends on n_shells and n_layers).
        """
        ...
    
    def band_structure(self, k_path, n_bands=10):
        """Compute bands along k-path.
        
        Returns
        -------
        energies : ndarray (n_k, n_bands)
        """
        energies = []
        for k in k_path:
            H = self.hamiltonian(k)
            evals = eigh(H, eigvals_only=True)
            energies.append(evals)
        return np.array(energies)
    
    def berry_curvature(self, k_mesh, band_index=0):
        """Compute Berry curvature Ω(k) on a discrete mesh.
        
        Uses Fukui's lattice gauge method for numerical stability.
        
        Returns
        -------
        omega : ndarray (n_kx, n_ky)
            Berry curvature at each k-point.
        chern : float
            Integrated Chern number (should be integer).
        """
        ...
    
    def flat_band_width(self):
        """Compute bandwidth W of the lowest moiré bands.
        
        At the magic angle θ ≈ 1.1°, W ≲ 10 meV.
        This validates the flat-band limit U/W → ∞.
        """
        ...
```

### 4.3 Folded Phonon Model

**Role:** Construct phonon dispersion in the moiré Brillouin zone from monolayer
graphene phonon branches + interlayer coupling.

**Physics:** The moiré folding takes graphene's phonon branches ω(q) and creates
replicas at q ± nK_M. Near the zone boundary, interlayer coupling opens gaps
and produces flat phonon bands. The resonance condition ω₀ = 2Δ*/ℏ picks out
specific phonon modes that mediate catalyzon formation.

```python
class MoirePhononModel:
    """Folded phonon dispersion in moiré Brillouin zone.
    
    Inputs (from Tier 0):
    - c_LA, c_TA: in-plane acoustic velocities
    - alpha_ZA: flexural mode coefficient
    - omega_S: interlayer shear frequency
    - K_M: moiré reciprocal lattice constant
    - interlayer_coupling: hybridization at zone boundary
    """
    
    def __init__(self, c_LA=21.4e3, c_TA=13.6e3, alpha_ZA=6.2e-7,
                 omega_S_cm=20.0, K_M=0.32e9, coupling_meV=0.5):
        ...
    
    def dispersion(self, q_path, n_replicas=3):
        """Compute folded phonon bands along q-path in moiré BZ."""
        ...
    
    def dos_at_frequency(self, omega):
        """Phonon density of states at frequency ω.
        
        Critical for the catalyzon resonance condition.
        """
        ...
    
    def casimir_mode_sum(self, geometry='fat_torus', r=8e-9, R=50e-9, n_modes=14):
        """Evaluate Casimir energy via mode summation on cavity geometry.
        
        This is the Tier 1 microscopic version of the locked vacuum
        E_vac = 0.1287 meV in the existing qvc/ code.
        """
        ...
```

### 4.4 kwant (Quantum Transport)

**Role:** Compute transport properties (conductance, Hall response) on finite-size
TEVC samples with realistic geometry.

**Use cases:**
- Landauer conductance through a moiré device
- Hall conductance vs. magnetic field (validate σ_xy = e²/2Nh staircase)
- Edge states and topological boundary modes
- Disorder effects on quantization

---

## 5. TIER 2: Effective Theory (Existing QVC Core)

This is the existing `qvc/` package. The new architecture connects it to
microscopically-derived parameters from Tier 1.

### 5.1 Parameter Bridge

The critical connection is a **parameter bridge** that translates Tier 1 outputs into
Tier 2 inputs:

```python
"""
qvc/parameter_bridge.py — Connect microscopic calculations to effective theory.
"""
from dataclasses import dataclass
from qvc.params import QVCParams


@dataclass
class MicroscopicInputs:
    """Parameters extracted from Tier 0/1 calculations."""
    
    # From BM model
    bandwidth_meV: float          # W: flat-band bandwidth
    effective_mass: float         # m*: flat-band effective mass (in m_e)
    berry_curvature_max: float    # max |Ω(k)| in BZ
    chern_number: int             # C: band Chern number
    
    # From phonon model
    c_LA_m_per_s: float          # longitudinal acoustic velocity
    c_TA_m_per_s: float          # transverse acoustic velocity
    omega_res_meV: float         # resonance frequency ω₀ = 2Δ/ℏ
    phonon_dos_at_res: float     # ρ_ph(ω₀): density of states at resonance
    
    # From LAMMPS relaxation
    domain_wall_width_nm: float  # width of SP domain walls
    corrugation_amplitude_A: float  # out-of-plane buckling
    aa_fraction: float           # fraction of cell that is AA-stacked
    
    # From Casimir mode sum
    casimir_energy_meV: float    # E_vac from microscopic mode sum


def bridge_to_qvc_params(micro: MicroscopicInputs) -> QVCParams:
    """Convert microscopic inputs to QVC effective theory parameters.
    
    The key identifications:
    - J_meV = bandwidth / 2  (exchange energy from Coulomb in flat band)
    - ℏc_s = gap × healing_length  (Bogoliubov velocity package)
    - E_vac from microscopic Casimir sum → validates/replaces locked value
    
    Returns a QVCParams instance with microscopically-grounded values.
    """
    J = micro.bandwidth_meV / 2  # Exchange ~ half bandwidth
    Delta0 = 2 * J               # Bare excitonic gap
    
    # Healing length from effective mass
    # a_H = ℏ / sqrt(2 m* Δ₀)
    import numpy as np
    hbar_eV_s = 6.582e-16
    m_star_kg = micro.effective_mass * 9.109e-31
    Delta0_J = Delta0 * 1.602e-22  # meV → J
    a_H_m = hbar_eV_s * 1.602e-19 / np.sqrt(2 * m_star_kg * Delta0_J)
    a_H_um = a_H_m * 1e6
    
    # Bogoliubov sound velocity
    c_s = Delta0_J / (hbar_eV_s * 1.602e-19) * a_H_m  # m/s
    hbar_cs = 6.582e-16 * c_s * 1e6 * 1e3  # meV·µm
    
    return QVCParams(
        J_meV=J,
        hbar_cs_meV_um=hbar_cs,
        N_layers=5,
        lambda_M_um=0.013,
        r_preprint_um=a_H_um,
        R_preprint_um=50e-3,  # 50 nm in µm
    )
```

### 5.2 Existing Modules (Preserved)

The current `qvc/` package remains the Tier 2 core:

| Module | Physics | Status |
|--------|---------|--------|
| `params.py` | Canonical parameter set | LOCKED |
| `gap_fold.py` | Amplitude map δ(α), saddle-node fold | LOCKED |
| `vacuum_freeze.py` | Option A* frozen geometry | LOCKED |
| `vacuum_convergence.py` | Bessel sum convergence | LOCKED |
| `bridge.py` | α–E_vac bridge | DERIVED |
| `coupling_matrix.py` | Democratic G = g₀(J-I) | DERIVED |
| `vortex_rg.py` | Kosterlitz flow equations | DERIVED |
| `hopf_integral.py` | Faddeev hopfion on L(N,1) | COMPUTED |
| `self_consistency.py` | Fixed-point iteration | COMPUTED |
| `retardation.py` | Causal propagator η_ret | LOCKED |

### 5.3 New Tier 2 Modules (To Build)

| Module | Physics | Dependencies |
|--------|---------|--------------|
| `dynamics/tqgl_solver.py` | Split-step Fourier TQGL | numpy.fft, scipy.integrate |
| `dynamics/dce_ode.py` | DCE self-limiting ODE | scipy.integrate.solve_ivp |
| `dynamics/lindblad.py` | Open-system dissipator | numpy |
| `dynamics/schwinger.py` | Pair production rate | existing params |
| `topology/berry_phase.py` | Berry phase from TQGL | numpy |
| `topology/hopfion_field.py` | 3D n-hat texture | numpy |
| `topology/hall.py` | σ_xy staircase | existing CS algebra |

---

## 6. TIER 3: Topology & Transport

### 6.1 Monte Carlo Engine

```python
"""
qvc/montecarlo/xy_model.py — 2D XY model Metropolis for vortex sector.

Purpose: Sample the phase correlator g(r) = ⟨e^{i[θ(r)-θ(0)]}⟩
at physical parameters (K, y) from the vortex-RG sector.

Validates:
- η = 1/(2πK) → 1/4 at KT transition
- Algebraic vs exponential decay
- Vortex density as function of temperature
"""
import numpy as np
from numba import njit


@njit
def metropolis_sweep(theta, K, L, beta=1.0):
    """Single Metropolis sweep over L×L XY lattice.
    
    Energy: H = -K Σ_{<ij>} cos(θ_i - θ_j)
    """
    for _ in range(L * L):
        i, j = np.random.randint(0, L), np.random.randint(0, L)
        dtheta = (np.random.random() - 0.5) * np.pi
        
        # Nearest neighbors (periodic BC)
        neighbors = [
            theta[(i+1) % L, j], theta[(i-1) % L, j],
            theta[i, (j+1) % L], theta[i, (j-1) % L]
        ]
        
        dE = 0.0
        for th_n in neighbors:
            dE += -K * (np.cos(theta[i,j] + dtheta - th_n) 
                       - np.cos(theta[i,j] - th_n))
        
        if dE < 0 or np.random.random() < np.exp(-beta * dE):
            theta[i,j] += dtheta
            theta[i,j] = theta[i,j] % (2 * np.pi)
    
    return theta


def phase_correlator(theta, L):
    """Compute g(r) = ⟨e^{i[θ(r)-θ(0)]}⟩ from equilibrated configuration."""
    ...
```

### 6.2 Berry Phase Monte Carlo

For the random-matrix comparison (R_mean = 0.908 from the paper):

```python
def berry_phase_ensemble(n_samples=10000, matrix_size=5):
    """Sample Berry phases from GUE/GOE random coupling matrices.
    
    Compares against the locked value R_mean = 0.908.
    """
    ...
```

---

## 7. TIER 4: Visualization

### 7.1 Publication Figures (matplotlib)

Style: colorblind-safe, vector PDF output, consistent font sizing.

```python
# qvc/viz/style.mplstyle
figure.figsize: 3.375, 2.5     # PRL single-column width
font.size: 8
axes.labelsize: 9
legend.fontsize: 7
lines.linewidth: 1.0
axes.prop_cycle: cycler('color', ['0C5DA5', 'FF2C00', '00B945', 'FF9500', '845B97'])
```

**Figure modules:**
- `material_structure.py`: Moiré lattice, BZ, stacking map
- `band_structure.py`: Flat bands, Berry curvature heatmap, DOS
- `phase_diagrams.py`: Gap map, competing-scale, T-window
- `dynamics.py`: TQGL snapshots, DCE trajectories, live gap
- `topology.py`: Hall staircase, hopfion field lines, Chern density

### 7.2 Interactive (HTML/Three.js)

Single-file HTML applications for exploration:

| App | Physics | Controls |
|-----|---------|----------|
| `tevc_explorer.html` | 3D layer stack + moiré | θ slider, layer toggle, zoom |
| `hopfion_viewer.html` | S³→S² field texture | rotation, cross-section plane |
| `gap_map_live.html` | Bifurcation diagram | α, λ sliders, branch highlighting |
| `phase_space.html` | (α, K, y) 3D | rotate, slice plane, trajectory overlay |
| `catalyzon_modes.html` | 5-mode hybridization | coupling g₀ slider, eigenvalue animation |
| `hall_staircase.html` | σ_xy vs B | N slider, step animation |

### 7.3 Structure Visualization (pymatgen + VESTA)

For publication-quality crystal structure images:

```python
from pymatgen.vis.structure_vtk import StructureVis

# Or export to VESTA format:
struct.to(filename="tevc.vasp")  # Open in VESTA for ray-traced rendering
```

VESTA produces ray-traced images of crystal structures that are standard in
materials science publications. pymatgen generates the structure; VESTA renders it.

---

## 8. Directory Structure (Full)

```
QVCCompute/
├── pyproject.toml                  # Package configuration
├── ARCHITECTURE.md                 # ← This document
│
├── qvc/                            # Tier 2: Effective theory (EXISTING)
│   ├── __init__.py
│   ├── params.py                   # Canonical parameters
│   ├── gap_fold.py                 # Amplitude map
│   ├── bridge.py                   # α–E_vac bridge
│   ├── vacuum_freeze.py            # Option A* freeze
│   ├── coupling_matrix.py          # Democratic G_ij
│   ├── vortex_rg.py                # Kosterlitz flow
│   ├── hopf_integral.py            # Faddeev hopfion
│   ├── self_consistency.py         # Fixed-point iteration
│   └── verify_all.py              # Lock validation
│
├── qvc/materials/                  # Tier 1: Structure & Hamiltonians (NEW)
│   ├── __init__.py
│   ├── lattice.py                  # Graphene lattice primitives
│   ├── moire.py                    # Moiré supercell builder (pymatgen)
│   ├── bistritzer_macdonald.py     # BM continuum model
│   ├── phonon_folding.py           # Folded phonon dispersion
│   ├── berry_curvature.py          # Ω(k) and Chern number
│   ├── relaxation.py              # LAMMPS input/output bridge
│   └── parameter_bridge.py         # Tier 0/1 → Tier 2 mapping
│
├── qvc/dynamics/                   # Tier 2 extension: PDE/ODE (NEW)
│   ├── __init__.py
│   ├── tqgl_solver.py             # Split-step Fourier TQGL
│   ├── dce_ode.py                  # DCE self-limiting dynamics
│   ├── lindblad.py                 # Open-system dissipation
│   └── schwinger.py                # Pair production rate
│
├── qvc/topology/                   # Tier 3: Topological computations (NEW)
│   ├── __init__.py
│   ├── hopfion_field.py            # 3D texture and Hopf charge
│   ├── chern_number.py             # Lattice gauge method
│   ├── hall_staircase.py           # σ_xy steps
│   └── berry_phase_stats.py        # Random-matrix Berry ensemble
│
├── qvc/montecarlo/                 # Tier 3: Statistical mechanics (NEW)
│   ├── __init__.py
│   ├── xy_model.py                 # 2D XY Metropolis
│   └── vortex_sampling.py          # Vortex configurations
│
├── qvc/viz/                        # Tier 4: Visualization (NEW)
│   ├── __init__.py
│   ├── publication/
│   │   ├── style.mplstyle
│   │   ├── material_structure.py
│   │   ├── band_structure.py
│   │   ├── phase_diagrams.py
│   │   ├── dynamics.py
│   │   └── topology.py
│   ├── interactive/
│   │   ├── tevc_explorer.html
│   │   ├── hopfion_viewer.html
│   │   ├── gap_map_live.html
│   │   ├── phase_space.html
│   │   └── catalyzon_modes.html
│   └── structure/
│       ├── vesta_export.py         # pymatgen → VESTA format
│       └── povray_render.py        # Ray-traced structure images
│
├── tier0/                          # Tier 0: Ab initio inputs (NEW)
│   ├── quantum_espresso/
│   │   ├── templates/              # QE input templates
│   │   ├── pseudopotentials/       # PAW PPs (or symlinks)
│   │   └── parse_results.py        # Output → JSON extraction
│   ├── lammps/
│   │   ├── templates/              # LAMMPS input scripts
│   │   ├── potentials/             # DRIP, KC, REBO files
│   │   └── parse_dump.py           # Relaxed structure extraction
│   └── wannier90/
│       ├── templates/              # Wannier90 input templates
│       └── extract_tb.py           # TB parameters → JSON
│
├── notebooks/                      # Interactive development
│   ├── 01_vacuum_convergence.ipynb # EXISTING
│   ├── 02_gij_wigner_seitz.ipynb   # EXISTING
│   ├── 03_vortex_rg_flow.ipynb     # EXISTING
│   ├── 04_hopf_integral.ipynb      # EXISTING
│   ├── 05_lindblad_tqgl.ipynb      # EXISTING
│   ├── 06_moire_structure.ipynb    # NEW: Build & visualize TEVC
│   ├── 07_band_structure.ipynb     # NEW: BM model flat bands
│   ├── 08_phonon_folding.ipynb     # NEW: Moiré phonon bands
│   ├── 09_lammps_relaxation.ipynb  # NEW: Structural relaxation
│   └── vacuum_bridge_lab.ipynb     # EXISTING
│
├── scripts/                        # Batch execution
│   ├── run_verification.py         # EXISTING
│   ├── generate_all_figures.py     # NEW: Master figure pipeline
│   ├── run_tier1_pipeline.py       # NEW: Structure → parameters
│   └── run_full_suite.py           # NEW: All tiers, all outputs
│
├── tests/                          # Validation
│   ├── test_locks.py               # EXISTING
│   ├── test_tqgl_v4.py            # EXISTING
│   ├── test_moire.py              # NEW: Structure generation
│   ├── test_bm_model.py           # NEW: Band structure
│   └── test_bridge.py             # NEW: Parameter bridge
│
├── data/                           # Cached results and references
│   ├── tier0_cache/                # Pre-computed DFT results (JSON)
│   ├── relaxed_structures/         # LAMMPS output structures
│   └── reference_values/           # Literature comparison data
│
└── output/                         # Generated figures and reports
    ├── figures/
    ├── verification_report.json    # EXISTING
    └── parameter_audit.json        # NEW: Full parameter provenance
```

---

## 9. Dependency Management

### 9.1 Core (always required)

```toml
[project]
dependencies = [
    "numpy>=1.24",
    "scipy>=1.10",
    "sympy>=1.12",
]
```

### 9.2 Materials (Tier 0/1)

```toml
[project.optional-dependencies]
materials = [
    "pymatgen>=2024.1",
    "ase>=3.22",          # Atomic Simulation Environment
    "phonopy>=2.20",       # Phonon calculations
    "spglib>=2.0",         # Symmetry detection
]
```

### 9.3 Visualization (Tier 4)

```toml
[project.optional-dependencies]
viz = [
    "matplotlib>=3.7",
    "plotly>=5.15",
    "kaleido>=0.2",        # Plotly static export
]
```

### 9.4 Computation (Tier 2/3)

```toml
[project.optional-dependencies]
compute = [
    "numba>=0.57",         # JIT for Monte Carlo inner loops
    "h5py>=3.8",           # Large dataset storage
]
```

### 9.5 Full suite

```toml
[project.optional-dependencies]
full = [
    "qvc[materials,viz,compute]",
    "jupyter>=1.0",
    "ipywidgets>=8.0",
]
```

### 9.6 External codes (not pip-installable)

| Code | Installation | Purpose |
|------|-------------|---------|
| Quantum ESPRESSO | system/HPC module | DFT + DFPT |
| LAMMPS | system/conda | Classical MD relaxation |
| Wannier90 | system/HPC module | TB extraction |
| VESTA | GUI application | Structure rendering |

---

## 10. Workflow: From Scratch to Full TEVC Simulation

### Phase 1: Structure (Weeks 1–2)

```
[pymatgen] Build 5-layer moiré structure
     ↓
[LAMMPS]  Relax atomic positions (DRIP + REBO)
     ↓
[pymatgen] Parse relaxed structure, identify domains
     ↓
[viz]     Render moiré pattern, layer stacking, domain map
```

**Deliverables:** Relaxed TEVC structure, domain wall width, corrugation profile.

### Phase 2: Electronic Structure (Weeks 3–4)

```
[BM model] Compute flat-band dispersion
     ↓
[Berry]   Calculate Ω(k), Chern number, Berry connection
     ↓
[DOS]     Flat-band density of states
     ↓
[viz]     Band structure plot, Berry curvature heatmap
```

**Deliverables:** Bandwidth W, Chern number C, validates flat-band regime.

### Phase 3: Phonon & Vacuum (Weeks 5–6)

```
[phonon]  Fold phonon branches into moiré BZ
     ↓
[Casimir] Mode sum on relaxed torus geometry
     ↓
[bridge]  Extract E_vac, compare to locked 0.1287 meV
     ↓
[viz]     Phonon dispersion, Casimir energy convergence
```

**Deliverables:** Microscopic E_vac, phonon DOS at resonance, validates vacuum lock.

### Phase 4: Dynamics & Transport (Weeks 7–8)

```
[TQGL]   Solve PDE with matched coefficients from Phases 1–3
     ↓
[DCE]    Couple to self-limiting ODE
     ↓
[MC]     Run XY model at physical K(α)
     ↓
[Hall]   Compute σ_xy staircase
     ↓
[viz]    Animations, phase diagrams, Hall plot
```

**Deliverables:** Full dynamical simulation, validates paper results end-to-end.

---

## 11. Validation Strategy

Every tier validates against the one below:

| Check | Expected | Source |
|-------|----------|--------|
| BM bandwidth at 1.1° | W ≲ 10 meV | Tier 1 vs. literature |
| Chern number | C = ±1 per valley | Tier 1 vs. exact |
| Phonon c_LA | ~21 km/s | Tier 1 vs. Tier 0 (QE) |
| E_vac from mode sum | 0.1287 meV | Tier 1 vs. Tier 2 lock |
| Gap Δ* from TQGL | 0.2324 meV | Tier 2 vs. lock |
| Phase coherence g^(1) | 0.8832 | Tier 2 vs. lock |
| η at primary branch | ≪ 1/4 | Tier 3 vs. C1 stiffness |
| Hall step | e²/(2Nh) | Tier 3 vs. CS algebra |

---

## 12. What Can Run Now vs. What Needs HPC

### Runs on your workstation (minutes):

- [x] Moiré structure generation (pymatgen)
- [x] BM model band structure (numpy eigensolve, ~seconds for 100×100 k-mesh)
- [x] Folded phonon model (parametric, ~seconds)
- [x] TQGL PDE on 48² grid (existing, ~seconds)
- [x] Vortex-RG flow integration (existing, ~seconds)
- [x] Berry curvature on BM bands (~minutes for 100×100 mesh)
- [x] XY Monte Carlo L=16 (existing, ~minutes)
- [x] All visualizations

### Needs modest computation (hours, single workstation):

- [ ] LAMMPS relaxation of 5-layer moiré (~50,000 atoms, ~1 hour)
- [ ] XY Monte Carlo L=64 for thermodynamic g(r) (~hours)
- [ ] TQGL on 128² or 256² grid (~hours)
- [ ] phonopy with LAMMPS forces (~hours)

### Needs HPC (days, cluster):

- [ ] QE DFT on commensurate moiré (~11,000 atoms per bilayer period)
- [ ] EPW electron-phonon interpolation
- [ ] Full 5-layer DFT (~55,000 atoms — at the frontier)
- [ ] Large-scale TQGL (512² or 3D) with noise

### Pragmatic approach:

For the 5-layer TEVC, **do NOT attempt full DFT** of the entire moiré cell.
Instead:

1. Use DFT on small cells (monolayer, bilayer at fixed stackings) to extract parameters
2. Use LAMMPS for full 5-layer relaxation (classical potentials, totally feasible)
3. Use BM continuum model for electronic structure (parameterized by DFT values)
4. Use parametric phonon model (calibrated against monolayer DFPT)

This gives you ab-initio-quality parameters at workstation cost.

---

## 13. Summary: The Tool Chain

```
                        STRUCTURE
                           │
              pymatgen ────┤──── LAMMPS
             (build)       │    (relax)
                           │
                        ELECTRONS
                           │
        Bistritzer-MacDonald ──── (QE/Wannier90 for params)
            (bands, Berry)
                           │
                        PHONONS
                           │
          phonopy/model ───┤──── (QE DFPT for calibration)
          (folded bands)   │
                           │
                     EFFECTIVE THEORY
                           │
              qvc/ package ┤ (gap, TQGL, DCE, vortex-RG)
                           │
                       OBSERVABLES
                           │
       Monte Carlo ────────┤──── topology (Hopf, Chern, Hall)
       (XY model, g(r))   │
                           │
                     VISUALIZATION
                           │
    matplotlib ──── Three.js ──── VESTA ──── Plotly
    (figures)    (interactive)   (3D struct) (dashboards)
```

The beauty of this architecture is that each piece is useful independently.
You don't need Quantum ESPRESSO running to use the BM model — just use
literature parameters. You don't need LAMMPS to do TQGL dynamics — just
use the existing locked parameters. But when you DO connect the pieces,
you get end-to-end parameter provenance from atoms to observables.
