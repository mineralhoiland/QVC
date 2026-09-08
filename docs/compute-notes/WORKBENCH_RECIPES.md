# TEVC / QVC–TDVT 3D Crystal & Materials Workbench — Recipes

The workbench (`qvc.workbench`) wires the existing QVCCompute solvers to 3D
fields, graphs and exporters. It adds **no new physics**: every number on
screen is read from `qvc.workbench.locks` (Option A\* freeze) and each viewport
is stamped with an epistemic layer:

| Layer | Meaning | Where it appears |
|---|---|---|
| **L1** | laboratory locks (fold, Spec(G), Hall step, E_vac, Δ\*) | lattice, dynamics, spectra viewports; lock table |
| **L2** | analogue geometry / hopfion texture diagnostics | fields viewport, optional analogue panel |
| **L3** | spacetime ontology — schematic text only | optional overlay note |

Withdrawn folklore numbers (0.244 meV, 92 %, 0.032 meV, 0.389 meV, 285 mK,
A_ZPF·ξ = 0.296, λ_min = −4g₀/−5g₀) are on a reject list; the lock builder
raises `DeprecatedNumberError` if any viewport would reproduce one.

---

## 1. Install

```bash
cd QVCCompute
python -m venv .venv && source .venv/bin/activate
pip install -e ".[workbench,materials]"     # panel, bokeh, plotly, kaleido, pyvista, pymatgen
```

## 2. Interactive workbench

```bash
python -m qvc.workbench                     # http://localhost:5006
python -m qvc.workbench --port 5010 --fields-grid 48 --gpe-grid 96 --gpe-t-ps 1.0
qvc-workbench --show                        # console script, opens a tab
```

Layout (four viewports + optional L2 row):

1. **Lattice (L1)** — magic-angle TEVC AA fragment from `build_tevc_stack`
   (θ, N sliders) or the Im-3m hopfion BCC cell with computed preimage fibres.
2. **Fields (L2 + L1 locks)** — Whitehead hopfion ρ_H isosurfaces with n̂
   glyphs restricted to the core; preimage-linking curves (Gauss linking → 1);
   Casimir torus with θ = Q_H φ branch cut (2D CS flux → 1/N) and the −E_vac
   envelope.
3. **Dynamics (L1)** — `qvc.tqgl` GPE↔DCE replay: live Δ(ψ) vs static fold
   Δ(A), A/A_c^fold, g⁽¹⁾, 80 ps DCE slaved to the fold; 2D |ψ|² frame and the
   same frame wrapped on the torus (an honest extrusion, not a 3D PDE).
4. **Spectra / barriers (L1)** — fold map (λ slider; λ ≥ λ_fold reports
   "fold-out, no physical branch"), Spec(G) democratic vs WS-overlap
   comparison, Mexican hat with U_eff(A) at fixed total energy, Arrhenius/WKB
   rate sketch, Hall staircase, BM Berry curvature.

Sidebar controls: θ, N, λ, g₀, A_ZPF/A_c, GPE frame. Layers: "Show L2
analogue-geometry panel", drain v₀/c_s (0 = pure fractional vortex, no sonic
horizon; ≥1 creates a diagnostic draining-bathtub horizon whose T_H is
*not* T\*), H-TOR1 toggle (parked — annotates only), L3 note.

Export menu writes to `output/workbench_exports/`:

| Choice | Files | Open with |
|---|---|---|
| CIF/POSCAR (VESTA) | `TEVC.cif`, `POSCAR`, `TEVC_fragment.cif`, `tevc.data`, `hopfion_crystal_BCC_Im3m.cif`, `hopfion_preimage_overlay.cif` | VESTA, OVITO (`tevc.data`) |
| CUBE + VTK | `hopfion_hopf_density.cube`, `hopfion_FN_energy_density.cube`, `hopfion_nz.cube`, `hopfion_casimir_envelope.cube`, `hopfion_fields.vti`, `hopfion_torus.vti` | VESTA (CUBE), ParaView (VTI) |
| GPE time series | `gpe_replay_2d.pvd` + `.vti` frames, `gpe_replay_torus.pvd` + `.vts` frames, `gpe_replay_timeseries.csv` | ParaView |
| Locks JSON | `workbench_locks.json` (locks + reject list) | any |

PNG: use the camera icon on any Plotly panel, or run the batch script.

## 3. Headless batch (build → solve → export)

```bash
python scripts/run_crystal_suite.py                 # → output/crystal_suite/, exit 1 if a gate fails
python scripts/run_crystal_suite.py --quick         # 32³ / 48² smoke run (~30 s)
python scripts/run_crystal_suite.py --html --out /tmp/suite   # + static workbench.html
```

Produces everything in the table above plus `fold_map.png`, `specG.png`,
`mexican_hat.png`, `rates.png`, `hall_staircase.png`, `dynamics_timeseries.png`,
`goldstone_L2.png`, `acoustic_L2.png`, `suite_report.json`, `acceptance.json`.

Acceptance gates (S9): locks match `docs/PHASE_TWO_PLAN.md`; Whitehead Q → 1
by Richardson extrapolation over 32/48/64 (the raw value at one grid is
reported, not gated); 2D CS flux = 1/N; VK bound; Hall step = e²/(2Nh);
λ_min = −g₀ (×N−1); fold α_c = 0.1818; GPE norm conserved; TA phonon has zero
TT overlap.

```bash
pytest tests/test_workbench.py -m "not slow"      # ~5 s
pytest tests/test_workbench.py                     # + GPE replay + Panel build
```

## 4. VESTA recipes

**Atomic TEVC** — `File ▸ Open ▸ TEVC.cif` (26 450 C at θ = 1.10°, 5 layers,
λ_M = 12.8 nm; unrelaxed analytic ±θ/2 stack). For a light view open
`TEVC_fragment.cif` (AA-site disc, r = 2 nm). `Objects ▸ Boundary` to extend;
`Edit ▸ Bonds` C–C 1.2–1.6 Å. The historical 13–22° CIFs in
`QVCCursor/TEVC_structures` are proxies, not the magic angle.

**Hopfion BCC** — open `hopfion_crystal_BCC_Im3m.cif`. Space group Im-3m
(229), Wyckoff 2a, a = ξ = 50 nm written as 50 display-Å. The species `Ho` is a
*display proxy* for a Q = 1 hopfion centre (not holmium). Overlay
`hopfion_preimage_overlay.cif` (two linked fibres as dummy sites; the fibre
linking is 1 while the L1 lens-space lock is Q_H = 1/N — different numbers).

**Volumes** — `File ▸ Open ▸ hopfion_hopf_density.cube`, then
`Properties ▸ Isosurfaces`: set the level to ~35 % of max (positive lobe, red)
and add a negative level (blue). The CUBE header records the Whitehead
normalisation `1/16π²` and the grid Q; the old
`TEVC_structures/cube/hopfion_charge_density.cube` header (Q ≈ 3.81 with a
1/4π² factor) is rejected. Repeat with `hopfion_FN_energy_density.cube`
(Ward-normalised; relaxed Q = 1 reference ≈ 1.22, VK bound 0.534|Q|^{3/4}) and
`hopfion_casimir_envelope.cube` (−E_vac = −0.1287 meV at the tube wall;
negative isosurfaces). `hopfion_nz.cube` at level 0 gives the n_z = 0 torus.

Vector overlay: VESTA renders vectors per site, not per voxel — use ParaView
(below) for n̂ glyphs; in VESTA, use the `nz` and `rho_H` volumes together.

## 5. ParaView recipes

**Static hopfion** — open `hopfion_fields.vti`. Arrays: `n` (vector), `rho_H`,
`energy_density`, `A`, `B`.
- *Contour* on `rho_H` (two levels ±0.35·max) → the linked-torus lobes.
- *Glyph* (Arrow) on `n`, mask by `rho_H > 0.03·max` via *Threshold* first,
  colour by `n_Z`.
- *Stream Tracer* on `B` (Berry curvature) seeded on a line through the core:
  the field lines are the Hopf fibres; two seeds link once.

**Casimir torus** — `hopfion_torus.vti`: `casimir_meV` (contour at −0.5·E_vac),
`rho_s`, `theta` (colour by `theta` with the *Phase*/*HSV* map to see the
Q_H φ branch cut).

**Dynamics** — open `gpe_replay_torus.pvd` (structured-grid torus coloured by
|ψ|²) or `gpe_replay_2d.pvd` (the actual 2D field). Use the animation toolbar;
`gpe_replay_timeseries.csv` (t, Δ_live, Δ_fold(A), A/A_c, g⁽¹⁾) can be loaded
with *Plot Data* for the live-gap graph. The torus is an extrusion of the 2D
solution along φ — it is not a 3D PDE (deferred; see plan §S5).

## 6. Honesty notes carried by the exporters

- Lattice: TEVC is the *unrelaxed* analytic stack (LAMMPS relaxation via
  `tevc.data` in OVITO/LAMMPS is the user's step).
- Fields: the hopfion n̂ is the stereographic Whitehead ansatz on the torus
  scale R (L2 geometry); the E_vac envelope and Q_H are the L1 locks.
- Spec(G): the democratic g₀(J − I) spectrum is the lock; the Wigner–Seitz
  overlap G^(H) is *compared* to it, never forced to J − I.
- Dynamics: live Δ(ψ) is the drive channel and may exceed Δ₀; the static fold
  Δ(A) → Δ\* = 0.2324 meV is the parallel channel. They are never stacked.
- L2 analogue panel: Unruh–Visser kinematics only. The pure fractional vortex
  has no sonic horizon (ρ_h inside the core); the drain slider fabricates one
  for illustration and its T_H is not the Schwinger lock T\*. H-TOR1 (torsion →
  G_ij) is parked; H-EIN1 (TA phonon as graviton) is refuted by the zero TT
  overlap.
- L3 overlay: schematic text. No Sakharov G, no gravitons, no cosmology.
