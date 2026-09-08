# TDVT GPU Particle Geometry Lab

GPU-optimized Three.js particle simulation (≥512K default, up to 1M) focused on TDVT topology and geometry — not mesh primitives.

## Run

```bash
cd sim-suites
npm install
npm run dev
```

Open http://localhost:5173

## Modes

| Mode | Description |
|------|-------------|
| Hopf S³ | Stereographic Hopf fibration time evolution |
| Hopfions | Linked hopfion tubes (two particle families) |
| TDVT lattice | BCC lattice of interlinked hopfion cores |
| Vortices | Particle lattice + multi-vortex Kelvin ripples |
| Catalyzon | Democratic Spec(G) swarm + fold melt |
| Horizon | Acoustic horizon flow past hopfion core |

## Architecture

- Positions/colors integrated in **vertex shaders** (no per-frame CPU particle loops)
- `THREE.Points` + additive blending + soft disc fragment shader
- React UI only for controls/readouts; WebGL owned by `SceneHost`
