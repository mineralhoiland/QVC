# Visualizations (not verification)

HTML/Three.js suites copied from `QVCCursor/sim-suites` and
`QVCCompute/web/tdvt3d`. They illustrate hopfion, cavity, and fold geometry.

Do **not** treat these pages as the numerical lock. Use
`../compute/scripts/run_verification.py` and `../locks/PARAMETERS.md`.

```bash
cd standalone
python3 -m http.server 8765
```

`node_modules/` is not shipped. From this directory: `npm install` if you
need the Vite app in `src/`.
