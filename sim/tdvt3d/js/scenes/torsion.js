// Tab 5 — Teleparallel torsion. ≥50k torus tracers, precessing tetrad glyphs. H-TOR1 parked. No FFT in loop.

import { LOCKS, fmt } from '../math/locks.js';
import { torsionProxy, parallelogramGap, H_TOR1_STATUS } from '../math/torsion.js';
import { makePointCloud, makeNmAxes, disposeObject, setPointStyle } from '../gfx/common.js';

const MIN_N = 50000;

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.fog = new THREE.FogExp2(0x05060a, 0.004);
  scene.add(makeNmAxes(THREE, 90));

  const state = {
    n: MIN_N,
    amp: 10,
    omega: LOCKS.omega0_star_ps_inv,
    size: 1.8,
    speed: 1,
    twist: 5,
    nGlyph: 16,
    glow: 1,
    htor1: false,
  };

  let tracers = null;
  let glyphs = null;
  let paraLine = null;
  let rhoH = null;
  let maxT = 0;
  let gapN = 0;
  let lastT = 0;
  const dummy = new THREE.Object3D();
  const gDir = new THREE.Vector3();
  const gQ = new THREE.Quaternion();
  const yUp = new THREE.Vector3(0, 1, 0);

  const R = LOCKS.R_nm;
  const r = LOCKS.r_nm;

  function rebuild() {
    if (tracers) { scene.remove(tracers); disposeObject(tracers); }
    if (glyphs) { scene.remove(glyphs); disposeObject(glyphs); }
    if (paraLine) { scene.remove(paraLine); disposeObject(paraLine); }

    const n = Math.max(MIN_N, state.n | 0);
    state.n = n;
    const pos = new Float32Array(n * 3);
    const col = new Float32Array(n * 3);
    rhoH = new Float64Array(n);
    tracers = makePointCloud(THREE, pos, col, { size: state.size, opacity: state.glow });
    scene.add(tracers);

    const nG = Math.max(4, state.nGlyph | 0);
    state.nGlyph = nG;
    const geo = new THREE.ConeGeometry(Math.max(1.1, state.size * 0.9), state.size * 4.2, 6);
    const mat = new THREE.MeshBasicMaterial({ color: 0xff8866, transparent: true, opacity: 0.55 + 0.4 * state.glow });
    glyphs = new THREE.InstancedMesh(geo, mat, nG);
    glyphs.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    scene.add(glyphs);

    paraLine = new THREE.Line(
      new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(new Float32Array(18), 3)),
      new THREE.LineBasicMaterial({ color: 0xffe08a, transparent: true, opacity: 0.9 }),
    );
    scene.add(paraLine);

    stamp(lastT);
  }

  function stamp(t) {
    if (!tracers || !glyphs) return;
    const n = state.n;
    const pos = tracers.geometry.getAttribute('position');
    const col = tracers.geometry.getAttribute('color');
    const twist = state.omega * state.speed * t;
    const amp = state.amp;
    const wind = state.twist;

    for (let i = 0; i < n; i++) {
      const u = i / n;
      const th = u * Math.PI * 2;
      const ph = wind * th + twist;
      const Rt = R + 0.12 * amp * Math.cos(2 * th + twist * 0.3);
      const rt = r + 0.18 * amp * Math.sin(th);
      const x = (Rt + rt * Math.cos(ph)) * Math.cos(th);
      const y = (Rt + rt * Math.cos(ph)) * Math.sin(th);
      const z = rt * Math.sin(ph);
      pos.setXYZ(i, x, y, z);
      const dens = 1 / (1 + ((Math.hypot(Math.hypot(x, y) - R, z) / r) ** 2));
      rhoH[i] = dens;
      col.setXYZ(i, 0.25 + 0.7 * dens, 0.35 * (1 - dens), 0.85);
    }
    pos.needsUpdate = true;
    col.needsUpdate = true;

    const proxy = torsionProxy(rhoH, { hTor1: false });
    let tmax = 0;
    for (let i = 0; i < proxy.T.length; i++) {
      const a = Math.abs(proxy.T[i]);
      if (a > tmax) tmax = a;
    }
    maxT = tmax;

    const origin = [R, 0, 0];
    const e1 = [0, 12, 0];
    const e2 = [0, 0, 12];
    const Taxial = 0.35 + 0.25 * Math.sin(twist);
    const pg = parallelogramGap(origin, e1, e2, Taxial, { amp: state.amp * 0.4 });
    gapN = Math.hypot(pg.gap[0], pg.gap[1], pg.gap[2]);
    const seq = [pg.a, pg.b, pg.c, pg.d, pg.c2, pg.a];
    const pa = paraLine.geometry.getAttribute('position');
    for (let i = 0; i < 6; i++) pa.setXYZ(i, seq[i][0], seq[i][1], seq[i][2]);
    pa.needsUpdate = true;

    const nG = state.nGlyph;
    for (let g = 0; g < nG; g++) {
      const th = (g / nG) * Math.PI * 2 + t * 0.12 * state.speed;
      const ph = 3 * th + twist;
      dummy.position.set(
        (R + r * 0.7 * Math.cos(ph)) * Math.cos(th),
        (R + r * 0.7 * Math.cos(ph)) * Math.sin(th),
        r * 0.7 * Math.sin(ph),
      );
      gDir.set(-Math.sin(th + twist * 0.4), Math.cos(th + twist * 0.4), 0.35 * Math.sin(ph)).normalize();
      dummy.quaternion.copy(gQ.setFromUnitVectors(yUp, gDir));
      dummy.updateMatrix();
      glyphs.setMatrixAt(g, dummy.matrix);
    }
    glyphs.instanceMatrix.needsUpdate = true;
    glyphs.material.opacity = 0.55 + 0.4 * state.glow;
  }

  function applyPreset(id) {
    if (id === 'tdvt') {
      Object.assign(state, { n: MIN_N, amp: 10, omega: LOCKS.omega0_star_ps_inv, size: 1.8, speed: 1, twist: 5, nGlyph: 16, glow: 1 });
    } else if (id === 'precision') {
      Object.assign(state, { n: 80000, amp: 8, omega: LOCKS.omega0_star_ps_inv * 0.7, size: 1.4, twist: 5, nGlyph: 24, glow: 0.9 });
    }
    rebuild();
  }

  rebuild();

  return {
    scene,
    viewScale: 140,
    get particleCount() { return state.n + state.nGlyph; },
    sliders() {
      return [
        { id: 'n', label: 'n particles', min: MIN_N, max: 120000, step: 5000, value: state.n, fmt: (v) => String(v | 0), oninput: (v) => { state.n = v | 0; rebuild(); } },
        { id: 'tw', label: 'twist', min: 1, max: 12, step: 1, value: state.twist, fmt: (v) => String(v | 0), oninput: (v) => { state.twist = v; stamp(lastT); } },
        { id: 'omega', label: 'ω (ps⁻¹)', min: 0.1, max: 2.4, step: 0.02, value: state.omega, fmt: (v) => v.toFixed(2), oninput: (v) => { state.omega = v; stamp(lastT); } },
        { id: 'gd', label: 'glyph density', min: 6, max: 48, step: 1, value: state.nGlyph, fmt: (v) => String(v | 0), oninput: (v) => { state.nGlyph = v | 0; rebuild(); } },
        { id: 'amp', label: 'amp', min: 0, max: 28, step: 0.5, value: state.amp, fmt: (v) => v.toFixed(1), oninput: (v) => { state.amp = v; stamp(lastT); } },
        { id: 'size', label: 'point size', min: 0.6, max: 5, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; setPointStyle(tracers, { size: v }); } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; setPointStyle(tracers, { opacity: v }); if (glyphs) glyphs.material.opacity = 0.55 + 0.4 * v; } },
        { id: 'speed', label: 'speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; stamp(lastT); } },
        {
          id: 'htor1', label: `H-TOR1 (${H_TOR1_STATUS})`, type: 'checkbox',
          value: false, parked: true,
          oninput: () => { state.htor1 = false; },
        },
      ];
    },
    numbers() {
      return [
        { label: 'max |T|', value: fmt(maxT, 4), kind: 'live' },
        { label: '□ gap', value: fmt(gapN, 3), kind: 'live' },
        { label: 'H-TOR1', value: H_TOR1_STATUS, kind: 'diag' },
        { label: 'feeds G_ij', value: 'no', kind: 'lock' },
        { label: 'Q_H', value: fmt(LOCKS.Q_H, 3), kind: 'lock' },
        { label: '(r, R)', value: `(${LOCKS.r_nm}, ${LOCKS.R_nm}) nm`, kind: 'lock' },
        { label: 'ω₀*', value: `${fmt(LOCKS.omega0_star_ps_inv, 3)} ps⁻¹`, kind: 'lock' },
        { label: 'FFT in loop', value: 'none', kind: 'diag' },
      ];
    },
    charts() {
      return [
        { label: 'maxT', value: maxT, color: '#ff8866', unit: '' },
        { label: 'gap', value: gapN, color: '#9dffb0', unit: '' },
      ];
    },
    banner() {
      return `L2/L3 · Axial torsion proxy T ∝ ρ_H on the (8,50) nm torus. ${state.n.toLocaleString()} tracers, glyphs precess. Parallelogram non-closure is Cartan torsion. H-TOR1 is ${H_TOR1_STATUS} — never fed into G_ij. No FFT in the loop.`;
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => applyPreset('tdvt') },
      precision: { label: 'Precision', apply: () => applyPreset('precision') },
    },
    reset() { lastT = 0; stamp(0); },
    update(_dt, t, { paused } = {}) {
      lastT = t;
      if (paused) return;
      stamp(t);
    },
    dispose() {},
  };
}
