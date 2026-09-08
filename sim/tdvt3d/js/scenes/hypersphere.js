// Tab 7 — Hyperspherical S³ flow (L2/L3). ≥200k particles, Hopf U(1), stereographic, Clifford torus, Hopf link.

import { LOCKS, fmt } from '../math/locks.js';
import { twoFibreLinking } from '../math/hopf.js';
import { HOPF_FLOW_VS, POINT_FS, makeLineLoop, makeNmAxes, disposeObject } from '../gfx/common.js';

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.fog = new THREE.FogExp2(0x05060a, 0.0035);
  scene.add(makeNmAxes(THREE, 55));

  const state = {
    stereo: true,
    omega: 0.55,
    showClifford: true,
    showLink: true,
    n: 220000,
    size: 2.0,
    speed: 1,
    glow: 1,
  };

  const uniforms = {
    uTime: { value: 0 },
    uOmega: { value: state.omega },
    uStereo: { value: 1 },
    uScale: { value: 42 },
    uSize: { value: state.size },
    uOpacity: { value: state.glow },
  };

  let pts = null;
  let lastT = 0;
  let nClifford = 0;

  function buildCloud() {
    if (pts) { scene.remove(pts); disposeObject(pts); pts = null; }
    const n = Math.max(40000, state.n | 0);
    state.n = n;
    const eta = new Float32Array(n);
    const psi = new Float32Array(n);
    const t0 = new Float32Array(n);
    const dummyPos = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      const u = (i + 0.5) / n;
      eta[i] = 0.5 * Math.acos(1 - 2 * ((i * 0.61803398875) % 1));
      if (i % 11 === 0) eta[i] = Math.PI / 4;
      psi[i] = Math.PI * 2 * ((i * 0.4142135623) % 1);
      t0[i] = Math.PI * 2 * u;
    }
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(dummyPos, 3));
    geo.setAttribute('aEta', new THREE.BufferAttribute(eta, 1));
    geo.setAttribute('aPsi', new THREE.BufferAttribute(psi, 1));
    geo.setAttribute('aT0', new THREE.BufferAttribute(t0, 1));
    const mat = new THREE.ShaderMaterial({
      uniforms,
      vertexShader: HOPF_FLOW_VS,
      fragmentShader: POINT_FS,
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    });
    pts = new THREE.Points(geo, mat);
    pts.frustumCulled = false;
    scene.add(pts);
    uniforms.uTime.value = lastT * state.speed;
  }

  buildCloud();

  const link = twoFibreLinking({
    nSeg: 320, stereo: true, scale: 42,
    eta1: Math.PI / 4, eta2: 0.7, psi1: 0, psi2: Math.PI / 2,
  });
  const f1 = makeLineLoop(THREE, link.c1, 0xffe08a);
  const f2 = makeLineLoop(THREE, link.c2, 0x7ee0ff);
  scene.add(f1, f2);

  const nC = 240;
  const cpack = new Float32Array(nC * nC * 3);
  let p = 0;
  for (let i = 0; i < nC; i += 4) {
    const ps = (i / nC) * Math.PI * 2;
    for (let j = 0; j < nC; j += 4) {
      const t = (j / nC) * Math.PI * 2;
      const c = Math.SQRT1_2;
      const ph1 = t + 0.5 * ps, ph2 = t - 0.5 * ps;
      const w0 = c * Math.cos(ph1), w1 = c * Math.sin(ph1);
      const w2 = c * Math.cos(ph2), w3 = c * Math.sin(ph2);
      const d = 1 + w3;
      cpack[p * 3] = w0 / d * 42;
      cpack[p * 3 + 1] = w1 / d * 42;
      cpack[p * 3 + 2] = w2 / d * 42;
      p++;
    }
  }
  nClifford = p;
  const clifford = new THREE.Points(
    new THREE.BufferGeometry().setAttribute('position', new THREE.BufferAttribute(cpack.slice(0, p * 3), 3)),
    new THREE.PointsMaterial({ color: 0xd48cff, size: 1.6, transparent: true, opacity: 0.55, sizeAttenuation: true }),
  );
  scene.add(clifford);

  function applyStyle() {
    uniforms.uSize.value = state.size;
    uniforms.uOpacity.value = state.glow;
    uniforms.uOmega.value = state.omega;
    uniforms.uStereo.value = state.stereo ? 1 : 0;
    clifford.visible = state.showClifford;
    f1.visible = state.showLink;
    f2.visible = state.showLink;
  }

  function applyPreset(id) {
    if (id === 'tdvt') {
      Object.assign(state, { n: 220000, omega: 0.55, stereo: true, size: 2.0, speed: 1, glow: 1, showClifford: true, showLink: true });
    } else if (id === 'precision') {
      Object.assign(state, { n: 280000, omega: 0.4, size: 1.6, speed: 0.8, glow: 0.9 });
    }
    buildCloud();
    applyStyle();
  }

  return {
    scene,
    viewScale: 90,
    get particleCount() { return state.n + nClifford; },
    sliders() {
      return [
        { id: 'n', label: 'n particles', min: 80000, max: 300000, step: 10000, value: state.n, fmt: (v) => String(v | 0), oninput: (v) => { state.n = v | 0; buildCloud(); } },
        { id: 'om', label: 'U(1) ω', min: 0, max: 1.6, step: 0.05, value: state.omega, fmt: (v) => v.toFixed(2), oninput: (v) => { state.omega = v; uniforms.uOmega.value = v; } },
        { id: 'st', label: 'stereographic', type: 'checkbox', value: state.stereo, oninput: (v) => { state.stereo = v; uniforms.uStereo.value = v ? 1 : 0; } },
        { id: 'size', label: 'point size', min: 0.6, max: 6, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; uniforms.uSize.value = v; } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; uniforms.uOpacity.value = v; } },
        { id: 'speed', label: 'speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; uniforms.uTime.value = lastT * v; } },
        { id: 'cl', label: 'Clifford torus', type: 'checkbox', value: state.showClifford, oninput: (v) => { state.showClifford = v; clifford.visible = v; } },
        { id: 'lk', label: 'Hopf link', type: 'checkbox', value: state.showLink, oninput: (v) => { state.showLink = v; f1.visible = v; f2.visible = v; } },
      ];
    },
    numbers() {
      return [
        { label: 'particles', value: state.n.toLocaleString(), kind: 'live' },
        { label: 'uTime', value: fmt(uniforms.uTime.value, 3), kind: 'live' },
        { label: 'flow', value: 'Hopf U(1) on fibres', kind: 'live' },
        { label: 'chart', value: state.stereo ? 'stereographic R³' : 'S³ embed (drop w3)', kind: 'live' },
        { label: 'Clifford torus', value: 'η = π/4', kind: 'live' },
        { label: 'Gauss lk (two fibres)', value: fmt(link.Lk, 4), kind: 'live' },
        { label: 'H-HS1', value: 'conjecture — schematic', kind: 'artistic' },
        { label: 'Q_H lock', value: fmt(LOCKS.Q_H, 3), kind: 'lock' },
      ];
    },
    charts() {
      return [
        { label: 'uTime', value: uniforms.uTime.value, color: '#d48cff', unit: '' },
        { label: 'ω', value: state.omega, color: '#7ee0ff', unit: '' },
      ];
    },
    banner() {
      return 'L2/L3 · Hyperspherical torsion-flow cartoon: particles stream along Hopf U(1) fibres. Magenta = Clifford torus. Gold/cyan = Hopf link. Stereographic toggle. H-HS1 (torsion ∝ mean curvature on the soliton surface) is a conjecture, not a lock.';
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => applyPreset('tdvt') },
      precision: { label: 'Precision', apply: () => applyPreset('precision') },
    },
    reset() { lastT = 0; uniforms.uTime.value = 0; },
    update(_dt, t, { paused } = {}) {
      lastT = t;
      if (paused) return;
      uniforms.uTime.value = t * state.speed;
    },
    dispose() {},
  };
}
