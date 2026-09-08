// Tab 3 — Living BCC Im-3m hopfion lattice. Goldstone u = A ê sin(k·R − ωt).

import { LOCKS, fmt, lensLocks } from '../math/locks.js';
import { goldstoneBranches, hEin1Overlap } from '../math/analogue.js';
import { makePointCloud, makeNmAxes, disposeObject, setPointStyle } from '../gfx/common.js';

const MODES = ['TA1', 'TA2', 'LA', 'magnon1', 'magnon2'];
const BRANCHES = goldstoneBranches();
const A_NM = LOCKS.a_nm;
const EIN = hEin1Overlap();

function omegaOf(kWant, mode) {
  const key = (mode === 'TA1' || mode === 'TA2') ? 'TA' : mode;
  const w = BRANCHES[key] || BRANCHES.TA;
  const k = BRANCHES.k;
  if (kWant <= k[0]) return w[0];
  const last = k.length - 1;
  if (kWant >= k[last]) return w[last];
  for (let i = 1; i <= last; i++) {
    if (kWant <= k[i]) {
      const t = (kWant - k[i - 1]) / (k[i] - k[i - 1] || 1);
      return w[i - 1] + t * (w[i] - w[i - 1]);
    }
  }
  return w[last];
}

function polarization(mode, kHat, out) {
  if (mode === 'LA') return out.copy(kHat);
  const up = Math.abs(kHat.y) < 0.92 ? { x: 0, y: 1, z: 0 } : { x: 1, y: 0, z: 0 };
  const tx = kHat.y * up.z - kHat.z * up.y;
  const ty = kHat.z * up.x - kHat.x * up.z;
  const tz = kHat.x * up.y - kHat.y * up.x;
  const tn = Math.hypot(tx, ty, tz) || 1;
  if (mode === 'TA2') {
    out.set(
      kHat.y * (tz / tn) - kHat.z * (ty / tn),
      kHat.z * (tx / tn) - kHat.x * (tz / tn),
      kHat.x * (ty / tn) - kHat.y * (tx / tn),
    );
    return out.normalize();
  }
  return out.set(tx / tn, ty / tn, tz / tn);
}

function bccSites(half) {
  const h = Math.max(1, Math.round(half));
  const sites = [];
  for (let i = -h; i <= h; i++) {
    for (let j = -h; j <= h; j++) {
      for (let k = -h; k <= h; k++) {
        sites.push(i * A_NM, j * A_NM, k * A_NM);
        sites.push((i + 0.5) * A_NM, (j + 0.5) * A_NM, (k + 0.5) * A_NM);
      }
    }
  }
  return sites;
}

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.fog = new THREE.FogExp2(0x05060a, 0.0035);

  const state = {
    N: LOCKS.N,
    half: 2,
    A: 2,
    kFrac: 0.35,
    mode: 'TA1',
    speed: 1,
    size: 2.4,
    coreSize: 2.4,
    fibrePer: 160,
    glow: 1,
    showFibres: true,
  };

  const dummy = new THREE.Object3D();
  const eHat = new THREE.Vector3(0, 1, 0);
  const kHat = new THREE.Vector3(1, 0, 0);
  const tmp = new THREE.Vector3();

  let axes = makeNmAxes(THREE, A_NM * 3);
  scene.add(axes);
  let cores = null;
  let fibres = null;
  let restCores = null;
  let restFibres = null;
  let omega = 0;
  let uMax = 0;
  let nFibres = 0;
  let lastT = 0;

  function rebuild() {
    if (cores) { scene.remove(cores); disposeObject(cores); cores = null; }
    if (fibres) { scene.remove(fibres); disposeObject(fibres); fibres = null; }
    if (axes) { scene.remove(axes); disposeObject(axes); }
    axes = makeNmAxes(THREE, A_NM * (state.half + 1.2));
    scene.add(axes);

    const packed = bccSites(state.half);
    const nCore = packed.length / 3;
    restCores = new Float32Array(packed);

    const geo = new THREE.SphereGeometry(Math.max(1.2, state.coreSize * 1.05), 8, 8);
    const mat = new THREE.MeshBasicMaterial({ color: 0x7ee0ff, transparent: true, opacity: 0.55 + 0.4 * state.glow });
    cores = new THREE.InstancedMesh(geo, mat, nCore);
    cores.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
    scene.add(cores);

    const per = Math.max(40, state.fibrePer | 0);
    nFibres = nCore * per;
    restFibres = new Float32Array(nFibres * 3);
    const col = new Float32Array(nFibres * 3);
    const R = LOCKS.r_nm;
    const r = 2.4;
    let p = 0;
    for (let c = 0; c < nCore; c++) {
      const cx = restCores[c * 3], cy = restCores[c * 3 + 1], cz = restCores[c * 3 + 2];
      for (let i = 0; i < per; i++) {
        const u = i / per;
        const th = u * Math.PI * 2 * 3;
        const ph = (i * 1.6180339887) * Math.PI * 2;
        restFibres[p * 3] = cx + (R + r * Math.cos(ph)) * Math.cos(th);
        restFibres[p * 3 + 1] = cy + (R + r * Math.cos(ph)) * Math.sin(th);
        restFibres[p * 3 + 2] = cz + r * Math.sin(ph);
        col[p * 3] = 0.85; col[p * 3 + 1] = 0.55 + 0.35 * Math.sin(th); col[p * 3 + 2] = 0.25;
        p++;
      }
    }
    const fpos = restFibres.slice();
    fibres = makePointCloud(THREE, fpos, col, { size: Math.max(1.2, state.size * 0.7), opacity: state.glow });
    fibres.visible = state.showFibres;
    scene.add(fibres);
    applyField(lastT);
  }

  function applyField(t) {
    if (!cores || !restCores) return;
    const kMag = state.kFrac * BRANCHES.kMax;
    omega = omegaOf(kMag, state.mode);
    polarization(state.mode, kHat, eHat);
    const A = state.A;
    const wt = omega * state.speed * t;
    let max = 0;
    const coreScale = Math.max(0.35, state.coreSize / 2.4);

    const nCore = restCores.length / 3;
    for (let i = 0; i < nCore; i++) {
      const x = restCores[i * 3], y = restCores[i * 3 + 1], z = restCores[i * 3 + 2];
      const u = A * Math.sin(kMag * x - wt);
      const mag = Math.abs(u);
      if (mag > max) max = mag;
      dummy.position.set(x + eHat.x * u, y + eHat.y * u, z + eHat.z * u);
      dummy.scale.setScalar(coreScale);
      dummy.rotation.set(0, 0, 0);
      dummy.updateMatrix();
      cores.setMatrixAt(i, dummy.matrix);
    }
    cores.instanceMatrix.needsUpdate = true;
    cores.material.opacity = 0.55 + 0.4 * state.glow;

    if (fibres && state.showFibres) {
      const pos = fibres.geometry.getAttribute('position');
      const nF = pos.count;
      for (let i = 0; i < nF; i++) {
        tmp.set(restFibres[i * 3], restFibres[i * 3 + 1], restFibres[i * 3 + 2]);
        const u = A * Math.sin(kMag * tmp.x - wt);
        pos.setXYZ(i, tmp.x + eHat.x * u, tmp.y + eHat.y * u, tmp.z + eHat.z * u);
      }
      pos.needsUpdate = true;
    }
    uMax = max;
  }

  function applyStyle() {
    if (fibres) setPointStyle(fibres, { size: Math.max(1.2, state.size * 0.7), opacity: state.glow });
    if (cores) cores.material.opacity = 0.55 + 0.4 * state.glow;
    if (fibres) fibres.visible = state.showFibres;
  }

  function applyPreset(id) {
    if (id === 'tdvt') {
      Object.assign(state, {
        N: LOCKS.N, half: 2, A: 2, mode: 'TA1', kFrac: 0.35,
        speed: 1, size: 2.4, coreSize: 2.4, fibrePer: 160, glow: 1, showFibres: true,
      });
    } else if (id === 'precision') {
      Object.assign(state, { half: 3, N: Math.max(state.N, LOCKS.N), size: 2.0, coreSize: 2.0, fibrePer: 90, glow: 0.85 });
    }
    rebuild();
  }

  rebuild();

  return {
    scene,
    get viewScale() { return A_NM * (state.half + 2); },
    get particleCount() { return (cores?.count || 0) + nFibres; },
    sliders() {
      return [
        { id: 'N', label: 'N (lens)', min: 3, max: 9, step: 1, value: state.N, fmt: (v) => String(v | 0), oninput: (v) => { state.N = v | 0; } },
        { id: 'half', label: 'half (cells)', min: 1, max: 4, step: 1, value: state.half, fmt: (v) => String(v | 0), oninput: (v) => { state.half = v | 0; rebuild(); } },
        { id: 'mode', label: 'Goldstone', type: 'select', value: state.mode,
          options: MODES.map((m) => ({ value: m, label: m })),
          oninput: (v) => { state.mode = MODES.includes(v) ? v : 'TA1'; applyField(lastT); } },
        { id: 'A', label: 'A (nm)', min: 0, max: 8, step: 0.1, value: state.A, fmt: (v) => v.toFixed(1), oninput: (v) => { state.A = v; applyField(lastT); } },
        { id: 'k', label: 'k / k_max', min: 0.05, max: 1, step: 0.01, value: state.kFrac, fmt: (v) => v.toFixed(2), oninput: (v) => { state.kFrac = v; applyField(lastT); } },
        { id: 'speed', label: 'speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; applyField(lastT); } },
        { id: 'fdens', label: 'fibre density', min: 40, max: 240, step: 10, value: state.fibrePer, fmt: (v) => String(v | 0), oninput: (v) => { state.fibrePer = v | 0; rebuild(); } },
        { id: 'core', label: 'core size', min: 1, max: 8, step: 0.1, value: state.coreSize, fmt: (v) => v.toFixed(1), oninput: (v) => { state.coreSize = v; applyField(lastT); } },
        { id: 'size', label: 'point size', min: 1, max: 8, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; applyStyle(); } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; applyStyle(); } },
        { id: 'fib', label: 'show fibres', type: 'checkbox', value: state.showFibres, oninput: (v) => { state.showFibres = v; applyStyle(); } },
      ];
    },
    numbers() {
      const L = lensLocks(state.N);
      return [
        { label: 'mode ω', value: fmt(omega, 4), kind: 'live' },
        { label: '|u|max', value: `${fmt(uMax, 3)} nm`, kind: 'live' },
        { label: 'a = ξ', value: `${LOCKS.a_nm} nm`, kind: 'lock' },
        { label: 'E_vac', value: `${fmt(LOCKS.E_vac_meV, 4)} meV`, kind: 'lock' },
        { label: 'Q_H = 1/N', value: fmt(L.Q_H, 3), kind: 'lock' },
        { label: 'θ_anyon', value: fmt(L.theta_anyon_rad, 4), kind: 'lock' },
        { label: 'Hall step', value: fmt(L.hall_step_e2_over_h, 3), kind: 'lock' },
        { label: 'H-EIN1 overlap', value: EIN.overlap_plus.toExponential(1), kind: 'live' },
        { label: 'TA vs graviton', value: 'helicity ±1, not TT', kind: 'diag' },
      ];
    },
    charts() {
      return [
        { label: 'ω', value: omega, color: '#7ee0ff', unit: '' },
        { label: '|u|max', value: uMax, color: '#ffe08a', unit: 'nm' },
      ];
    },
    banner() {
      return `L2 · BCC Im-3m a=${LOCKS.a_nm} nm. Goldstone ${state.mode}  u=A ê sin(k·R−ωt)  A=${fmt(state.A, 2)} nm  ω=${fmt(omega, 4)}  |u|max=${fmt(uMax, 3)} nm. TA is helicity ±1, not a graviton.`;
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => applyPreset('tdvt') },
      precision: { label: 'Precision', apply: () => applyPreset('precision') },
    },
    reset() { lastT = 0; applyField(0); },
    update(_dt, t, { paused } = {}) {
      lastT = t;
      if (paused) return;
      applyField(t);
    },
    dispose() {},
  };
}
