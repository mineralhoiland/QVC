// Tab 6 — Fold drive λ(t), δ(t), Mexican-hat particles, moving operating-point marker.

import { LOCKS, fmt } from '../math/locks.js';
import { foldCurve, specG, mexicanHat, SYK_NOTE, solveFoldBranches } from '../math/spectra.js';
import { makePointCloud, makeNmAxes, makeLine, disposeObject, setPointStyle } from '../gfx/common.js';

const LAMBDA_STAR = LOCKS.lambda_star;
const SCALE = 36;

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.add(makeNmAxes(THREE, 50));

  const state = {
    lambda0: LAMBDA_STAR,
    amp: 0.06,
    omega: 0.55,
    N: LOCKS.N,
    size: 2.2,
    speed: 1,
    glow: 1,
    g0: LOCKS.g0_star_meV,
    nHat: 48 * 72,
  };

  let hat = null;
  let fountain = null;
  let foldLine = null;
  let marker = null;
  let lambda = state.lambda0;
  let delta = LOCKS.delta_star;
  let Delta = LOCKS.Delta_star_meV;
  let foldOut = false;
  let fountPhase = null;
  let fountRho = null;
  let lastT = 0;

  function foldPacked() {
    const fc = foldCurve(LOCKS.ell, 140);
    const pts = [];
    for (let i = 0; i < fc.phys.length; i++) {
      if (!Number.isFinite(fc.phys[i])) continue;
      pts.push((fc.alphaOver[i] - 0.5) * SCALE * 1.6, fc.phys[i] * SCALE, 8);
    }
    for (let i = 0; i < fc.unst.length; i++) {
      if (!Number.isFinite(fc.unst[i])) continue;
      pts.push((fc.alphaOver[i] - 0.5) * SCALE * 1.6, fc.unst[i] * SCALE, 4);
    }
    return new Float32Array(pts);
  }

  function hatRes() {
    const nHat = Math.max(800, state.nHat | 0);
    state.nHat = nHat;
    const nR = Math.max(24, Math.round(Math.sqrt(nHat / 1.5)));
    const nPhi = Math.max(24, Math.round(nHat / nR));
    return { nR, nPhi };
  }

  function rebuildHat(DeltaMeV) {
    if (hat) { scene.remove(hat); disposeObject(hat); }
    if (fountain) { scene.remove(fountain); disposeObject(fountain); }
    const { nR, nPhi } = hatRes();
    const mesh = mexicanHat(DeltaMeV, { nR, nPhi });
    const n = mesh.pos.length / 3;
    const pos = new Float32Array(n * 3);
    for (let i = 0; i < n; i++) {
      pos[i * 3] = mesh.pos[i * 3] * SCALE * 0.45;
      pos[i * 3 + 1] = mesh.pos[i * 3 + 1] * SCALE * 0.45;
      pos[i * 3 + 2] = mesh.pos[i * 3 + 2] * SCALE * 2.2;
    }
    hat = makePointCloud(THREE, pos, mesh.col, { size: state.size * 0.7, opacity: state.glow });
    scene.add(hat);

    const nf = Math.max(800, Math.round(state.nHat * 0.7));
    const fp = new Float32Array(nf * 3);
    const fc = new Float32Array(nf * 3);
    fountPhase = new Float32Array(nf);
    fountRho = new Float32Array(nf);
    const rhoWell = mesh.rhoMin || 0.4;
    for (let i = 0; i < nf; i++) {
      fountPhase[i] = Math.random() * Math.PI * 2;
      fountRho[i] = rhoWell * (0.7 + 0.5 * Math.random());
      fc[i * 3] = 0.95; fc[i * 3 + 1] = 0.85; fc[i * 3 + 2] = 0.35;
    }
    fountain = makePointCloud(THREE, fp, fc, { size: Math.max(1.4, state.size), opacity: state.glow });
    scene.add(fountain);
  }

  function rebuildStatic() {
    if (foldLine) { scene.remove(foldLine); disposeObject(foldLine); }
    if (marker) { scene.remove(marker); disposeObject(marker); }
    foldLine = makeLine(THREE, foldPacked(), 0x7ee0ff);
    scene.add(foldLine);
    marker = new THREE.Mesh(
      new THREE.SphereGeometry(2.4, 14, 14),
      new THREE.MeshBasicMaterial({ color: 0xffe08a }),
    );
    scene.add(marker);
    rebuildHat(LOCKS.Delta_star_meV);
    drive(lastT);
  }

  function drive(t) {
    lambda = state.lambda0 + state.amp * Math.sin(state.omega * state.speed * t);
    const alpha = LOCKS.alpha_star * (lambda / LAMBDA_STAR);
    const br = solveFoldBranches(alpha, LOCKS.ell);
    foldOut = br.phys == null;
    delta = br.phys == null ? (br.delta_c || LOCKS.delta_c) : br.phys;
    Delta = LOCKS.Delta0_meV * delta;

    const aOver = alpha / LOCKS.alpha_c;
    marker.position.set((aOver - 0.5) * SCALE * 1.6, delta * SCALE, 10);
    marker.material.color.set(foldOut ? 0xff6b7a : 0xffe08a);
    marker.scale.setScalar(1 + 0.12 * Math.sin(t * 4 * state.speed));

    if (fountain && fountPhase) {
      const pos = fountain.geometry.getAttribute('position');
      const rho0 = Math.sqrt(Math.max(0, Delta / LOCKS.Delta0_meV)) * 0.9 || 0.35;
      for (let i = 0; i < fountPhase.length; i++) {
        const ph = fountPhase[i] + t * state.speed * 0.7;
        const rho = (fountRho[i] * 0.35 + rho0 * 0.65) * SCALE * 0.45;
        pos.setXYZ(i, rho * Math.cos(ph), rho * Math.sin(ph), -Delta * SCALE * 0.15);
      }
      pos.needsUpdate = true;
    }
  }

  function applyStyle() {
    if (hat) setPointStyle(hat, { size: state.size * 0.7, opacity: state.glow });
    if (fountain) setPointStyle(fountain, { size: Math.max(1.4, state.size), opacity: state.glow });
  }

  function applyPreset(id) {
    if (id === 'tdvt') {
      Object.assign(state, {
        lambda0: LAMBDA_STAR, N: LOCKS.N, amp: 0.06, omega: 0.55, size: 2.2, speed: 1,
        glow: 1, g0: LOCKS.g0_star_meV, nHat: 48 * 72,
      });
    } else if (id === 'precision') {
      Object.assign(state, { lambda0: LAMBDA_STAR, N: LOCKS.N, amp: 0.02, omega: 0.3, size: 1.8, glow: 0.9, nHat: 64 * 80 });
    }
    rebuildStatic();
  }

  rebuildStatic();

  return {
    scene,
    viewScale: 70,
    get particleCount() {
      const h = hat?.geometry?.getAttribute('position')?.count || 0;
      const f = fountain?.geometry?.getAttribute('position')?.count || 0;
      return h + f;
    },
    sliders() {
      return [
        { id: 'lam', label: 'λ drive', min: 0.45, max: 1.05, step: 0.0001, value: state.lambda0, fmt: (v) => v.toFixed(4), oninput: (v) => { state.lambda0 = v; drive(lastT); } },
        { id: 'amp', label: 'amp', min: 0, max: 0.2, step: 0.002, value: state.amp, fmt: (v) => v.toFixed(3), oninput: (v) => { state.amp = v; drive(lastT); } },
        { id: 'om', label: 'Ω', min: 0.05, max: 3, step: 0.05, value: state.omega, fmt: (v) => v.toFixed(2), oninput: (v) => { state.omega = v; drive(lastT); } },
        { id: 'N', label: 'N', min: 3, max: 9, step: 1, value: state.N, fmt: (v) => String(v | 0), oninput: (v) => { state.N = v | 0; } },
        { id: 'g0', label: 'g₀ (diag)', min: 0.05, max: 0.2, step: 0.001, value: state.g0, fmt: (v) => v.toFixed(4), oninput: (v) => { state.g0 = v; } },
        { id: 'nh', label: 'n̂ particles', min: 1500, max: 12000, step: 300, value: state.nHat, fmt: (v) => String(v | 0), oninput: (v) => { state.nHat = v | 0; rebuildHat(Delta); drive(lastT); } },
        { id: 'size', label: 'point size', min: 0.6, max: 6, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; applyStyle(); } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; applyStyle(); } },
        { id: 'speed', label: 'speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; drive(lastT); } },
      ];
    },
    numbers() {
      const g = specG(state.N, state.g0);
      return [
        { label: 'λ(t)', value: fmt(lambda, 4), kind: 'live' },
        { label: 'λ*', value: fmt(LAMBDA_STAR, 4), kind: 'lock' },
        { label: 'δ(t)', value: foldOut ? 'fold-out' : fmt(delta, 4), kind: foldOut ? 'artistic' : 'live' },
        { label: 'Δ(t)', value: `${fmt(Delta, 4)} meV`, kind: 'live' },
        { label: 'Δ*', value: `${fmt(LOCKS.Delta_star_meV, 4)} meV`, kind: 'lock' },
        { label: 'λ_max Spec(G)', value: `${fmt(g.lambda_max, 4)} meV`, kind: 'diag' },
        { label: 'g₀*', value: `${fmt(LOCKS.g0_star_meV, 4)} meV`, kind: 'lock' },
        { label: 'T*', value: `${LOCKS.T_star_mK} mK`, kind: 'lock' },
        { label: 'a = ξ', value: `${LOCKS.a_nm} nm`, kind: 'lock' },
        { label: '(r, R)', value: `(${LOCKS.r_nm}, ${LOCKS.R_nm}) nm`, kind: 'lock' },
        { label: 'SYK', value: 'analogy only', kind: 'artistic' },
      ];
    },
    charts() {
      return [
        { label: 'δ', value: delta, color: '#d48cff', unit: '' },
        { label: 'Δ', value: Delta, color: '#ffe08a', unit: 'meV' },
      ];
    },
    banner() {
      const flash = foldOut ? '  FOLD-OUT (α ≥ α_c).' : '';
      return `L1 · Drive λ(t) about λ*=${fmt(LAMBDA_STAR, 4)}, N=${state.N}. δ(t) from the fold; Δ=Δ₀ δ. Hat particles are an artistic fountain in V(ρ).${flash} ${SYK_NOTE}`;
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => applyPreset('tdvt') },
      precision: { label: 'Precision', apply: () => applyPreset('precision') },
    },
    reset() { lastT = 0; drive(0); },
    update(_dt, t, { paused } = {}) {
      lastT = t;
      if (paused) return;
      drive(t);
    },
    dispose() {},
  };
}
