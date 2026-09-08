// Tab 4 — Acoustic horizon (L2). Unruh–Visser g_tt = −(c_s²−v²). Drain slider. ≥100k tracers.

import { LOCKS, fmt } from '../math/locks.js';
import { acousticMetric, hEin1Overlap } from '../math/analogue.js';
import { hsvToRgb } from '../math/hopf.js';
import { makePointCloud, makeNmAxes, disposeObject, setPointStyle } from '../gfx/common.js';

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.fog = new THREE.FogExp2(0x05060a, 0.006);
  scene.add(makeNmAxes(THREE, 80));

  const state = {
    drain: 1.2,
    n: 120000,
    size: 1.8,
    speed: 1,
    noise: 0.15,
    glow: 1,
  };
  let metric = acousticMetric({ drainV0OverCs: state.drain });

  let pos = null, col = null, phi = null, rho = null, z0 = null;
  let tracers = null;
  let meanRho = 0;
  let meanMach = 0;
  const rgb = [0, 0, 0];
  const ein = hEin1Overlap();

  const horizonRing = new THREE.Mesh(
    new THREE.TorusGeometry(60, 1.2, 8, 96),
    new THREE.MeshBasicMaterial({ color: 0xffc857, transparent: true, opacity: 0.9 }),
  );
  horizonRing.rotation.x = Math.PI / 2;
  scene.add(horizonRing);

  const coreRing = new THREE.Mesh(
    new THREE.TorusGeometry(LOCKS.r_nm, 0.6, 8, 64),
    new THREE.MeshBasicMaterial({ color: 0x6688aa, transparent: true, opacity: 0.5 }),
  );
  coreRing.rotation.x = Math.PI / 2;
  scene.add(coreRing);

  function alloc() {
    const n = Math.max(20000, state.n | 0);
    state.n = n;
    pos = new Float32Array(n * 3);
    col = new Float32Array(n * 3);
    phi = new Float32Array(n);
    rho = new Float32Array(n);
    z0 = new Float32Array(n);
  }

  function seed() {
    const n = state.n;
    for (let i = 0; i < n; i++) {
      rho[i] = 8 + Math.pow(Math.random(), 0.55) * 110;
      phi[i] = Math.random() * Math.PI * 2;
      z0[i] = (Math.random() - 0.5) * 14;
    }
  }

  function paint(t = 0) {
    const cs = metric.cs;
    const hm = metric.hm;
    const Rd = LOCKS.R_nm * 1e-3;
    const n = state.n;
    let sRho = 0, sMach = 0;
    for (let i = 0; i < n; i++) {
      const rnm = rho[i];
      const rum = Math.max(rnm * 1e-3, 1e-5);
      const vphi = hm * LOCKS.Q_H / rum;
      const vr = -state.drain * cs * (Rd / rum);
      const speed = Math.hypot(vphi, vr);
      const mach = speed / cs;
      const z = z0[i] + state.noise * 5 * Math.sin(phi[i] * 2.7 + i * 0.013 + t);
      pos[i * 3] = rnm * Math.cos(phi[i]);
      pos[i * 3 + 1] = rnm * Math.sin(phi[i]);
      pos[i * 3 + 2] = z;
      const hue = mach > 1 ? 0.02 : 0.55 - 0.25 * Math.min(1, mach);
      hsvToRgb(hue, 0.75, 0.35 + 0.55 * Math.min(1.2, mach), rgb);
      col[i * 3] = rgb[0]; col[i * 3 + 1] = rgb[1]; col[i * 3 + 2] = rgb[2];
      sRho += rnm;
      sMach += mach;
    }
    meanRho = sRho / n;
    meanMach = sMach / n;
  }

  function rebuild() {
    if (tracers) { scene.remove(tracers); disposeObject(tracers); tracers = null; }
    alloc();
    seed();
    paint(0);
    tracers = makePointCloud(THREE, pos, col, { size: state.size, opacity: state.glow });
    scene.add(tracers);
    applyDrain(state.drain, true);
  }

  function applyDrain(d, skipPaint = false) {
    state.drain = d;
    metric = acousticMetric({ drainV0OverCs: d });
    const rh = Number.isFinite(metric.rho_h_nm) ? metric.rho_h_nm : 0.5;
    horizonRing.scale.setScalar(Math.max(0.02, rh / 60));
    horizonRing.visible = d > 0 && metric.exteriorHorizon;
    if (!skipPaint) {
      paint(0);
      if (tracers) {
        tracers.geometry.attributes.position.needsUpdate = true;
        tracers.geometry.attributes.color.needsUpdate = true;
      }
    }
  }

  rebuild();

  return {
    scene,
    viewScale: 130,
    get particleCount() { return state.n; },
    sliders() {
      return [
        { id: 'dr', label: 'drain v₀/c_s', min: 0, max: 2.0, step: 0.05, value: state.drain, fmt: (v) => v.toFixed(2), oninput: (v) => applyDrain(v) },
        { id: 'n', label: 'n tracers', min: 40000, max: 180000, step: 10000, value: state.n, fmt: (v) => String(v | 0), oninput: (v) => { state.n = v | 0; rebuild(); } },
        { id: 'size', label: 'point size', min: 0.6, max: 5, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; setPointStyle(tracers, { size: v }); } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; setPointStyle(tracers, { opacity: v }); } },
        { id: 'speed', label: 'speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; } },
        { id: 'noise', label: 'noise', min: 0, max: 1.2, step: 0.02, value: state.noise, fmt: (v) => v.toFixed(2), oninput: (v) => { state.noise = v; } },
      ];
    },
    numbers() {
      const noH = state.drain <= 0 || !metric.exteriorHorizon;
      return [
        { label: 'ρ_h', value: Number.isFinite(metric.rho_h_nm) ? `${fmt(metric.rho_h_nm, 2)} nm` : 'nan', kind: 'live' },
        { label: 'κ', value: Number.isFinite(metric.kappa_ps_inv) ? `${fmt(metric.kappa_ps_inv, 3)} ps⁻¹` : 'nan', kind: 'live' },
        { label: 'T_H', value: Number.isFinite(metric.T_H_mK) ? `${fmt(metric.T_H_mK, 3)} mK` : 'nan', kind: 'diag' },
        { label: 'T* (Schwinger lock)', value: `${LOCKS.T_star_mK} mK`, kind: 'lock' },
        { label: '⟨ρ⟩', value: `${fmt(meanRho, 1)} nm`, kind: 'live' },
        { label: '⟨Mach⟩', value: fmt(meanMach, 3), kind: 'live' },
        { label: 'horizon', value: noH ? (state.drain <= 0 ? 'none (drain=0)' : 'inside core — no exterior horizon') : 'sonic surface |v|=c_s', kind: 'live' },
        { label: 'g_tt', value: '−(c_s² − v²)', kind: 'live' },
        { label: 'H-EIN1 overlap', value: ein.overlap_plus.toExponential(1), kind: 'live' },
        { label: 'note', value: 'T_H is diagnostic, not T*', kind: 'diag' },
      ];
    },
    charts() {
      return [
        { label: '⟨ρ⟩', value: meanRho, color: '#7ee0ff', unit: 'nm' },
        { label: '⟨Mach⟩', value: meanMach, color: '#ffc857', unit: '' },
      ];
    },
    banner() {
      return state.drain <= 0
        ? 'L2 · drain = 0: no exterior acoustic horizon. Vortex-only sonic radius sits inside the 8 nm healing core. T_H is Unruh–Visser kinematics, not the locked T*=118.7 mK.'
        : 'L2 · Unruh–Visser acoustic metric. Gold ring = sonic surface |v|=c_s when an exterior horizon exists. T_H vs T*=118.7 mK is a diagnostic (H-HAWK1), not a lock.';
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => { Object.assign(state, { drain: 1.2, n: 120000, size: 1.8, speed: 1, noise: 0.15, glow: 1 }); rebuild(); } },
      precision: { label: 'Precision', apply: () => { Object.assign(state, { drain: 1.0, n: 160000, size: 1.4, speed: 0.7, noise: 0.08, glow: 0.9 }); rebuild(); } },
    },
    reset() { seed(); paint(0); if (tracers) { tracers.geometry.attributes.position.needsUpdate = true; tracers.geometry.attributes.color.needsUpdate = true; } },
    update(dt, t, { paused } = {}) {
      if (paused || !tracers) return;
      const cs = metric.cs, hm = metric.hm;
      const Rd = LOCKS.R_nm * 1e-3;
      const n = state.n;
      const sp = state.speed;
      let sRho = 0, sMach = 0;
      for (let i = 0; i < n; i++) {
        const rum = Math.max(rho[i] * 1e-3, 1e-5);
        const vphi = hm * LOCKS.Q_H / rum;
        const vr = -state.drain * cs * (Rd / rum);
        phi[i] += (vphi / rum) * dt * 0.35 * sp;
        if (state.noise > 0) phi[i] += (Math.random() - 0.5) * state.noise * dt * 0.9;
        rho[i] = Math.max(4, rho[i] + vr * 1e3 * dt * 0.08 * sp);
        if (rho[i] > 130) rho[i] = 8 + Math.random() * 20;
        const z = z0[i] + state.noise * 5 * Math.sin(phi[i] * 2.7 + i * 0.013 + t);
        pos[i * 3] = rho[i] * Math.cos(phi[i]);
        pos[i * 3 + 1] = rho[i] * Math.sin(phi[i]);
        pos[i * 3 + 2] = z;
        const speedV = Math.hypot(vphi, vr);
        const mach = speedV / cs;
        sRho += rho[i];
        sMach += mach;
      }
      meanRho = sRho / n;
      meanMach = sMach / n;
      tracers.geometry.attributes.position.needsUpdate = true;
    },
    dispose() {},
  };
}
