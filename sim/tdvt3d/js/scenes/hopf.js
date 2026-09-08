// Tab 1 — Hopf fibration / fibre bundle (L2). ≥100k points, color by S² base, Gauss linking.

import { LOCKS, fmt, lensLocks } from '../math/locks.js';
import { generateFibres, twoFibreLinking, hopfFibrePoint, stereoFromS3, s3EmbedXYZ, nColor } from '../math/hopf.js';
import { HOPF_FLOW_VS, POINT_FS, makeLineLoop, makeNmAxes, makePointCloud, disposeObject } from '../gfx/common.js';

const ETA1 = 0.45, ETA2 = 1.05, PSI1 = 0.3, PSI2 = 2.6;

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.fog = new THREE.FogExp2(0x05060a, 0.004);
  scene.add(makeNmAxes(THREE, 50));

  const state = {
    nFibres: 220,
    nSeg: 512,
    stereo: true,
    branchCut: false,
    scale: 40,
    N: LOCKS.N,
    speed: 1,
    size: 2.1,
    glow: 1,
    showFibres: true,
  };

  const uniforms = {
    uTime: { value: 0 },
    uOmega: { value: 1 },
    uStereo: { value: 1 },
    uScale: { value: state.scale },
    uSize: { value: state.size },
    uOpacity: { value: state.glow },
  };

  let cloud = null;
  let fibreA = null, fibreB = null;
  let cutMarks = null;
  let particleCount = 0;
  let Lk = NaN;
  let lastT = 0;

  function writeLoop(line, eta, psi, phase) {
    if (!line) return;
    const pos = line.geometry.getAttribute('position');
    const n = pos.count;
    const project = state.stereo ? stereoFromS3 : s3EmbedXYZ;
    for (let s = 0; s < n; s++) {
      const h = hopfFibrePoint(eta, psi, Math.PI * 2 * s / n + phase);
      const r = project(h.w0, h.w1, h.w2, h.w3);
      pos.setXYZ(s, r.x * state.scale, r.y * state.scale, r.z * state.scale);
    }
    pos.needsUpdate = true;
  }

  function stampHighlights(t) {
    const phase = t * state.speed;
    writeLoop(fibreA, ETA1, PSI1, phase);
    writeLoop(fibreB, ETA2, PSI2, phase);
    if (cutMarks && fibreA) {
      const src = fibreA.geometry.getAttribute('position');
      const dst = cutMarks.geometry.getAttribute('position');
      const nCut = dst.count;
      const nSrc = src.count;
      for (let i = 0; i < nCut; i++) {
        const idx = Math.floor((i / nCut) * (nSrc - 1));
        dst.setXYZ(i, src.getX(idx), src.getY(idx), src.getZ(idx));
      }
      dst.needsUpdate = true;
    }
  }

  function applyUniforms() {
    uniforms.uStereo.value = state.stereo ? 1 : 0;
    uniforms.uScale.value = state.scale;
    uniforms.uSize.value = state.size;
    uniforms.uOpacity.value = state.glow;
    uniforms.uTime.value = lastT * state.speed;
    if (fibreA) fibreA.visible = state.showFibres;
    if (fibreB) fibreB.visible = state.showFibres;
  }

  function rebuild() {
    if (cloud) { scene.remove(cloud); disposeObject(cloud); cloud = null; }
    if (fibreA) { scene.remove(fibreA); disposeObject(fibreA); fibreA = null; }
    if (fibreB) { scene.remove(fibreB); disposeObject(fibreB); fibreB = null; }
    if (cutMarks) { scene.remove(cutMarks); disposeObject(cutMarks); cutMarks = null; }

    const pack = generateFibres({
      nFibres: state.nFibres,
      nSeg: state.nSeg,
      stereo: state.stereo,
      scale: state.scale,
    });
    const geo = new THREE.BufferGeometry();
    geo.setAttribute('position', new THREE.BufferAttribute(pack.pos, 3));
    geo.setAttribute('aEta', new THREE.BufferAttribute(pack.eta, 1));
    geo.setAttribute('aPsi', new THREE.BufferAttribute(pack.psi, 1));
    geo.setAttribute('aT0', new THREE.BufferAttribute(pack.t0, 1));
    const mat = new THREE.ShaderMaterial({
      uniforms,
      vertexShader: HOPF_FLOW_VS,
      fragmentShader: POINT_FS,
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
    });
    cloud = new THREE.Points(geo, mat);
    cloud.frustumCulled = false;
    scene.add(cloud);
    particleCount = pack.nPts;

    const link = twoFibreLinking({
      nSeg: 288, stereo: state.stereo, scale: state.scale,
      eta1: ETA1, eta2: ETA2, psi1: PSI1, psi2: PSI2,
    });
    Lk = link.Lk;
    fibreA = makeLineLoop(THREE, link.c1, 0xffe08a);
    fibreB = makeLineLoop(THREE, link.c2, 0x7ee0ff);
    fibreA.visible = state.showFibres;
    fibreB.visible = state.showFibres;
    scene.add(fibreA, fibreB);

    if (state.branchCut) {
      const nCut = state.N;
      const cpos = new Float32Array(nCut * 3);
      const ccol = new Float32Array(nCut * 3);
      const rgb = [0, 0, 0];
      for (let i = 0; i < nCut; i++) {
        const tt = (i / nCut) * Math.PI * 2;
        nColor(Math.cos(tt), Math.sin(tt), 0, rgb);
        ccol[i * 3] = rgb[0]; ccol[i * 3 + 1] = rgb[1]; ccol[i * 3 + 2] = rgb[2];
      }
      cutMarks = makePointCloud(THREE, cpos, ccol, { size: 8, opacity: state.glow });
      scene.add(cutMarks);
    }
    applyUniforms();
    stampHighlights(lastT);
  }

  function applyPreset(id) {
    if (id === 'tdvt') {
      Object.assign(state, {
        nFibres: 220, nSeg: 512, stereo: true, branchCut: false, scale: 40,
        N: LOCKS.N, speed: 1, size: 2.1, glow: 1, showFibres: true,
      });
    } else if (id === 'precision') {
      Object.assign(state, { nFibres: 300, nSeg: 640, size: 1.8, glow: 0.9, speed: 0.8 });
    }
    rebuild();
  }

  rebuild();

  return {
    scene,
    viewScale: 80,
    get particleCount() { return particleCount; },
    sliders() {
      return [
        { id: 'nf', label: 'fibres', min: 80, max: 360, step: 10, value: state.nFibres, fmt: (v) => String(v | 0), oninput: (v) => { state.nFibres = v | 0; rebuild(); } },
        { id: 'ns', label: 'pts / fibre', min: 128, max: 768, step: 64, value: state.nSeg, fmt: (v) => String(v | 0), oninput: (v) => { state.nSeg = v | 0; rebuild(); } },
        { id: 'st', label: 'stereographic', type: 'checkbox', value: state.stereo, oninput: (v) => { state.stereo = v; applyUniforms(); stampHighlights(lastT); } },
        { id: 'sp', label: 'phase speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; applyUniforms(); } },
        { id: 'sc', label: 'scale', min: 16, max: 80, step: 1, value: state.scale, fmt: (v) => String(v | 0), oninput: (v) => { state.scale = v; applyUniforms(); stampHighlights(lastT); } },
        { id: 'N', label: 'N (Q_H=1/N)', min: 3, max: 9, step: 1, value: state.N, fmt: (v) => String(v | 0), oninput: (v) => { state.N = v | 0; if (state.branchCut) rebuild(); } },
        { id: 'bc', label: 'Q_H cut', type: 'checkbox', value: state.branchCut, oninput: (v) => { state.branchCut = v; rebuild(); } },
        { id: 'size', label: 'point size', min: 0.6, max: 8, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; uniforms.uSize.value = v; } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; uniforms.uOpacity.value = v; } },
        { id: 'fib', label: 'show fibres', type: 'checkbox', value: state.showFibres, oninput: (v) => { state.showFibres = v; applyUniforms(); } },
      ];
    },
    numbers() {
      const L = lensLocks(state.N);
      return [
        { label: 'Gauss linking (live)', value: fmt(Lk, 4), kind: 'live' },
        { label: 'expected lk', value: '±1', kind: 'lock' },
        { label: 'particles', value: particleCount.toLocaleString(), kind: 'live' },
        { label: 'Q_H = 1/N', value: fmt(L.Q_H, 3), kind: 'lock' },
        { label: 'phase', value: fmt((lastT * state.speed) % (Math.PI * 2), 3), kind: 'live' },
        { label: 'color', value: 'S² base (HSV)', kind: 'live' },
        { label: 'branch cut', value: state.branchCut ? `Z_${state.N} on fibre` : 'off', kind: 'artistic' },
      ];
    },
    charts() {
      const phase = (lastT * state.speed) % (Math.PI * 2);
      return [
        { label: 'phase', value: phase, color: '#7ee0ff', unit: '' },
        { label: 'Lk', value: Lk, color: '#ffe08a', unit: '' },
      ];
    },
    banner() {
      return `L2 · Hopf bundle π: S³ → S². Colour = base point on S². Two highlighted fibres: Gauss linking (live) → ±1. Phase speed ${fmt(state.speed, 2)}. Branch-cut toggle is the lens-space identification, not a Whitehead integral.`;
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => applyPreset('tdvt') },
      precision: { label: 'Precision', apply: () => applyPreset('precision') },
    },
    reset() { lastT = 0; uniforms.uTime.value = 0; stampHighlights(0); },
    update(_dt, t, { paused } = {}) {
      lastT = t;
      if (paused) return;
      uniforms.uTime.value = t * state.speed;
      stampHighlights(t);
    },
    dispose() {},
  };
}
