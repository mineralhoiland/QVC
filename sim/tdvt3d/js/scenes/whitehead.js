// Tab 2 — Whitehead hopfion vortex core (L2). n-field, ρ_H cloud, preimages, torus (8,50) nm.

import { LOCKS, fmt } from '../math/locks.js';
import { nColor } from '../math/hopf.js';
import {
  computeHopfion, richardsonQ, visualizationSamples, preimageLinking,
  torusThetaGrid, WARD_VK, WARD_Q1_MIN,
} from '../math/whitehead.js';
import { makePointCloud, makeLineLoop, makeNmAxes, makeTorusWire, disposeObject, setPointStyle } from '../gfx/common.js';

export function createScene({ THREE }) {
  const scene = new THREE.Scene();
  scene.background = new THREE.Color(0x05060a);
  scene.add(makeNmAxes(THREE, 90));
  scene.add(makeTorusWire(THREE, LOCKS.R_nm, LOCKS.r_nm, 0x4a6a88));

  const state = {
    showGlyphs: true,
    showPre: true,
    breathe: true,
    advect: true,
    isoFrac: 0.18,
    size: 2.6,
    speed: 1,
    glow: 1,
    nAdv: 2400,
  };

  const g16 = computeHopfion(16);
  const g24 = computeHopfion(24);
  const rich = richardsonQ(g16.wh.Q, g16.wh.dx, g24.wh.Q, g24.wh.dx);
  const vk = g24.vk;
  const flux = torusThetaGrid(64, LOCKS.R_nm, LOCKS.r_nm, LOCKS.Q_H).flux;
  const pre = preimageLinking(g24.tex);

  let cloud = null;
  let nCloud = 0;
  const glyphGroup = new THREE.Group();
  scene.add(glyphGroup);

  const gDir = new THREE.Vector3();
  const gQ = new THREE.Quaternion();
  const yUp = new THREE.Vector3(0, 1, 0);
  const dummy = new THREE.Object3D();
  const rgb = [0, 0, 0];
  let nGlyphs = 0;
  let arrows = null;

  function rebuildCloud() {
    if (cloud) { scene.remove(cloud); disposeObject(cloud); cloud = null; }
    const vis = visualizationSamples(g24, { glyphStride: 3, cloudFrac: state.isoFrac });
    nCloud = vis.cloud.length / 7;
    const cpos = new Float32Array(nCloud * 3);
    const ccol = new Float32Array(nCloud * 3);
    const csz = new Float32Array(nCloud);
    for (let i = 0; i < nCloud; i++) {
      cpos[i * 3] = vis.cloud[i * 7];
      cpos[i * 3 + 1] = vis.cloud[i * 7 + 1];
      cpos[i * 3 + 2] = vis.cloud[i * 7 + 2];
      ccol[i * 3] = vis.cloud[i * 7 + 3];
      ccol[i * 3 + 1] = vis.cloud[i * 7 + 4];
      ccol[i * 3 + 2] = vis.cloud[i * 7 + 5];
      csz[i] = 2.2 + 6 * Math.abs(vis.cloud[i * 7 + 6]) / (g24.wh.rhoMax || 1e-12);
    }
    cloud = makePointCloud(THREE, cpos, ccol, { size: state.size, sizes: csz, opacity: state.glow });
    scene.add(cloud);

    if (arrows) { glyphGroup.remove(arrows); disposeObject(arrows); arrows = null; }
    nGlyphs = vis.glyphs.length;
    const arrowGeo = new THREE.ConeGeometry(1.4, 5.5, 5);
    const arrowMat = new THREE.MeshBasicMaterial({ color: 0xaad4ff, transparent: true, opacity: 0.35 + 0.4 * state.glow });
    arrows = new THREE.InstancedMesh(arrowGeo, arrowMat, nGlyphs);
    for (let i = 0; i < nGlyphs; i++) {
      const g = vis.glyphs[i];
      dummy.position.set(g.p[0], g.p[1], g.p[2]);
      gDir.set(g.n[0], g.n[1], g.n[2]).normalize();
      dummy.quaternion.copy(gQ.setFromUnitVectors(yUp, gDir));
      dummy.updateMatrix();
      arrows.setMatrixAt(i, dummy.matrix);
      nColor(g.n[0], g.n[1], g.n[2], rgb);
      arrows.setColorAt(i, new THREE.Color(rgb[0], rgb[1], rgb[2]));
    }
    arrows.instanceMatrix.needsUpdate = true;
    if (arrows.instanceColor) arrows.instanceColor.needsUpdate = true;
    glyphGroup.add(arrows);
    glyphGroup.visible = state.showGlyphs;
  }

  rebuildCloud();

  let preA = null, preB = null;
  if (pre.a.n > 8) {
    preA = makeLineLoop(THREE, pre.a.packed, 0xffe08a);
    scene.add(preA);
  }
  if (pre.b.n > 8) {
    preB = makeLineLoop(THREE, pre.b.packed, 0x7ee0ff);
    scene.add(preB);
  }

  let apos = null, acol = null, aPhase = null, aRho = null, aZ = null, adv = null;

  function seedAdv() {
    const nAdv = Math.max(400, state.nAdv | 0);
    state.nAdv = nAdv;
    apos = new Float32Array(nAdv * 3);
    acol = new Float32Array(nAdv * 3);
    aPhase = new Float32Array(nAdv);
    aRho = new Float32Array(nAdv);
    aZ = new Float32Array(nAdv);
    for (let i = 0; i < nAdv; i++) {
      const phi = Math.random() * Math.PI * 2;
      const chi = Math.random() * Math.PI * 2;
      const rt = LOCKS.r_nm * (0.15 + 0.85 * Math.random());
      aPhase[i] = phi;
      aRho[i] = LOCKS.R_nm + rt * Math.cos(chi);
      aZ[i] = rt * Math.sin(chi);
      nColor(Math.cos(LOCKS.Q_H * phi), Math.sin(LOCKS.Q_H * phi), 0, rgb);
      acol[i * 3] = rgb[0]; acol[i * 3 + 1] = rgb[1]; acol[i * 3 + 2] = rgb[2];
      apos[i * 3] = aRho[i] * Math.cos(phi);
      apos[i * 3 + 1] = aRho[i] * Math.sin(phi);
      apos[i * 3 + 2] = aZ[i];
    }
  }

  function rebuildAdv() {
    if (adv) { scene.remove(adv); disposeObject(adv); adv = null; }
    seedAdv();
    adv = makePointCloud(THREE, apos, acol, { size: Math.max(1.4, state.size * 0.75), opacity: state.glow });
    adv.visible = state.advect;
    scene.add(adv);
  }
  rebuildAdv();

  function setPre(on) {
    if (preA) preA.visible = on;
    if (preB) preB.visible = on;
  }
  setPre(state.showPre);

  let lastT = 0;
  let breatheScale = 1;

  function applyStyle() {
    if (cloud) setPointStyle(cloud, { size: state.size, opacity: state.glow });
    if (adv) setPointStyle(adv, { size: Math.max(1.4, state.size * 0.75), opacity: state.glow });
    if (arrows) arrows.material.opacity = 0.35 + 0.4 * state.glow;
    glyphGroup.visible = state.showGlyphs;
    if (adv) adv.visible = state.advect;
    setPre(state.showPre);
  }

  function stampAdv(dt) {
    if (!state.advect || !adv) return;
    const nAdv = state.nAdv;
    const w = LOCKS.Q_H * 0.55 * state.speed;
    for (let i = 0; i < nAdv; i++) {
      aPhase[i] += w * dt / Math.max(0.2, aRho[i] / LOCKS.R_nm);
      const phi = aPhase[i];
      apos[i * 3] = aRho[i] * Math.cos(phi);
      apos[i * 3 + 1] = aRho[i] * Math.sin(phi);
      apos[i * 3 + 2] = aZ[i];
      nColor(Math.cos(LOCKS.Q_H * phi), Math.sin(LOCKS.Q_H * phi), 0, rgb);
      acol[i * 3] = rgb[0]; acol[i * 3 + 1] = rgb[1]; acol[i * 3 + 2] = rgb[2];
    }
    adv.geometry.attributes.position.needsUpdate = true;
    adv.geometry.attributes.color.needsUpdate = true;
  }

  function applyPreset(id) {
    if (id === 'tdvt') {
      Object.assign(state, {
        showGlyphs: true, showPre: true, breathe: true, advect: true,
        isoFrac: 0.18, size: 2.6, speed: 1, glow: 1, nAdv: 2400,
      });
    } else if (id === 'precision') {
      Object.assign(state, { isoFrac: 0.12, size: 2.0, speed: 0.7, glow: 0.9, nAdv: 3600, breathe: false });
    }
    rebuildCloud();
    rebuildAdv();
    applyStyle();
  }

  return {
    scene,
    viewScale: 140,
    get particleCount() { return nCloud + state.nAdv + nGlyphs; },
    sliders() {
      return [
        { id: 'cf', label: 'cloud frac', min: 0.06, max: 0.45, step: 0.01, value: state.isoFrac, fmt: (v) => v.toFixed(2), oninput: (v) => { state.isoFrac = v; rebuildCloud(); applyStyle(); } },
        { id: 'na', label: 'advectors', min: 800, max: 8000, step: 200, value: state.nAdv, fmt: (v) => String(v | 0), oninput: (v) => { state.nAdv = v | 0; rebuildAdv(); } },
        { id: 'gl', label: 'n̂ glyphs', type: 'checkbox', value: state.showGlyphs, oninput: (v) => { state.showGlyphs = v; glyphGroup.visible = v; } },
        { id: 'pr', label: 'preimages', type: 'checkbox', value: state.showPre, oninput: (v) => { state.showPre = v; setPre(v); } },
        { id: 'ad', label: 'advection v∝∇θ', type: 'checkbox', value: state.advect, oninput: (v) => { state.advect = v; if (adv) adv.visible = v; } },
        { id: 'br', label: 'breathing (artistic)', type: 'checkbox', value: state.breathe, oninput: (v) => { state.breathe = v; if (!v && cloud) cloud.scale.setScalar(1); } },
        { id: 'size', label: 'point size', min: 0.8, max: 7, step: 0.1, value: state.size, fmt: (v) => v.toFixed(1), oninput: (v) => { state.size = v; applyStyle(); } },
        { id: 'glow', label: 'glow', min: 0.15, max: 1.4, step: 0.05, value: state.glow, fmt: (v) => v.toFixed(2), oninput: (v) => { state.glow = v; applyStyle(); } },
        { id: 'speed', label: 'speed', min: 0.1, max: 4, step: 0.05, value: state.speed, fmt: (v) => v.toFixed(2), oninput: (v) => { state.speed = v; } },
      ];
    },
    numbers() {
      return [
        { label: 'Whitehead Q 24³ (live)', value: fmt(g24.wh.Q, 4), kind: 'live' },
        { label: 'Whitehead Q 16³ (live)', value: fmt(g16.wh.Q, 4), kind: 'live' },
        { label: 'Richardson 16/24', value: fmt(rich, 4), kind: 'live' },
        { label: '|n|−1 max', value: g24.tex.normDev.toExponential(2), kind: 'live' },
        { label: 'FN Ward E 24³', value: fmt(g24.fn.E, 4), kind: 'live' },
        { label: `VK 0.534|Q|^{3/4}`, value: fmt(vk.bound, 4), kind: 'live' },
        { label: 'relaxed Q=1 ref', value: String(WARD_Q1_MIN), kind: 'lock' },
        { label: 'VK C', value: fmt(WARD_VK, 4), kind: 'lock' },
        { label: 'grid preimage lk', value: Number.isFinite(pre.Lk) ? fmt(pre.Lk, 3) : 'n/a', kind: 'live' },
        { label: 'CS-flux Q_H', value: fmt(flux.Q_H_flux, 4), kind: 'live' },
        { label: 'Q_H lock', value: fmt(LOCKS.Q_H, 3), kind: 'lock' },
        { label: '(r, R)', value: `(${LOCKS.r_nm}, ${LOCKS.R_nm}) nm`, kind: 'lock' },
        { label: 'fibre-phase / breathe', value: 'artistic', kind: 'artistic' },
      ];
    },
    charts() {
      return [
        { label: 'breathe', value: breatheScale, color: '#d48cff', unit: '' },
        { label: 'Q24', value: g24.wh.Q, color: '#7ee0ff', unit: '' },
      ];
    },
    banner() {
      return 'L2 · Stereographic Hopf ansatz (unrelaxed). Whitehead Q = (1/16π²)∫A·B, distinct from Q_H=1/N. Advection v∝∇θ with θ=Q_H φ is kinematics. Breathing / fibre-phase pulse is artistic, not a dynamical solver.';
    },
    presets: {
      tdvt: { label: 'TDVT', apply: () => applyPreset('tdvt') },
      precision: { label: 'Precision', apply: () => applyPreset('precision') },
    },
    reset() {
      lastT = 0;
      rebuildAdv();
      if (cloud) cloud.scale.setScalar(1);
      breatheScale = 1;
    },
    update(dt, t, { paused } = {}) {
      lastT = t;
      if (paused) return;
      if (state.breathe && cloud) {
        breatheScale = 1 + 0.035 * Math.sin(t * 1.4 * state.speed);
        cloud.scale.setScalar(breatheScale);
      } else {
        breatheScale = 1;
        if (cloud) cloud.scale.setScalar(1);
      }
      stampAdv(dt);
    },
    dispose() {},
  };
}
