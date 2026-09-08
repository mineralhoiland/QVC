import * as THREE from 'three';
import { OrbitControls } from '../vendor/jsm/controls/OrbitControls.js';
import { EffectComposer } from '../vendor/jsm/postprocessing/EffectComposer.js';
import { RenderPass } from '../vendor/jsm/postprocessing/RenderPass.js';
import { UnrealBloomPass } from '../vendor/jsm/postprocessing/UnrealBloomPass.js';
import { createScene as hopf } from './scenes/hopf.js';
import { createScene as whitehead } from './scenes/whitehead.js';
import { createScene as lattice } from './scenes/lattice.js';
import { createScene as horizon } from './scenes/horizon.js';
import { createScene as torsion } from './scenes/torsion.js';
import { createScene as spectra } from './scenes/spectra.js';
import { createScene as hypersphere } from './scenes/hypersphere.js';
import { applyPreset, disposeObject } from './gfx/common.js';
import { createHUD, fmtSelftestHTML } from './ui/hud.js';
import { isSelftestQuery, runSelftest } from './math/verify.js';

const FACTORIES = { hopf, whitehead, lattice, horizon, torsion, spectra, hypersphere };
const HIST_MAX = 300;

const canvas = document.getElementById('c');
const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: false });
renderer.setPixelRatio(Math.min(2, window.devicePixelRatio || 1));
renderer.setSize(window.innerWidth, window.innerHeight);
renderer.outputColorSpace = THREE.SRGBColorSpace;

const camera = new THREE.PerspectiveCamera(50, window.innerWidth / window.innerHeight, 0.1, 4000);
const controls = new OrbitControls(camera, canvas);
controls.enableDamping = true;
controls.dampingFactor = 0.06;
controls.autoRotateSpeed = 1.4;

const composer = new EffectComposer(renderer);
const renderPass = new RenderPass(new THREE.Scene(), camera);
composer.addPass(renderPass);
const bloom = new UnrealBloomPass(
  new THREE.Vector2(window.innerWidth, window.innerHeight),
  0.4, 0.6, 0.8,
);
composer.addPass(bloom);

let current = null;
let currentId = null;
let autoOrbit = false;
let playing = true;
let simTime = 0;
let hudAcc = 0;
let frames = 0;
let fps = 0;
const history = new Map();
const clock = new THREE.Clock();
const _offset = new THREE.Vector3();
const _axis = new THREE.Vector3(0, 1, 0);

function bannerText(tab) {
  if (!tab) return '';
  return typeof tab.banner === 'function' ? tab.banner() : tab.banner;
}

function refreshHud() {
  if (!current) return;
  hud.setNumbers(current.numbers());
  hud.setBanner(bannerText(current));
  hud.setStats({
    fps,
    particles: current.particleCount,
    pixelRatio: renderer.getPixelRatio(),
  });
}

function bindSliders(tab) {
  const defs = tab.sliders().map((d) => {
    const orig = d.oninput;
    return {
      ...d,
      oninput: (v) => {
        orig?.(v);
        refreshHud();
      },
    };
  });
  hud.setSliders(defs);
}

function resetSim() {
  simTime = 0;
  history.clear();
  if (!current) return;
  current.reset?.();
  current.update(0, 0, { paused: !playing });
  refreshHud();
  pushCharts();
}

function bindSimPresets(tab) {
  const p = tab.presets || {};
  const items = [];
  for (const id of ['tdvt', 'precision']) {
    const slot = p[id];
    if (!slot) continue;
    items.push({
      id,
      label: slot.label || (id === 'tdvt' ? 'TDVT' : 'Precision'),
      onClick: () => {
        simTime = 0;
        slot.apply?.();
        current.reset?.();
        bindSliders(tab);
        history.clear();
        current.update(0, 0, { paused: !playing });
        refreshHud();
        pushCharts();
      },
    });
  }
  items.push({
    id: 'reset',
    label: 'Reset t',
    onClick: resetSim,
  });
  hud.setSimPresets(items);
}

function pushCharts() {
  if (!current?.charts) {
    hud.setCharts([]);
    return;
  }
  const series = current.charts() || [];
  const out = [];
  for (const s of series) {
    if (!history.has(s.label)) history.set(s.label, []);
    const arr = history.get(s.label);
    if (playing) {
      arr.push(s.value);
      if (arr.length > HIST_MAX) arr.shift();
    }
    out.push({ label: s.label, data: arr, color: s.color, unit: s.unit });
  }
  hud.setCharts(out);
}

function switchTab(id) {
  const factory = FACTORIES[id];
  if (!factory) return;

  if (current) {
    try { current.dispose(); } catch (_) { /* empty scene dispose */ }
    disposeObject(current.scene);
  }

  current = factory({ THREE });
  currentId = id;
  simTime = 0;
  renderPass.scene = current.scene;
  applyPreset(camera, controls, 'iso', current.viewScale);

  history.clear();
  hud.setTab(id);
  bindSliders(current);
  bindSimPresets(current);
  hud.setDots([]);
  refreshHud();
  pushCharts();
}

const hud = createHUD({
  onTab: switchTab,
  onPreset: (name) => {
    if (current) applyPreset(camera, controls, name, current.viewScale);
  },
  onOrbit: (on) => { autoOrbit = on; },
  onPlay: (on) => { playing = on; },
});

window.addEventListener('keydown', (e) => {
  if (e.code !== 'Space') return;
  const tag = e.target && e.target.tagName;
  if (tag === 'INPUT' || tag === 'TEXTAREA' || tag === 'SELECT' || e.target.isContentEditable) return;
  e.preventDefault();
  playing = !playing;
  hud.setPlaying(playing);
});

function onResize() {
  const w = window.innerWidth;
  const h = window.innerHeight;
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
  renderer.setSize(w, h);
  composer.setSize(w, h);
  bloom.setSize(w, h);
}

window.addEventListener('resize', onResize);

function tick() {
  requestAnimationFrame(tick);
  const dt = Math.min(clock.getDelta(), 0.05);
  if (document.hidden) return;

  if (playing) simTime += dt;
  if (current) current.update(dt, simTime, { paused: !playing });

  if (autoOrbit) {
    const target = controls.target;
    _offset.copy(camera.position).sub(target);
    _offset.applyAxisAngle(_axis, dt * 0.35);
    camera.position.copy(target).add(_offset);
    camera.lookAt(target);
  }
  controls.update();
  composer.render();

  if (current) pushCharts();

  frames += 1;
  hudAcc += dt;
  if (hudAcc >= 1) {
    fps = frames / hudAcc;
    frames = 0;
    hudAcc = 0;
    refreshHud();
  }
}

switchTab('hopf');
tick();

if (isSelftestQuery()) {
  hud.setSelftest(fmtSelftestHTML(null, 'Running live compute…'), true);
  runSelftest({
    onProgress: (_rows, msg) => {
      hud.setSelftest(fmtSelftestHTML(null, msg || 'Running live compute…'), true);
    },
  }).then((result) => {
    hud.setSelftest(fmtSelftestHTML(result), true);
    hud.setDots(result.rows.map((r) => ({ pass: r.pass, label: r.label })));
  }).catch((err) => {
    hud.setSelftest(fmtSelftestHTML(null, `Self-test failed: ${err.message || err}`), true);
    console.error(err);
  });
}
