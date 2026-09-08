import * as THREE from "three";
import { particleFrag, particleVert } from "../shaders/particleShaders";

export const DEFAULT_COUNT = 524_288; // 512K
export const MAX_COUNT = 1_048_576; // 1M

export type ParticleUniforms = {
  mode: number;
  speed: number;
  Q: number;
  scale: number;
  twist: number;
  lattice: number;
  link: number;
  flow: number;
  alpha: number;
  N: number;
  g0: number;
  pointSize: number;
  spread: number;
  glow: number;
  opacity: number;
};

const defaultUniforms: ParticleUniforms = {
  mode: 0,
  speed: 1,
  Q: 1,
  scale: 1.4,
  twist: 1,
  lattice: 1.6,
  link: 0.85,
  flow: 1.25,
  alpha: 0.45,
  N: 5,
  g0: 0.21,
  pointSize: 0.85,
  spread: 1,
  glow: 0.4,
  opacity: 0.55,
};

export class ParticleSystem {
  readonly points: THREE.Points;
  readonly material: THREE.ShaderMaterial;
  private geometry: THREE.BufferGeometry;
  private capacity: number;
  private count: number;
  private uniforms: ParticleUniforms;

  constructor(initialCount = DEFAULT_COUNT) {
    this.capacity = MAX_COUNT;
    this.count = Math.min(MAX_COUNT, Math.max(16_384, initialCount));
    this.uniforms = { ...defaultUniforms };
    this.geometry = new THREE.BufferGeometry();
    this.material = new THREE.ShaderMaterial({
      vertexShader: particleVert,
      fragmentShader: particleFrag,
      transparent: true,
      depthWrite: false,
      blending: THREE.AdditiveBlending,
      uniforms: {
        uTime: { value: 0 },
        uSpeed: { value: this.uniforms.speed },
        uMode: { value: this.uniforms.mode },
        uQ: { value: this.uniforms.Q },
        uScale: { value: this.uniforms.scale },
        uTwist: { value: this.uniforms.twist },
        uLattice: { value: this.uniforms.lattice },
        uLink: { value: this.uniforms.link },
        uFlow: { value: this.uniforms.flow },
        uAlpha: { value: this.uniforms.alpha },
        uN: { value: this.uniforms.N },
        uG0: { value: this.uniforms.g0 },
        uPointSize: { value: this.uniforms.pointSize },
        uSpread: { value: this.uniforms.spread },
        uGlow: { value: this.uniforms.glow },
        uOpacity: { value: this.uniforms.opacity },
      },
    });
    // Pre-allocate 1M seeds so the count slider can go to full capacity without realloc.
    this.allocate(MAX_COUNT);
    this.geometry.setDrawRange(0, this.count);
    this.points = new THREE.Points(this.geometry, this.material);
    this.points.frustumCulled = false;
  }

  private allocate(n: number) {
    const seed = new Float32Array(n * 4);
    const rand = new Float32Array(n);
    // Stratified-ish filling for denser fiber coverage
    const golden = 0.61803398875;
    for (let i = 0; i < n; i++) {
      const f = (i + 0.5) / n;
      const eta = Math.acos(1 - 2 * ((i * golden) % 1));
      const phi = TAU * ((i * golden * 1.324) % 1);
      const psi0 = TAU * ((i * golden * 2.618) % 1);
      const sid = (i % 54) + f;
      seed[i * 4] = eta;
      seed[i * 4 + 1] = phi;
      seed[i * 4 + 2] = psi0;
      seed[i * 4 + 3] = sid;
      rand[i] = (Math.sin(i * 12.9898) * 43758.5453) % 1;
      if (rand[i] < 0) rand[i] += 1;
    }
    this.geometry.setAttribute("aSeed", new THREE.BufferAttribute(seed, 4));
    this.geometry.setAttribute("aRand", new THREE.BufferAttribute(rand, 1));
    // Dummy position attribute required by Three.js
    this.geometry.setAttribute("position", new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    this.geometry.setDrawRange(0, n);
    this.capacity = n;
    this.count = n;
  }

  setCount(n: number) {
    const next = Math.min(this.capacity, Math.max(16_384, Math.floor(n)));
    this.count = next;
    this.geometry.setDrawRange(0, next);
  }

  getCount() {
    return this.count;
  }

  getCapacity() {
    return this.capacity;
  }

  setUniforms(partial: Partial<ParticleUniforms>) {
    Object.assign(this.uniforms, partial);
    const u = this.material.uniforms;
    u.uSpeed.value = this.uniforms.speed;
    u.uMode.value = this.uniforms.mode;
    u.uQ.value = this.uniforms.Q;
    u.uScale.value = this.uniforms.scale;
    u.uTwist.value = this.uniforms.twist;
    u.uLattice.value = this.uniforms.lattice;
    u.uLink.value = this.uniforms.link;
    u.uFlow.value = this.uniforms.flow;
    u.uAlpha.value = this.uniforms.alpha;
    u.uN.value = this.uniforms.N;
    u.uG0.value = this.uniforms.g0;
    u.uPointSize.value = this.uniforms.pointSize;
    u.uSpread.value = this.uniforms.spread;
    u.uGlow.value = this.uniforms.glow;
    u.uOpacity.value = this.uniforms.opacity;
  }

  getUniforms() {
    return { ...this.uniforms };
  }

  update(time: number, running: boolean) {
    if (running) this.material.uniforms.uTime.value = time;
  }

  dispose() {
    this.geometry.dispose();
    this.material.dispose();
  }
}

const TAU = Math.PI * 2;
