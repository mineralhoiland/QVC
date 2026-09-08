import * as THREE from "three";
import { OrbitControls } from "three/examples/jsm/controls/OrbitControls.js";
import { DEFAULT_COUNT, ParticleSystem, type ParticleUniforms } from "./ParticleSystem";

export class SceneHost {
  readonly renderer: THREE.WebGLRenderer;
  readonly scene: THREE.Scene;
  readonly camera: THREE.PerspectiveCamera;
  readonly controls: OrbitControls;
  readonly particles: ParticleSystem;
  private raf = 0;
  private running = true;
  private clock = new THREE.Clock();
  private mounted = false;
  private onFps?: (fps: number, count: number) => void;
  private frames = 0;
  private fpsT = 0;

  constructor(canvas: HTMLCanvasElement, count = DEFAULT_COUNT) {
    this.renderer = new THREE.WebGLRenderer({
      canvas,
      antialias: false,
      powerPreference: "high-performance",
      alpha: false,
    });
    this.renderer.setPixelRatio(Math.min(window.devicePixelRatio, 2));
    this.renderer.setClearColor(0x03060c, 1);
    this.renderer.sortObjects = false;

    this.scene = new THREE.Scene();
    // subtle vignette-like fog
    this.scene.fog = new THREE.FogExp2(0x03060c, 0.028);

    this.camera = new THREE.PerspectiveCamera(55, 1, 0.05, 200);
    this.camera.position.set(3.8, 2.2, 4.6);

    this.controls = new OrbitControls(this.camera, canvas);
    this.controls.enableDamping = true;
    this.controls.dampingFactor = 0.06;
    this.controls.minDistance = 0.8;
    this.controls.maxDistance = 40;

    this.particles = new ParticleSystem(count);
    this.scene.add(this.particles.points);

    // faint reference sphere (topology ambient)
    const shell = new THREE.Mesh(
      new THREE.SphereGeometry(3.2, 32, 24),
      new THREE.MeshBasicMaterial({
        color: 0x114455,
        wireframe: true,
        transparent: true,
        opacity: 0.06,
      }),
    );
    this.scene.add(shell);
  }

  mount(container: HTMLElement) {
    const resize = () => {
      const w = container.clientWidth;
      const h = container.clientHeight;
      this.camera.aspect = w / Math.max(h, 1);
      this.camera.updateProjectionMatrix();
      this.renderer.setSize(w, h, false);
    };
    resize();
    const ro = new ResizeObserver(resize);
    ro.observe(container);
    this.mounted = true;
    const loop = () => {
      this.raf = requestAnimationFrame(loop);
      const dt = this.clock.getDelta();
      const t = this.clock.elapsedTime;
      this.controls.update();
      this.particles.update(t, this.running);
      this.renderer.render(this.scene, this.camera);
      this.frames++;
      this.fpsT += dt;
      if (this.fpsT >= 0.5) {
        const fps = this.frames / this.fpsT;
        this.onFps?.(fps, this.particles.getCount());
        this.frames = 0;
        this.fpsT = 0;
      }
    };
    loop();
    return () => {
      cancelAnimationFrame(this.raf);
      ro.disconnect();
      this.dispose();
    };
  }

  setRunning(v: boolean) {
    this.running = v;
  }

  setUniforms(u: Partial<ParticleUniforms>) {
    this.particles.setUniforms(u);
  }

  setCount(n: number) {
    this.particles.setCount(n);
  }

  setFpsCallback(cb: (fps: number, count: number) => void) {
    this.onFps = cb;
  }

  dispose() {
    if (!this.mounted) return;
    this.mounted = false;
    cancelAnimationFrame(this.raf);
    this.particles.dispose();
    this.controls.dispose();
    this.renderer.dispose();
  }
}
