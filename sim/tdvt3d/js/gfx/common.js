// Shared Three.js helpers: glow points, axes in nm, disposal, camera presets.

export const POINT_VS = /* glsl */ `
uniform float uSize;
attribute float aSize;
varying vec3 vColor;
void main() {
  vColor = color;
  vec4 mv = modelViewMatrix * vec4(position, 1.0);
  float ps = aSize * uSize * (220.0 / max(1.0, -mv.z));
  gl_PointSize = clamp(ps, 1.0, 18.0);
  gl_Position = projectionMatrix * mv;
}
`;

export const POINT_FS = /* glsl */ `
uniform float uOpacity;
varying vec3 vColor;
void main() {
  vec2 p = gl_PointCoord * 2.0 - 1.0;
  float r2 = dot(p, p);
  if (r2 > 1.0) discard;
  float a = exp(-3.2 * r2) * uOpacity;
  gl_FragColor = vec4(vColor * a, a);
}
`;

export const HOPF_FLOW_VS = /* glsl */ `
uniform float uTime;
uniform float uOmega;
uniform float uStereo;
uniform float uScale;
uniform float uSize;
attribute float aEta;
attribute float aPsi;
attribute float aT0;
varying vec3 vColor;

vec3 hsv(float h, float s, float v) {
  vec3 c = clamp(abs(mod(h * 6.0 + vec3(0.0, 4.0, 2.0), 6.0) - 3.0) - 1.0, 0.0, 1.0);
  return v * mix(vec3(1.0), c, s);
}

void main() {
  float t = aT0 + uTime * uOmega;
  float c = cos(aEta), s = sin(aEta);
  float ph1 = t + aPsi * 0.5;
  float ph2 = t - aPsi * 0.5;
  vec4 w = vec4(c * cos(ph1), c * sin(ph1), s * cos(ph2), s * sin(ph2));
  vec3 p;
  if (uStereo > 0.5) {
    float d = max(1e-4, 1.0 + w.w);
    p = w.xyz / d * uScale;
  } else {
    p = w.xyz * uScale;
  }
  float nx = sin(2.0 * aEta) * cos(aPsi);
  float ny = sin(2.0 * aEta) * sin(aPsi);
  float nz = cos(2.0 * aEta);
  float hue = atan(ny, nx) / 6.28318530718 + 0.5;
  vColor = hsv(hue, 0.75, 0.45 + 0.55 * (0.5 * (nz + 1.0)));
  vec4 mv = modelViewMatrix * vec4(p, 1.0);
  gl_PointSize = clamp(uSize * (200.0 / max(1.0, -mv.z)), 1.0, 16.0);
  gl_Position = projectionMatrix * mv;
}
`;

export function makeGlowMaterial(THREE, { size = 2.4, opacity = 1 } = {}) {
  return new THREE.ShaderMaterial({
    uniforms: {
      uSize: { value: size },
      uOpacity: { value: opacity },
    },
    vertexShader: POINT_VS,
    fragmentShader: POINT_FS,
    vertexColors: true,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
  });
}

export function makePointCloud(THREE, positions, colors, { size = 2.4, sizes = null, opacity = 1 } = {}) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
  const n = positions.length / 3;
  const sz = new Float32Array(n);
  if (sizes) {
    const den = size || 1;
    for (let i = 0; i < n; i++) sz[i] = sizes[i] / den;
  } else {
    sz.fill(1);
  }
  geo.setAttribute('aSize', new THREE.BufferAttribute(sz, 1));
  const mat = makeGlowMaterial(THREE, { size, opacity });
  const pts = new THREE.Points(geo, mat);
  pts.frustumCulled = false;
  pts.userData.baseSize = size;
  return pts;
}

export function setPointStyle(pts, { size, opacity } = {}) {
  const u = pts?.material?.uniforms;
  if (!u) return;
  if (size != null && u.uSize) u.uSize.value = size;
  if (opacity != null && u.uOpacity) u.uOpacity.value = opacity;
}

export function makeLineLoop(THREE, packed, color, { linewidth = 1 } = {}) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(packed, 3));
  const mat = new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.95 });
  const line = new THREE.LineLoop(geo, mat);
  line.userData.linewidth = linewidth;
  return line;
}

export function makeLine(THREE, packed, color) {
  const geo = new THREE.BufferGeometry();
  geo.setAttribute('position', new THREE.BufferAttribute(packed, 3));
  const mat = new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.9 });
  return new THREE.Line(geo, mat);
}

export function makeNmAxes(THREE, size = 80) {
  const g = new THREE.Group();
  g.add(new THREE.AxesHelper(size));
  const mk = (text, color, pos) => {
    const c = document.createElement('canvas');
    c.width = 256; c.height = 64;
    const ctx = c.getContext('2d');
    ctx.clearRect(0, 0, 256, 64);
    ctx.fillStyle = color;
    ctx.font = '28px ui-monospace, SFMono-Regular, Menlo, monospace';
    ctx.fillText(text, 8, 40);
    const tex = new THREE.CanvasTexture(c);
    const mat = new THREE.SpriteMaterial({ map: tex, transparent: true, depthWrite: false });
    const spr = new THREE.Sprite(mat);
    spr.position.copy(pos);
    spr.scale.set(size * 0.45, size * 0.12, 1);
    g.add(spr);
  };
  mk('x (nm)', '#ff6b6b', new THREE.Vector3(size * 1.05, 0, 0));
  mk('y (nm)', '#7dffb3', new THREE.Vector3(0, size * 1.05, 0));
  mk('z (nm)', '#6bb8ff', new THREE.Vector3(0, 0, size * 1.05));
  return g;
}

export function makeTorusWire(THREE, R, r, color = 0x6688aa) {
  const geo = new THREE.TorusGeometry(R, r, 12, 64);
  const mat = new THREE.MeshBasicMaterial({
    color, wireframe: true, transparent: true, opacity: 0.28,
  });
  const mesh = new THREE.Mesh(geo, mat);
  mesh.rotation.x = Math.PI / 2;
  return mesh;
}

export function disposeObject(obj) {
  obj.traverse((o) => {
    if (o.geometry) o.geometry.dispose();
    const mats = o.material ? (Array.isArray(o.material) ? o.material : [o.material]) : [];
    for (const m of mats) {
      if (!m) continue;
      for (const k of Object.keys(m)) {
        const v = m[k];
        if (v && v.isTexture) v.dispose();
      }
      m.dispose();
    }
  });
}

export const CAMERA_PRESETS = {
  iso: (scale) => ({ pos: [1.15 * scale, 0.85 * scale, 1.15 * scale], target: [0, 0, 0] }),
  top: (scale) => ({ pos: [0, 1.7 * scale, 0.02 * scale], target: [0, 0, 0] }),
  core: (scale) => ({ pos: [0.35 * scale, 0.22 * scale, 0.35 * scale], target: [0, 0, 0] }),
};

export function applyPreset(camera, controls, name, scale) {
  const p = CAMERA_PRESETS[name](scale);
  camera.position.set(...p.pos);
  controls.target.set(...p.target);
  controls.update();
}
