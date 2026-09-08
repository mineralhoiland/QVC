/** GPU particle shaders — positions/colors computed entirely on GPU. */

export const particleVert = /* glsl */ `
precision highp float;

attribute vec4 aSeed; // eta, phi, psi0, site/fiber id
attribute float aRand;

uniform float uTime;
uniform float uSpeed;
uniform int uMode;       // 0 hopf fiber, 1 hopfion dyn, 2 lattice, 3 vortex, 4 catalyzon, 5 horizon
uniform float uQ;        // Hopf charge / winding
uniform float uScale;
uniform float uTwist;
uniform float uLattice;  // lattice spacing
uniform float uLink;     // interlink strength
uniform float uFlow;     // v_s / c_s
uniform float uAlpha;    // fold / ZPF proxy
uniform float uN;        // mode count (catalyzon)
uniform float uG0;
uniform float uPointSize;
uniform float uSpread;

varying vec3 vColor;
varying float vAlpha;

#define PI 3.141592653589793
#define TAU 6.283185307179586

vec3 hsv2rgb(vec3 c) {
  vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
  vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
  return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

// Stereographic Hopf map S^3 → R^3 for fiber coords (η, φ, ψ)
vec3 hopfStereo(float eta, float phi, float psi) {
  float ce = cos(eta * 0.5);
  float se = sin(eta * 0.5);
  float x = ce * cos(psi);
  float y = ce * sin(psi);
  float z = se * cos(phi);
  float w = se * sin(phi);
  // Project from near (0,0,0,1); soft-clip pole so density stays finite
  float d = max(1.05 - w, 0.12);
  vec3 p = vec3(x, y, z) / d;
  float L = length(p);
  // Soft radial compress keeps structure readable at 500K+
  if (L > 3.5) p *= 3.5 * (1.0 + 0.35 * log(L / 3.5)) / L;
  return p;
}

// Toroidal hopfion density shell (particle lives near |m|=1 core tube)
vec3 hopfionTube(float u, float v, float R, float r, float Q, float t) {
  float ang = u * TAU * Q + t;
  float tube = v * TAU + t * 0.37;
  float rr = r * (0.85 + 0.15 * sin(3.0 * tube + ang));
  float cx = (R + rr * cos(tube)) * cos(ang);
  float cy = rr * sin(tube);
  float cz = (R + rr * cos(tube)) * sin(ang);
  return vec3(cx, cy, cz);
}

vec3 bccSite(float id, float a) {
  float n = floor(id + 0.5);
  float ix = mod(n, 3.0) - 1.0;
  float iy = mod(floor(n / 3.0), 3.0) - 1.0;
  float iz = mod(floor(n / 9.0), 3.0) - 1.0;
  vec3 c = vec3(ix, iy, iz) * a;
  // BCC body-center offset for odd parity sites
  if (mod(n, 2.0) > 0.5) c += vec3(0.5 * a);
  return c;
}

void main() {
  float eta = aSeed.x;
  float phi = aSeed.y;
  float psi0 = aSeed.z;
  float sid = aSeed.w;
  float t = uTime * uSpeed;

  vec3 pos = vec3(0.0);
  vec3 col = vec3(1.0);
  float alpha = 0.85;

  if (uMode == 0) {
    // —— Hopf fibration time evolution ——
    float psi = psi0 + t * (1.0 + 0.35 * uQ);
    // mild base-sphere drift
    float eta2 = clamp(eta + 0.04 * sin(t * 0.2 + phi), 0.02, PI - 0.02);
    float phi2 = phi + t * 0.07 * uTwist;
    pos = hopfStereo(eta2, phi2, psi) * uScale * 0.55;
    // color by base S² (η,φ) — fiber identity
    float hue = fract(phi2 / TAU + 0.08 * eta2 / PI);
    float sat = 0.85 + 0.15 * sin(psi);
    float val = 0.65 + 0.35 * (1.0 - eta2 / PI);
    col = hsv2rgb(vec3(hue, sat, val));
    alpha = 0.22 + 0.18 * (0.5 + 0.5 * cos(psi));
  }
  else if (uMode == 1) {
    // —— Dynamic linked hopfions ——
    float family = step(0.5, fract(sid * 0.173));
    float Q = max(1.0, uQ);
    float R = 1.15;
    float r = 0.38;
    float u = eta / PI;
    float v = phi / TAU;
    float phase = t * (0.9 + 0.4 * family) + psi0;
    // relative motion / linking
    float linkOff = uLink * 0.55 * family;
    vec3 p1 = hopfionTube(u + linkOff, v, R, r, Q, phase);
    // second hopfion rotated / phase-shifted for linking
    vec3 p2 = hopfionTube(v, u, R * 0.92, r * 1.05, Q, phase * 1.13 + PI * 0.5);
    p2 = vec3(p2.x, p2.z, -p2.y); // reorient for Hopf link
    pos = mix(p1, p2, family) * uScale * 0.7;
    // slight orbital dance of centers
    pos += vec3(sin(t * 0.3 + family * PI), 0.15 * cos(t * 0.5), cos(t * 0.3)) * 0.2 * uLink;
    float hue = family > 0.5 ? 0.92 : 0.52; // magenta / cyan families
    hue += 0.08 * sin(phase * Q);
    col = hsv2rgb(vec3(fract(hue), 0.92, 0.95));
    alpha = 0.28;
  }
  else if (uMode == 2) {
    // —— Lattice spacetime: interlinked hopfions on BCC ——
    float a = uLattice;
    vec3 center = bccSite(mod(sid, 27.0), a);
    float Q = max(1.0, uQ);
    float u = fract(eta / PI + t * 0.02);
    float v = fract(phi / TAU + psi0);
    float localT = t + sid * 0.11;
    vec3 local = hopfionTube(u, v, 0.42, 0.16, Q, localT);
    // neighbor linking filaments
    vec3 nOff = normalize(center + vec3(0.001)) * (0.15 * uLink);
    float fiberMix = step(0.72, aRand);
    pos = (center + local + nOff * fiberMix * sin(localT + v * TAU)) * uScale * 0.55;
    float hue = fract(length(center) * 0.15 + sid * 0.037 + 0.55);
    col = hsv2rgb(vec3(hue, 0.9, 0.95));
    alpha = 0.18 + 0.16 * fiberMix;
  }
  else if (uMode == 3) {
    // —— Vortex mechanics on particle lattice ——
    float nx = floor(eta / PI * 48.0);
    float ny = floor(phi / TAU * 48.0);
    float nz = floor(fract(sid) * 16.0);
    vec3 grid = (vec3(nx, ny, nz) - vec3(24.0, 24.0, 8.0)) * (0.085 * uLattice);
    // vortex cores along z
    float nV = max(1.0, floor(uQ));
    vec2 core = vec2(0.0);
    float minD = 100.0;
    for (int k = 0; k < 6; k++) {
      if (float(k) >= nV) break;
      float ang = float(k) * TAU / nV + t * 0.15;
      vec2 c = 0.9 * vec2(cos(ang), sin(ang));
      float d = length(grid.xy - c);
      if (d < minD) { minD = d; core = c; }
    }
    float circ = atan(grid.y - core.y, grid.x - core.x);
    float rad = max(minD, 0.02);
    // azimuthal velocity + Kelvin wave ripple
    float vTheta = (1.0 / rad) * uTwist;
    grid.xy += vec2(-sin(circ), cos(circ)) * vTheta * 0.02 * sin(t * 2.0 + grid.z * 4.0);
    grid.z += 0.08 * sin(circ * uQ + t * 3.0 - rad * 5.0);
    pos = grid * uScale;
    float hue = fract(0.08 + circ / TAU + 0.12 * (1.0 - exp(-rad)));
    col = hsv2rgb(vec3(hue, 0.95, mix(1.0, 0.55, exp(-rad * 2.0))));
    alpha = mix(0.45, 0.12, smoothstep(0.0, 2.5, rad));
  }
  else if (uMode == 4) {
    // —— Catalyzon Spec(G): democratic λ_min=-g0, λ_max=(N-1)g0 ——
    float N = max(2.0, uN);
    float mode = mod(floor(sid), N);
    float isBright = step(N - 1.5, mode); // collective
    float lam = mix(-uG0, (N - 1.0) * uG0, isBright);
    float ang = psi0 * TAU + t * (0.4 + 0.25 * abs(lam)) + mode;
    float rad = 0.6 + 0.35 * mode / N + 0.15 * sin(t + eta * 3.0);
    // fold melting: as α rises, dark modes scatter
    float melt = smoothstep(0.3, 0.95, uAlpha);
    rad += melt * aRand * 1.2 * (1.0 - isBright);
    pos = vec3(cos(ang) * rad, lam * 0.55 + 0.2 * sin(ang * 2.0 + t), sin(ang) * rad);
    pos *= uScale * 0.85;
    float hue = mix(0.58, 0.08, isBright); // teal dark / amber bright
    hue += 0.05 * mode / N;
    col = hsv2rgb(vec3(fract(hue), 0.9, 0.9));
    alpha = mix(0.16, 0.4, isBright);
  }
  else {
    // —— Acoustic horizon flow past hopfion core ——
    float R = 1.0;
    float r = 0.32;
    float u = eta / PI;
    float v = phi / TAU;
    vec3 tube = hopfionTube(u, v, R, r, max(1.0, uQ), t * 0.2);
    // radial inflow / wind
    float rad = length(tube.xz);
    float rH = 1.0 / max(uFlow, 0.15); // horizon proxy ~ cs/vs scaled
    vec3 flow = -normalize(tube + vec3(0.001)) * (0.35 * uFlow);
    float nearH = exp(-abs(rad - rH) * 6.0);
    pos = (tube + flow * (t * 0.15 + aRand) + vec3(0.0, nearH * 0.25, 0.0)) * uScale * 0.7;
    // wrap particles that fell in
    if (length(pos) < 0.15) pos *= 8.0;
    float hue = mix(0.5, 0.02, nearH);
    col = hsv2rgb(vec3(hue, 0.92, 0.7 + 0.3 * nearH));
    alpha = 0.14 + 0.28 * nearH;
  }

  // soft depth cue
  float dist = length(pos);
  pos += normalize(pos + vec3(0.001)) * uSpread * (aRand - 0.5) * 0.05;

  vColor = col;
  vAlpha = alpha * smoothstep(10.0, 0.4, dist);

  vec4 mv = modelViewMatrix * vec4(pos, 1.0);
  gl_Position = projectionMatrix * mv;
  // Smaller sprites at 500K+ so fibers remain resolved
  float psz = uPointSize * (180.0 / max(2.0, -mv.z));
  gl_PointSize = clamp(psz, 0.4, 4.5);
}
`;

export const particleFrag = /* glsl */ `
precision highp float;

varying vec3 vColor;
varying float vAlpha;

uniform float uGlow;
uniform float uOpacity;

void main() {
  vec2 uv = gl_PointCoord * 2.0 - 1.0;
  float d = dot(uv, uv);
  if (d > 1.0) discard;
  float core = exp(-d * 4.5);
  float halo = exp(-d * 1.6) * uGlow * 0.35;
  float a = (core * 0.85 + halo) * vAlpha * uOpacity;
  if (a < 0.02) discard;
  vec3 col = vColor * (0.55 + 0.7 * core);
  gl_FragColor = vec4(col, a);
}
`;
