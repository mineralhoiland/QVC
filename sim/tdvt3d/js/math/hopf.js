// Hopf fibration π: S³ → S², stereographic R³, Gauss linking of fibres.
// Linking of two distinct fibres is the geometric origin of the Hopf invariant (Paper C § hopf).

const TWO_PI = Math.PI * 2;

export function hsvToRgb(h, s, v, out, o = 0) {
  const i = Math.floor(h * 6);
  const f = h * 6 - i;
  const p = v * (1 - s);
  const q = v * (1 - f * s);
  const t = v * (1 - (1 - f) * s);
  let r, g, b;
  switch (i % 6) {
    case 0: r = v; g = t; b = p; break;
    case 1: r = q; g = v; b = p; break;
    case 2: r = p; g = v; b = t; break;
    case 3: r = p; g = q; b = v; break;
    case 4: r = t; g = p; b = v; break;
    default: r = v; g = p; b = q;
  }
  out[o] = r; out[o + 1] = g; out[o + 2] = b;
}

export function nColor(nx, ny, nz, out, o = 0) {
  const hue = (Math.atan2(ny, nx) / TWO_PI + 1) % 1;
  const sat = 0.55 + 0.4 * Math.sqrt(nx * nx + ny * ny);
  const val = 0.45 + 0.55 * (0.5 * (nz + 1));
  hsvToRgb(hue, Math.min(1, sat), val, out, o);
}

/** Hopf map π(z1, z2) ∈ S². */
export function hopfMap(w0, w1, w2, w3) {
  const nx = 2 * (w0 * w2 + w1 * w3);
  const ny = 2 * (w1 * w2 - w0 * w3);
  const nz = (w0 * w0 + w1 * w1) - (w2 * w2 + w3 * w3);
  return { nx, ny, nz };
}

/**
 * Point on the Hopf fibre over n(η, ψ).
 * φ1 − φ2 = ψ, (φ1 + φ2)/2 = t  ⇒  fibre coordinate t ∈ [0, 2π).
 */
export function hopfFibrePoint(eta, psi, t) {
  const c = Math.cos(eta), s = Math.sin(eta);
  const ph1 = t + 0.5 * psi;
  const ph2 = t - 0.5 * psi;
  const w0 = c * Math.cos(ph1);
  const w1 = c * Math.sin(ph1);
  const w2 = s * Math.cos(ph2);
  const w3 = s * Math.sin(ph2);
  const n = hopfMap(w0, w1, w2, w3);
  return { w0, w1, w2, w3, ...n };
}

/** Stereographic S³ → R³ from the pole w3 = −1, inverse of the hopfion_texture chart. */
export function stereoFromS3(w0, w1, w2, w3) {
  const d = 1 + w3;
  if (Math.abs(d) < 1e-12) return { x: w0 * 1e6, y: w1 * 1e6, z: w2 * 1e6 };
  return { x: w0 / d, y: w1 / d, z: w2 / d };
}

export function s3EmbedXYZ(w0, w1, w2, _w3) {
  return { x: w0, y: w1, z: w2 };
}

/** Clifford torus |z1| = |z2| = 1/√2, i.e. η = π/4. */
export function isCliffordEta(eta, tol = 1e-6) {
  return Math.abs(eta - Math.PI / 4) < tol;
}

/**
 * Sample nFibres × nSeg closed Hopf fibres.
 * Returns packed Float32Arrays: pos (3), color (3), nBase (3), eta, psi.
 */
export function generateFibres({
  nFibres = 220,
  nSeg = 480,
  stereo = true,
  scale = 40,
  cliffordBias = 0.12,
} = {}) {
  const nPts = nFibres * nSeg;
  const pos = new Float32Array(nPts * 3);
  const col = new Float32Array(nPts * 3);
  const nBase = new Float32Array(nPts * 3);
  const etaA = new Float32Array(nPts);
  const psiA = new Float32Array(nPts);
  const tA = new Float32Array(nPts);
  const rgb = [0, 0, 0];
  const project = stereo ? stereoFromS3 : s3EmbedXYZ;
  let p = 0;
  for (let f = 0; f < nFibres; f++) {
    const u = (f + 0.5) / nFibres;
    // slight oversampling of the Clifford torus (η = π/4)
    let eta = Math.acos(1 - 2 * u) * 0.5; // η ∈ (0, π/2) from z-uniform on S²
    if (cliffordBias > 0 && Math.abs(u - 0.5) < cliffordBias) eta = Math.PI / 4;
    const psi = TWO_PI * ((f * 0.61803398875) % 1);
    for (let s = 0; s < nSeg; s++) {
      const t = TWO_PI * (s / nSeg);
      const h = hopfFibrePoint(eta, psi, t);
      const r = project(h.w0, h.w1, h.w2, h.w3);
      pos[p * 3] = r.x * scale;
      pos[p * 3 + 1] = r.y * scale;
      pos[p * 3 + 2] = r.z * scale;
      nColor(h.nx, h.ny, h.nz, rgb);
      col[p * 3] = rgb[0]; col[p * 3 + 1] = rgb[1]; col[p * 3 + 2] = rgb[2];
      nBase[p * 3] = h.nx; nBase[p * 3 + 1] = h.ny; nBase[p * 3 + 2] = h.nz;
      etaA[p] = eta; psiA[p] = psi; tA[p] = t;
      p++;
    }
  }
  return { pos, col, nBase, eta: etaA, psi: psiA, t0: tA, nPts, nFibres, nSeg, scale };
}

/** Packed polyline of one fibre (closed). */
export function fibrePolyline(eta, psi, nSeg, stereo, scale) {
  const pts = new Float32Array(nSeg * 3);
  const project = stereo ? stereoFromS3 : s3EmbedXYZ;
  for (let s = 0; s < nSeg; s++) {
    const h = hopfFibrePoint(eta, psi, TWO_PI * s / nSeg);
    const r = project(h.w0, h.w1, h.w2, h.w3);
    pts[s * 3] = r.x * scale;
    pts[s * 3 + 1] = r.y * scale;
    pts[s * 3 + 2] = r.z * scale;
  }
  return pts;
}

/**
 * Gauss double integral for two closed polylines (N×3 packed Float32/Float64).
 * Distinct Hopf fibres → ±1.
 */
export function gaussLinking(c1, c2, n1, n2) {
  if (n1 < 3 || n2 < 3) return NaN;
  let acc = 0;
  for (let i = 0; i < n1; i++) {
    const i1 = ((i + 1) % n1) * 3, i0 = i * 3;
    const ax = c1[i0], ay = c1[i0 + 1], az = c1[i0 + 2];
    const dx1 = c1[i1] - ax, dy1 = c1[i1 + 1] - ay, dz1 = c1[i1 + 2] - az;
    const mx1 = ax + 0.5 * dx1, my1 = ay + 0.5 * dy1, mz1 = az + 0.5 * dz1;
    for (let j = 0; j < n2; j++) {
      const j1 = ((j + 1) % n2) * 3, j0 = j * 3;
      const bx = c2[j0], by = c2[j0 + 1], bz = c2[j0 + 2];
      const dx2 = c2[j1] - bx, dy2 = c2[j1 + 1] - by, dz2 = c2[j1 + 2] - bz;
      const mx2 = bx + 0.5 * dx2, my2 = by + 0.5 * dy2, mz2 = bz + 0.5 * dz2;
      const rx = mx1 - mx2, ry = my1 - my2, rz = mz1 - mz2;
      const r2 = rx * rx + ry * ry + rz * rz;
      const rn = Math.pow(r2, 1.5);
      if (rn < 1e-18) continue;
      const cx = dy1 * dz2 - dz1 * dy2;
      const cy = dz1 * dx2 - dx1 * dz2;
      const cz = dx1 * dy2 - dy1 * dx2;
      acc += (cx * rx + cy * ry + cz * rz) / rn;
    }
  }
  return acc / (4 * Math.PI);
}

export function twoFibreLinking({ nSeg = 256, stereo = true, scale = 40, eta1 = 0.4, eta2 = 1.1, psi1 = 0.2, psi2 = 2.4 } = {}) {
  const c1 = fibrePolyline(eta1, psi1, nSeg, stereo, scale);
  const c2 = fibrePolyline(eta2, psi2, nSeg, stereo, scale);
  return { Lk: gaussLinking(c1, c2, nSeg, nSeg), c1, c2, nSeg };
}
