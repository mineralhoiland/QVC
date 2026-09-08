// Stereographic Hopf hopfion, Whitehead Q = (1/16π²)∫A·B, Ward FN energy, CS flux, preimages.
// Port of qvc/viz/fields3d.py (windowed texture, Coulomb-gauge FFT inversion).

import { fft3d, angularFreqs } from './fft.js';
import { gaussLinking, nColor } from './hopf.js';
import { LOCKS } from './locks.js';

export const WARD_VK = (3 / 16) ** (3 / 8); // ≈ 0.5338
export const WARD_Q1_MIN = 1.22;
const PI2 = Math.PI * Math.PI;

export function makeGrid(n, Lnm) {
  const x = new Float64Array(n);
  const dx = (2 * Lnm) / n; // linspace(-L, L, n, endpoint=False)
  for (let i = 0; i < n; i++) x[i] = -Lnm + i * dx;
  return { n, Lnm, dx, x };
}

function wrap(i, n) { return (i + n) % n; }
export function idx3(ix, iy, iz, n) { return (ix * n + iy) * n + iz; }

function pgrad(f, dx, n, axis, out) {
  const n2 = n * n;
  for (let ix = 0; ix < n; ix++) {
    for (let iy = 0; iy < n; iy++) {
      for (let iz = 0; iz < n; iz++) {
        let a, b;
        if (axis === 0) {
          a = f[idx3(wrap(ix + 1, n), iy, iz, n)];
          b = f[idx3(wrap(ix - 1, n), iy, iz, n)];
        } else if (axis === 1) {
          a = f[idx3(ix, wrap(iy + 1, n), iz, n)];
          b = f[idx3(ix, wrap(iy - 1, n), iz, n)];
        } else {
          a = f[idx3(ix, iy, wrap(iz + 1, n), n)];
          b = f[idx3(ix, iy, wrap(iz - 1, n), n)];
        }
        out[idx3(ix, iy, iz, n)] = (a - b) / (2 * dx);
      }
    }
  }
}

/** Windowed stereographic Hopf map. Returns {nx,ny,nz, n, dx, x, Lnm, coreNm}. */
export function hopfionTexture(nGrid, {
  coreNm = 50,
  Lnm = null,
  RcutOverCore = 2.5,
  windowWidthOverCore = 0.4,
  winding = 1,
} = {}) {
  const L = Lnm == null ? 4 * coreNm : Lnm;
  const g = makeGrid(nGrid, L);
  const N = nGrid, N3 = N * N * N;
  const nx = new Float64Array(N3), ny = new Float64Array(N3), nz = new Float64Array(N3);
  let maxDev = 0;
  for (let ix = 0; ix < N; ix++) {
    const X = g.x[ix] / coreNm;
    for (let iy = 0; iy < N; iy++) {
      const Y = g.x[iy] / coreNm;
      for (let iz = 0; iz < N; iz++) {
        const Z = g.x[iz] / coreNm;
        const r2 = X * X + Y * Y + Z * Z;
        const r = Math.sqrt(r2);
        const den = 1 + r2;
        let w0 = 2 * X / den, w1 = 2 * Y / den, w2 = 2 * Z / den, w3 = (1 - r2) / den;
        if (winding !== 1) {
          const a1 = Math.hypot(w0, w1);
          const ph = Math.atan2(w1, w0) * winding;
          w0 = a1 * Math.cos(ph); w1 = a1 * Math.sin(ph);
        }
        const z1r = w0, z1i = w1, z2r = w2, z2i = w3;
        // z1 conj(z2)
        const zr = z1r * z2r + z1i * z2i;
        const zi = z1i * z2r - z1r * z2i;
        let nxx = 2 * zr, nyy = 2 * zi;
        const a1s = z1r * z1r + z1i * z1i, a2s = z2r * z2r + z2i * z2i;
        let nzz = a1s - a2s;
        const t = 0.5 * (1 + Math.tanh((r - RcutOverCore) / windowWidthOverCore));
        nxx *= (1 - t); nyy *= (1 - t); nzz = (1 - t) * nzz - t;
        const nn = Math.hypot(nxx, nyy, nzz) || 1e-15;
        const i = idx3(ix, iy, iz, N);
        nx[i] = nxx / nn; ny[i] = nyy / nn; nz[i] = nzz / nn;
        const dev = Math.abs(nn / nn - 1); // after normalize, check unit
        const nd = Math.abs(Math.hypot(nx[i], ny[i], nz[i]) - 1);
        if (nd > maxDev) maxDev = nd;
      }
    }
  }
  return { nx, ny, nz, n: N, dx: g.dx, x: g.x, Lnm: L, coreNm, normDev: maxDev };
}

function berryB(tex) {
  const { nx, ny, nz, n: N, dx } = tex;
  const N3 = N * N * N;
  const dnx = [new Float64Array(N3), new Float64Array(N3), new Float64Array(N3)];
  const dny = [new Float64Array(N3), new Float64Array(N3), new Float64Array(N3)];
  const dnz = [new Float64Array(N3), new Float64Array(N3), new Float64Array(N3)];
  pgrad(nx, dx, N, 0, dnx[0]); pgrad(nx, dx, N, 1, dnx[1]); pgrad(nx, dx, N, 2, dnx[2]);
  pgrad(ny, dx, N, 0, dny[0]); pgrad(ny, dx, N, 1, dny[1]); pgrad(ny, dx, N, 2, dny[2]);
  pgrad(nz, dx, N, 0, dnz[0]); pgrad(nz, dx, N, 1, dnz[1]); pgrad(nz, dx, N, 2, dnz[2]);
  const Bx = new Float64Array(N3), By = new Float64Array(N3), Bz = new Float64Array(N3);
  const F = (i, j, k) => {
    const cx = dny[i][k] * dnz[j][k] - dnz[i][k] * dny[j][k];
    const cy = dnz[i][k] * dnx[j][k] - dnx[i][k] * dnz[j][k];
    const cz = dnx[i][k] * dny[j][k] - dny[i][k] * dnx[j][k];
    return nx[k] * cx + ny[k] * cy + nz[k] * cz;
  };
  for (let k = 0; k < N3; k++) {
    Bx[k] = F(1, 2, k);
    By[k] = F(2, 0, k);
    Bz[k] = F(0, 1, k);
  }
  return { Bx, By, Bz, dnx, dny, dnz };
}

function fftField(src, N, inverse) {
  const re = new Float64Array(src);
  const im = new Float64Array(src.length);
  fft3d(re, im, N, inverse);
  return { re, im };
}

function coulombA(Bx, By, Bz, dx, N) {
  const kx = angularFreqs(N, dx);
  const Bhx = fftField(Bx, N, false);
  const Bhy = fftField(By, N, false);
  const Bhz = fftField(Bz, N, false);
  const N3 = N * N * N;
  const Ahx = { re: new Float64Array(N3), im: new Float64Array(N3) };
  const Ahy = { re: new Float64Array(N3), im: new Float64Array(N3) };
  const Ahz = { re: new Float64Array(N3), im: new Float64Array(N3) };
  const Bcx = { re: new Float64Array(N3), im: new Float64Array(N3) };
  const Bcy = { re: new Float64Array(N3), im: new Float64Array(N3) };
  const Bcz = { re: new Float64Array(N3), im: new Float64Array(N3) };
  for (let ix = 0; ix < N; ix++) {
    const KX = kx[ix];
    for (let iy = 0; iy < N; iy++) {
      const KY = kx[iy];
      for (let iz = 0; iz < N; iz++) {
        const KZ = kx[iz];
        const i = idx3(ix, iy, iz, N);
        const K2 = KX * KX + KY * KY + KZ * KZ;
        const K2s = K2 > 0 ? K2 : 1;
        const kdotBr = KX * Bhx.re[i] + KY * Bhy.re[i] + KZ * Bhz.re[i];
        const kdotBi = KX * Bhx.im[i] + KY * Bhy.im[i] + KZ * Bhz.im[i];
        let bcxr = Bhx.re[i] - KX * kdotBr / K2s;
        let bcxi = Bhx.im[i] - KX * kdotBi / K2s;
        let bcyr = Bhy.re[i] - KY * kdotBr / K2s;
        let bcyi = Bhy.im[i] - KY * kdotBi / K2s;
        let bczr = Bhz.re[i] - KZ * kdotBr / K2s;
        let bczi = Bhz.im[i] - KZ * kdotBi / K2s;
        if (ix === 0 && iy === 0 && iz === 0) {
          bcxr = bcxi = bcyr = bcyi = bczr = bczi = 0;
        }
        Bcx.re[i] = bcxr; Bcx.im[i] = bcxi;
        Bcy.re[i] = bcyr; Bcy.im[i] = bcyi;
        Bcz.re[i] = bczr; Bcz.im[i] = bczi;
        // (k × Bc)
        const cxr = KY * bczr - KZ * bcyr, cxi = KY * bczi - KZ * bcyi;
        const cyr = KZ * bcxr - KX * bczr, cyi = KZ * bcxi - KX * bczi;
        const czr = KX * bcyr - KY * bcxr, czi = KX * bcyi - KY * bcxi;
        // A(k) = -i (k×B)/k²  →  re = Im(c)/k², im = -Re(c)/k²
        if (ix === 0 && iy === 0 && iz === 0) {
          Ahx.re[i] = Ahx.im[i] = Ahy.re[i] = Ahy.im[i] = Ahz.re[i] = Ahz.im[i] = 0;
        } else {
          Ahx.re[i] = cxi / K2s; Ahx.im[i] = -cxr / K2s;
          Ahy.re[i] = cyi / K2s; Ahy.im[i] = -cyr / K2s;
          Ahz.re[i] = czi / K2s; Ahz.im[i] = -czr / K2s;
        }
      }
    }
  }
  fft3d(Ahx.re, Ahx.im, N, true);
  fft3d(Ahy.re, Ahy.im, N, true);
  fft3d(Ahz.re, Ahz.im, N, true);
  fft3d(Bcx.re, Bcx.im, N, true);
  fft3d(Bcy.re, Bcy.im, N, true);
  fft3d(Bcz.re, Bcz.im, N, true);
  return {
    Ax: Ahx.re, Ay: Ahy.re, Az: Ahz.re,
    Bx: Bcx.re, By: Bcy.re, Bz: Bcz.re,
  };
}

export function whiteheadCharge(tex) {
  const { n: N, dx } = tex;
  const N3 = N * N * N;
  const { Bx: BxRaw, By: ByRaw, Bz: BzRaw, dnx, dny, dnz } = berryB(tex);
  const A = coulombA(BxRaw, ByRaw, BzRaw, dx, N);
  const rhoH = new Float64Array(N3);
  const inv16 = 1 / (16 * PI2);
  let raw = 0;
  let rhoMax = 0;
  const dV = dx * dx * dx;
  for (let i = 0; i < N3; i++) {
    const adotb = A.Ax[i] * A.Bx[i] + A.Ay[i] * A.By[i] + A.Az[i] * A.Bz[i];
    raw += adotb * dV;
    rhoH[i] = adotb * inv16;
    const a = Math.abs(rhoH[i]);
    if (a > rhoMax) rhoMax = a;
  }
  const Q = raw * inv16;
  return {
    Q,
    Q_4pi2: raw / (4 * PI2),
    raw,
    rhoH,
    rhoMax,
    A, Braw: { Bx: BxRaw, By: ByRaw, Bz: BzRaw },
    dnx, dny, dnz,
    normDev: tex.normDev,
    n: N,
    dx,
  };
}

export function faddeevNiemiEnergy(tex, berry, { lengthUnitNm = null } = {}) {
  const s = lengthUnitNm == null ? tex.coreNm : lengthUnitNm;
  const { n: N, dx } = tex;
  const N3 = N * N * N;
  const { dnx, dny, dnz } = berry;
  const { Bx, By, Bz } = berry.Braw;
  let E2s = 0, E4s = 0;
  const dVc = (dx / s) ** 3;
  const inv32 = 1 / (32 * PI2);
  const eps = new Float64Array(N3);
  for (let i = 0; i < N3; i++) {
    let g2 = 0;
    for (let ax = 0; ax < 3; ax++) {
      g2 += dnx[ax][i] ** 2 + dny[ax][i] ** 2 + dnz[ax][i] ** 2;
    }
    const BB = Bx[i] * Bx[i] + By[i] * By[i] + Bz[i] * Bz[i];
    const F2 = 2 * BB;
    const dens = (g2 * s * s + 0.5 * F2 * s ** 4) * inv32;
    eps[i] = dens;
    E2s += g2 * s * s * dVc;
    E4s += 0.5 * F2 * s ** 4 * dVc;
  }
  const E2 = E2s * inv32;
  const E4 = E4s * inv32;
  return { E: E2 + E4, E2, E4, eps, VK: WARD_VK, relaxed: WARD_Q1_MIN };
}

export function vkCheck(E, Q) {
  const bound = WARD_VK * Math.abs(Q) ** 0.75;
  return { E, bound, ratio: bound ? E / bound : Infinity, ok: E >= bound - 1e-9 };
}

export function richardsonQ(Q1, h1, Q2, h2) {
  return (Q2 * h1 * h1 - Q1 * h2 * h2) / (h1 * h1 - h2 * h2);
}

/** Live Whitehead + FN on one grid (used by HUD and verify). */
export function computeHopfion(nGrid, opts = {}) {
  const tex = hopfionTexture(nGrid, opts);
  const wh = whiteheadCharge(tex);
  const fn = faddeevNiemiEnergy(tex, wh, { lengthUnitNm: tex.coreNm });
  const vk = vkCheck(fn.E, wh.Q);
  return { tex, wh, fn, vk };
}

export function csFlux2d(thetaSlice, x, Rnm, nPts = 720) {
  const n = x.length;
  const th = new Float64Array(nPts);
  for (let k = 0; k < nPts; k++) {
    const ang = 2 * Math.PI * k / nPts;
    const px = Rnm * Math.cos(ang);
    const py = Rnm * Math.sin(ang);
    let ix = 0, iy = 0;
    // x is uniform
    const dx = x[1] - x[0];
    ix = Math.round((px - x[0]) / dx);
    iy = Math.round((py - x[0]) / dx);
    ix = Math.max(0, Math.min(n - 1, ix));
    iy = Math.max(0, Math.min(n - 1, iy));
    th[k] = thetaSlice[ix * n + iy];
  }
  const dth = new Float64Array(nPts);
  let cut = 0, cutAbs = 0;
  for (let k = 0; k < nPts; k++) {
    const a = th[k], b = th[(k + 1) % nPts];
    dth[k] = b - a;
    const ad = Math.abs(dth[k]);
    if (ad > cutAbs) { cutAbs = ad; cut = k; }
  }
  let smooth = 0;
  for (let k = 0; k < nPts; k++) if (k !== cut) smooth += dth[k];
  return {
    Q_H_flux: smooth / (2 * Math.PI),
    Q_H_from_cut: -dth[cut] / (2 * Math.PI),
    n_samples: nPts,
  };
}

export function torusThetaGrid(n, Rnm, rnm, QH) {
  const L = Rnm + 4 * rnm;
  const g = makeGrid(n, L);
  const N = n;
  const theta = new Float64Array(N * N); // z = mid plane, layout ix*N+iy
  const mid = (N / 2) | 0;
  const z = g.x[mid];
  for (let ix = 0; ix < N; ix++) {
    for (let iy = 0; iy < N; iy++) {
      const phi = Math.atan2(g.x[iy], g.x[ix]);
      theta[ix * N + iy] = QH * phi;
    }
  }
  return { theta, x: g.x, z, flux: csFlux2d(theta, g.x, Rnm) };
}

export function preimagePoints(tex, p, { tolRad = 0.18, maxPts = 800 } = {}) {
  const pn = Math.hypot(p[0], p[1], p[2]) || 1;
  const px = p[0] / pn, py = p[1] / pn, pz = p[2] / pn;
  const { nx, ny, nz, n: N, x } = tex;
  const pts = [];
  for (let ix = 0; ix < N; ix++) {
    for (let iy = 0; iy < N; iy++) {
      for (let iz = 0; iz < N; iz++) {
        const i = idx3(ix, iy, iz, N);
        const c = Math.max(-1, Math.min(1, nx[i] * px + ny[i] * py + nz[i] * pz));
        if (Math.acos(c) < tolRad) pts.push([x[ix], x[iy], x[iz]]);
      }
    }
  }
  if (pts.length < 8) return { packed: new Float32Array(0), n: 0 };
  const c0 = [0, 0, 0];
  for (const q of pts) { c0[0] += q[0]; c0[1] += q[1]; c0[2] += q[2]; }
  c0[0] /= pts.length; c0[1] /= pts.length; c0[2] /= pts.length;
  // PCA plane via covariance 3×3
  let cxx = 0, cxy = 0, cxz = 0, cyy = 0, cyz = 0, czz = 0;
  for (const q of pts) {
    const dx = q[0] - c0[0], dy = q[1] - c0[1], dz = q[2] - c0[2];
    cxx += dx * dx; cxy += dx * dy; cxz += dx * dz;
    cyy += dy * dy; cyz += dy * dz; czz += dz * dz;
  }
  // power iteration for first two eigenvectors (good enough for a circle)
  const ev = (Mapply, v) => {
    let x = v[0], y = v[1], z = v[2];
    for (let k = 0; k < 24; k++) {
      const r = Mapply(x, y, z);
      const n = Math.hypot(r[0], r[1], r[2]) || 1;
      x = r[0] / n; y = r[1] / n; z = r[2] / n;
    }
    return [x, y, z];
  };
  const M = (x, y, z) => [cxx * x + cxy * y + cxz * z, cxy * x + cyy * y + cyz * z, cxz * x + cyz * y + czz * z];
  const e1 = ev(M, [1, 0.2, 0.1]);
  const M2 = (x, y, z) => {
    const r = M(x, y, z);
    const d = e1[0] * x + e1[1] * y + e1[2] * z;
    return [r[0] - 1e6 * d * e1[0], r[1] - 1e6 * d * e1[1], r[2] - 1e6 * d * e1[2]];
  };
  const e2 = ev(M2, [0.1, 1, 0.2]);
  pts.sort((a, b) => {
    const ax = (a[0] - c0[0]) * e1[0] + (a[1] - c0[1]) * e1[1] + (a[2] - c0[2]) * e1[2];
    const ay = (a[0] - c0[0]) * e2[0] + (a[1] - c0[1]) * e2[1] + (a[2] - c0[2]) * e2[2];
    const bx = (b[0] - c0[0]) * e1[0] + (b[1] - c0[1]) * e1[1] + (b[2] - c0[2]) * e1[2];
    const by = (b[0] - c0[0]) * e2[0] + (b[1] - c0[1]) * e2[1] + (b[2] - c0[2]) * e2[2];
    return Math.atan2(ay, ax) - Math.atan2(by, bx);
  });
  const take = Math.min(maxPts, pts.length);
  const packed = new Float32Array(take * 3);
  for (let i = 0; i < take; i++) {
    const q = pts[Math.round(i * (pts.length - 1) / Math.max(1, take - 1))];
    packed[i * 3] = q[0]; packed[i * 3 + 1] = q[1]; packed[i * 3 + 2] = q[2];
  }
  return { packed, n: take };
}

export function preimageLinking(tex) {
  const a = preimagePoints(tex, [0, 0, 1]);
  const b = preimagePoints(tex, [1, 0, 0]);
  const Lk = gaussLinking(a.packed, b.packed, a.n, b.n);
  return { a, b, Lk };
}

/** Subsampled n̂ glyphs + ρ_H cloud for the Whitehead scene. */
export function visualizationSamples(bundle, { glyphStride = 3, cloudFrac = 0.18 } = {}) {
  const { tex, wh } = bundle;
  const { nx, ny, nz, n: N, x } = tex;
  const glyphs = [];
  for (let ix = 0; ix < N; ix += glyphStride) {
    for (let iy = 0; iy < N; iy += glyphStride) {
      for (let iz = 0; iz < N; iz += glyphStride) {
        const i = idx3(ix, iy, iz, N);
        glyphs.push({
          p: [x[ix], x[iy], x[iz]],
          n: [nx[i], ny[i], nz[i]],
        });
      }
    }
  }
  const thresh = wh.rhoMax * cloudFrac;
  const cloud = [];
  const rgb = [0, 0, 0];
  for (let ix = 0; ix < N; ix++) {
    for (let iy = 0; iy < N; iy++) {
      for (let iz = 0; iz < N; iz++) {
        const i = idx3(ix, iy, iz, N);
        if (Math.abs(wh.rhoH[i]) < thresh) continue;
        nColor(nx[i], ny[i], nz[i], rgb);
        cloud.push(x[ix], x[iy], x[iz], rgb[0], rgb[1], rgb[2], wh.rhoH[i]);
      }
    }
  }
  return { glyphs, cloud };
}

export function defaultLocksGeometry() {
  return { Rnm: LOCKS.R_nm, rnm: LOCKS.r_nm, QH: LOCKS.Q_H };
}

/** Single-point windowed stereographic Hopf n-field (no FFT). */
export function hopfionNAt(Xnm, Ynm, Znm, {
  coreNm = LOCKS.a_nm,
  RcutOverCore = 2.5,
  windowWidthOverCore = 0.4,
} = {}) {
  const X = Xnm / coreNm, Y = Ynm / coreNm, Z = Znm / coreNm;
  const r2 = X * X + Y * Y + Z * Z;
  const r = Math.sqrt(r2);
  const den = 1 + r2;
  const w0 = 2 * X / den, w1 = 2 * Y / den, w2 = 2 * Z / den, w3 = (1 - r2) / den;
  const zr = w0 * w2 + w1 * w3;
  const zi = w1 * w2 - w0 * w3;
  let nxx = 2 * zr, nyy = 2 * zi;
  const a1s = w0 * w0 + w1 * w1, a2s = w2 * w2 + w3 * w3;
  let nzz = a1s - a2s;
  const t = 0.5 * (1 + Math.tanh((r - RcutOverCore) / windowWidthOverCore));
  nxx *= (1 - t); nyy *= (1 - t); nzz = (1 - t) * nzz - t;
  const nn = Math.hypot(nxx, nyy, nzz) || 1e-15;
  return { nx: nxx / nn, ny: nyy / nn, nz: nzz / nn };
}

/** Cheap analytic Hopf density proxy ∝ 1/(1+r²)⁴. Not a live FFT Whitehead integrand. */
export function analyticRhoH(Xnm, Ynm, Znm, coreNm = LOCKS.a_nm) {
  const X = Xnm / coreNm, Y = Ynm / coreNm, Z = Znm / coreNm;
  const d = 1 + X * X + Y * Y + Z * Z;
  return 16 / (d * d * d * d);
}
