// Fold map, democratic Spec(G), Mexican-hat potential, FN E vs Q / VK.
// L1 locks. SYK is recorded in Paper C as a chaos analogy only — not identified with the lattice.

import { LOCKS, lensLocks } from './locks.js';
import { WARD_VK, WARD_Q1_MIN } from './whitehead.js';

export function foldCondition(ell) {
  let lo = 1e-12, hi = 1;
  const f = (d) => d * (1 - Math.log(d) - ell) - 1;
  for (let i = 0; i < 120; i++) {
    const mid = 0.5 * (lo + hi);
    if (f(mid) > 0) hi = mid; else lo = mid;
  }
  const dc = 0.5 * (lo + hi);
  return { delta_c: dc, alpha_c: dc };
}

export function foldResidual(delta, alpha, ell) {
  if (delta <= 0) return NaN;
  return delta - (1 + alpha * (Math.log(delta) + ell));
}

export function solveFoldBranches(alpha, ell, nScan = 4000) {
  const { delta_c, alpha_c } = foldCondition(ell);
  if (alpha <= 0) return { phys: 1, unst: null, alpha_c, delta_c };
  if (alpha >= alpha_c * (1 - 1e-12)) return { phys: null, unst: null, alpha_c, delta_c };
  const roots = [];
  const dLo = 1e-6, dHi = 1;
  let prevD = dLo, prevF = foldResidual(prevD, alpha, ell);
  for (let i = 1; i <= nScan; i++) {
    const d = dLo + (dHi - dLo) * i / nScan;
    const f = foldResidual(d, alpha, ell);
    if (prevF === 0) roots.push(prevD);
    else if (prevF * f < 0) {
      const t = -prevF / (f - prevF);
      roots.push(prevD + t * (d - prevD));
    }
    prevD = d; prevF = f;
  }
  roots.sort((a, b) => b - a);
  return { phys: roots[0] ?? null, unst: roots[1] ?? null, alpha_c, delta_c };
}

export function foldCurve(ell = LOCKS.ell, n = 120) {
  const { delta_c, alpha_c } = foldCondition(ell);
  const alphaOver = new Float64Array(n);
  const phys = new Float64Array(n);
  const unst = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    const a = alpha_c * 0.995 * i / Math.max(n - 1, 1);
    const br = solveFoldBranches(a, ell);
    alphaOver[i] = a / alpha_c;
    phys[i] = br.phys == null ? NaN : br.phys;
    unst[i] = br.unst == null ? NaN : br.unst;
  }
  return {
    alphaOver, phys, unst, alpha_c, delta_c, ell,
    operating: { alpha_over_ac: LOCKS.alpha_star_over_ac, delta: LOCKS.delta_star },
  };
}

export function democraticG(N, g0) {
  const G = Array.from({ length: N }, (_, i) =>
    Array.from({ length: N }, (_, j) => (i === j ? 0 : g0)));
  return G;
}

/** Exact Spec(G) for G = g0 (J − I). Also Jacobi-rotate a numerical check. */
export function specG(N = LOCKS.N, g0 = LOCKS.g0_star_meV) {
  const exact = {
    lambda_max: (N - 1) * g0,
    lambda_max_mult: 1,
    lambda_min: -g0,
    lambda_min_mult: N - 1,
  };
  const G = democraticG(N, g0);
  const ev = jacobiEigenvalues(G);
  ev.sort((a, b) => b - a);
  return { ...exact, numerical: ev, N, g0, lens: lensLocks(N) };
}

function jacobiEigenvalues(A, tol = 1e-12, maxSweeps = 40) {
  const n = A.length;
  const M = A.map((row) => row.slice());
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let max = 0, p = 0, q = 1;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      const a = Math.abs(M[i][j]);
      if (a > max) { max = a; p = i; q = j; }
    }
    if (max < tol) break;
    const app = M[p][p], aqq = M[q][q], apq = M[p][q];
    const tau = (aqq - app) / (2 * apq);
    const t = Math.sign(tau) / (Math.abs(tau) + Math.hypot(1, tau));
    const c = 1 / Math.hypot(1, t), s = t * c;
    M[p][p] = app - t * apq;
    M[q][q] = aqq + t * apq;
    M[p][q] = M[q][p] = 0;
    for (let k = 0; k < n; k++) {
      if (k === p || k === q) continue;
      const mkp = M[k][p], mkq = M[k][q];
      M[k][p] = M[p][k] = c * mkp - s * mkq;
      M[k][q] = M[q][k] = s * mkp + c * mkq;
    }
  }
  return M.map((_, i) => M[i][i]);
}

export function mexicanHat(DeltaMeV, { Delta0 = LOCKS.Delta0_meV, n0 = 1, nR = 80, nPhi = 96 } = {}) {
  const a = -DeltaMeV;
  const b = Delta0 / n0;
  const rhoMin = a < 0 ? Math.sqrt(-a / b) : 0;
  const depth = a < 0 ? a * a / (2 * b) : 0;
  const rmax = 1.7 * Math.max(rhoMin, 1);
  const pos = new Float32Array(nR * nPhi * 3);
  const col = new Float32Array(nR * nPhi * 3);
  let p = 0;
  for (let ir = 0; ir < nR; ir++) {
    const rho = rmax * ir / (nR - 1);
    const V = a * rho * rho + 0.5 * b * rho ** 4;
    for (let ip = 0; ip < nPhi; ip++) {
      const ph = 2 * Math.PI * ip / nPhi;
      pos[p * 3] = rho * Math.cos(ph);
      pos[p * 3 + 1] = rho * Math.sin(ph);
      pos[p * 3 + 2] = V;
      const t = (V - (-depth)) / Math.max(1e-9, depth + Math.abs(a) * rmax * rmax);
      col[p * 3] = 0.2 + 0.8 * Math.max(0, 1 - t);
      col[p * 3 + 1] = 0.35;
      col[p * 3 + 2] = 0.55 + 0.4 * Math.min(1, Math.max(0, t));
      p++;
    }
  }
  return { pos, col, nR, nPhi, rhoMin, depth, a, b, DeltaMeV, rmax };
}

export function vkCurve(Qmax = 4, n = 80) {
  const Q = new Float64Array(n);
  const E = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    Q[i] = Qmax * i / (n - 1);
    E[i] = WARD_VK * Math.abs(Q[i]) ** 0.75;
  }
  return { Q, E, C: WARD_VK, relaxedQ1: WARD_Q1_MIN };
}

export const SYK_NOTE =
  'Paper C § syk records a chaos analogy with SYK and explicitly does not identify the hopfion lattice with SYK. ' +
  'No hopfion-lattice OTOC or Lyapunov exponent is computed. Not a lock.';
