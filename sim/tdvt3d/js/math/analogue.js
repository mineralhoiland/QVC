// L2 analogue geometry: Unruh–Visser acoustic metric, T_H diagnostic, Goldstones, H-EIN1.
// Port of qvc/workbench/analogue.py. Nothing here is promoted to L1 or L3.

import { LOCKS } from './locks.js';

export function hbarOverMstarUm2Ps() {
  const ke = LOCKS.Delta0_meV * (LOCKS.r_nm * 1e-3) ** 2; // Δ0 r² with r in µm
  return 2 * ke / LOCKS.hbar_meV_ps;
}

export function soundSpeedUmPs() {
  return LOCKS.hbar_cs_meV_um / LOCKS.hbar_meV_ps;
}

/**
 * Unruh–Visser: ds² = −(c_s² − v²) dt² − 2 v·dx dt + dx².
 * v = (ħ/m*) Q_H φ̂/ρ + drain v_r = −v0 (R/ρ).
 * drain = 0 → sonic radius sits inside the healing core (report as no exterior horizon).
 */
export function acousticMetric({
  nGrid = 96,
  Lum = 0.12,
  QH = LOCKS.Q_H,
  drainV0OverCs = 0,
  drainRum = LOCKS.R_nm * 1e-3,
} = {}) {
  const cs = soundSpeedUmPs();
  const hm = hbarOverMstarUm2Ps();
  const Rd = drainRum;
  const x = new Float64Array(nGrid);
  const dx = (2 * Lum) / (nGrid - 1);
  for (let i = 0; i < nGrid; i++) x[i] = -Lum + i * dx;
  const gtt = new Float64Array(nGrid * nGrid);
  const mach = new Float64Array(nGrid * nGrid);
  const vx = new Float64Array(nGrid * nGrid);
  const vy = new Float64Array(nGrid * nGrid);
  const rhoMin = 0.25 * dx;
  for (let iy = 0; iy < nGrid; iy++) {
    for (let ix = 0; ix < nGrid; ix++) {
      const X = x[ix], Y = x[iy];
      const rho = Math.max(Math.hypot(X, Y), rhoMin);
      const vphi = hm * QH / rho;
      const vr = -drainV0OverCs * cs * (Rd / rho);
      const vxi = vr * X / rho - vphi * Y / rho;
      const vyi = vr * Y / rho + vphi * X / rho;
      const v2 = vxi * vxi + vyi * vyi;
      const i = iy * nGrid + ix;
      vx[i] = vxi; vy[i] = vyi;
      gtt[i] = -(cs * cs - v2);
      mach[i] = Math.sqrt(v2) / cs;
    }
  }
  const rhoHVortex = hm * QH / cs;
  let rhoH, kappa;
  if (drainV0OverCs > 0) {
    const nLine = 4000;
    const r0 = 0.2 * dx, r1 = Lum;
    let found = false;
    let rPrev = r0;
    const fOf = (r) => {
      const vr = drainV0OverCs * cs * Rd / r;
      const vp = hm * QH / r;
      return Math.hypot(vr, vp) - cs;
    };
    let fPrev = fOf(r0);
    for (let k = 1; k < nLine; k++) {
      const r = r0 + (r1 - r0) * k / (nLine - 1);
      const f = fOf(r);
      if (fPrev === 0 || fPrev * f < 0) {
        const t = f - fPrev !== 0 ? -fPrev / (f - fPrev) : 0;
        rhoH = rPrev + t * (r - rPrev);
        // κ ≈ |d|v_r|/dρ| at horizon
        const dr = Math.max(1e-8, r - rPrev);
        const vr0 = drainV0OverCs * cs * Rd / Math.max(rhoH, 1e-12);
        const vr1 = drainV0OverCs * cs * Rd / (rhoH + dr);
        kappa = Math.abs((vr1 - vr0) / dr);
        found = true;
      }
      rPrev = r; fPrev = f;
    }
    if (!found) { rhoH = NaN; kappa = NaN; }
  } else {
    rhoH = rhoHVortex;
    kappa = rhoH > 0 ? hm * QH / (rhoH * rhoH) : NaN;
  }
  const T_H_mK = Number.isFinite(kappa)
    ? 1e3 * LOCKS.hbar_meV_ps * kappa / (2 * Math.PI * LOCKS.kB_meV_per_K)
    : NaN;
  const rCoreUm = LOCKS.r_nm * 1e-3;
  const horizonInsideCore = rhoH < rCoreUm;
  const exteriorHorizon = drainV0OverCs > 0 && Number.isFinite(rhoH) && !horizonInsideCore;
  return {
    x, gtt, mach, vx, vy, nGrid, Lum, cs, hm,
    rho_h_um: rhoH,
    rho_h_nm: rhoH * 1e3,
    rho_h_vortex_nm: rhoHVortex * 1e3,
    kappa_ps_inv: kappa,
    T_H_mK,
    T_star_mK: LOCKS.T_star_mK,
    drain: drainV0OverCs,
    horizonInsideCore,
    exteriorHorizon,
    note: 'Analogue kinematics only (T_H = κ/2π). Not Einstein gravity; T* is a Schwinger-channel lock, not T_H.',
  };
}

export function goldstoneBranches({
  cLA = 21400, cTA = 13600, aNm = 50, torsionZeta = 0.08,
  magnonGap = 0.05, cMagnon = 8000, nK = 200,
} = {}) {
  const conv = LOCKS.hbar_meV_ps * 1e-3;
  const kMax = Math.PI / aNm;
  const k = new Float64Array(nK);
  const LA = new Float64Array(nK), TA = new Float64Array(nK);
  const m1 = new Float64Array(nK), m2 = new Float64Array(nK);
  for (let i = 0; i < nK; i++) {
    k[i] = kMax * i / (nK - 1);
    LA[i] = conv * cLA * k[i];
    TA[i] = conv * cTA * k[i] * (1 - torsionZeta * (k[i] * aNm / Math.PI) ** 2);
    m1[i] = Math.hypot(magnonGap, conv * cMagnon * k[i]);
    m2[i] = Math.hypot(1.3 * magnonGap, conv * cMagnon * k[i]);
  }
  return { k, LA, TA, TA2: TA, magnon1: m1, magnon2: m2, kMax };
}

/** Frobenius overlap of a TA (helicity ±1) strain with TT helicity-±2 polarisations → 0. */
export function hEin1Overlap(kDir = [0, 0, 1], eDir = [1, 0, 0]) {
  const kn = Math.hypot(...kDir) || 1;
  const k = kDir.map((v) => v / kn);
  let e = eDir.slice();
  const ke = k[0] * e[0] + k[1] * e[1] + k[2] * e[2];
  e = e.map((v, i) => v - k[i] * ke);
  const en = Math.hypot(...e) || 1;
  e = e.map((v) => v / en);
  const strain = Array.from({ length: 3 }, (_, i) =>
    Array.from({ length: 3 }, (_, j) => 0.5 * (k[i] * e[j] + e[i] * k[j])));
  const f = [
    k[1] * e[2] - k[2] * e[1],
    k[2] * e[0] - k[0] * e[2],
    k[0] * e[1] - k[1] * e[0],
  ];
  const ePlus = Array.from({ length: 3 }, (_, i) =>
    Array.from({ length: 3 }, (_, j) => e[i] * e[j] - f[i] * f[j]));
  const eCross = Array.from({ length: 3 }, (_, i) =>
    Array.from({ length: 3 }, (_, j) => e[i] * f[j] + f[i] * e[j]));
  const fro = (A, B) => {
    let s = 0, nA = 0, nB = 0;
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
      s += A[i][j] * B[i][j];
      nA += A[i][j] * A[i][j];
      nB += B[i][j] * B[i][j];
    }
    return s / (Math.sqrt(nA) * Math.sqrt(nB));
  };
  const ovPlus = fro(strain, ePlus);
  const ovCross = fro(strain, eCross);
  return {
    overlap_plus: ovPlus,
    overlap_cross: ovCross,
    isGravitonLike: Math.max(Math.abs(ovPlus), Math.abs(ovCross)) > 1e-8,
  };
}
