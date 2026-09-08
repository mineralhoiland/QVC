// L1 laboratory locks (Option A*) — mirrored from qvc/workbench/locks.py (get_locks()).
// Every number here is a LOCK: it is displayed, never re-fitted. data/locks.json is the
// machine-generated copy produced by the Python solver chain; this module is the JS mirror
// used synchronously by the scenes and by verify.js.

export const LAYER = {
  L1: { id: 'L1', title: 'L1 · laboratory lock (QVC / TEVC)', short: 'lab lock' },
  L2: { id: 'L2', title: 'L2 · analogue geometry (TDVT) — diagnosed, not Einstein', short: 'analogue geometry' },
  L3: { id: 'L3', title: 'L3 · spacetime ontology — schematic, not locked', short: 'ontology (schematic)' },
};

export const LOCKS = Object.freeze({
  N: 5,
  r_nm: 8.0,
  R_nm: 50.0,
  a_nm: 50.0, // BCC lattice constant a = xi = R
  hbar_cs_meV_um: 0.07,
  Delta0_meV: 0.6,
  x_aspect: 1.0053096491487339,
  F_bessel: 0.5806751779874608,
  E_vac_meV: 0.12870089825559464,
  ell: -2.7958108397939263,
  alpha_c: 0.18179560182195154,
  delta_c: 0.18179560182195154,
  Delta_c_meV: 0.10907736109317091,
  max_static_reduction_pct: 81.82043981780484,
  lambda_fold: 0.8475260279578454,
  lambda_star: 0.7627734251620608,
  alpha_star: 0.16361604163975638,
  alpha_star_over_ac: 0.9,
  Delta_star_meV: 0.23244395506499999,
  delta_star: 0.387406591775,
  gap_reduction_star_pct: 61.25934082249999,
  g0_star_meV: 0.09816962498385381,
  lambda_max_meV: 0.39267849993541526,
  lambda_min_meV: -0.09816962498385381,
  lambda_min_multiplicity: 4,
  Delta_cat_meV: 0.5018303750161461,
  Q_H: 0.2,
  theta_anyon_rad: Math.PI / 5, // locks.py: anyonic_exchange_angle(N) = pi/N
  hall_step_e2_over_h: 0.1,
  hall_step_siemens: 3.87404614e-6,
  hall_base_e2_over_h: 5.0,
  chern_simons_level_k: 10,
  T_star_mK: 118.7,
  E_pair_meV: 0.1963,
  omega0_star_ps_inv: 0.7062890688274581,
  nu0_star_THz: 0.11240939655566184,
  hbar_meV_ps: 0.6582119569,
  kB_meV_per_K: 0.08617333262,
  E2_OVER_H_SIEMENS: 3.87404614e-5,
});

// Reject list (docs/PHASE_TWO_PLAN.md "DEPRECATED (do NOT use)"). Displaying one of these throws.
export const DEPRECATED = [
  { key: 'E_Cas_0p244_meV', value: 0.244, unit: 'meV', why: 'withdrawn Casimir number; Option A* gives E_vac = 0.1287 meV' },
  { key: 'gap_reduction_92pct', value: 92.0, unit: '%', why: 'lies below the fold; max static reduction is 81.8%' },
  { key: 'Delta_1ps_0p032_meV', value: 0.032, unit: 'meV', why: 'scheduled live gap from the withdrawn solver' },
  { key: 'Delta_QVC_0p389_meV', value: 0.389, unit: 'meV', why: '35% static loop at E_Cas = 0.244 meV' },
  { key: 'T_star_285_290_mK', value: 285.0, unit: 'mK', why: 'superseded by T* = 118.7 mK' },
  { key: 'A_ZPF_xi_0p296', value: 0.296, unit: 'dimensionless', why: 'A_ZPF*xi is not the coupling constant' },
  { key: 'lambda_min_minus_4g0', value: -4.0, unit: 'g0', why: 'democratic G has lambda_min = -g0 (mult N-1)' },
  { key: 'lambda_min_minus_5g0', value: -5.0, unit: 'g0', why: 'see lambda_min_minus_4g0' },
];

export class DeprecatedNumberError extends Error {}

export function assertNotDeprecated(value, unit, relTol = 1e-3) {
  for (const d of DEPRECATED) {
    if (d.unit !== unit) continue;
    if (Math.abs(value - d.value) <= Math.max(relTol * Math.abs(d.value), 1e-12)) {
      throw new DeprecatedNumberError(`${value} ${unit} matches withdrawn number '${d.key}': ${d.why}`);
    }
  }
  return true;
}

// Layer-dependent locks for a different layer number N (the lens-space identities scale exactly).
export function lensLocks(N) {
  return {
    N,
    Q_H: 1 / N,
    theta_anyon_rad: Math.PI / N,
    theta_paper_k2N_rad: Math.PI / (2 * N), // Paper C Sec. anyons quotes k_CS = 2N -> theta = pi/(2N)
    hall_step_e2_over_h: 1 / (2 * N),
    hall_step_siemens: LOCKS.E2_OVER_H_SIEMENS / (2 * N),
    chern_simons_level_k: 2 * N,
    lambda_max_meV: (N - 1) * LOCKS.g0_star_meV,
    lambda_min_meV: -LOCKS.g0_star_meV,
    lambda_min_multiplicity: N - 1,
  };
}

export function fmt(v, digits = 4) {
  if (!Number.isFinite(v)) return 'nan';
  const a = Math.abs(v);
  if (a !== 0 && (a < 1e-3 || a >= 1e5)) return v.toExponential(Math.max(1, digits - 1));
  return v.toFixed(digits);
}
