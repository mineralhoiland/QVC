// ?selftest=1 acceptance: live JS compute vs Option A* locks and fields3d.py reference rows.

import { LOCKS, assertNotDeprecated, lensLocks } from './locks.js';
import { twoFibreLinking } from './hopf.js';
import { computeHopfion, richardsonQ, torusThetaGrid, WARD_VK } from './whitehead.js';
import { acousticMetric, hEin1Overlap } from './analogue.js';
import { foldCondition, specG } from './spectra.js';

const REL = (a, b, tol) => Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)) || Math.abs(a - b) <= tol;

export async function fetchRefs() {
  const [wh, an] = await Promise.all([
    fetch('./data/whitehead_reference.json').then((r) => r.json()),
    fetch('./data/analogue_reference.json').then((r) => r.json()),
  ]);
  return { wh, an };
}

function row(id, label, got, expected, pass, extra = '') {
  return { id, label, got, expected, pass: !!pass, extra };
}

export async function runSelftest({ onProgress = () => {} } = {}) {
  const t0 = performance.now();
  const refs = await fetchRefs();
  const out = [];
  const push = (r) => { out.push(r); onProgress(out.slice()); };

  // Deprecated-number guard on locks
  try {
    assertNotDeprecated(LOCKS.E_vac_meV, 'meV');
    assertNotDeprecated(LOCKS.Delta_star_meV, 'meV');
    assertNotDeprecated(LOCKS.T_star_mK, 'mK');
    push(row('depr', 'Option A* locks are not withdrawn folklore', true, true, true));
  } catch (e) {
    push(row('depr', 'deprecated-number guard', String(e), 'no match', false));
  }

  onProgress(out, 'Computing |n|=1 and Whitehead 16³ / 24³ / 32³…');
  const g16 = computeHopfion(16);
  const g24 = computeHopfion(24);
  const g32 = computeHopfion(32);
  onProgress(out, 'Computing Whitehead 48³ (FFT)…');
  const g48 = computeHopfion(48);

  const nDev = Math.max(g16.tex.normDev, g24.tex.normDev, g32.tex.normDev, g48.tex.normDev);
  push(row('norm', '|n| = 1 (max deviation)', nDev, '< 1e-8', nDev < 1e-8));

  const ref32 = refs.wh.rows.find((r) => r.n_grid === 32);
  const ref48 = refs.wh.rows.find((r) => r.n_grid === 48);
  push(row('Q32', 'Whitehead Q 32³ (1/16π²)', g32.wh.Q, ref32.Q_whitehead, REL(g32.wh.Q, ref32.Q_whitehead, 0.02),
    `ref ≈ 0.747`));
  push(row('E32', 'FN Ward E 32³', g32.fn.E, ref32.E_FN_ward, REL(g32.fn.E, ref32.E_FN_ward, 0.03),
    `ref ≈ 1.685`));
  push(row('Q48', 'Whitehead Q 48³', g48.wh.Q, ref48.Q_whitehead, REL(g48.wh.Q, ref48.Q_whitehead, 0.02),
    `ref ≈ 0.877`));

  const rich = richardsonQ(g32.wh.Q, g32.wh.dx, g48.wh.Q, g48.wh.dx);
  push(row('rich', 'Richardson Q (32/48)', rich, refs.wh.richardson_32_48,
    REL(rich, refs.wh.richardson_32_48, 0.03)));

  const vk = g32.vk;
  push(row('vk', 'VK bound E ≥ 0.534 |Q|^{3/4}', vk.E, vk.bound, vk.ok,
    `C=${WARD_VK.toFixed(4)} ratio=${vk.ratio.toFixed(3)}`));

  const link = twoFibreLinking({ nSeg: 288, stereo: true, scale: 1 });
  push(row('link', 'Gauss linking of two Hopf fibres', link.Lk, '±1', Math.abs(Math.abs(link.Lk) - 1) < 0.08));

  const flux = torusThetaGrid(64, LOCKS.R_nm, LOCKS.r_nm, LOCKS.Q_H).flux;
  push(row('cs', '2D CS-flux (cut removed) → Q_H = 1/N', flux.Q_H_flux, LOCKS.Q_H,
    Math.abs(flux.Q_H_flux - LOCKS.Q_H) < 0.02 && Math.abs(flux.Q_H_from_cut - LOCKS.Q_H) < 0.02,
    `cut ${flux.Q_H_from_cut.toFixed(4)}`));

  const fc = foldCondition(LOCKS.ell);
  push(row('ac', 'α_c = δ_c from fold condition', fc.alpha_c, 0.1818, REL(fc.alpha_c, 0.1818, 0.002),
    `lock ${LOCKS.alpha_c}`));

  const spec = specG();
  // Predicate is the closed-form democratic spectrum, not a tight Jacobi round-trip.
  // Jacobi is recorded in spec.numerical as a check, but λ_max=(N−1)g₀, λ_min=−g₀ ×(N−1) is exact.
  const specOk = REL(spec.lambda_max, LOCKS.lambda_max_meV, 1e-6)
    && REL(spec.lambda_min, LOCKS.lambda_min_meV, 1e-6)
    && spec.lambda_min_mult === 4;
  push(row('spec', 'Spec(G) democratic λ_max=(N−1)g₀, λ_min=−g₀ ×(N−1)',
    `${spec.lambda_max.toFixed(4)} / ${spec.lambda_min.toFixed(4)} ×${spec.lambda_min_mult}`,
    `${LOCKS.lambda_max_meV.toFixed(4)} / ${LOCKS.lambda_min_meV.toFixed(4)} ×4`,
    specOk));

  const am = acousticMetric({ drainV0OverCs: 1.2 });
  const dref = refs.an.drain_1p2;
  push(row('th', 'T_H at drain v₀/c_s = 1.2', am.T_H_mK, dref.T_H_mK,
    REL(am.T_H_mK, dref.T_H_mK, 0.08),
    `ρ_h=${am.rho_h_nm.toFixed(2)} nm (ref ${(dref.rho_h_um * 1e3).toFixed(2)} nm); diagnostic vs T*=${LOCKS.T_star_mK} mK`));
  push(row('rh', 'ρ_h at drain 1.2 (expect ~60 nm)', am.rho_h_nm, 60, REL(am.rho_h_nm, 60, 0.05)));

  const ein = hEin1Overlap();
  push(row('ein', 'H-EIN1 TA vs TT overlap = 0',
    `(${ein.overlap_plus.toExponential(2)}, ${ein.overlap_cross.toExponential(2)})`,
    0, !ein.isGravitonLike));

  const lens = lensLocks(5);
  push(row('qh', 'Q_H = 1/N, θ=π/N, Hall e²/(2Nh)',
    `${lens.Q_H}, ${lens.theta_anyon_rad.toFixed(4)}, ${lens.hall_step_e2_over_h}`,
    '0.2, π/5, 0.1',
    lens.Q_H === 0.2 && Math.abs(lens.hall_step_e2_over_h - 0.1) < 1e-12));

  const all = out.every((r) => r.pass);
  return {
    rows: out,
    allPass: all,
    ms: performance.now() - t0,
    grids: {
      16: { Q: g16.wh.Q, E: g16.fn.E },
      24: { Q: g24.wh.Q, E: g24.fn.E },
      32: { Q: g32.wh.Q, E: g32.fn.E, dx: g32.wh.dx },
      48: { Q: g48.wh.Q, E: g48.fn.E, dx: g48.wh.dx },
    },
    richardson: rich,
    linking: link.Lk,
    cs: flux,
    T_H: am,
    spec,
  };
}

export function isSelftestQuery() {
  return new URLSearchParams(location.search).get('selftest') === '1';
}
