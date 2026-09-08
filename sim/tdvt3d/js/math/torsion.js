// Axial torsion proxy T ∝ J_Hopf ∝ ρ_H. H-TOR1 is parked: never fed into G_ij.
// Parallelogram non-closure is the Cartan definition of torsion (TEGR / Paper C § torsion).

/**
 * Normalised axial torsion from Hopf density. Toggle only changes the annotation.
 */
export function torsionProxy(rhoH, { hTor1 = false } = {}) {
  let scale = 0;
  for (let i = 0; i < rhoH.length; i++) {
    const a = Math.abs(rhoH[i]);
    if (a > scale) scale = a;
  }
  scale = scale || 1;
  const T = new Float64Array(rhoH.length);
  for (let i = 0; i < rhoH.length; i++) T[i] = rhoH[i] / scale;
  return {
    T,
    scale,
    hTor1,
    feedsGij: false,
    note: hTor1
      ? 'H-TOR1 toggle is ON: this remains a parked experiment; Spec(G) shown elsewhere is unchanged.'
      : 'H-TOR1 parked (toggle off).',
  };
}

/**
 * Infinitesimal parallelogram spanned by e1, e2 fails to close by T(e1,e2) ∝ ρ_H (e1 × e2).
 * Returns the two paths and the closure gap vector.
 */
export function parallelogramGap(origin, e1, e2, Taxial, { amp = 8 } = {}) {
  const gap = [
    amp * Taxial * (e1[1] * e2[2] - e1[2] * e2[1]),
    amp * Taxial * (e1[2] * e2[0] - e1[0] * e2[2]),
    amp * Taxial * (e1[0] * e2[1] - e1[1] * e2[0]),
  ];
  const a = origin;
  const b = [a[0] + e1[0], a[1] + e1[1], a[2] + e1[2]];
  const c = [b[0] + e2[0], b[1] + e2[1], b[2] + e2[2]];
  const d = [a[0] + e2[0], a[1] + e2[1], a[2] + e2[2]];
  const c2 = [d[0] + e1[0] + gap[0], d[1] + e1[1] + gap[1], d[2] + e1[2] + gap[2]];
  return { a, b, c, d, c2, gap };
}

export const H_TOR1_STATUS = 'PARKED';
