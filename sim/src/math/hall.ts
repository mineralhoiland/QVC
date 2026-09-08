/** Hall / anyon helpers: Q_H = 1/N, θ = π/(2N), δσ = e²/(2Nh) */

export function fractionalQH(N: number): number {
  return 1 / N;
}

export function exchangeTheta(N: number): number {
  // k = 2N → θ = π/k = π/(2N)
  return Math.PI / (2 * N);
}

/** Conductance step in units of e²/h */
export function hallStepUnits(N: number): number {
  return 1 / (2 * N);
}

export function hallStaircase(N: number, steps = 10): { nu: number; sigma: number }[] {
  const h = hallStepUnits(N);
  const out = [];
  for (let k = 0; k <= steps; k++) {
    out.push({ nu: k / steps, sigma: k * h });
  }
  return out;
}

export function csInvariant(k: number, N: number): number {
  // CS[A_k] = k²/N mod 1
  return ((k * k) / N) % 1;
}
