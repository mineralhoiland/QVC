/** Saddle-node fold map: δ = 1 + α(ln δ + ℓ) */

export function foldResidual(delta: number, alpha: number, ell: number): number {
  if (delta <= 0) return Infinity;
  return delta - (1 + alpha * (Math.log(delta) + ell));
}

export function foldCritical(ell: number): { alphaC: number; deltaC: number } {
  // δ_c (1 - ln δ_c - ℓ) = 1, α_c = δ_c
  let lo = 1e-4;
  let hi = 1;
  for (let i = 0; i < 60; i++) {
    const mid = 0.5 * (lo + hi);
    const f = mid * (1 - Math.log(mid) - ell) - 1;
    if (f > 0) hi = mid;
    else lo = mid;
  }
  const deltaC = 0.5 * (lo + hi);
  return { alphaC: deltaC, deltaC };
}

export function solveFoldBranches(
  alpha: number,
  ell: number,
): { phys: number | null; unst: number | null } {
  const { alphaC } = foldCritical(ell);
  if (alpha <= 0) return { phys: 1, unst: null };
  if (alpha >= alphaC * 0.999999) return { phys: null, unst: null };

  const roots: number[] = [];
  const ds = 4000;
  let prev = foldResidual(1e-4, alpha, ell);
  for (let i = 1; i <= ds; i++) {
    const d = 1e-4 + ((1 - 1e-4) * i) / ds;
    const f = foldResidual(d, alpha, ell);
    if (prev === 0) roots.push(d);
    else if (prev * f < 0) {
      const t = -prev / (f - prev);
      roots.push(d - (1 - 1e-4) / ds + t * ((1 - 1e-4) / ds));
    }
    prev = f;
  }
  roots.sort((a, b) => b - a);
  return { phys: roots[0] ?? null, unst: roots[1] ?? null };
}

export function foldCurve(ell: number, n = 80) {
  const { alphaC, deltaC } = foldCritical(ell);
  const alphas: number[] = [];
  const phys: (number | null)[] = [];
  const unst: (number | null)[] = [];
  for (let i = 0; i < n; i++) {
    const a = (alphaC * 0.995 * i) / (n - 1);
    const br = solveFoldBranches(a, ell);
    alphas.push(a / alphaC);
    phys.push(br.phys);
    unst.push(br.unst);
  }
  return { alphas, phys, unst, alphaC, deltaC };
}
