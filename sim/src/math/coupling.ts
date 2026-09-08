/** Democratic Spec(G): G = g0(J-I). λ_max=(N-1)g0, λ_min=-g0 (mult N-1). */

export function democraticMatrix(N: number, g0: number): number[][] {
  const G = Array.from({ length: N }, () => Array(N).fill(0));
  for (let i = 0; i < N; i++) {
    for (let j = 0; j < N; j++) {
      G[i][j] = i === j ? 0 : g0;
    }
  }
  return G;
}

export function democraticSpectrum(N: number, g0: number) {
  return {
    lambdaMax: (N - 1) * g0,
    lambdaMin: -g0,
    multiplicityMin: N - 1,
    eigenvalues: [...Array(N - 1).fill(-g0), (N - 1) * g0],
  };
}

/** Jacobi-style Helmert dark basis (sum-zero). */
export function darkBasis(N: number): number[][] {
  const basis: number[][] = [];
  for (let k = 1; k < N; k++) {
    const v = Array(N).fill(0);
    for (let i = 0; i < k; i++) v[i] = 1;
    v[k] = -k;
    const norm = Math.hypot(...v);
    basis.push(v.map((x) => x / norm));
  }
  return basis;
}
