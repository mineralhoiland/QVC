/** Modified Bessel K0 approximations + vacuum factor F(x)=Σ K0(nx). */

export function besselK0(z: number): number {
  if (z <= 0) return Infinity;
  if (z < 0.5) {
    // Series: K0(z) ≈ -ln(z/2)γ - ln(z/2) I0 + …
    const z2 = (z / 2) * (z / 2);
    let i0 = 1;
    let term = 1;
    for (let k = 1; k < 20; k++) {
      term *= z2 / (k * k);
      i0 += term;
    }
    const gammaE = 0.5772156649;
    return -(Math.log(z / 2) + gammaE) * i0;
  }
  // Uniform asymptotic for moderate/large z
  const inv = 1 / z;
  const poly =
    1 -
    0.125 * inv +
    (9 / 128) * inv * inv -
    (225 / 3072) * inv ** 3 +
    (11025 / 98304) * inv ** 4;
  return Math.sqrt(Math.PI / (2 * z)) * Math.exp(-z) * poly;
}

/** Locked spectral factor F(x)=Σ_{n=1}^{nmax} K0(n x) */
export function vacuumFactorF(x: number, nMax = 50): number {
  let s = 0;
  for (let n = 1; n <= nMax; n++) s += besselK0(n * x);
  return s;
}

/** Tier-C locked form: E = (1/4π²)(ℏcs/r) F(2π r/R) */
export function eVacLocked(rNm: number, RNm: number, hbarCs = 0.07): number {
  const r = rNm / 1000; // µm
  const x = (2 * Math.PI * rNm) / RNm;
  const F = vacuumFactorF(x);
  return (1 / (4 * Math.PI ** 2)) * (hbarCs / r) * F;
}

/** Older dimensional ansatz (audit-failed for 0.244 on fat torus) */
export function eVacOld(rNm: number, RNm: number, hbarCs = 0.07): number {
  const r = rNm / 1000;
  const x = (2 * Math.PI * rNm) / RNm;
  const F = vacuumFactorF(x);
  return (1 / (4 * Math.PI)) * (hbarCs / r) * (RNm / rNm) * F;
}
