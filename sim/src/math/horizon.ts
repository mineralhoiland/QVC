/** Toy acoustic horizon + twin-channel rates */

export function hawkingT(kappa: number, hbarOverKb = 1): number {
  return (hbarOverKb * kappa) / (2 * Math.PI);
}

export function schwingerRate(omega0: number, aEffOverAc: number): number {
  if (aEffOverAc <= 0) return 0;
  return omega0 * Math.exp(-Math.PI / aEffOverAc);
}

export function thermalRate(omega0: number, Epair: number, T: number, kB = 0.0862): number {
  // Epair in meV, T in K, kB meV/K
  if (T <= 0) return 0;
  return omega0 * Math.exp(-Epair / (kB * T));
}

export function crossoverT(EpairMeV: number, aEffOverAc: number, kB = 0.0862): number {
  const expo = Math.PI / aEffOverAc;
  return EpairMeV / (expo * kB);
}
