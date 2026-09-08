// Complex FFT utilities: radix-2 iterative for power-of-two sizes, cached DFT matrix otherwise.
// 3D transforms operate on separate Float64Array re/im with layout idx = (ix*N + iy)*N + iz.

const dftCache = new Map();

function isPow2(n) { return n > 0 && (n & (n - 1)) === 0; }

function fftPow2(re, im, inverse) {
  const n = re.length;
  // bit reversal
  for (let i = 1, j = 0; i < n; i++) {
    let bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) {
      let t = re[i]; re[i] = re[j]; re[j] = t;
      t = im[i]; im[i] = im[j]; im[j] = t;
    }
  }
  for (let len = 2; len <= n; len <<= 1) {
    const ang = (2 * Math.PI / len) * (inverse ? 1 : -1);
    const wr = Math.cos(ang), wi = Math.sin(ang);
    for (let i = 0; i < n; i += len) {
      let cr = 1, ci = 0;
      const half = len >> 1;
      for (let j = 0; j < half; j++) {
        const ur = re[i + j], ui = im[i + j];
        const xr = re[i + j + half], xi = im[i + j + half];
        const vr = xr * cr - xi * ci, vi = xr * ci + xi * cr;
        re[i + j] = ur + vr; im[i + j] = ui + vi;
        re[i + j + half] = ur - vr; im[i + j + half] = ui - vi;
        const ncr = cr * wr - ci * wi; ci = cr * wi + ci * wr; cr = ncr;
      }
    }
  }
  if (inverse) { for (let i = 0; i < n; i++) { re[i] /= n; im[i] /= n; } }
}

function dftMatrix(n, inverse) {
  const key = n * (inverse ? -1 : 1);
  if (dftCache.has(key)) return dftCache.get(key);
  const cr = new Float64Array(n * n), ci = new Float64Array(n * n);
  const s = inverse ? 1 : -1;
  for (let k = 0; k < n; k++) for (let j = 0; j < n; j++) {
    const ang = s * 2 * Math.PI * ((k * j) % n) / n;
    cr[k * n + j] = Math.cos(ang); ci[k * n + j] = Math.sin(ang);
  }
  const m = { cr, ci };
  dftCache.set(key, m);
  return m;
}

function dftGeneral(re, im, inverse) {
  const n = re.length;
  const { cr, ci } = dftMatrix(n, inverse);
  const outR = new Float64Array(n), outI = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    let sr = 0, si = 0;
    const row = k * n;
    for (let j = 0; j < n; j++) {
      const a = cr[row + j], b = ci[row + j];
      sr += re[j] * a - im[j] * b; si += re[j] * b + im[j] * a;
    }
    outR[k] = sr; outI[k] = si;
  }
  const inv = inverse ? 1 / n : 1;
  for (let k = 0; k < n; k++) { re[k] = outR[k] * inv; im[k] = outI[k] * inv; }
}

export function fft1d(re, im, inverse = false) {
  if (isPow2(re.length)) fftPow2(re, im, inverse); else dftGeneral(re, im, inverse);
}

/** In-place 3D FFT on cube of side N (re, im Float64Array length N^3). */
export function fft3d(re, im, N, inverse = false) {
  const lr = new Float64Array(N), li = new Float64Array(N);
  const N2 = N * N;
  // axis 2 (stride 1)
  for (let a = 0; a < N2; a++) {
    const base = a * N;
    for (let i = 0; i < N; i++) { lr[i] = re[base + i]; li[i] = im[base + i]; }
    fft1d(lr, li, inverse);
    for (let i = 0; i < N; i++) { re[base + i] = lr[i]; im[base + i] = li[i]; }
  }
  // axis 1 (stride N)
  for (let ix = 0; ix < N; ix++) for (let iz = 0; iz < N; iz++) {
    const base = ix * N2 + iz;
    for (let i = 0; i < N; i++) { lr[i] = re[base + i * N]; li[i] = im[base + i * N]; }
    fft1d(lr, li, inverse);
    for (let i = 0; i < N; i++) { re[base + i * N] = lr[i]; im[base + i * N] = li[i]; }
  }
  // axis 0 (stride N^2)
  for (let a = 0; a < N2; a++) {
    for (let i = 0; i < N; i++) { lr[i] = re[a + i * N2]; li[i] = im[a + i * N2]; }
    fft1d(lr, li, inverse);
    for (let i = 0; i < N; i++) { re[a + i * N2] = lr[i]; im[a + i * N2] = li[i]; }
  }
}

/** numpy.fft.fftfreq(N, d) * 2*pi */
export function angularFreqs(N, dx) {
  const k = new Float64Array(N);
  const half = Math.floor((N - 1) / 2) + 1;
  for (let i = 0; i < half; i++) k[i] = i;
  for (let i = half; i < N; i++) k[i] = i - N;
  const f = 2 * Math.PI / (N * dx);
  for (let i = 0; i < N; i++) k[i] *= f;
  return k;
}
