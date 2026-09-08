import { useMemo } from "react";
import { democraticSpectrum } from "../math/coupling";
import { foldCritical, foldCurve } from "../math/fold";
import { MODES, useSimStore } from "../store";
import { LineChart } from "./LineChart";

export function Analysis() {
  const p = useSimStore((s) => s.params);
  const fps = useSimStore((s) => s.fps);
  const meta = MODES[p.mode];

  const spec = useMemo(() => democraticSpectrum(p.N, p.g0), [p.N, p.g0]);
  const fold = useMemo(() => foldCritical(-2.791), []);
  const foldSeries = useMemo(() => {
    const c = foldCurve(-2.791, 60);
    return [
      { xs: c.alphas, ys: c.phys, color: "#22d3ee" },
      { xs: c.alphas, ys: c.unst, color: "#f472b6" },
    ];
  }, []);

  const hopfSeries = useMemo(() => {
    const xs = Array.from({ length: 40 }, (_, i) => 1 + i * 0.05);
    return [
      {
        xs,
        ys: xs.map((q) => Math.pow(q, 0.75)),
        color: "#fbbf24",
      },
      {
        xs,
        ys: xs.map((q) => 0.85 * Math.pow(q, 0.75)),
        color: "#34d399",
      },
    ];
  }, []);

  const horizonTH = (1 / (2 * Math.PI)) * (p.flow * 1.2);
  const rH = 1 / Math.max(p.flow, 0.15);

  return (
    <>
      <div className="readout">
        <div>
          <span className="k">Mode </span>
          {meta.name}
        </div>
        <div>
          <span className="k">Particles </span>
          {p.count.toLocaleString()}
        </div>
        <div>
          <span className="k">FPS </span>
          {fps.toFixed(0)}
        </div>
        <div>
          <span className="k">GPU path </span>
          vertex shader (no CPU integrate)
        </div>
      </div>

      {p.mode <= 2 && (
        <>
          <div className="readout">
            <div>
              <span className="k">Hopf charge Q </span>
              {p.Q.toFixed(0)}
            </div>
            <div>
              <span className="k">VK proxy E ≥ C|Q|³⁄⁴ </span>
              {Math.pow(Math.abs(p.Q), 0.75).toFixed(3)}
            </div>
            <div>
              <span className="k">Link strength </span>
              {p.link.toFixed(2)}
            </div>
            <div>
              <span className="k">Map </span>
              Hopf S³ → ℝ³ (stereo)
            </div>
          </div>
          <LineChart series={hopfSeries} xLabel="|Q|" yLabel="E_min" />
        </>
      )}

      {p.mode === 3 && (
        <div className="readout">
          <div>
            <span className="k">Vortex count ~ Q </span>
            {Math.max(1, Math.floor(p.Q))}
          </div>
          <div>
            <span className="k">Circulation twist </span>
            {p.twist.toFixed(2)}
          </div>
          <div>
            <span className="k">Lattice a </span>
            {p.lattice.toFixed(2)}
          </div>
          <p className="hint">Azimuthal particle drift ∼ 1/r around cores; ripples mimic Kelvin waves.</p>
        </div>
      )}

      {p.mode === 4 && (
        <>
          <div className="readout">
            <div>
              <span className="k">λ_max = (N−1)g₀ </span>
              {spec.lambdaMax.toFixed(4)}
            </div>
            <div>
              <span className="k">λ_min = −g₀ (×{spec.multiplicityMin}) </span>
              {spec.lambdaMin.toFixed(4)}
            </div>
            <div>
              <span className="k">α / α_c (fold) </span>
              {(p.alpha / fold.alphaC).toFixed(3)}
            </div>
            <div>
              <span className="k">α_c (ℓ=−2.791) </span>
              {fold.alphaC.toFixed(4)}
            </div>
          </div>
          <LineChart series={foldSeries} xLabel="α/α_c" yLabel="δ" />
          <p className="hint">Amber = collective bright; teal = dark multiplet. Raising α melts dark modes (fold).</p>
        </>
      )}

      {p.mode === 5 && (
        <div className="readout">
          <div>
            <span className="k">v/c proxy (flow) </span>
            {p.flow.toFixed(3)}
          </div>
          <div>
            <span className="k">Horizon radius r_H </span>
            {rH.toFixed(3)}
          </div>
          <div>
            <span className="k">T_H ∼ κ/2π </span>
            {horizonTH.toFixed(3)}
          </div>
          <p className="hint">Warm particles mark the acoustic horizon shell around the hopfion core.</p>
        </div>
      )}

      <p className="hint">{meta.blurb}</p>
    </>
  );
}
