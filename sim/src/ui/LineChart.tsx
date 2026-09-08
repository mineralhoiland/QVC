import { useMemo } from "react";

export type Series = {
  xs: number[];
  ys: (number | null)[];
  color: string;
};

type Props = {
  series: Series[];
  width?: number;
  height?: number;
  xLabel?: string;
  yLabel?: string;
};

export function LineChart({ series, width = 268, height = 140, xLabel, yLabel }: Props) {
  const pad = { l: 34, r: 6, t: 10, b: 24 };
  const { paths, yTicks } = useMemo(() => {
    let xmin = Infinity,
      xmax = -Infinity,
      ymin = Infinity,
      ymax = -Infinity;
    for (const s of series) {
      for (let i = 0; i < s.xs.length; i++) {
        const y = s.ys[i];
        if (y == null || Number.isNaN(y)) continue;
        xmin = Math.min(xmin, s.xs[i]);
        xmax = Math.max(xmax, s.xs[i]);
        ymin = Math.min(ymin, y);
        ymax = Math.max(ymax, y);
      }
    }
    if (!Number.isFinite(xmin)) {
      xmin = 0;
      xmax = 1;
      ymin = 0;
      ymax = 1;
    }
    if (ymax === ymin) {
      ymax += 1;
      ymin -= 1;
    }
    const iw = width - pad.l - pad.r;
    const ih = height - pad.t - pad.b;
    const sx = (x: number) => pad.l + ((x - xmin) / (xmax - xmin || 1)) * iw;
    const sy = (y: number) => pad.t + (1 - (y - ymin) / (ymax - ymin || 1)) * ih;
    const paths = series.map((s) => {
      let d = "";
      let started = false;
      for (let i = 0; i < s.xs.length; i++) {
        const y = s.ys[i];
        if (y == null || Number.isNaN(y)) {
          started = false;
          continue;
        }
        d += `${started ? "L" : "M"}${sx(s.xs[i]).toFixed(1)},${sy(y).toFixed(1)} `;
        started = true;
      }
      return { d, color: s.color };
    });
    const yTicks = [ymin, 0.5 * (ymin + ymax), ymax];
    return { paths, yTicks };
  }, [series, width, height]);

  return (
    <svg className="chart" viewBox={`0 0 ${width} ${height}`} width="100%" height={height}>
      {yTicks.map((y, i) => (
        <text key={i} x={3} y={pad.t + ((2 - i) / 2) * (height - pad.t - pad.b) + 3} fill="#7dd3fc" fontSize="9">
          {y.toFixed(2)}
        </text>
      ))}
      {paths.map((p, i) => (
        <path key={i} d={p.d} fill="none" stroke={p.color} strokeWidth="1.8" />
      ))}
      {xLabel && (
        <text x={width / 2} y={height - 3} textAnchor="middle" fill="#94a3b8" fontSize="9">
          {xLabel}
        </text>
      )}
      {yLabel && (
        <text x={12} y={12} fill="#94a3b8" fontSize="9">
          {yLabel}
        </text>
      )}
    </svg>
  );
}
