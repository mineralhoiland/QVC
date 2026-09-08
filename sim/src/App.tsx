import { MODES, useSimStore, type ModeId } from "./store";
import { Analysis } from "./ui/Analysis";
import { CanvasView } from "./ui/CanvasView";
import { Slider } from "./ui/Slider";

export default function App() {
  const params = useSimStore((s) => s.params);
  const setMode = useSimStore((s) => s.setMode);
  const setParam = useSimStore((s) => s.setParam);
  const reset = useSimStore((s) => s.reset);
  const fps = useSimStore((s) => s.fps);
  const meta = MODES[params.mode];

  return (
    <div className="app">
      <header className="topbar">
        <div className="brand">
          <span className="brand-mark">TDVT</span>
          <span className="brand-sub">GPU Particle Geometry Lab</span>
        </div>
        <nav className="suite-nav">
          {MODES.map((m) => (
            <button
              key={m.id}
              type="button"
              className={`suite-btn${params.mode === m.id ? " active" : ""}`}
              onClick={() => setMode(m.id as ModeId)}
              title={m.name}
            >
              {m.short}
            </button>
          ))}
        </nav>
        <span className={`badge ${meta.badge}`}>{meta.badge}</span>
        <span className="fps-pill">
          {params.count >= 1e6 ? `${(params.count / 1e6).toFixed(2)}M` : `${Math.round(params.count / 1000)}K`}
          {" · "}
          {fps.toFixed(0)} fps
        </span>
      </header>

      <div className="title-strip">
        <h1>{meta.name}</h1>
        <p>{meta.blurb}</p>
      </div>

      <div className="workspace">
        <aside className="panel">
          <h3>Interactive controls</h3>

          <Slider
            label="Particle count"
            value={params.count}
            min={65536}
            max={1048576}
            step={16384}
            format={(v) => v.toLocaleString()}
            onChange={(v) => setParam("count", v)}
          />
          <Slider label="Speed" value={params.speed} min={0} max={3} onChange={(v) => setParam("speed", v)} />
          <Slider
            label="Hopf charge Q"
            value={params.Q}
            min={1}
            max={5}
            step={1}
            format={(v) => String(v)}
            onChange={(v) => setParam("Q", v)}
          />
          <Slider label="Scale" value={params.scale} min={0.4} max={3} onChange={(v) => setParam("scale", v)} />
          <Slider label="Twist / winding" value={params.twist} min={0} max={3} onChange={(v) => setParam("twist", v)} />
          <Slider label="Link strength" value={params.link} min={0} max={1.5} onChange={(v) => setParam("link", v)} />
          <Slider label="Lattice a" value={params.lattice} min={0.6} max={3} onChange={(v) => setParam("lattice", v)} />

          {params.mode === 4 && (
            <>
              <Slider label="α (fold / ZPF)" value={params.alpha} min={0} max={1} onChange={(v) => setParam("alpha", v)} />
              <Slider
                label="N modes"
                value={params.N}
                min={2}
                max={8}
                step={1}
                format={(v) => String(v)}
                onChange={(v) => setParam("N", v)}
              />
              <Slider label="g₀" value={params.g0} min={0.05} max={0.6} onChange={(v) => setParam("g0", v)} />
            </>
          )}

          {params.mode === 5 && (
            <Slider label="Flow v/c" value={params.flow} min={0.3} max={2.5} onChange={(v) => setParam("flow", v)} />
          )}

          <h3>Look</h3>
          <Slider label="Point size" value={params.pointSize} min={0.4} max={3} onChange={(v) => setParam("pointSize", v)} />
          <Slider label="Glow" value={params.glow} min={0} max={1.5} onChange={(v) => setParam("glow", v)} />
          <Slider label="Opacity" value={params.opacity} min={0.2} max={1} onChange={(v) => setParam("opacity", v)} />

          <div className="btn-row">
            <button
              type="button"
              className={`btn${params.running ? " primary" : ""}`}
              onClick={() => setParam("running", !params.running)}
            >
              {params.running ? "Pause" : "Run"}
            </button>
            <button type="button" className="btn" onClick={reset}>
              Reset
            </button>
            <button type="button" className="btn" onClick={() => setParam("count", 524288)}>
              512K
            </button>
            <button type="button" className="btn" onClick={() => setParam("count", 1048576)}>
              1M
            </button>
          </div>
        </aside>

        <CanvasView />

        <aside className="panel right">
          <h3>Mathematical analysis</h3>
          <Analysis />
        </aside>
      </div>

      <footer className="footer">
        All motion in vertex shaders · additive particles · locked Spec(G) & Hopf stereo map · orbit drag to explore
      </footer>
    </div>
  );
}
