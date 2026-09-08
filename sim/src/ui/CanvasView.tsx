import { useEffect, useRef } from "react";
import { SceneHost } from "../engine/SceneHost";
import { useSimStore } from "../store";

export function CanvasView() {
  const wrapRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const hostRef = useRef<SceneHost | null>(null);
  const params = useSimStore((s) => s.params);
  const setFps = useSimStore((s) => s.setFps);

  // Mount once
  useEffect(() => {
    const wrap = wrapRef.current;
    const canvas = canvasRef.current;
    if (!wrap || !canvas) return;
    const host = new SceneHost(canvas, useSimStore.getState().params.count);
    hostRef.current = host;
    host.setFpsCallback((fps) => setFps(fps));
    const dispose = host.mount(wrap);
    return () => {
      dispose();
      hostRef.current = null;
    };
  }, [setFps]);

  // Sync uniforms / count / running
  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;
    host.setRunning(params.running);
    host.setCount(params.count);
    // Density compensation so 512K/1M stay structured, not whiteout
    const dens = Math.sqrt(65536 / Math.max(params.count, 65536));
    host.setUniforms({
      mode: params.mode,
      speed: params.speed,
      Q: params.Q,
      scale: params.scale,
      twist: params.twist,
      lattice: params.lattice,
      link: params.link,
      flow: params.flow,
      alpha: params.alpha,
      N: params.N,
      g0: params.g0,
      pointSize: params.pointSize * (0.65 + 0.35 * dens),
      glow: params.glow,
      opacity: params.opacity * dens,
    });
  }, [params]);

  return (
    <div className="viewport" ref={wrapRef}>
      <canvas ref={canvasRef} className="gpu-canvas" />
      <div className="viewport-badge">GPU particles · additive · shader integrate</div>
    </div>
  );
}
