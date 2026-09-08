import { create } from "zustand";

export type ModeId = 0 | 1 | 2 | 3 | 4 | 5;

export type ModeMeta = {
  id: ModeId;
  name: string;
  short: string;
  badge: "Locked" | "Provisional" | "Ideation";
  blurb: string;
};

export const MODES: ModeMeta[] = [
  {
    id: 0,
    name: "Hopf Fibration Evolution",
    short: "Hopf S³",
    badge: "Locked",
    blurb: "Stereographic Hopf map S³→ℝ³. Particles ride fibers; color encodes base S².",
  },
  {
    id: 1,
    name: "Linked Hopfion Dynamics",
    short: "Hopfions",
    badge: "Locked",
    blurb: "Two particle families on linked toroidal tubes — Hopf-link topology in motion.",
  },
  {
    id: 2,
    name: "Lattice Spacetime",
    short: "TDVT lattice",
    badge: "Locked",
    blurb: "BCC lattice of hopfion cores with interlink filaments — discrete TDVT geometry.",
  },
  {
    id: 3,
    name: "Vortex Mechanics",
    short: "Vortices",
    badge: "Provisional",
    blurb: "Particle lattice + multi-vortex cores, azimuthal flow and Kelvin-like ripples.",
  },
  {
    id: 4,
    name: "Catalyzon Spec(G)",
    short: "Catalyzon",
    badge: "Locked",
    blurb: "Democratic G=g₀(J−I): λ_min=−g₀ (×N−1), λ_max=(N−1)g₀. Fold melt via α.",
  },
  {
    id: 5,
    name: "Acoustic Horizon Flow",
    short: "Horizon",
    badge: "Provisional",
    blurb: "Particles stream past a hopfion core; glow marks the |v|∼cₛ acoustic horizon.",
  },
];

export type Params = {
  mode: ModeId;
  count: number;
  speed: number;
  Q: number;
  scale: number;
  twist: number;
  lattice: number;
  link: number;
  flow: number;
  alpha: number;
  N: number;
  g0: number;
  pointSize: number;
  glow: number;
  opacity: number;
  running: boolean;
};

const defaults: Params = {
  mode: 0,
  count: 524_288,
  speed: 1,
  Q: 1,
  scale: 1.4,
  twist: 1,
  lattice: 1.6,
  link: 0.85,
  flow: 1.25,
  alpha: 0.45,
  N: 5,
  g0: 0.21,
  pointSize: 0.85,
  glow: 0.4,
  opacity: 0.55,
  running: true,
};

type Store = {
  params: Params;
  fps: number;
  setMode: (m: ModeId) => void;
  setParam: <K extends keyof Params>(k: K, v: Params[K]) => void;
  setFps: (f: number) => void;
  reset: () => void;
};

export const useSimStore = create<Store>((set) => ({
  params: { ...defaults },
  fps: 0,
  setMode: (m) => set((s) => ({ params: { ...s.params, mode: m } })),
  setParam: (k, v) => set((s) => ({ params: { ...s.params, [k]: v } })),
  setFps: (f) => set({ fps: f }),
  reset: () => set({ params: { ...defaults } }),
}));
