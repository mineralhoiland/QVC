"""S8 — Panel four-viewport workbench.

    python -m qvc.workbench            # serves on http://localhost:5006
    python -m qvc.workbench --port 5010

Viewports
---------
1. Lattice   — magic-angle TEVC AA fragment (build_tevc_stack) or Im-3m hopfion BCC
2. Fields    — hopfion ρ_H isosurfaces + n glyphs, preimage linking, Casimir torus
3. Dynamics  — qvc.tqgl GPE↔DCE replay: live Δ(ψ), fold Δ(A), g⁽¹⁾, A/A_c, 3D torus frame
4. Spectra   — fold map, Spec(G) (democratic vs WS overlap), Mexican hat / U_eff(A), rates, Hall, Berry

Every viewport carries its L1/L2/L3 stamp. The optional L2 panel (acoustic
metric, torsion proxy, Goldstones) is off by default. L3 is a text overlay only.
All numbers come from :mod:`qvc.workbench.locks` and the existing solvers.
"""

from __future__ import annotations

import json
import time
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from qvc.workbench import analogue, barriers, catalyzon, dynamics, topology
from qvc.workbench.locks import DEPRECATED, Layer, get_locks, stamp

EXPORT_DIR = Path(__file__).resolve().parents[2] / "output" / "workbench_exports"

LAYER_COLOR = {"L1": "#1b7f3b", "L2": "#b8860b", "L3": "#8b0000"}


def _layer_badge(layer: str, text: str) -> str:
    c = LAYER_COLOR.get(layer, "#444")
    return (
        f'<div style="border-left:6px solid {c};padding:4px 8px;margin:2px 0 6px 0;background:#f7f7f7;'
        f'font-size:12px"><b style="color:{c}">{layer}</b> &nbsp; {text}</div>'
    )


# --------------------------------------------------------------------------
# Cached heavy computations
# --------------------------------------------------------------------------


@lru_cache(maxsize=4)
def _fields(n_grid: int):
    from qvc.viz.fields3d import compute_static_fields

    L = get_locks()
    return compute_static_fields(n_grid=n_grid, R_nm=L.R_nm, r_nm=L.r_nm, E_vac_meV=L.E_vac_meV, Q_H=L.Q_H)


@lru_cache(maxsize=4)
def _preimages(n_grid: int):
    f = _fields(n_grid)
    return topology.preimage_panel(f["n"], f["grid"], n_fibres=4)


@lru_cache(maxsize=6)
def _geometry(theta_deg: float, n_layers: int, frag_nm: float):
    from qvc.viz.structure import build_tevc_geometry

    return build_tevc_geometry(theta_deg, n_layers, fragment_radius_nm=frag_nm)


@lru_cache(maxsize=2)
def _replay(n_grid: int, t_end_ps: float, n_frames: int):
    res = dynamics.run_replay(n_grid=n_grid, t_end_ps=t_end_ps, n_frames=n_frames)
    frames = dynamics.replay_frames_3d(res, n_phi=90, n_chi=30)
    return res, frames


@lru_cache(maxsize=2)
def _slow_dce():
    return dynamics.run_slow_dce()


@lru_cache(maxsize=4)
def _berry(theta_deg: float, k_mesh: int):
    return topology.berry_panel(theta_deg=theta_deg, k_mesh=k_mesh, n_shells=3 if k_mesh <= 16 else 2)


@lru_cache(maxsize=2)
def _ws():
    return catalyzon.ws_overlap_comparison(n_grid=32)


# --------------------------------------------------------------------------
# Plotly figure builders
# --------------------------------------------------------------------------


def fig_lattice(theta_deg: float, n_layers: int, mode: str, frag_nm: float = 2.0, n_grid_fields: int = 40):
    import plotly.graph_objects as go

    L = get_locks()
    fig = go.Figure()
    if mode.startswith("TEVC"):
        g = _geometry(theta_deg, n_layers, frag_nm)
        fr = g.fragment
        pos = fr["cart_A"] / 10.0  # nm
        fig.add_trace(
            go.Scatter3d(
                x=pos[:, 0], y=pos[:, 1], z=pos[:, 2] * 3.0,  # z exaggerated ×3 for layer visibility
                mode="markers",
                marker=dict(size=2.2, color=fr["layer_index"], colorscale="Viridis", showscale=True, colorbar=dict(title="layer")),
                name="C",
            )
        )
        title = f"TEVC fragment: θ={theta_deg:.2f}°, {n_layers} layers, λ_M={g.lambda_M_A/10:.2f} nm, cell {g.n_atoms} C; fragment {fr['n_atoms']} atoms (AA site, r={frag_nm} nm, z×3)"
    else:
        a = 50.0  # nm = ξ display lattice
        pts = np.array([[i, j, k] for i in (0, 1) for j in (0, 1) for k in (0, 1)], float) * a
        fig.add_trace(go.Scatter3d(x=pts[:, 0], y=pts[:, 1], z=pts[:, 2], mode="markers", marker=dict(size=9, color="#1f77b4"), name="Hopf centre (2a corners)"))
        fig.add_trace(go.Scatter3d(x=[a / 2], y=[a / 2], z=[a / 2], mode="markers", marker=dict(size=9, color="#ff7f0e"), name="Hopf centre (body)"))
        # edges
        for i, p in enumerate(pts):
            for q in pts[i + 1 :]:
                if np.sum(np.abs(p - q) > 0) == 1:
                    fig.add_trace(go.Scatter3d(x=[p[0], q[0]], y=[p[1], q[1]], z=[p[2], q[2]], mode="lines", line=dict(color="#999", width=2), showlegend=False))
        pp = _preimages(n_grid_fields)
        cols = ["#d62728", "#2ca02c", "#9467bd", "#17becf"]
        for k, c in enumerate(pp["curves_nm"]):
            if len(c) < 3:
                continue
            cc = np.vstack([c, c[:1]]) * (0.3 * a / max(1e-9, np.max(np.abs(c)))) + a / 2
            fig.add_trace(go.Scatter3d(x=cc[:, 0], y=cc[:, 1], z=cc[:, 2], mode="lines", line=dict(color=cols[k % 4], width=6), name=f"fibre n⁻¹(p{k+1})"))
        title = f"Hopfion BCC (Im-3m, a=ξ=50 nm) with computed preimage fibres — pairwise linking {pp['mean_abs_linking']:.2f} [L2 schematic]"
    fig.update_layout(title=dict(text=title, font=dict(size=12)), margin=dict(l=0, r=0, t=40, b=0), scene=dict(aspectmode="data"), height=430, legend=dict(font=dict(size=10)))
    return fig


def fig_fields(mode: str, n_grid: int, iso_frac: float = 0.35, glyph_stride: int = 2):
    import plotly.graph_objects as go

    L = get_locks()
    f = _fields(n_grid)
    g = f["grid"]
    fig = go.Figure()
    if mode.startswith("Hopf"):
        rho = f["rho_H"]
        m = float(np.max(np.abs(rho)))
        fig.add_trace(
            go.Isosurface(
                x=g.X.ravel(), y=g.Y.ravel(), z=g.Z.ravel(), value=rho.ravel(),
                isomin=iso_frac * m, isomax=m, surface_count=2, opacity=0.45, colorscale="Reds", showscale=False, caps=dict(x_show=False, y_show=False, z_show=False), name="ρ_H>0",
            )
        )
        if float(rho.min()) < -1e-12:
            fig.add_trace(
                go.Isosurface(
                    x=g.X.ravel(), y=g.Y.ravel(), z=g.Z.ravel(), value=rho.ravel(),
                    isomin=float(rho.min()), isomax=-iso_frac * m, surface_count=1, opacity=0.35, colorscale="Blues", showscale=False, caps=dict(x_show=False, y_show=False, z_show=False), name="ρ_H<0",
                )
            )
        # n̂ glyphs only where the Hopf density is non-negligible (the core torus);
        # |n̂| = 1 everywhere so an unrestricted cone field just fills the box.
        s = slice(None, None, glyph_stride)
        n = f["n"]
        keep = np.abs(rho[s, s, s]) > 0.03 * m
        Xs, Ys, Zs = g.X[s, s, s][keep], g.Y[s, s, s][keep], g.Z[s, s, s][keep]
        ns = n[s, s, s][keep]
        # scale each cone by the local |ρ_H| so the glyphs read as "where the texture winds"
        amp = (np.abs(rho[s, s, s][keep]) / m) ** 0.5
        if Xs.size:
            fig.add_trace(
                go.Cone(
                    x=Xs, y=Ys, z=Zs,
                    u=ns[:, 0] * amp, v=ns[:, 1] * amp, w=ns[:, 2] * amp,
                    sizemode="scaled", sizeref=1.2, colorscale="Viridis", showscale=False, opacity=0.85, name="n̂ (core)",
                    hovertemplate="n̂ glyph<extra></extra>",
                )
            )
        title = (
            f"Hopf density ρ_H = A·B/16π² and n̂ glyphs — Whitehead Q = {f['Q_whitehead']:.3f} at {n_grid}³ "
            f"(Richardson → 1.00; old 1/4π² header would read {f['Q_if_4pi2_normalisation']:.2f}, rejected)"
        )
    elif mode.startswith("Preimage"):
        pp = _preimages(n_grid)
        cols = ["#d62728", "#2ca02c", "#9467bd", "#17becf"]
        for k, c in enumerate(pp["curves_nm"]):
            if len(c) < 3:
                continue
            cc = np.vstack([c, c[:1]])
            fig.add_trace(go.Scatter3d(x=cc[:, 0], y=cc[:, 1], z=cc[:, 2], mode="lines", line=dict(color=cols[k % 4], width=7), name=f"n⁻¹(p{k+1}) p={np.round(pp['points_S2'][k],1)}"))
        title = f"Preimage fibres of n̂ (Whitehead hopfion): pairwise Gauss linking = {pp['mean_abs_linking']:.2f} (expect 1). Q_H = 1/N = {L.Q_H:.2f} is a different (lens-space) number"
    else:
        # Casimir torus: surface coloured by the ZPF envelope value at the tube wall, plus phase θ
        R, r = L.R_nm, L.r_nm
        phi = np.linspace(0, 2 * np.pi, 120)
        chi = np.linspace(0, 2 * np.pi, 40)
        PHI, CHI = np.meshgrid(phi, chi, indexing="ij")
        X = (R + r * np.cos(CHI)) * np.cos(PHI)
        Y = (R + r * np.cos(CHI)) * np.sin(PHI)
        Z = r * np.sin(CHI)
        theta = (L.Q_H * PHI) % (2 * np.pi)
        fig.add_trace(go.Surface(x=X, y=Y, z=Z, surfacecolor=theta, colorscale="Phase", cmin=0, cmax=2 * np.pi, colorbar=dict(title="θ = Q_H φ"), name="torus"))
        # cross-section of the envelope in the z=0 plane
        gt = f["grid_torus"]
        mid = gt.n // 2
        fig.add_trace(go.Surface(x=gt.X[:, :, mid], y=gt.Y[:, :, mid], z=np.full_like(gt.X[:, :, mid], -2.5 * r), surfacecolor=f["casimir_envelope_meV"][:, :, mid], colorscale="Viridis", showscale=True, colorbar=dict(title="V_ZPF (meV)", x=1.15), opacity=0.9, name="−E_vac envelope (z=0)"))
        title = f"Fat torus (r,R)=({r:.0f},{R:.0f}) nm: phase θ = Q_H φ (branch cut = lens-space identification, 2D CS-flux {f['Q_H_flux_2d']:.3f}) and ZPF envelope −E_vac = −{L.E_vac_meV:.4f} meV"
    scene: dict[str, Any] = dict(aspectmode="data")
    if mode.startswith("Hopf"):
        # the box is ±4·core so the vacuum window fits the Whitehead integral; zoom the camera to the torus
        lim = 1.6 * f["core_nm"]
        scene = dict(
            aspectmode="cube",
            xaxis=dict(range=[-lim, lim], title="x (nm)"),
            yaxis=dict(range=[-lim, lim], title="y (nm)"),
            zaxis=dict(range=[-lim, lim], title="z (nm)"),
        )
    fig.update_layout(title=dict(text=title, font=dict(size=11)), margin=dict(l=0, r=0, t=40, b=0), scene=scene, height=430, legend=dict(font=dict(size=10)))
    return fig


def fig_dynamics_timeseries(res: dict[str, Any], slow: dict[str, Any] | None):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    L = get_locks()
    fig = make_subplots(rows=1, cols=3, subplot_titles=("live Δ(ψ) vs fold Δ(A) [meV]", "A(t)/A_c^fold and g⁽¹⁾", "80 ps DCE slaved to fold"))
    t = res["t_ps"]
    fig.add_trace(go.Scatter(x=t, y=res["Delta_live_meV"], name="Δ_live(ψ) — drive channel", line=dict(color="#d62728")), 1, 1)
    fig.add_trace(go.Scatter(x=t, y=res["Delta_fold_of_A_meV"], name="Δ_fold(A(t)) — static branch", line=dict(color="#1f77b4", dash="dash")), 1, 1)
    fig.add_hline(y=L.Delta_star_meV, line=dict(color="#1f77b4", dash="dot"), annotation_text="Δ* = 0.2324", row=1, col=1)
    fig.add_hline(y=L.tqgl.Delta_c_meV, line=dict(color="#7f7f7f", dash="dot"), annotation_text="Δ_c", row=1, col=1)
    fig.add_hline(y=L.Delta0_meV, line=dict(color="#2ca02c", dash="dot"), annotation_text="Δ₀", row=1, col=1)
    fig.add_trace(go.Scatter(x=t, y=res["A_over_Ac"], name="A/A_c^fold", line=dict(color="#9467bd")), 1, 2)
    fig.add_trace(go.Scatter(x=t, y=res["g1"], name="g⁽¹⁾", line=dict(color="#8c564b")), 1, 2)
    if slow is not None:
        for key, col in (("resonant", "#e377c2"), ("detuning", "#17becf")):
            d = slow[key]
            fig.add_trace(go.Scatter(x=d["t_ps"], y=d["A_um_inv"] / slow["A_c_fold_um_inv"], name=f"A/A_c ({key}) → {d['A_star_over_Ac']:.2f}", line=dict(color=col)), 1, 3)
        fig.add_hline(y=1.0, line=dict(color="#7f7f7f", dash="dot"), annotation_text="fold", row=1, col=3)
    fig.update_xaxes(title_text="t (ps)")
    fig.update_layout(height=300, margin=dict(l=10, r=10, t=40, b=30), legend=dict(orientation="h", font=dict(size=10), y=-0.25))
    return fig


def fig_torus_frame(frames: list[dict[str, Any]], k: int):
    import plotly.graph_objects as go

    fr = frames[min(k, len(frames) - 1)]
    fig = go.Figure(go.Surface(x=fr["X_um"] * 1e3, y=fr["Y_um"] * 1e3, z=fr["Z_um"] * 1e3, surfacecolor=fr["color"], colorscale="Inferno", colorbar=dict(title="|ψ|²")))
    fig.update_layout(title=dict(text=f"|ψ|² at t = {fr['t_ps']:.3f} ps wrapped on the torus — {dynamics.CAPTION_3D}", font=dict(size=11)), margin=dict(l=0, r=0, t=40, b=0), scene=dict(aspectmode="data", xaxis_title="nm", yaxis_title="nm", zaxis_title="nm"), height=330)
    return fig


def fig_density_2d(frames: list[dict[str, Any]], k: int, x_um: np.ndarray):
    import plotly.graph_objects as go

    fr = frames[min(k, len(frames) - 1)]
    fig = go.Figure(go.Heatmap(x=x_um * 1e3, y=x_um * 1e3, z=fr["density_2d"], colorscale="Inferno", colorbar=dict(title="|ψ|²")))
    fig.update_layout(title=dict(text=f"2D GPE |ψ|²(x,y) at t = {fr['t_ps']:.3f} ps (the actual evolved field)", font=dict(size=11)), margin=dict(l=10, r=10, t=40, b=30), height=330, xaxis_title="x (nm)", yaxis_title="y (nm)", yaxis=dict(scaleanchor="x"))
    return fig


def fig_fold(lam: float):
    import plotly.graph_objects as go

    L = get_locks()
    fm = barriers.fold_map()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=fm["alpha_over_ac"], y=fm["delta_phys"], name="physical branch", line=dict(color="#1f77b4", width=3)))
    fig.add_trace(go.Scatter(x=fm["alpha_over_ac"], y=fm["delta_unst"], name="unstable branch", line=dict(color="#d62728", dash="dash")))
    fig.add_trace(go.Scatter(x=[1.0], y=[fm["delta_c"]], mode="markers+text", text=["δ_c=α_c=0.1818"], textposition="bottom right", marker=dict(size=10, color="black"), name="fold"))
    op = fm["operating_point"]
    fig.add_trace(go.Scatter(x=[op["alpha_over_ac"]], y=[op["delta"]], mode="markers+text", text=[f"λ*={L.lambda_star:.3f}: Δ*={L.Delta_star_meV:.4f} meV"], textposition="top right", marker=dict(size=10, color="#2ca02c"), name="λ* operating point"))
    a_lam = lam * L.E_vac_meV / L.Delta0_meV / L.tqgl.alpha_c
    from qvc.gap_fold import solve_fold_branches

    br = solve_fold_branches(lam * L.E_vac_meV / L.Delta0_meV, L.tqgl.ell)
    if br["phys"] is not None:
        fig.add_trace(go.Scatter(x=[a_lam], y=[br["phys"]], mode="markers", marker=dict(size=12, symbol="diamond", color="#ff7f0e"), name=f"slider λ={lam:.3f}: δ={br['phys']:.3f} (Δ={br['phys']*L.Delta0_meV:.3f} meV)"))
    else:
        fig.add_vline(x=a_lam, line=dict(color="#ff7f0e", dash="dot"), annotation_text=f"λ={lam:.2f}: fold-out (no physical branch)")
    fig.update_layout(title=dict(text=fm["caption"], font=dict(size=11)), xaxis_title="α/α_c", yaxis_title="δ = Δ/Δ₀", height=330, margin=dict(l=10, r=10, t=40, b=30), legend=dict(font=dict(size=9), orientation="h", y=-0.3))
    return fig


def fig_specG(N: int, g0: float):
    import plotly.graph_objects as go

    L = get_locks()
    sw = catalyzon.spec_sweep(N_values=(N,))
    d = catalyzon.democratic_panel(N=N, g0_meV=g0)
    ws = _ws()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=sw["g0_meV"], y=sw["by_N"][N]["lambda_max"], name=f"λ_max=(N−1)g₀, N={N}", line=dict(color="#d62728")))
    fig.add_trace(go.Scatter(x=sw["g0_meV"], y=sw["by_N"][N]["lambda_min"], name="λ_min = −g₀ (×N−1, dark/catalyzon)", line=dict(color="#1f77b4")))
    fig.add_trace(go.Scatter(x=[g0] * N, y=d["eigenvalues_numeric"], mode="markers", marker=dict(size=10, color="black", symbol="x"), name=f"numeric Spec(G) at g₀={g0:.4f}"))
    fig.add_vline(x=L.g0_star_meV, line=dict(dash="dot", color="#2ca02c"), annotation_text=f"g₀*={L.g0_star_meV:.4f}")
    # WS overlap (dimensionless λ/g0_eff scaled to current g0) — compared, not forced
    for row in ws["rows"]:
        if abs(row["xi_over_lambda_M"] - 3.85) < 0.01 or row["xi_over_lambda_M"] == 1.0:
            fig.add_trace(go.Scatter(x=[g0 * 1.02] * len(row["eigenvalues_over_g0eff"]), y=np.asarray(row["eigenvalues_over_g0eff"]) * g0, mode="markers", marker=dict(size=8, symbol="diamond-open", color="#ff7f0e" if row["xi_over_lambda_M"] == 1.0 else "#9467bd"), name=f"WS overlap ξ_K={row['xi_nm']:.0f} nm (×g₀, rel‖·‖_F to J−I = {row['rel_frobenius_to_democratic']:.2f})"))
    fig.update_layout(title=dict(text=stamp(Layer.L1, f"Spec(G): democratic lock vs WS overlap G^(H) (compared, not forced). Δ_cat = Δ_ex − g₀ = {d['Delta_cat_meV']:.4f} meV ({d['reduction_pct']:.1f}%, parallel channel)"), font=dict(size=11)), xaxis_title="g₀ (meV)", yaxis_title="eigenvalue (meV)", height=330, margin=dict(l=10, r=10, t=40, b=30), legend=dict(font=dict(size=9), orientation="h", y=-0.3))
    return fig


def fig_mexican(A_over_Ac: float):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    L = get_locks()
    mh = barriers.mexican_hat_tracks()
    A = A_over_Ac * L.tqgl.A_c_fold_um_inv
    u = barriers.u_eff_of_A(A)
    fig = make_subplots(rows=1, cols=2, specs=[[{"type": "xy"}, {"type": "surface"}]], subplot_titles=("V(ρ) tracks and U_eff(ρ; A) at fixed energy", "3D Mexican hat (catalyzon Δ_cat)"))
    for key, col in (("bare_Delta0", "#2ca02c"), ("catalyzon_Delta_cat", "#d62728"), ("static_fold_Delta_star", "#1f77b4")):
        tr = mh["tracks"][key]
        fig.add_trace(go.Scatter(x=tr["rho"], y=tr["V_meV"], name=f"{key.replace('_',' ')} (Δ={tr['Delta_meV']:.3f}, B={tr['well_depth_meV']:.3f})", line=dict(color=col)), 1, 1)
    if not u["folded_out"]:
        fig.add_trace(go.Scatter(x=u["rho"], y=u["V_meV_fixed_energy"], name=f"U_eff(A/A_c={A_over_Ac:.2f}): Δ(A)={u['Delta_track_meV']:.3f}, offset pinned {u['energy_offset_meV']:+.3f} meV", line=dict(color="#ff7f0e", dash="dash", width=3)), 1, 1)
    surf = barriers.mexican_hat_surface(L.Delta_cat_meV)
    fig.add_trace(go.Surface(x=surf["X"], y=surf["Y"], z=surf["Z"], colorscale="Viridis", showscale=False), 1, 2)
    fig.update_xaxes(title_text="ρ = |Ψ|", row=1, col=1)
    fig.update_yaxes(title_text="V (meV)", row=1, col=1)
    fig.update_layout(title=dict(text=stamp(Layer.L1, f"B_cat/B₀ = (Δ_cat/Δ₀)² = {mh['B_cat_over_B0']:.3f}; fold channel B*/B₀ = {mh['B_fold_over_B0']:.3f} — parallel, never stacked"), font=dict(size=11)), height=330, margin=dict(l=10, r=10, t=40, b=30), legend=dict(font=dict(size=9), orientation="h", y=-0.3))
    return fig


def fig_rates():
    import plotly.graph_objects as go

    rk = barriers.rate_catalysis_sketch()
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=rk["T_K"], y=rk["arrhenius_ratio_catalyzon"], name="Arrhenius exp(g₀*/k_BT) (catalyzon)", line=dict(color="#d62728")))
    fig.add_trace(go.Scatter(x=rk["T_K"], y=rk["arrhenius_ratio_fold_channel"], name="Arrhenius exp((Δ₀−Δ*)/k_BT) (fold channel)", line=dict(color="#1f77b4")))
    w = rk["wkb"]
    fig.add_hline(y=w["ratio_catalyzon"], line=dict(dash="dot", color="#d62728"), annotation_text=f"WKB ratio (cat) {w['ratio_catalyzon']:.2f}")
    fig.add_hline(y=w["ratio_fold"], line=dict(dash="dot", color="#1f77b4"), annotation_text=f"WKB ratio (fold) {w['ratio_fold']:.2f}")
    fig.update_yaxes(type="log", title_text="Γ_QVC / Γ_bare")
    fig.update_xaxes(title_text="T (K)")
    fig.update_layout(title=dict(text=rk["caption"], font=dict(size=11)), height=330, margin=dict(l=10, r=10, t=40, b=30), legend=dict(font=dict(size=9), orientation="h", y=-0.3))
    return fig


def fig_hall(N: int):
    import plotly.graph_objects as go

    h = topology.hall_panel(N_values=tuple(sorted({3, 4, 5, 6, N})))
    fig = go.Figure()
    for n_, v in h["by_N"].items():
        fig.add_trace(go.Scatter(x=h["B_over_B0"], y=v["sigma_xy_e2h"], name=f"N={n_}: step e²/{2*n_}h", line=dict(width=3 if n_ == N else 1.2)))
    fig.update_layout(title=dict(text=h["caption"], font=dict(size=11)), xaxis_title="B / B₀", yaxis_title="σ_xy (e²/h)", height=330, margin=dict(l=10, r=10, t=40, b=30), legend=dict(font=dict(size=9), orientation="h", y=-0.3))
    return fig


def fig_berry(theta_deg: float, k_mesh: int):
    import plotly.graph_objects as go

    b = _berry(theta_deg, k_mesh)
    fig = go.Figure(go.Heatmap(z=b["omega"], colorscale="RdBu", colorbar=dict(title="Ω(k)")))
    fig.update_layout(title=dict(text=b["caption"], font=dict(size=11)), xaxis_title="k₂ index (moiré BZ)", yaxis_title="k₁ index", height=330, margin=dict(l=10, r=10, t=40, b=30))
    return fig


def fig_acoustic(drain: float):
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

    am = analogue.acoustic_metric_slice(drain_v0_over_cs=drain)
    x = am["x_um"] * 1e3
    fig = make_subplots(rows=1, cols=2, subplot_titles=("Mach |v|/c_s with sonic surface", "g_tt = −(c_s² − v²)"))
    fig.add_trace(go.Heatmap(x=x, y=x, z=np.clip(am["mach"], 0, 3), colorscale="Viridis", colorbar=dict(title="Mach", x=0.45)), 1, 1)
    fig.add_trace(go.Contour(x=x, y=x, z=am["mach"], contours=dict(start=1, end=1, size=1, coloring="none"), line=dict(color="white", width=2), showscale=False, name="sonic"), 1, 1)
    fig.add_trace(go.Heatmap(x=x, y=x, z=np.clip(am["g_tt"], -3 * am["c_s_um_ps"] ** 2, 3 * am["c_s_um_ps"] ** 2), colorscale="RdBu", colorbar=dict(title="g_tt")), 1, 2)
    fig.update_layout(title=dict(text=am["caption"], font=dict(size=11)), height=330, margin=dict(l=10, r=10, t=40, b=30))
    return fig


def fig_goldstone():
    import plotly.graph_objects as go

    g = analogue.goldstone_cartoon()
    ein = analogue.h_ein1_overlap()
    fig = go.Figure()
    for name, w in g["branches"].items():
        fig.add_trace(go.Scatter(x=g["k_nm_inv"], y=w, name=name))
    fig.update_layout(title=dict(text=g["caption"] + f"<br>{ein['caption']}", font=dict(size=10)), xaxis_title="k (nm⁻¹)", yaxis_title="ħω (meV)", height=330, margin=dict(l=10, r=10, t=60, b=30), legend=dict(font=dict(size=9)))
    return fig


# --------------------------------------------------------------------------
# Exports
# --------------------------------------------------------------------------


def do_export(kind: str, *, theta_deg: float, n_layers: int, n_grid_fields: int, replay=None) -> str:
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if kind == "CIF/POSCAR (VESTA)":
        from qvc.viz.structure import export_tevc, hopfion_bcc_cif, preimage_overlay_cif

        g = _geometry(theta_deg, n_layers, 2.0)
        r = export_tevc(EXPORT_DIR, g)
        hopfion_bcc_cif(EXPORT_DIR / "hopfion_crystal_BCC_Im3m.cif")
        pp = _preimages(n_grid_fields)
        preimage_overlay_cif(EXPORT_DIR / "hopfion_preimage_overlay.cif", pp["curves_nm"][:2])
        return f"wrote {len(r['files'])+2} files to {EXPORT_DIR} ({time.time()-t0:.1f}s): TEVC.cif ({r['n_atoms']} C), POSCAR, tevc.data, TEVC_fragment.cif, hopfion_crystal_BCC_Im3m.cif, hopfion_preimage_overlay.cif"
    if kind == "CUBE + VTK (VESTA/ParaView)":
        from qvc.viz.fields3d import export_static_fields

        f = _fields(n_grid_fields)
        files = export_static_fields(f, EXPORT_DIR)
        return f"wrote {len(files)} volume files to {EXPORT_DIR} ({time.time()-t0:.1f}s); Whitehead Q={f['Q_whitehead']:.3f} at {n_grid_fields}³"
    if kind == "GPE time series (.pvd)":
        res, _ = replay
        ex = dynamics.export_replay_paraview(res, EXPORT_DIR)
        return f"wrote {ex['n_frames']} frames → {ex['pvd_torus']} and {ex['pvd_2d']} ({time.time()-t0:.1f}s)"
    if kind == "Locks JSON":
        L = get_locks()
        p = EXPORT_DIR / "workbench_locks.json"
        p.write_text(json.dumps({"locks": L.as_dict(), "deprecated_reject_list": DEPRECATED}, indent=2, default=str))
        return f"wrote {p}"
    return "unknown export"


# --------------------------------------------------------------------------
# App
# --------------------------------------------------------------------------


def build_app(*, fields_grid: int = 40, gpe_grid: int = 64, gpe_t_ps: float = 0.5, gpe_frames: int = 12):
    import panel as pn

    pn.extension("plotly", sizing_mode="stretch_width")
    L = get_locks()
    t = L.tqgl

    # ---------------- controls ----------------
    theta = pn.widgets.FloatSlider(name="θ (°) — lattice / Berry", start=0.9, end=1.4, step=0.05, value=1.1)
    N = pn.widgets.IntSlider(name="N layers — Spec(G), Hall, Q_H=1/N", start=3, end=7, value=L.N)
    lam = pn.widgets.FloatSlider(name="λ participation (fold α = λE_vac/Δ₀)", start=0.0, end=1.2, step=0.01, value=round(L.lambda_star, 4))
    g0 = pn.widgets.FloatSlider(name="g₀ (meV) — Spec(G)", start=0.0, end=0.3, step=0.002, value=round(L.g0_star_meV, 4))
    A_frac = pn.widgets.FloatSlider(name="A_ZPF / A_c^fold — U_eff reshaping", start=0.0, end=0.99, step=0.01, value=round(t.A_star_um_inv / t.A_c_fold_um_inv, 2))
    frame = pn.widgets.IntSlider(name="GPE frame", start=0, end=gpe_frames, value=gpe_frames)
    lattice_mode = pn.widgets.RadioButtonGroup(name="Lattice", options=["TEVC AA fragment", "Hopfion BCC + fibres"], value="TEVC AA fragment")
    field_mode = pn.widgets.RadioButtonGroup(name="Fields", options=["Hopf density + n̂", "Preimage linking", "Casimir torus"], value="Hopf density + n̂")
    spectra_tab = pn.widgets.RadioButtonGroup(name="Spectra", options=["Fold", "Spec(G)", "Mexican hat", "Rates", "Hall", "Berry"], value="Fold")
    show_l2 = pn.widgets.Checkbox(name="Show L2 analogue-geometry panel (tagged)", value=False)
    drain = pn.widgets.FloatSlider(name="L2: drain v₀/c_s (0 = pure fractional vortex)", start=0.0, end=2.0, step=0.1, value=0.0)
    h_tor1 = pn.widgets.Checkbox(name="L2: H-TOR1 toggle (parked — annotates only)", value=False)
    show_l3 = pn.widgets.Checkbox(name="Show L3 schematic overlay note", value=False)
    export_kind = pn.widgets.Select(name="Export", options=["CIF/POSCAR (VESTA)", "CUBE + VTK (VESTA/ParaView)", "GPE time series (.pvd)", "Locks JSON"])
    export_btn = pn.widgets.Button(name="Export", button_type="primary")
    export_status = pn.pane.Markdown("", styles={"font-size": "11px"})

    # ---------------- lock table ----------------
    rows = "".join(f"<tr><td><b>{s}</b></td><td>{v}</td><td style='color:{LAYER_COLOR[l]}'>{l}</td></tr>" for s, v, l in L.table())
    lock_html = pn.pane.HTML(
        f"<details open><summary><b>S0 frozen locks (Option A*)</b> — deprecated: "
        + ", ".join(f"{d['value']} {d['unit']}" for d in DEPRECATED.values())
        + f"</summary><table style='font-size:11px'>{rows}</table></details>"
    )

    # ---------------- viewports ----------------
    replay = _replay(gpe_grid, gpe_t_ps, gpe_frames)
    res, frames = replay
    frame.end = len(frames) - 1
    frame.value = len(frames) - 1
    slow = _slow_dce()

    lattice_pane = pn.pane.Plotly(pn.bind(lambda th, nl, m: fig_lattice(th, nl, m, n_grid_fields=fields_grid), theta, N, lattice_mode), height=440)
    fields_pane = pn.pane.Plotly(pn.bind(lambda m: fig_fields(m, fields_grid), field_mode), height=440)
    dyn_ts = pn.pane.Plotly(fig_dynamics_timeseries(res, slow), height=310)
    dyn_3d = pn.pane.Plotly(pn.bind(lambda k: fig_torus_frame(frames, k), frame), height=340)
    dyn_2d = pn.pane.Plotly(pn.bind(lambda k: fig_density_2d(frames, k, res["x_um"]), frame), height=340)

    def spectra_fig(tab, lam_v, N_v, g0_v, A_v, th):
        if tab == "Fold":
            return fig_fold(lam_v)
        if tab == "Spec(G)":
            return fig_specG(int(N_v), float(g0_v))
        if tab == "Mexican hat":
            return fig_mexican(float(A_v))
        if tab == "Rates":
            return fig_rates()
        if tab == "Hall":
            return fig_hall(int(N_v))
        return fig_berry(float(th), 12)

    spectra_pane = pn.pane.Plotly(pn.bind(spectra_fig, spectra_tab, lam, N, g0, A_frac, theta), height=340)

    l2_panel = pn.Column(
        pn.pane.HTML(_layer_badge("L2", Layer.L2.title + " — acoustic metric, torsion proxy, Goldstones. Nothing here feeds L1 locks.")),
        pn.pane.Plotly(pn.bind(fig_acoustic, drain), height=340),
        pn.pane.Plotly(fig_goldstone(), height=340),
        pn.pane.Markdown(pn.bind(lambda tog: "**" + analogue.torsion_proxy(_fields(fields_grid)["rho_H"], h_tor1_toggle=tog)["note"] + "**", h_tor1)),
        visible=False,
    )

    def _toggle_l2(event):
        l2_panel.visible = bool(event.new)

    show_l2.param.watch(_toggle_l2, "value")
    l3_note = pn.pane.HTML(pn.bind(lambda s: _layer_badge("L3", analogue.L3_SCHEMATIC_NOTE) if s else "", show_l3))

    def _export(event):
        export_status.object = "exporting…"
        try:
            export_status.object = do_export(export_kind.value, theta_deg=float(theta.value), n_layers=int(N.value), n_grid_fields=fields_grid, replay=replay)
        except Exception as exc:  # pragma: no cover
            export_status.object = f"export failed: {exc}"

    export_btn.on_click(_export)

    sidebar = pn.Column(
        pn.pane.Markdown("### Controls"),
        theta, N, lam, g0, A_frac, frame,
        pn.pane.Markdown("### Layers"),
        show_l2, drain, h_tor1, show_l3,
        pn.pane.Markdown("### Export"),
        export_kind, export_btn, export_status,
        pn.pane.Markdown(f"<small>GPE replay: qvc.tqgl {gpe_grid}² × {gpe_t_ps} ps; fields {fields_grid}³ (Whitehead Q at this grid = {_fields(fields_grid)['Q_whitehead']:.3f}, Richardson → 1.00). Units: Å atomic, nm mesoscale, µm/ps/meV GPE.</small>"),
        width=330,
    )

    vp1 = pn.Column(pn.pane.HTML(_layer_badge("L1", "Lattice — magic-angle TEVC from build_tevc_stack (13–22° CIFs are proxies); BCC hopfion view is L2 schematic")), lattice_mode, lattice_pane)
    vp2 = pn.Column(pn.pane.HTML(_layer_badge("L2", "Fields — Whitehead hopfion n̂, ρ_H, FN energy (L2 geometry); Casimir envelope uses the L1 lock E_vac; Q_H=1/N is the L1 lens-space number")), field_mode, fields_pane)
    vp3 = pn.Column(pn.pane.HTML(_layer_badge("L1", "Dynamics — " + dynamics.CAPTION_CHANNELS)), dyn_ts, pn.Row(dyn_2d, dyn_3d))
    vp4 = pn.Column(pn.pane.HTML(_layer_badge("L1", "Spectra / barriers — fold, Spec(G) democratic vs WS overlap, Mexican hat at fixed energy, rate sketch, Hall staircase, BM Berry curvature")), spectra_tab, spectra_pane)

    main = pn.Column(
        lock_html,
        l3_note,
        pn.Row(vp1, vp2),
        pn.Row(vp3, vp4),
        l2_panel,
    )
    template = pn.template.FastListTemplate(
        title="TEVC / QVC–TDVT crystal & materials workbench — QVCCompute (VESTA/ParaView export)",
        sidebar=[sidebar],
        main=[main],
        accent_base_color="#1b7f3b",
        header_background="#1b7f3b",
    )
    return template


def serve(port: int = 5006, show: bool = False, **kw):
    import panel as pn

    app = build_app(**kw)
    pn.serve(app, port=port, show=show, title="TEVC QVC-TDVT workbench")
