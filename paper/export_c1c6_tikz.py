#!/usr/bin/env python3
"""Export C1–C6 summary.json into pgfplots .dat files for QVC_arXiv_v4."""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

ROOT = Path(__file__).resolve().parent
SUMMARY = Path(
    "/Users/mineralhoiland/Code/QVCCursor/notebooks/bkt_suite/output/summary.json"
)
DATA = ROOT / "data"
PI = math.pi


def integrate_kt(K_init: float, y_init: float, ell_ir: float = 8.0, n_eval: int = 161):
    u0 = 1.0 / max(K_init, 1e-12)
    y0 = max(float(y_init), 1e-18)

    def rhs(_ell, z):
        u, y = z
        u = max(float(u), 1e-12)
        K = 1.0 / u
        return [4.0 * PI**3 * y * y, (2.0 - PI * K) * y]

    def plasma(_ell, z):
        return 0.4 - z[1]

    plasma.terminal = True
    plasma.direction = -1.0

    sol = solve_ivp(
        rhs,
        (0.0, ell_ir),
        [u0, y0],
        rtol=1e-7,
        atol=1e-10,
        dense_output=True,
        events=plasma,
        max_step=0.05,
    )
    t_end = min(float(sol.t[-1]), ell_ir)
    ell = np.linspace(0.0, t_end, n_eval)
    uv = sol.sol(ell)
    K = 1.0 / np.maximum(uv[0], 1e-12)
    y = np.maximum(uv[1], 0.0)
    return ell, K, y


def trim_flow(ell, K, y, ycut: float = 1e-10):
    y = np.maximum(np.asarray(y, dtype=float), 1e-16)
    ell = np.asarray(ell, dtype=float)
    K = np.asarray(K, dtype=float)
    below = np.where(y < ycut)[0]
    if below.size > 2:
        cut = int(below[2]) + 1
        return ell[:cut], K[:cut], y[:cut]
    return ell, K, y


def dump(path: Path, header: str, rows: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(header if header.endswith("\n") else header + "\n")
        np.savetxt(f, rows)


def main() -> None:
    bundle = json.loads(SUMMARY.read_text())
    c1 = bundle["C1"]
    c3 = bundle["C3"]
    c5 = bundle["C5"]
    c6 = bundle["C6"]
    kn = bundle["key_numbers"]

    prim = c1["primary"]["rows"]
    conv = c1["variants"]["conv_only"]["rows"]
    n = min(len(prim), len(conv))
    competing = np.column_stack(
        [
            [r["alpha_over_ac"] for r in prim[:n]],
            [r["K_Tstar"] for r in prim[:n]],
            [r["K_Tstar"] for r in conv[:n]],
            [r["kappa_meV"] for r in prim[:n]],
            [r["kappa_meV"] for r in conv[:n]],
            [r["delta"] for r in prim[:n]],
        ]
    )
    dump(
        DATA / "competing_scales.dat",
        "# alpha_over_ac  K_primary  K_conv  kappa_primary  kappa_conv  delta",
        competing,
    )

    conv_scan = c3["scan_Tstar_conv_only"]["rows"]
    plasma_flags = np.array([float(r["plasma"]) for r in conv_scan[:n]])
    dump(
        DATA / "competing_scales_ir.dat",
        "# alpha_over_ac  plasma_conv",
        np.column_stack(
            [np.array([r["alpha_over_ac"] for r in conv_scan[:n]]), plasma_flags]
        ),
    )

    for st in c3["flows_Tstar"]:
        ell, K, y = integrate_kt(st["K_bare"], st["y_bare"])
        ell, K, y = trim_flow(ell, K, y)
        rows = np.column_stack([ell, K, y, 1.0 / np.maximum(K, 1e-12)])
        dump(DATA / f"vortex_rg_flow_{st['name']}.dat", "# ell  K  y  Kinv", rows)

    plasma_row = next(r for r in conv_scan if r["phase"] == "plasma")
    ell, K, y = integrate_kt(plasma_row["K_bare"], plasma_row["y_bare"])
    y = np.maximum(y, 1e-16)
    dump(
        DATA / "vortex_rg_flow_conv_plasma.dat",
        f"# ell  K  y  Kinv  # conv-only IR plasma at alpha/ac={plasma_row['alpha_over_ac']:.4f}",
        np.column_stack([ell, K, y, 1.0 / np.maximum(K, 1e-12)]),
    )

    teff_rows = []
    for rec in sorted(c3["T_eff_grid"], key=lambda r: r["T_K"]):
        p_a = rec["primary_alpha_kt_over_ac"]
        c_a = rec["conv_only_alpha_kt_over_ac"]
        teff_rows.append(
            [
                rec["T_K"],
                1.0 if rec["primary_kt_before_fold"] else 0.0,
                1.0 if p_a is None else p_a,
                1.0 if rec["conv_only_kt_before_fold"] else 0.0,
                1.0 if c_a is None else c_a,
            ]
        )
    dump(
        DATA / "Teff_scan.dat",
        "# T_K  primary_kt  primary_alpha_or_fold  conv_kt  conv_alpha_or_fold",
        np.asarray(teff_rows),
    )

    star = min(prim, key=lambda r: abs(r["alpha_over_ac"] - 0.90))
    K_star = float(star["K_Tstar"])
    eta_star = 1.0 / (2.0 * PI * K_star)
    r = np.logspace(0.0, math.log10(64.0), 64)
    g_star = r ** (-eta_star)
    g_kt = r ** (-0.25)
    g_plasma = np.exp(-(r - 1.0) / 8.0)
    dump(
        DATA / "g_r_scaling.dat",
        f"# r  g_star  g_kt  g_plasma  # K(lambda*)={K_star:.4f} eta={eta_star:.6f}",
        np.column_stack([r, g_star, g_kt, g_plasma]),
    )

    eta_rows = []
    for rec in prim:
        K = float(rec["K_Tstar"])
        eta = 1.0 / (2.0 * PI * K)
        eta_rows.append([rec["alpha_over_ac"], K, eta, eta, eta])
    dump(
        DATA / "eta_eff.dat",
        "# alpha_over_ac  K_primary  eta_sw  eta_sw  eta_sw",
        np.asarray(eta_rows),
    )

    ahns = c6.get("scan_GPE_rows") or []
    if ahns:
        dump(
            DATA / "ahns_ratio.dat",
            "# alpha_over_ac  kappa  ratio",
            np.array(
                [[r["alpha_over_ac"], r["kappa_meV"], r["ratio"]] for r in ahns],
                dtype=float,
            ),
        )

    cur = c5["curves"]
    dump(
        DATA / "delta_vs_tA.dat",
        "# t_A  delta_phys  delta_sqrt  delta_bkt  delta_c",
        np.column_stack(
            [
                np.asarray(cur["t_A"], dtype=float),
                np.asarray(cur["delta_phys"], dtype=float),
                np.asarray(cur["delta_sqrt_local"], dtype=float),
                np.asarray(cur["delta_bkt_toy"], dtype=float),
                np.full(len(cur["t_A"]), kn["delta_c"]),
            ]
        ),
    )

    summary = DATA / "vortex_rg_summary.txt"
    summary.write_text(
        "QVC vortex-RG C1–C6 summary (frozen α, Option A*)\n"
        f"kappa_UV={kn['kappa_uv_meV']:.4f} meV  kappa_fold={kn['kappa_fold_meV']:.4f} meV\n"
        f"K(0,T*)={kn['K_uv_Tstar']:.3f}  K(alpha_c,T*)={kn['K_fold_Tstar']:.3f}  K_c=2/pi={kn['K_c_bosonic']:.4f}\n"
        f"E_core/kappa={kn['E_core_over_kappa']:.4f}  a_c=xi={kn['a_c_um']:.4f} um\n"
        f"primary T* KT-before-fold={kn['kt_before_fold_Tstar_primary']}\n"
        f"conv-only T* KT-before-fold={kn['kt_before_fold_Tstar_conv_only']}  "
        f"IR plasma alpha/ac={c3['scan_Tstar_conv_only']['alpha_kt_ir_over_ac']:.3f}\n"
        f"C5 free slope={kn['C5_fold_slope_free']:.3f}  dlle={kn['C5_delta_ll_sqrt_minus_essential']:.1f}\n"
        f"C6 jump_vs_crossover={kn['C6_jump_vs_crossover']}\n"
        "No T_BKT^bare input. No diagnostic beta_alpha. CS truncated in default scan.\n"
    )
    print("wrote", DATA)
    print("K_star", K_star, "eta_star", eta_star)
    print("conv plasma station", plasma_row["alpha_over_ac"], plasma_row["K_bare"], plasma_row["y_bare"])


if __name__ == "__main__":
    main()
