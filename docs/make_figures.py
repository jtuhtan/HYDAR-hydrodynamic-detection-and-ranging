"""Regenerate the static README figures in Jeff Tuhtan's paper figure style.

Run from the repository root:   python docs/make_figures.py
Writes docs/img/figNN_*.png (600 dpi) plus .pdf and .svg (editable text).
Sizes assume a 170 mm page text width (A4, 20 mm margins).
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "docs"))
import paperfig as pf  # noqa: E402
from hydar import DropObject, Grid, SensorArray, Solver, locate  # noqa: E402

OUT = ROOT / "docs" / "img"
pf.page(170)
font = pf.use(grid=False)
print(f"font: {font}")
C = pf.C
FIELD_CMAP = "RdBu_r"


def field(ax, sim, eta, vlim):
    return ax.imshow(eta.T * 1e3, origin="lower", extent=sim.grid.extent, cmap=FIELD_CMAP,
                     vmin=-vlim * 1e3, vmax=vlim * 1e3, interpolation="bilinear", rasterized=True)


# ------------------------------------------------------------------ fig 1
def fig_open_boundary():
    def run(n_pml, snaps_at=()):
        grid = Grid(lx=20.0, ly=12.0, dx=0.05, n_pml=n_pml)
        sim = Solver(grid, depth=1.0)
        X, Y = np.meshgrid(grid.x, grid.y, indexing="ij")
        sim.set_eta(0.01 * np.exp(-((X - 7) ** 2 + (Y - 6) ** 2) / 0.4 ** 2))
        e0, t, e, snaps = sim.energy(), [0.0], [1.0], {}
        targets = list(snaps_at)
        if targets and targets[0] == 0:
            snaps[0.0] = sim.eta[grid.interior].copy()
            targets.pop(0)
        while sim.t < 8.0:
            sim.step()
            t.append(sim.t)
            e.append(sim.energy() / e0)
            if targets and sim.t >= targets[0]:
                snaps[targets.pop(0)] = sim.eta[grid.interior].copy()
        return sim, np.array(t), np.array(e), snaps

    sim, t, e_pml, snaps = run(40, (1.5, 3.0, 4.5))
    _, t2, e_rigid, _ = run(0)

    fig, axes = pf.figure("full", height_mm=46, ncols=4,
                          gridspec_kw={"width_ratios": [1, 1, 1, 1.05]})
    for ax, (ts, eta) in zip(axes[:3], snaps.items()):
        m = field(ax, sim, eta, 2e-3)
        ax.set_aspect("equal")
        ax.set_xticks([0, 10, 20])
        ax.set_yticks([0, 6, 12])
        ax.text(0.97, 0.95, f"t = {ts:g} s", transform=ax.transAxes, ha="right", va="top",
                fontsize=7.5)
        ax.set_xlabel("x (m)")
    axes[0].set_ylabel("y (m)")
    for ax in axes[1:3]:
        ax.tick_params(labelleft=False)
    pf.colorbar(fig, m, list(axes[:3]), "η (mm)")

    ax = axes[3]
    ax.grid(True)
    ax.semilogy(t2, e_rigid, color=C["dark"], lw=1.0, label="Rigid walls")
    ax.semilogy(t, np.maximum(e_pml, 1e-12), color=C["blue"], lw=1.2, label="PML walls")
    ax.set(xlabel="Time (s)", ylabel="Wave energy, E/E$_0$ (–)", xlim=(0, 8), ylim=(1e-6, 30))
    ax.set_yticks([1e-6, 1e-4, 1e-2, 1])
    ax.set_xticks([0, 2, 4, 6, 8])
    ax.text(7.9, 1.7, "Rigid walls", color=C["dark"], ha="right", va="bottom", fontsize=7.5)
    ax.text(5.6, 2e-3, "PML walls", color=C["blue"], ha="left", va="bottom", fontsize=7.5)
    pf.panel_labels(axes[:3], loc="inside")
    pf.panel_labels([axes[3]], start=3)
    pf.save(fig, "fig01_open_boundary", outdir=OUT)


# ------------------------------------------------------------------ fig 2
def fig_drop():
    D, depth = 0.2, 0.5
    runs = {}
    for name, rho in (("Floating, 600 kg/m$^3$", 600.0), ("Sinking, 2500 kg/m$^3$", 2500.0)):
        grid = Grid(lx=8.0, ly=8.0, dx=D / 10, n_pml=40)
        sim = Solver(grid, depth=depth)
        ball = DropObject(diameter=D, density=rho, x=4.0, y=4.0, z0=0.05)
        sim.add(ball)
        out = sim.run(1.2, gauges=[(5.0, 4.0)])
        snap = sim.eta[grid.interior].copy()
        out2 = sim.run(5.0, gauges=[(5.0, 4.0)])
        g_t = np.r_[out["t"], out2["t"]]
        g_e = np.r_[out["gauges"][:, 0], out2["gauges"][:, 0]]
        runs[name] = (sim, ball, snap, g_t, g_e)

    fig = plt.figure(figsize=(pf.width_mm("full") * pf.MM, 128 * pf.MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.45, 1])
    f1, f2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    lz = fig.add_subplot(gs[1, 0])
    le = fig.add_subplot(gs[1, 1], sharex=lz)
    colors = [C["blue"], C["red"]]
    vlim = max(np.abs(r[2]).max() for r in runs.values())
    for ax, (name, (sim, ball, snap, _, _)), col in zip((f1, f2), runs.items(), colors):
        m = field(ax, sim, snap, vlim)
        ax.set_aspect("equal")
        ax.set(xlabel="x (m)", xticks=[0, 2, 4, 6, 8], yticks=[0, 2, 4, 6, 8])
        ax.text(0.97, 0.95, f"{name.split(',')[0]}, t = 1.2 s", transform=ax.transAxes,
                ha="right", va="top", fontsize=7.5, color=col)
        ax.plot(5.0, 4.0, marker="v", color=C["ink"], ms=4, ls="none")
    f1.set_ylabel("y (m)")
    f2.tick_params(labelleft=False)
    pf.colorbar(fig, m, [f1, f2], "η (mm)")

    lz.grid(True)
    le.grid(True)
    lines = []
    for (name, (sim, ball, _, g_t, g_e)), col in zip(runs.items(), colors):
        h = ball.history
        lines.append(lz.plot(h["t"], np.array(h["z"]) - ball.R, color=col, lw=1.1,
                             label=name.split(",")[0])[0])
        le.plot(g_t, g_e * 1e3, color=col, lw=1.0, label=name.split(",")[0])
    ball = runs["Floating, 600 kg/m$^3$"][1]
    pf.ref_line(lz, y=-ball.equilibrium_draft, color=C["blue"])
    lz.text(4.95, -ball.equilibrium_draft - 0.03, "Archimedes draft", color=C["blue"],
            ha="right", va="top", fontsize=7.5)
    pf.ref_line(lz, y=-depth, text="Bottom", color=C["dark"], where="right", ls="-")
    lz.set(xlabel="Time (s)", ylabel="Sphere bottom, z (m)", xlim=(0, 5), ylim=(-0.56, 0.08))
    le.set(xlabel="Time (s)", ylabel="η at 1 m from drop (mm)", ylim=(-6, 14))
    pf.legend(le, loc="upper right")
    pf.panel_labels([f1, f2], loc="inside")
    pf.panel_labels([lz, le], start=2)
    pf.save(fig, "fig02_drop_sphere", outdir=OUT)


# ------------------------------------------------------------------ fig 3 + 4
def hydar_scene():
    grid = Grid(lx=30.0, ly=30.0, dx=0.05, n_pml=40)
    sim = Solver(grid, depth=1.0)
    arrays = [SensorArray("A", 6.0, 6.0, z=-0.2), SensorArray("B", 24.0, 6.0, z=-0.5),
              SensorArray("C", 24.0, 24.0, z=-0.8), SensorArray("D", 6.0, 24.0, z=-0.5)]
    for a in arrays:
        sim.add(a)
    ball = DropObject(diameter=0.5, density=2500.0, x=11.0, y=17.0)
    sim.add(ball)
    sim.run(3.0)
    snap = sim.eta[grid.interior].copy()
    sim.run(9.0)
    return sim, arrays, ball, snap, locate(arrays, c=sim.c, t0=0.0)


def fig_timeseries(arrays, res):
    fig, axes = pf.stack(4, width="full", row_mm=24)
    for ax, a, f in zip(axes, arrays, res["fixes"]):
        ax.grid(True)
        p = a.dynamic_pressure
        for k in range(a.n_sensors):
            ax.plot(a.time, p[k], color=C["dark"], lw=0.6, alpha=0.85)
        first = np.nanmin(f.arrivals)
        last = np.nanmax(f.arrivals)
        pf.shade_span(ax, first, last)
        ax.text(0.99, 0.92, f"Array {a.name}, z = {a.z:g} m".replace("-", "\u2212"), transform=ax.transAxes,
                ha="right", va="top", fontsize=7.5)
        ax.set(ylabel="p$_{dyn}$ (Pa)", ylim=(-90, 150), yticks=[-50, 0, 50, 100])
    axes[0].annotate("arrivals across the array", xy=(np.mean([np.nanmin(res["fixes"][0].arrivals),
                     np.nanmax(res["fixes"][0].arrivals)]), 1),
                     xycoords=("data", "axes fraction"), xytext=(0, 3), textcoords="offset points",
                     ha="center", va="bottom", color=C["cyan"], fontsize=7.5)
    axes[-1].set(xlabel="Time since drop (s)", xlim=(1.5, 7.5))
    pf.panel_labels(axes)
    pf.save(fig, "fig03_hydar_timeseries", outdir=OUT)


def fig_map(sim, arrays, ball, snap, res):
    fig, ax = pf.figure("half", height_mm=86)
    m = field(ax, sim, snap, np.abs(snap).max())
    for f in res["fixes"]:
        ax.plot([f.center[0], f.position[0]], [f.center[1], f.position[1]], color=C["ink"],
                lw=0.8, ls="--", zorder=3)
    for a in arrays:
        ax.plot(*a.positions.T, "s", color=C["ink"], ms=1.8, zorder=4)
        left = a.x < 15
        ax.annotate(a.name, (a.positions[:, 0].min() if left else a.positions[:, 0].max(), a.y),
                    xytext=(-5 if left else 5, 0), textcoords="offset points",
                    ha="right" if left else "left", va="center", fontsize=7.5)
    ax.plot(ball.x, ball.y, "o", mfc="none", mec=C["ink"], ms=7, mew=0.9, label="True drop", zorder=5)
    err = 100 * np.linalg.norm(res["position"] - [ball.x, ball.y])
    ax.plot(*res["position"], "+", color=C["red"], ms=8, mew=1.4,
            label=f"HYDAR fix, error {err:.0f} cm", zorder=6)
    ax.plot([], [], color=C["ink"], lw=0.8, ls="--", label="Bearing and range per array")
    ax.set(xlabel="x (m)", ylabel="y (m)", xlim=(0, 30), ylim=(0, 30),
           xticks=[0, 10, 20, 30], yticks=[0, 10, 20, 30])
    ax.set_aspect("equal")
    pf.colorbar(fig, m, ax, "η (mm)")
    fig.legend(*ax.get_legend_handles_labels(), loc="outside upper center", ncol=1,
               frameon=False)
    pf.save(fig, "fig04_hydar_map", outdir=OUT)


if __name__ == "__main__":
    which = sys.argv[1:] or ["1", "2", "3"]
    if "1" in which:
        fig_open_boundary()
    if "2" in which:
        fig_drop()
    if "3" in which:
        sim, arrays, ball, snap, res = hydar_scene()
        fig_timeseries(arrays, res)
        fig_map(sim, arrays, ball, snap, res)
