"""Plotting helpers (matplotlib)."""

from __future__ import annotations

import numpy as np


def plot_timeseries(arrays, dynamic: bool = True, fixes=None, axes=None):
    """One panel per array with the pressure record of every sensor.

    dynamic : plot the wave-induced pressure (True) or total gauge pressure.
    fixes   : optional ArrayFix list; picked arrival times are marked.
    """
    import matplotlib.pyplot as plt

    if axes is None:
        fig, axes = plt.subplots(len(arrays), 1, figsize=(9, 2.4 * len(arrays)),
                                 sharex=True, squeeze=False)
        axes = axes[:, 0]
    for n, (a, ax) in enumerate(zip(arrays, axes)):
        p = a.dynamic_pressure if dynamic else a.pressure
        for k in range(a.n_sensors):
            ax.plot(a.time, p[k], lw=0.9, label=f"s{k}" if n == 0 else None)
            if fixes is not None and np.isfinite(fixes[n].arrivals[k]):
                ax.axvline(fixes[n].arrivals[k], color="k", lw=0.4, alpha=0.4)
        ax.set_ylabel("p_dyn [Pa]" if dynamic else "p [Pa]")
        ax.set_title(f"{a.name}  (x={a.x:g} m, y={a.y:g} m, z={a.z:g} m, {a.fs:g} Hz)",
                     fontsize=9, loc="left")
        ax.grid(alpha=0.3)
    axes[-1].set_xlabel("time [s]")
    if len(arrays):
        axes[0].legend(ncol=9, fontsize=7, loc="upper right")
    return axes


def plot_map(solver, arrays=(), drops=(), result=None, eta=None, ax=None,
             vlim=None):
    """Plan view: surface elevation, arrays, true drop points and estimates."""
    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots(figsize=(7, 6))
    eta = solver.eta[solver.grid.interior] if eta is None else eta
    v = vlim or (np.abs(eta).max() or 1.0)
    im = ax.imshow(eta.T, origin="lower", extent=solver.grid.extent,
                   cmap="RdBu_r", vmin=-v, vmax=v)
    plt.colorbar(im, ax=ax, label="eta [m]", shrink=0.8)
    for a in arrays:
        ax.plot(*a.positions.T, "k.", ms=4)
        ax.annotate(a.name, a.center, xytext=(5, 5), textcoords="offset points", fontsize=8)
    for d in drops:
        ax.plot(d.x, d.y, "*", color="gold", mec="k", ms=14, label="true drop")
    if result is not None:
        for f in result["fixes"]:
            L = f.range if f.range is not None else 0.5 * max(solver.grid.lx, solver.grid.ly)
            end = f.center + L * f.direction
            ax.plot([f.center[0], end[0]], [f.center[1], end[1]], "g--", lw=1)
            if f.position is not None:
                ax.plot(*f.position, "go", ms=5, mfc="none")
        ax.plot(*result["position"], "mx", ms=12, mew=2, label="HYDAR estimate")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_aspect("equal")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        uniq = dict(zip(labels, handles))
        ax.legend(uniq.values(), uniq.keys(), fontsize=8, loc="lower right")
    return ax
