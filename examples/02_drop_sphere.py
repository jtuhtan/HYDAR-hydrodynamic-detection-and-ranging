"""Drop a floating and a sinking sphere and compare their motion and waves.

Run:  python examples/02_drop_sphere.py   ->  drop_sphere.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hydar import DropObject, Grid, Solver, recommend_grid

D, depth = 0.2, 0.5
print("suggested grid:", recommend_grid(D, depth))

fig, axes = plt.subplots(2, 2, figsize=(11, 7))
for col, rho_s in enumerate([600.0, 2500.0]):
    grid = Grid(lx=8.0, ly=8.0, dx=D / 10, n_pml=40)
    sim = Solver(grid, depth=depth)
    ball = DropObject(diameter=D, density=rho_s, x=4.0, y=4.0, z0=0.05)
    sim.add(ball)
    out = sim.run(1.2, gauges=[(5.0, 4.0), (6.0, 4.0)])
    snap = sim.eta[grid.interior].copy()
    sim.run(5.0)

    ax = axes[0, col]
    v = np.abs(snap).max()
    im = ax.imshow(snap.T, origin="lower", extent=grid.extent, cmap="RdBu_r", vmin=-v, vmax=v)
    plt.colorbar(im, ax=ax, label="eta [m]")
    ax.set_title(f"density {rho_s:g} kg/m3: surface at t = 1.2 s")

    ax = axes[1, col]
    h = ball.history
    ax.plot(h["t"], np.array(h["z"]) - ball.R, label="sphere bottom z [m]")
    ax.axhline(0, color="k", lw=0.5)
    if ball.equilibrium_draft is not None:
        ax.axhline(-ball.equilibrium_draft, ls=":", color="C0", label="Archimedes draft")
    ax.plot(out["t"], 20 * out["gauges"][:, 0], label="20 x eta at 1 m [m]")
    ax.set_xlabel("time [s]")
    ax.legend(fontsize=8)
    state = "floating" if ball.density < 1000 else ("on bottom" if ball.on_bottom else "sinking")
    print(f"density {rho_s:g}: final state {state}, bottom at z = {ball.z - ball.R:.3f} m")

plt.tight_layout()
plt.savefig("drop_sphere.png", dpi=110)
