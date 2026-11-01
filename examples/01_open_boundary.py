"""A Gaussian hump of water spreads out and leaves through the absorbing walls.

Run:  python examples/01_open_boundary.py   ->  open_boundary.png
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from hydar import Grid, Solver

grid = Grid(lx=20.0, ly=12.0, dx=0.05, n_pml=40)
sim = Solver(grid, depth=1.0)
X, Y = np.meshgrid(grid.x, grid.y, indexing="ij")
sim.set_eta(0.01 * np.exp(-((X - 7) ** 2 + (Y - 6) ** 2) / 0.4 ** 2))
e0 = sim.energy()

times = [0.0, 1.5, 3.0, 4.5, 8.0]
fig, axes = plt.subplots(1, len(times), figsize=(16, 2.8))
for t, ax in zip(times, axes):
    sim.run(t)
    eta = sim.eta[grid.interior]
    ax.imshow(eta.T, origin="lower", extent=grid.extent, cmap="RdBu_r", vmin=-2e-3, vmax=2e-3)
    ax.set_title(f"t = {sim.t:.1f} s, E/E0 = {sim.energy() / e0:.1e}", fontsize=9)
    ax.set_xlabel("x [m]")
axes[0].set_ylabel("y [m]")
plt.tight_layout()
plt.savefig("open_boundary.png", dpi=110)
print(f"energy left in the domain after {sim.t:.1f} s: {sim.energy() / e0:.2e} of the initial value")
