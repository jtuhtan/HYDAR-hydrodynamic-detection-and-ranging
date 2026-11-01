"""The domain walls must absorb outgoing waves (open-water behaviour)."""

import numpy as np

from hydar import Grid, Solver


def _residual_energy(n_pml):
    grid = Grid(lx=20.0, ly=12.0, dx=0.1, n_pml=n_pml)
    sim = Solver(grid, depth=1.0)
    X, Y = np.meshgrid(grid.x, grid.y, indexing="ij")
    sim.set_eta(0.01 * np.exp(-((X - 10) ** 2 + (Y - 6) ** 2) / 0.5 ** 2))
    e0 = sim.energy()
    t_cross = np.hypot(grid.lx, grid.ly) / sim.c
    while sim.t < 2.0 * t_cross:
        sim.step()
    return sim.energy() / e0


def test_pml_absorbs_outgoing_waves():
    assert _residual_energy(n_pml=30) < 1e-4


def test_rigid_walls_keep_energy():
    # Sanity check of the test itself: without the PML the walls reflect.
    assert _residual_energy(n_pml=0) > 0.3
