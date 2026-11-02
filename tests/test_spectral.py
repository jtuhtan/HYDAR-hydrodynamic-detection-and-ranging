"""Dispersive solver: exact linear dispersion and absorbing sponge."""

import numpy as np
import pytest

from hydar import Grid, SpectralSolver


def test_standing_wave_has_exact_dispersion():
    grid = Grid(lx=16.0, ly=4.0, dx=0.1, n_pml=0)
    sim = SpectralSolver(grid, depth=2.0, dt=0.01)
    lx = grid.nx * grid.dx
    k = 2 * np.pi * 4 / lx                         # 4 wavelengths across the periodic domain
    X, _ = np.meshgrid(grid.x, grid.y, indexing="ij")
    sim.set_eta(0.01 * np.cos(k * X))
    omega = np.sqrt(9.81 * k * np.tanh(k * 2.0))
    T = 2 * np.pi / omega
    sim.run(2.5 * T)                               # half-period phase: eta flips sign
    expected = 0.01 * np.cos(omega * sim.t) * np.cos(k * X)
    assert np.abs(sim.eta - expected).max() < 2e-4


def test_sponge_absorbs_outgoing_waves():
    grid = Grid(lx=30.0, ly=30.0, dx=0.2, n_pml=50)
    sim = SpectralSolver(grid, depth=3.0)
    X, Y = np.meshgrid(grid.x, grid.y, indexing="ij")
    sim.set_eta(0.01 * np.exp(-((X - 15) ** 2 + (Y - 15) ** 2) / 1.0))
    e0 = sim.energy(interior_only=False)
    sim.run(1.0)
    assert sim.energy(interior_only=False) / e0 == pytest.approx(1.0, abs=1e-3)
    sim.run(40.0)
    assert sim.energy() / e0 < 1e-3
