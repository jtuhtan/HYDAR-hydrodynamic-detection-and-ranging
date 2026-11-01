import numpy as np
import pytest

from hydar import DropObject, Grid, Solver, cap_volume


def _sim(density, depth=0.5, t_end=2.0, D=0.1):
    grid = Grid(lx=4.0, ly=4.0, dx=D / 10, n_pml=30)
    sim = Solver(grid, depth=depth)
    ball = DropObject(diameter=D, density=density, x=2.0, y=2.0)
    sim.add(ball)
    sim.run(t_end)
    return sim, ball


def test_cap_volume_limits():
    R = 0.3
    assert cap_volume(0, R) == 0
    assert cap_volume(2 * R, R) == pytest.approx(4 / 3 * np.pi * R ** 3)
    assert cap_volume(R, R) == pytest.approx(2 / 3 * np.pi * R ** 3)


def test_displaced_volume_enters_the_water():
    # Before waves reach the PML, the volume of eta equals V_sub.
    sim, ball = _sim(density=500.0, t_end=0.3)
    assert sim.volume() == pytest.approx(ball.V_sub, rel=0.02)


def test_light_sphere_floats_at_archimedes_draft():
    sim, ball = _sim(density=500.0, t_end=8.0)
    draft = ball.R - ball.z
    assert draft == pytest.approx(ball.equilibrium_draft, rel=0.05)
    assert abs(ball.w) < 1e-2


def test_heavy_sphere_sinks_to_bottom():
    sim, ball = _sim(density=2500.0, depth=0.5, t_end=3.0)
    assert ball.on_bottom
    assert ball.V_sub == pytest.approx(ball.V)


def test_waves_radiate_outward():
    grid = Grid(lx=4.0, ly=4.0, dx=0.01, n_pml=30)
    sim = Solver(grid, depth=0.5)
    sim.add(DropObject(diameter=0.1, density=2500.0, x=2.0, y=2.0))
    out = sim.run(0.8, gauges=[(2.5, 2.0), (3.0, 2.0)])
    arrival = [out["t"][np.argmax(np.abs(out["gauges"][:, k]) > 1e-5)] for k in range(2)]
    # front travels at about c = sqrt(g h) between the two gauges
    c_est = 0.5 / (arrival[1] - arrival[0])
    assert c_est == pytest.approx(np.sqrt(9.81 * 0.5), rel=0.15)
