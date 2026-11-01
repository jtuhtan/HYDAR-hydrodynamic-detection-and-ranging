"""End-to-end: sensors record the drop and HYDAR recovers its location."""

import numpy as np
import pytest

from hydar import DropObject, Grid, SensorArray, Solver, locate, wavenumber


def test_dispersion_relation():
    h, g = 2.0, 9.81
    k = wavenumber(np.array([0.5, 2.0, 10.0]), h, g)
    assert np.allclose(g * k * np.tanh(k * h), [0.25, 4.0, 100.0])


@pytest.fixture(scope="module")
def scene():
    grid = Grid(lx=16.0, ly=16.0, dx=0.05, n_pml=30)
    sim = Solver(grid, depth=1.0)
    arrays = [SensorArray("A", 3.0, 3.0, z=-0.3),
              SensorArray("B", 13.0, 3.0, z=-0.6),
              SensorArray("C", 8.0, 13.0, z=-0.5)]
    for a in arrays:
        sim.add(a)
    ball = DropObject(diameter=0.5, density=2500.0, x=6.0, y=8.5)
    sim.add(ball)
    sim.run(5.0)
    return sim, arrays, ball


def test_sampling_rate(scene):
    _, arrays, _ = scene
    a = arrays[0]
    assert np.allclose(np.diff(a.time), 0.01)
    assert a.pressure.shape == (9, a.time.size)
    # static head rho*g*|z| before the waves arrive
    assert a.pressure[0, 0] == pytest.approx(1000 * 9.81 * 0.3)


def test_bearings(scene):
    sim, arrays, ball = scene
    res = locate(arrays, c=sim.c, t0=0.0)
    for f in res["fixes"]:
        true = np.degrees(np.arctan2(ball.y - f.center[1], ball.x - f.center[0]))
        assert abs((f.bearing_deg - true + 180) % 360 - 180) < 2.0


def test_ranges_and_fused_position(scene):
    sim, arrays, ball = scene
    res = locate(arrays, c=sim.c, t0=0.0)
    for f in res["fixes"]:
        true = np.hypot(ball.x - f.center[0], ball.y - f.center[1])
        assert f.range == pytest.approx(true, abs=0.5)
    assert np.linalg.norm(res["position"] - [ball.x, ball.y]) < 0.3


def test_unknown_drop_time(scene):
    sim, arrays, ball = scene
    res = locate(arrays, c=sim.c)
    assert np.linalg.norm(res["position"] - [ball.x, ball.y]) < 0.3
