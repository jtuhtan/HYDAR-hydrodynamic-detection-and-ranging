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


def test_pick_arrivals_empty_window_gives_no_picks():
    from hydar import pick_arrivals
    t = np.arange(0, 5, 0.01)
    assert np.isnan(pick_arrivals(t, np.ones((3, t.size)), t_start=10.0)).all()


def test_parallel_bearings_do_not_intersect():
    from hydar import ArrayFix, intersect_bearings
    d = np.array([1.0, 0.0])
    fixes = [ArrayFix("a", np.array([0.0, y]), None, 0.0, 0.0, d, 1.0, None, None) for y in (0.0, 1.0)]
    with pytest.raises(ValueError, match="parallel"):
        intersect_bearings(fixes)


def test_sensor_noise_is_not_filtered_by_depth():
    # 2 Pa of sensor noise stays 2 Pa below 2.8 m of water (the Airy filter must not apply to it)
    t = np.arange(0, 200, 0.01)
    a = SensorArray("H", 0.0, 0.0, z=-2.8, pressure_model="airy", noise_pa=2.0, noise_seed=3)
    a.load(t, np.zeros((9, t.size)), depth=3.0)
    p = a.dynamic_pressure
    assert p.std() == pytest.approx(2.0, rel=0.02)
    assert np.array_equal(p, a.dynamic_pressure)          # same seed, same record
    quiet = SensorArray("H", 0.0, 0.0, z=-2.8, pressure_model="airy")
    quiet.load(t, np.zeros((9, t.size)), depth=3.0)
    assert not quiet.dynamic_pressure.any()


def test_custom_layout_records_in_a_solver():
    layout = np.array([[-1.0, 0.0], [0.3, 0.2], [0.9, -0.6], [0.0, 1.1]])
    sim = Solver(Grid(lx=8.0, ly=8.0, dx=0.05, n_pml=20), depth=1.0)
    a = SensorArray("L", 4.0, 4.0, layout=layout, rotation_deg=90.0)
    sim.add(a)
    sim.add(DropObject(diameter=0.4, density=2500.0, x=2.0, y=4.0))
    sim.run(1.5)
    assert a.n_sensors == 4 and a.eta.shape == (4, a.time.size)
    assert np.allclose(a.positions[0], [4.0, 3.0])        # (-1, 0) rotated by 90 deg
    assert np.abs(a.eta).max() > 1e-4
