"""Series 60 hull geometry and HYDAR boat tracking (Mach cone and Kelvin wake)."""

import numpy as np
import pytest

from hydar import Boat, Grid, SensorArray, Solver, SpectralSolver, track_boat, wake_arrival
from hydar.boat import SERIES60, Series60Hull


def test_series60_form_coefficients():
    hull = Series60Hull(length=4.5)
    assert hull.beam == pytest.approx(4.5 / 7.5)
    assert hull.draft == pytest.approx(hull.beam / 2.5)
    co = hull.coefficients()
    for key in ("CB", "CP", "CM", "CWP", "LCB"):
        assert co[key] == pytest.approx(SERIES60[key], abs=0.003), key


def test_boat_adds_no_net_water():
    grid = Grid(lx=20.0, ly=10.0, dx=0.05, n_pml=30)
    sim = Solver(grid, depth=1.0)
    boat = Boat(speed=2.0, x=6.0, y=5.0, ramp=0.5)
    sim.add(boat)
    sim.run(0.6)          # fully ramped, waves still inside the domain
    # displaced water raised around the hull equals the hull volume
    assert sim.volume() == pytest.approx(boat.volume, rel=0.03)


def test_kelvin_wake_arrival_geometry():
    # the Kelvin cusp line (19.47 deg) passes a point when the boat is 2*sqrt(2)*d past it
    t = wake_arrival(np.array([[10.0, 3.0]]), (0.0, 0.0), 0.0, 2.0)
    assert t[0] == pytest.approx((10 + 3 * 2 * np.sqrt(2)) / 2)


def test_wake_arrival_geometry():
    # a point abeam of the track is reached when the boat is d*sqrt(M^2-1) past it
    t = wake_arrival(np.array([[10.0, 3.0]]), (0.0, 0.0), 0.0, 6.0, 3.0)
    assert t[0] == pytest.approx((10 + 3 * np.sqrt(3)) / 6)


@pytest.mark.slow
def test_track_supercritical_boat():
    grid = Grid(lx=30.0, ly=20.0, dx=0.05, n_pml=60)
    sim = Solver(grid, depth=1.0)
    arrays = [SensorArray("A", 8.0, 3.0), SensorArray("B", 22.0, 3.5),
              SensorArray("C", 9.0, 17.0), SensorArray("D", 23.0, 17.0)]
    for a in arrays:
        sim.add(a)
    boat = Boat(speed=5.0, x=-2.5, y=9.0, heading_deg=5.0)
    sim.add(boat)
    sim.run(9.0)
    tr = track_boat(arrays, c=sim.c, wake="mach")
    assert tr.speed == pytest.approx(5.0, rel=0.03)
    assert tr.heading_deg == pytest.approx(5.0, abs=1.0)
    bow = boat.position(4.0) + 0.5 * boat.length * boat.e
    assert np.linalg.norm(tr.position(4.0) - bow) < 0.5


@pytest.mark.slow
def test_track_boat_from_kelvin_wake():
    grid = Grid(lx=90.0, ly=40.0, dx=0.15, n_pml=60)
    sim = SpectralSolver(grid, depth=3.0)
    arrays = [SensorArray("A", 50.0, 10.0, z=-0.4, pressure_model="solver"),
              SensorArray("B", 62.0, 11.0, z=-0.4, pressure_model="solver"),
              SensorArray("C", 52.0, 30.0, z=-0.4, pressure_model="solver")]
    for a in arrays:
        sim.add(a)
    boat = Boat(speed=3.0, x=-7.0, y=20.0, heading_deg=0.0, ramp=4.0)
    sim.add(boat)
    sim.run(36.0)
    tr = track_boat(arrays)
    assert tr.speed == pytest.approx(3.0, rel=0.05)
    assert tr.heading_deg == pytest.approx(0.0, abs=2.0) or tr.heading_deg == pytest.approx(360.0, abs=2.0)
    assert np.linalg.norm(tr.position(20.0) - boat.position(20.0)) < 2.5
