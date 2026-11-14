"""Case study 3: Kelvin wake theory and vessel localization."""

import numpy as np
import pytest

from hydar import SensorArray
from hydar.kelvin import (THETA_CUSP, KelvinShip, branch_frequencies, fit_chirp, locate_vessel,
                          record_passage, wave_direction, wind_sea)


def test_branches_meet_at_cusp_and_approach_transverse_frequency():
    U, d = 10.0, 1000.0
    f_t = 9.81 / (2 * np.pi * U)
    t_c = 2 * np.sqrt(2) * d / U
    ft, fd = branch_frequencies([t_c, 50 * t_c], U, d, 0.0)
    assert ft[0] == pytest.approx(fd[0], rel=1e-6)
    assert ft[0] == pytest.approx(np.sqrt(1.5) * f_t, rel=1e-6)
    assert ft[1] == pytest.approx(f_t, rel=0.01)
    assert np.degrees(THETA_CUSP) == pytest.approx(35.26, abs=0.01)


def _plane_wave_array(theta_deg=-50.0, f0=0.15, t_end=120.0):
    a = SensorArray("H", 0.0, 0.0, spacing=2.5, rotation_deg=22.5)
    t = np.arange(0, t_end, 0.01)
    k = (2 * np.pi * f0) ** 2 / 9.81
    n = np.array([np.cos(np.radians(theta_deg)), np.sin(np.radians(theta_deg))])
    a.load(t, 0.05 * np.cos(k * (a.positions @ n)[:, None] - 2 * np.pi * f0 * t), depth=3.0)
    return a, n, 2 * np.pi * f0 / k


def test_wave_direction_of_a_plane_wave():
    a, n, c = _plane_wave_array()
    n_est, c_est = wave_direction(a, 30.0, 90.0)
    assert np.degrees(np.arccos(np.clip(n_est @ n, -1, 1))) < 1.0
    assert c_est == pytest.approx(c, rel=0.08)


def test_analysis_outside_the_record_is_a_clear_error():
    a, _, _ = _plane_wave_array()
    with pytest.raises(RuntimeError, match="too few samples"):
        wave_direction(a, 500.0, 540.0)
    with pytest.raises(RuntimeError, match="ridge points"):
        fit_chirp(np.array([]), 100.0)


def test_wind_sea_has_the_requested_height():
    t = np.arange(0, 1800, 0.25)
    for seed in range(3):
        eta = wind_sea(np.zeros((1, 2)), t, hs=0.3, f_peak=0.4, seed=seed)
        assert 4 * eta.std() == pytest.approx(0.3, rel=0.05)
    # pressure head at depth is smaller than the surface elevation, same components
    head = wind_sea(np.zeros((1, 2)), t, hs=0.3, f_peak=0.4, seed=0, depth=3.0, z=-2.8)
    assert 0.1 < head.std() / wind_sea(np.zeros((1, 2)), t, hs=0.3, f_peak=0.4, seed=0).std() < 0.9


@pytest.mark.slow
def test_locate_passing_ferry():
    arr = SensorArray("H", 0.0, 0.0, z=-2.8, spacing=2.5, rotation_deg=22.5, pressure_model="airy")
    course = -69.0
    e = np.array([np.cos(np.radians(course)), np.sin(np.radians(course))])
    n = np.array([-e[1], e[0]])
    p0 = -2000.0 * n - 11.0 * 200.0 * e
    ship = KelvinShip(speed=11.0, course_deg=course, x0=p0[0], y0=p0[1])
    record_passage([arr], ship, t_end=1300.0, depth=3.0,
                   wind={"hs": 0.15, "f_peak": 0.4, "mean_dir_deg": 112.5})
    fix = locate_vessel(arr)
    assert fix.speed == pytest.approx(11.0, rel=0.03)
    assert fix.distance == pytest.approx(2000.0, rel=0.15)
    assert abs((fix.course_deg - course + 180) % 360 - 180) < 5
    assert np.linalg.norm(fix.position(fix.t_cusp) - ship.position(fix.t_cusp)) < 600
