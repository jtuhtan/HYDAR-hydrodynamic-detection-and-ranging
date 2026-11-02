"""Case study 3: Kelvin wake theory and vessel localization."""

import numpy as np
import pytest

from hydar import SensorArray
from hydar.kelvin import THETA_CUSP, KelvinShip, branch_frequencies, locate_vessel, record_passage


def test_branches_meet_at_cusp_and_approach_transverse_frequency():
    U, d = 10.0, 1000.0
    f_t = 9.81 / (2 * np.pi * U)
    t_c = 2 * np.sqrt(2) * d / U
    ft, fd = branch_frequencies([t_c, 50 * t_c], U, d, 0.0)
    assert ft[0] == pytest.approx(fd[0], rel=1e-6)
    assert ft[0] == pytest.approx(np.sqrt(1.5) * f_t, rel=1e-6)
    assert ft[1] == pytest.approx(f_t, rel=0.01)
    assert np.degrees(THETA_CUSP) == pytest.approx(35.26, abs=0.01)


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
