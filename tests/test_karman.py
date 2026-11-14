"""Case study 5: a fish locates a rock from its Karman vortex street."""

import numpy as np
import pytest

from hydar.karman import (Fish, HemisphereRock, KarmanStreet, RiverTurbulence, Survey,
                          TrapezoidalChannel, fish_sensor_positions, lamb_oseen_pressure,
                          lateral_survey, locate_rock)

_trapz = getattr(np, "trapezoid", None) or np.trapz


def test_trapezoidal_channel_geometry_and_manning():
    ch = TrapezoidalChannel(bed_width=4.0, side_slope=2.0, depth=0.25, slope=0.002, manning_n=0.035)
    assert ch.top_width == pytest.approx(5.0)
    assert ch.area == pytest.approx(1.125)
    assert ch.wetted_perimeter == pytest.approx(4 + 0.5 * np.sqrt(5))
    R = ch.area / ch.wetted_perimeter
    assert ch.mean_velocity == pytest.approx(R ** (2 / 3) * np.sqrt(0.002) / 0.035)
    assert ch.local_depth([0.0, 2.0, 2.25, 2.5, 3.0]) == pytest.approx([0.25, 0.25, 0.125, 0.0, 0.0])
    # the strip velocities carry the whole discharge
    y = np.linspace(-2.5, 2.5, 4001)
    assert _trapz(ch.local_velocity(y) * ch.local_depth(y), y) == pytest.approx(ch.discharge, rel=1e-3)


def test_rock_half_out_of_the_water():
    ch = TrapezoidalChannel(depth=0.25)
    rock = HemisphereRock.from_protrusion(ch, 0.0, 0.3, protrusion=0.5)
    assert rock.radius == pytest.approx(0.5)                       # half its height above water
    assert rock.waterline_diameter(0.25) == pytest.approx(2 * np.sqrt(0.5 ** 2 - 0.25 ** 2))
    # depth-averaged diameter of a hemisphere of radius 2h: (sqrt 3 + 2 pi / 3) h
    assert rock.shedding_diameter(0.25) == pytest.approx((np.sqrt(3) + 2 * np.pi / 3) * 0.25, rel=1e-4)


def test_fish_has_five_sensors_on_head_and_middle_thirds_and_none_on_the_tail():
    s = fish_sensor_positions(0.25)
    assert s.size == 10
    assert np.sum(s < 0.25 / 3) == 5
    assert np.sum((s > 0.25 / 3) & (s < 0.5 / 1.5)) == 5
    assert np.all(s < 2 * 0.25 / 3)                                # tail third is bare
    assert np.diff(s) == pytest.approx(np.full(9, 0.25 / 15))
    fish = Fish(1.0, 0.2)
    arr = fish.array()
    assert arr.n_sensors == 10
    assert arr.positions[:, 1] == pytest.approx(np.full(10, 0.2))
    assert arr.positions[:, 0] == pytest.approx(1.0 + s)


def test_lamb_oseen_pressure_is_cyclostrophic():
    G, rc, rho = 0.7, 0.2, 1000.0
    assert lamb_oseen_pressure(0.0, G, rc) == pytest.approx(-rho * G ** 2 * np.log(2) / (4 * np.pi ** 2 * rc ** 2),
                                                            rel=1e-6)
    r = np.linspace(0.05, 1.0, 2001)
    u = G / (2 * np.pi * r) * (1 - np.exp(-(r / rc) ** 2))
    dpdr = np.gradient(lamb_oseen_pressure(r, G, rc), r)
    assert dpdr[5:-5] == pytest.approx((rho * u ** 2 / r)[5:-5], rel=1e-3)


def _street(y=0.0, **kw):
    ch = TrapezoidalChannel()
    return KarmanStreet(ch, HemisphereRock.from_protrusion(ch, 0.0, y), **kw)


def test_street_follows_strouhal_and_von_karman():
    st = _street()
    assert st.f == pytest.approx(0.2 * st.U / st.D)
    assert st.Uc == pytest.approx(0.88 * st.U)
    assert st.row_spacing / st.wavelength == pytest.approx(0.281)
    # von Karman: self-induced speed of an infinite street
    assert st.gamma0 / (2 * st.wavelength) * np.tanh(np.pi * 0.281) == pytest.approx(st.U - st.Uc)
    assert st.stability_number < 0.2          # shallow wake, still a vortex street


def test_street_pattern_is_mirror_symmetric_about_the_rock():
    st = _street(0.0, n_images=0)
    y = np.linspace(-1.5, 1.5, 61)
    A1, A2 = st.harmonic_profile(6.0, y)
    assert A1 == pytest.approx(A1[::-1], rel=1e-6, abs=1e-6)
    assert A1[30] < 0.02 * A1.max()            # the rows cancel at f on the centreline
    assert A2[30] > 3 * A1[30]                 # ... where only 2 f is left
    assert abs(y[np.argmax(A1)]) == pytest.approx(0.5 * st.row_spacing, abs=0.15)   # near each row


def test_street_is_carried_downstream_at_the_convection_speed():
    st = _street(n_images=0)
    t = np.arange(0, 60, 0.05)
    p = st.pressure([6.0, 6.5], [0.5, 0.5], t)
    lags = np.arange(0.5, 2.0, 0.01)
    corr = [np.corrcoef(p[0, :-200], np.interp(t[:-200] + L, t, p[1]))[0, 1] for L in lags]
    assert 0.5 / lags[np.argmax(corr)] == pytest.approx(st.Uc, rel=0.05)


def test_turbulence_level():
    tu = RiverTurbulence(0.5, intensity=0.1, seed=3)
    assert tu.rms == pytest.approx(1000 * 0.05 ** 2, rel=1e-6)
    p = tu.pressure([0.0, 0.5], np.arange(0, 100, 0.01))
    # frozen field: the downstream point sees the same signal 1 s later
    assert p[1, 100:] == pytest.approx(p[0, :-100], abs=1e-9)


def test_locate_rock_without_noise():
    st = _street(0.4, phase=3.0)
    sv = lateral_survey(st, 6.0 - 0.25 / 3, np.round(np.arange(-1.8, 1.81, 0.2), 2), noise_pa=0.0)
    fix = locate_rock(sv)
    assert fix.detected and fix.upstream
    assert fix.f == pytest.approx(st.f, rel=0.005)
    assert fix.Uc == pytest.approx(st.Uc, rel=0.02)
    assert fix.diameter == pytest.approx(st.D, rel=0.03)
    assert fix.y == pytest.approx(0.4, abs=0.03)
    assert fix.distance == pytest.approx(6.0, abs=0.25)


def test_no_street_no_detection():
    rng = np.random.default_rng(0)
    ys = np.round(np.arange(-1.8, 1.81, 0.2), 2)
    t = np.arange(len(ys))[:, None] * 60 + np.arange(6000)[None, :] / 100
    sv = Survey(5.9, ys, t, rng.normal(0, 2.0, (len(ys), 10, 6000)),
                fish_sensor_positions() - fish_sensor_positions().mean(), 100.0, 0.25)
    assert not locate_rock(sv).detected
