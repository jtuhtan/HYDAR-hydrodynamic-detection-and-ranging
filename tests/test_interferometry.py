"""Case study 4: correlative interferometry, array design and track fit."""

import numpy as np
import pytest

from hydar import SensorArray
from hydar.interferometry import (Observations, ci_scan, design_array, dispersion, evaluate_layout,
                                  fit_track, grid_layout, peak_sidelobe)

G = 9.81


def test_dispersion_limits():
    k, c, cg = dispersion(0.2)                       # deep water
    assert c[0] == pytest.approx(G / (2 * np.pi * 0.2))
    assert cg[0] == pytest.approx(c[0] / 2)
    k, c, cg = dispersion(0.01, depth=2.0)           # shallow-water limit
    assert c[0] == pytest.approx(np.sqrt(G * 2.0), rel=1e-3)
    assert cg[0] == pytest.approx(c[0], rel=1e-3)


def test_ci_scan_finds_plane_wave_direction():
    a = SensorArray("H", 0.0, 0.0, spacing=2.5, fs=10.0)
    t = np.arange(0, 240, 0.1)
    th0 = np.radians(-50.0)
    f0 = 0.15
    k = dispersion(f0)[0][0]
    eta = 0.01 * np.cos(k * (a.positions @ [np.cos(th0), np.sin(th0)])[:, None] - 2 * np.pi * f0 * t)
    eta += np.random.default_rng(0).normal(0, 0.001, eta.shape)
    a.load(t, eta, depth=3.0)
    T, f, TH, R, P = ci_scan(a, (0.1, 0.2), window=120, hop=60)
    i = np.argmax(P[0])
    err = (np.degrees(TH[0, i]) - np.degrees(th0) + 180) % 360 - 180
    assert abs(err) < 1.0
    assert R[0, i] > 0.99


def test_regular_grid_aliases_and_design_avoids_it():
    # a 1 m grid sees waves shorter than 2 m: grating lobes (sidelobe ~ 1)
    base = grid_layout([-1, 0, 1], [-1, 0, 1])
    k_short = dispersion(1.2, depth=3.0)[0][0]
    assert peak_sidelobe(base, k_short) > 0.95
    d = design_array((0.45, 1.2), depth=3.0, max_aperture=2.83, min_spacing=0.25, maxiter=40)
    assert d.aperture <= 2.83 * 1.01
    assert d.peak_sidelobe < 0.7
    crlb_base, _ = evaluate_layout(base, (0.45, 1.2), 3.0)
    assert d.mean_crlb_deg < crlb_base


def _synthetic_observations(V, P0, r, n_obs=60, seed=0):
    """Observations that obey equations (1) and (2) exactly (deep water)."""
    rng = np.random.default_rng(seed)
    U = np.linalg.norm(V)
    e = V / U
    rows = []
    while len(rows) < n_obs:
        th = rng.uniform(-np.radians(80), np.radians(80))    # wave angle to the course
        n = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]]) @ e
        c = U * np.cos(th)
        f = G / (2 * np.pi * c)
        cg = c / 2
        te, tau = np.linalg.solve(np.column_stack([V, cg * n]), r - P0)
        if tau <= 0:
            continue
        rows.append((te + tau, f, np.arctan2(n[1], n[0]), c, cg))
    t, f, theta, c, cg = map(np.array, zip(*rows))
    N = len(t)
    return Observations(t, f, theta, np.ones(N), np.ones(N), c, cg, np.tile(r, (N, 1)), np.zeros(N, int))


def test_fit_track_recovers_velocity_and_position_with_outliers():
    V = 11.0 * np.array([np.cos(np.radians(-69)), np.sin(np.radians(-69))])
    P0 = np.array([-3500.0, 2200.0])
    obs = _synthetic_observations(V, P0, np.array([0.0, 0.0]))
    # add 30 % wind-sea-like outliers: one direction, random frequencies, all times
    rng = np.random.default_rng(1)
    m = int(0.3 * len(obs))
    fw = rng.uniform(0.25, 0.6, m)
    cw = G / (2 * np.pi * fw)
    out = Observations(rng.uniform(obs.t.min(), obs.t.max(), m), fw,
                       np.radians(112.5 + 20 * rng.standard_normal(m)), np.ones(m), np.ones(m),
                       cw, cw / 2, np.zeros((m, 2)), np.zeros(m, int))
    tr = fit_track(Observations.concat([obs, out]), tol_t=60.0)
    assert tr.speed == pytest.approx(11.0, rel=1e-3)
    assert tr.heading_deg == pytest.approx(-69.0, abs=0.1)
    assert np.linalg.norm(tr.p0 - P0) < 5.0
    assert tr.n_obs == len(obs)


@pytest.mark.slow
def test_track_ferry_without_kelvin_theory():
    from hydar.interferometry import track_vessel
    from hydar.kelvin import KelvinShip, record_passage
    arr = SensorArray("H", 0.0, 0.0, z=-2.8, spacing=2.5, rotation_deg=22.5, pressure_model="airy")
    course = -69.0
    e = np.array([np.cos(np.radians(course)), np.sin(np.radians(course))])
    n = np.array([-e[1], e[0]])
    p0 = -2500.0 * n - 11.0 * 300.0 * e
    ship = KelvinShip(speed=11.0, course_deg=course, x0=p0[0], y0=p0[1])
    record_passage([arr], ship, t_end=1800.0, depth=3.0,
                   wind={"hs": 0.15, "f_peak": 0.4, "mean_dir_deg": 112.5})
    tr = track_vessel([arr], band=(0.07, 0.3), depth=np.inf, window=120.0, hop=20.0)
    assert tr.speed == pytest.approx(11.0, rel=0.01)
    assert abs((tr.heading_deg - course + 180) % 360 - 180) < 1.0
    tt = np.linspace(0, 1800, 61)
    assert np.sqrt(np.mean(np.sum((tr.position(tt) - ship.position(tt)) ** 2, axis=1))) < 300
