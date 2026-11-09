"""Case study 4 benchmark: correlative interferometry vs the Kelvin methods.

Monte Carlo comparison of accuracy (bias), precision (spread) and robustness
to wind noise, for both earlier vessel case studies:

* Ferry (case study 3): one seabed array 2.5 km from an 11 m/s ship.
  Methods: locate_vessel (Kelvin chirp + cusp angle) and track_vessel
  (interferometry). Geometries: the Rätsep et al. grid (2.5 m spacing), the
  optimal grid with the same footprint, and the optimal grid in 20 m.
* Boat (case study 2): four arrays 8 to 18 m from a 3 m/s, 4.5 m Series 60
  boat in the dispersive solver. Methods: track_boat (Kelvin cusp line) and
  track_vessel. Geometries: the 1 m grid and the optimal grid with the same
  footprint.

For each wind-sea height Hs, N_SEEDS realizations of the wind sea (random
phases, frequencies and directions) and of the sensor noise are added to the
same wake. Wind is added linearly to the recorded wake (linear theory).

Run:  python examples/07_benchmark_interferometry.py   (about 30 min on 4 cores)
Writes docs/data/ci_benchmark.csv and prints the summary tables.
"""

import csv
import os
import sys
import warnings
from multiprocessing import Pool
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline

from hydar import Boat, Grid, SensorArray, SpectralSolver, track_boat
from hydar.boat import Series60Hull
from hydar.interferometry import design_array, track_vessel, wind_sea_head
from hydar.kelvin import KelvinShip, locate_vessel, wind_sea

warnings.filterwarnings("ignore")
G = 9.81
N_SEEDS = int(os.environ.get("N_SEEDS", 20))
OUT = Path(__file__).resolve().parents[1] / "docs" / "data" / "ci_benchmark.csv"

# ------------------------------------------------------------------ ferry (case 3)
F_DEPTH, F_BAND, F_TEND, F_FSM = 3.0, (0.07, 0.3), 1800.0, 4.0
F_HS = [0.0, 0.15, 0.3, 0.45, 0.6]
course = 90.0 - 159.0
_e = np.array([np.cos(np.radians(course)), np.sin(np.radians(course))])
_p0 = -2500.0 * np.array([-_e[1], _e[0]]) - 11.0 * 300.0 * _e
SHIP = KelvinShip(speed=11.0, course_deg=course, x0=_p0[0], y0=_p0[1], amplitude=0.25)
F_TM = np.arange(0.0, F_TEND, 1 / F_FSM)
F_TT = np.linspace(0, F_TEND, 61)

# ------------------------------------------------------------------ boat (case 2)
B_DEPTH, B_BAND = 3.0, (0.45, 1.2)
B_HS, B_FPEAK = [0.0, 0.025, 0.05, 0.1, 0.15, 0.2], 0.8
B_SPEC = [("A", 50.0, 14.0, -0.3), ("B", 75.0, 16.0, -0.5), ("C", 55.0, 36.0, -0.5), ("D", 80.0, 40.0, -0.8)]
B_TT = np.linspace(15, 44, 30)
BOAT = dict(speed=3.0, x=-8.0, y=22.0, heading_deg=4.0)

STATE = {}


def ferry_array(layout):
    kw = {"layout": layout} if layout is not None else {"spacing": 2.5}
    return SensorArray("H", 0.0, 0.0, z=-(F_DEPTH - 0.2), rotation_deg=22.5, pressure_model="airy", **kw)


def boat_true(t):
    e = np.array([np.cos(np.radians(BOAT["heading_deg"])), np.sin(np.radians(BOAT["heading_deg"]))])
    return np.array([BOAT["x"], BOAT["y"]]) + np.multiply.outer(np.asarray(t, float) * BOAT["speed"], e)


def score(speed, heading, pos, true_speed, true_heading, true_pos, tt):
    err = np.linalg.norm(pos(tt) - true_pos(tt), axis=1)
    return (speed - true_speed, (heading - true_heading + 180) % 360 - 180,
            float(np.sqrt(np.mean(err ** 2))))


def run_ferry(job):
    geo, hs, seed = job
    layout, wake = STATE["ferry"][geo]
    arr = ferry_array(layout)
    eta_m = wake.copy()
    if hs > 0:
        eta_m += wind_sea(arr.positions, F_TM, hs=hs, f_peak=0.4, mean_dir_deg=112.5, seed=seed)
    tf = np.arange(0.0, F_TM[-1], 1 / arr.fs)
    eta = CubicSpline(F_TM, eta_m, axis=1)(tf)
    eta += np.random.default_rng(seed + 1000).normal(0, 2.0 / (1000 * G), eta.shape)
    arr.load(tf, eta, depth=F_DEPTH)
    rows = []
    for method in ("kelvin", "interferometry"):
        try:
            if method == "kelvin":
                f = locate_vessel(arr)
                r = score(f.speed, f.course_deg, f.position, 11.0, course, SHIP.position, F_TT)
            else:
                f = track_vessel([arr], F_BAND, depth=np.inf, window=120.0, hop=20.0)
                r = score(f.speed, f.heading_deg, f.position, 11.0, course, SHIP.position, F_TT)
            rows.append(("ferry", geo, method, hs, seed, 1, *r))
        except Exception:
            rows.append(("ferry", geo, method, hs, seed, 0, np.nan, np.nan, np.nan))
    return rows


def run_boat(job):
    geo, hs, seed = job
    layout, recs, t = STATE["boat"][geo]
    rng = np.random.default_rng(seed + 1000)
    arrays = []
    for (name, x, y, z), eta0 in zip(B_SPEC, recs):
        a = SensorArray(name, x, y, z=z, layout=layout)
        eta = eta0.copy()
        if hs > 0:
            eta += wind_sea_head(a.positions, t, B_DEPTH, z, hs, B_FPEAK, mean_dir_deg=112.5, seed=seed)
        eta += rng.normal(0, 2.0 / (1000 * G), eta.shape)
        a.load(t, eta, depth=B_DEPTH)
        arrays.append(a)
    rows = []
    for method in ("kelvin", "interferometry"):
        try:
            if method == "kelvin":
                f = track_boat(arrays)
            else:
                f = track_vessel(arrays, B_BAND, depth=B_DEPTH, window=6.0, hop=0.5, time_stamp="centroid")
            r = score(f.speed, f.heading_deg, f.position, BOAT["speed"], BOAT["heading_deg"], boat_true, B_TT)
            rows.append(("boat", geo, method, hs, seed, 1, *r))
        except Exception:
            rows.append(("boat", geo, method, hs, seed, 0, np.nan, np.nan, np.nan))
    return rows


def boat_records(layouts):
    """One dispersive-solver run with every geometry's arrays (sensors are passive)."""
    grid = Grid(lx=120.0, ly=50.0, dx=0.125, n_pml=80)
    sim = SpectralSolver(grid, depth=B_DEPTH, dt=0.01)
    arrays = {g: [SensorArray(n, x, y, z=z, pressure_model="solver", layout=lay) for n, x, y, z in B_SPEC]
              for g, lay in layouts.items()}
    for arrs in arrays.values():
        for a in arrs:
            sim.add(a)
    sim.add(Boat(hull=Series60Hull(length=4.5), ramp=4.0, **BOAT))
    sim.run(44.0)
    return {g: (layouts[g], [a.eta for a in arrs], arrs[0].time) for g, arrs in arrays.items()}


def _dispatch(job):
    fn, args = job
    return fn(args)


def summarize(rows):
    print("\n| Case | Geometry | Method | Hs (m) | fixes | speed bias (m/s) | speed SD | "
          "course bias (deg) | course SD | track RMS median | track RMS 90th pct |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    keys = sorted({(r[0], r[1], r[2], r[3]) for r in rows}, key=lambda k: (k[0], k[1], k[3], k[2]))
    for case, geo, method, hs in keys:
        R = [r for r in rows if r[:4] == (case, geo, method, hs)]
        ok = np.array([r[6:] for r in R if r[5]], float)
        unit = "m"
        if len(ok) == 0:
            print(f"| {case} | {geo} | {method} | {hs} | 0/{len(R)} | | | | | | |")
            continue
        fmt = "{:.0f}" if case == "ferry" else "{:.2f}"
        print(f"| {case} | {geo} | {method} | {hs} | {len(ok)}/{len(R)} | {ok[:, 0].mean():+.3f} | "
              f"{ok[:, 0].std():.3f} | {ok[:, 1].mean():+.2f} | {ok[:, 1].std():.2f} | "
              f"{fmt.format(np.median(ok[:, 2]))} {unit} | {fmt.format(np.percentile(ok[:, 2], 90))} {unit} |")


if __name__ == "__main__":
    which = sys.argv[1:] or ["ferry", "boat"]
    print("designing arrays ...")
    STATE["ferry"], STATE["boat"] = {}, {}
    if "ferry" in which:
        f_lay = {"Ratsep 2.5 m grid": None,
                 "optimal, same footprint": design_array(F_BAND, max_aperture=7.08, min_spacing=0.5).offsets,
                 "optimal, 20 m footprint": design_array(F_BAND, max_aperture=20.0, min_spacing=0.5).offsets}
        print("synthesizing the ferry wake for each geometry ...")
        for g, lay in f_lay.items():
            STATE["ferry"][g] = (lay, SHIP.wake(ferry_array(lay).positions, F_TM))
    if "boat" in which:
        b_lay = {"1 m grid": None,
                 "optimal, same footprint": design_array(B_BAND, depth=B_DEPTH, max_aperture=2.83,
                                                         min_spacing=0.25).offsets}
        print("running the dispersive solver for the boat (about 6 min) ...")
        STATE["boat"] = boat_records(b_lay)
    jobs = []
    if "ferry" in which:
        jobs += [(run_ferry, (g, hs, s)) for g in STATE["ferry"] for hs in F_HS for s in range(N_SEEDS)]
    if "boat" in which:
        jobs += [(run_boat, (g, hs, s)) for g in STATE["boat"] for hs in B_HS for s in range(N_SEEDS)]
    print(f"{len(jobs)} Monte Carlo runs ...")
    with Pool() as pool:
        res = pool.map(_dispatch, jobs, chunksize=1)
    rows = [r for rr in res for r in rr]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "geometry", "method", "hs_m", "seed", "ok", "speed_err_mps",
                    "course_err_deg", "track_rms_m"])
        w.writerows(rows)
    print(f"wrote {OUT}")
    summarize(rows)
