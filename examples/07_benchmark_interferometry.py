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
same wake. Wind is added linearly to the recorded wake (linear theory). Every
geometry and method sees the same seeds, so runs can be compared in pairs.

Run:  python examples/07_benchmark_interferometry.py            (100 seeds, about 2 h on 4 cores)
      N_SEEDS=5 python examples/07_benchmark_interferometry.py ferry      (a quick partial run)

The full run writes docs/data/ci_benchmark.csv. A partial run (one case only,
or another N_SEEDS) writes ci_benchmark_partial.csv in the current directory
instead, so it cannot replace the full data set. Both print the summary and
the paired comparisons.
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
N_SEEDS_FULL = 100
N_SEEDS = int(os.environ.get("N_SEEDS", N_SEEDS_FULL))
OUT = Path(__file__).resolve().parents[1] / "docs" / "data" / "ci_benchmark.csv"
OUT_PARTIAL = Path("ci_benchmark_partial.csv")
NOISE_PA = 2.0                # white sensor noise [Pa]

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

STATE = {}                    # wakes and solver records, filled in the parent process


def _init_worker(state):
    """Give each worker process the shared records (works with fork, spawn and forkserver)."""
    STATE.update(state)


def _attempt(fn):
    """Run one method. Returns (ok, reason, result).

    A RuntimeError is the method saying it found no fix. Any other exception is
    a defect; it is recorded as such instead of being counted silently.
    """
    try:
        return 1, "", fn()
    except RuntimeError as exc:
        return 0, f"no fix: {exc}", (np.nan, np.nan, np.nan)
    except Exception as exc:  # noqa: BLE001  (keep the run alive, but name the error)
        return 0, f"error: {type(exc).__name__}: {exc}", (np.nan, np.nan, np.nan)


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
    arr.load(tf, CubicSpline(F_TM, eta_m, axis=1)(tf), depth=F_DEPTH)
    arr.noise_pa, arr.noise_seed = NOISE_PA, seed + 1000

    def kelvin():
        f = locate_vessel(arr)
        return score(f.speed, f.course_deg, f.position, 11.0, course, SHIP.position, F_TT)

    def interferometry():
        f = track_vessel([arr], F_BAND, depth=np.inf, window=120.0, hop=20.0)
        return score(f.speed, f.heading_deg, f.position, 11.0, course, SHIP.position, F_TT)

    rows = []
    for method, fn in (("kelvin", kelvin), ("interferometry", interferometry)):
        ok, reason, r = _attempt(fn)
        rows.append(("ferry", geo, method, hs, seed, ok, *r, reason))
    return rows


def run_boat(job):
    geo, hs, seed = job
    layout, recs, t = STATE["boat"][geo]
    arrays = []
    for i, ((name, x, y, z), eta0) in enumerate(zip(B_SPEC, recs)):
        a = SensorArray(name, x, y, z=z, layout=layout, noise_pa=NOISE_PA, noise_seed=1000 * (seed + 1) + i)
        eta = eta0.copy()
        if hs > 0:
            eta += wind_sea_head(a.positions, t, B_DEPTH, z, hs, B_FPEAK, mean_dir_deg=112.5, seed=seed)
        a.load(t, eta, depth=B_DEPTH)
        arrays.append(a)

    def kelvin():
        return track_boat(arrays)

    def interferometry():
        return track_vessel(arrays, B_BAND, depth=B_DEPTH, window=6.0, hop=0.5, time_stamp="centroid")

    rows = []
    for method, fn in (("kelvin", kelvin), ("interferometry", interferometry)):
        def scored(fn=fn):
            f = fn()
            return score(f.speed, f.heading_deg, f.position, BOAT["speed"], BOAT["heading_deg"], boat_true, B_TT)
        ok, reason, r = _attempt(scored)
        rows.append(("boat", geo, method, hs, seed, ok, *r, reason))
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
        ok = np.array([r[6:9] for r in R if r[5]], float)
        unit = "m"
        if len(ok) == 0:
            print(f"| {case} | {geo} | {method} | {hs} | 0/{len(R)} | | | | | | |")
            continue
        fmt = "{:.0f}" if case == "ferry" else "{:.2f}"
        print(f"| {case} | {geo} | {method} | {hs} | {len(ok)}/{len(R)} | {ok[:, 0].mean():+.3f} | "
              f"{ok[:, 0].std():.3f} | {ok[:, 1].mean():+.2f} | {ok[:, 1].std():.2f} | "
              f"{fmt.format(np.median(ok[:, 2]))} {unit} | {fmt.format(np.percentile(ok[:, 2], 90))} {unit} |")


def paired(rows, case, method, geo_a, geo_b):
    """Seed-by-seed comparison of the track error of two geometries (same method)."""
    from scipy.stats import binomtest, wilcoxon
    err = {(r[1], r[3], r[4]): r[8] for r in rows if r[0] == case and r[2] == method and r[5]}
    print(f"\n{case}, {method}: '{geo_a}' against '{geo_b}' (paired by seed, track RMS)")
    print("| Hs (m) | pairs | first is better | sign test p | median difference (95 % CI) | Wilcoxon p |")
    print("|---|---|---|---|---|---|")
    rng = np.random.default_rng(0)
    pooled = []
    for hs in sorted({k[1] for k in err}):
        seeds = sorted({k[2] for k in err if k[:2] == (geo_a, hs)} & {k[2] for k in err if k[:2] == (geo_b, hs)})
        if len(seeds) < 6:
            print(f"| {hs} | {len(seeds)} | | | | |")
            continue
        d = np.array([err[(geo_a, hs, s)] - err[(geo_b, hs, s)] for s in seeds])
        wins = int((d < 0).sum())
        boot = np.median(rng.choice(d, (2000, len(d))), axis=1)
        lo, hi = np.percentile(boot, [2.5, 97.5])
        unit = "{:+.0f} ({:+.0f} to {:+.0f}) m" if case == "ferry" else "{:+.2f} ({:+.2f} to {:+.2f}) m"
        print(f"| {hs} | {len(d)} | {wins} | {binomtest(wins, len(d)).pvalue:.3f} | "
              f"{unit.format(np.median(d), lo, hi)} | {wilcoxon(d).pvalue if d.any() else 1.0:.3f} |")
        if hs > 0:
            pooled.append(d)
    if pooled:
        d = np.concatenate(pooled)
        wins = int((d < 0).sum())
        print(f"| all Hs > 0 | {len(d)} | {wins} | {binomtest(wins, len(d)).pvalue:.3f} | | {wilcoxon(d).pvalue:.3f} |")


def failures(rows):
    """Count the recorded reasons for runs without a fix; defects are listed separately."""
    from collections import Counter
    c = Counter()
    for r in rows:
        if not r[5]:
            kind = "no fix" if r[9].startswith("no fix") else r[9]
            c[(r[0], r[2], kind)] += 1
    print("\nRuns without a fix:" + ("  none" if not c else ""))
    for (case, method, kind), n in sorted(c.items()):
        print(f"  {case:5s} {method:14s} {n:4d}  {kind}")
    if any(k[2] != "no fix" for k in c):
        print("  WARNING: rows marked 'error' are defects, not detection failures.")


if __name__ == "__main__":
    which = sys.argv[1:] or ["ferry", "boat"]
    full = set(which) == {"ferry", "boat"} and N_SEEDS == N_SEEDS_FULL
    out = OUT if full else OUT_PARTIAL
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
    with Pool(initializer=_init_worker, initargs=(STATE,)) as pool:
        res = pool.map(_dispatch, jobs, chunksize=1)
    rows = [r for rr in res for r in rr]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "geometry", "method", "hs_m", "seed", "ok", "speed_err_mps",
                    "course_err_deg", "track_rms_m", "reason"])
        w.writerows(rows)
    print(f"wrote {out}" + ("" if full else "  (partial run; docs/data/ci_benchmark.csv is unchanged)"))
    summarize(rows)
    failures(rows)
    if "ferry" in which:
        paired(rows, "ferry", "interferometry", "optimal, 20 m footprint", "Ratsep 2.5 m grid")
        paired(rows, "ferry", "interferometry", "optimal, same footprint", "Ratsep 2.5 m grid")
    if "boat" in which:
        paired(rows, "boat", "interferometry", "optimal, same footprint", "1 m grid")
