"""Time the hot paths of HYDAR (see repository_review/2026-11-16-performance-review.md).

Run:  python examples/10_speed.py                    (with Numba, if installed)
      HYDAR_NO_NUMBA=1 python examples/10_speed.py   (NumPy only)

Each line is the median of a few repeats after one warm-up call, so Numba's
compile time (cached on disk after the first run) is not counted.
"""

import time

import numpy as np

from hydar import Boat, Grid, Solver, SpectralSolver
from hydar.boat import Series60Hull
from hydar.interferometry import peak_sidelobe, grid_layout
from hydar.kelvin import KelvinShip


def timed(fn, repeat=3):
    fn()
    ts = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        ts.append(time.perf_counter() - t0)
    return float(np.median(ts))


def shallow_water():
    sim = Solver(Grid(lx=4.0, ly=4.0, dx=0.01, n_pml=30), depth=0.5)
    X, Y = np.meshgrid(sim.grid.x, sim.grid.y, indexing="ij")
    sim.set_eta(0.01 * np.exp(-((X - 2) ** 2 + (Y - 1.7) ** 2) / 0.02))
    return lambda: [sim.step() for _ in range(50)], 50, sim.ex.shape


def spectral_boat(precision="double"):
    kw = {} if precision == "double" else {"precision": precision}
    sim = SpectralSolver(Grid(lx=120.0, ly=50.0, dx=0.125, n_pml=80), depth=3.0, dt=0.01, **kw)
    sim.add(Boat(hull=Series60Hull(length=4.5), ramp=4.0, speed=3.0, x=-8.0, y=22.0, heading_deg=4.0))
    return lambda: [sim.step() for _ in range(10)], 10, sim.eta.shape


def fish_job():
    import importlib.util
    from pathlib import Path
    path = Path(__file__).with_name("09_benchmark_fish_rock.py")
    spec = importlib.util.spec_from_file_location("benchmark_fish_rock", path)
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    return lambda: bench.run(("distance", "nominal", 6.0, 0.10, 0)), 1, "survey + locate_rock"


def sidelobe():
    off = grid_layout([0, 1.1, 2.5], [0, 0.4, 2.0])
    return lambda: [peak_sidelobe(off, 1.3, n_src=12) for _ in range(20)], 20, "3 x 3 grid"


def kelvin():
    ship = KelvinShip()
    t = np.arange(0.0, 1800.0, 2.0)
    return lambda: ship.wake(np.array([[0.0, 0.0]]), t), 1, f"1 point, {t.size} times"


import hydar
try:
    from hydar._kernels import HAVE_NUMBA
except ImportError:            # v1.3.0 and older have no compiled kernels
    HAVE_NUMBA = False
print(f"hydar {hydar.__version__} from {hydar.__file__}; Numba {'on' if HAVE_NUMBA else 'off'}")
for name, make, unit in (("shallow-water Solver", shallow_water, "step"),
                         ("SpectralSolver with boat", spectral_boat, "step"),
                         ("  same, single precision", lambda: spectral_boat("single"), "step"),
                         ("benchmark 09, one run", fish_job, "run"),
                         ("peak_sidelobe", sidelobe, "call"),
                         ("KelvinShip.wake", kelvin, "call")):
    fn, n, size = make()
    print(f"{name:26s} {1e3 * timed(fn) / n:9.2f} ms per {unit}   ({size})")
