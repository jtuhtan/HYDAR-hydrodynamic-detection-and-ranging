"""Case study 5 benchmark: how well does the fish locate the rock?

Monte Carlo runs of the lateral survey of example 08 (19 stations 0.2 m apart,
60 s each) with a random lateral rock position (|y| <= 0.8 m), a random
shedding phase and random turbulence and sensor noise for every seed.

* distance : the fish at 3, 6, 9 and 12 m behind the rock, turbulence
             intensity 0, 5, 10 and 15 %. The fish's wake model equals the
             true one, so this measures the effect of noise and of distance.
* mismatch : 6 m behind the rock, 10 % turbulence, while the true wake
             differs from what the fish assumes: Strouhal number, convection
             ratio, eddy viscosity of the cores and roughness of the bed.
             This measures how much the answer rests on the wake model.

Run:  python examples/09_benchmark_fish_rock.py              (25 seeds, about 25 min on 4 cores)
      N_SEEDS=3 python examples/09_benchmark_fish_rock.py    (a quick partial run)

The full run writes docs/data/fish_rock_benchmark.csv; a run with another
N_SEEDS writes fish_rock_benchmark_partial.csv in the current directory.
"""

import csv
import os
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from hydar.karman import (HemisphereRock, KarmanStreet, RiverTurbulence, TrapezoidalChannel,
                          lateral_survey, locate_rock)

N_SEEDS_FULL = 25
N_SEEDS = int(os.environ.get("N_SEEDS", N_SEEDS_FULL))
OUT = Path(__file__).resolve().parents[1] / "docs" / "data" / "fish_rock_benchmark.csv"
OUT_PARTIAL = Path("fish_rock_benchmark_partial.csv")
STATIONS = np.round(np.arange(-1.8, 1.81, 0.2), 2)
FISH_LENGTH = 0.25
DISTANCES = [3.0, 6.0, 9.0, 12.0]
TURBULENCE = [0.0, 0.05, 0.10, 0.15]
# true wake (street keywords, channel keywords) that differs from the fish's assumptions
MISMATCH = {
    "nominal": ({}, {}),
    "St 0.17": ({"strouhal": 0.17}, {}),
    "St 0.23": ({"strouhal": 0.23}, {}),
    "kappa 0.80": ({"convection": 0.80}, {}),
    "kappa 0.95": ({"convection": 0.95}, {}),
    "nu_t x0.5": ({"eddy_viscosity": 0.005}, {}),
    "nu_t x2": ({"eddy_viscosity": 0.02}, {}),
    "Manning n 0.025": ({}, {"manning_n": 0.025}),
    "Manning n 0.050": ({}, {"manning_n": 0.05}),
}
FIELDS = ["case", "variant", "distance_m", "turbulence", "seed", "y_true", "f_true", "Uc_true",
          "D_true", "detected", "f_est", "Uc_est", "D_est", "distance_est", "y_est", "snr",
          "err_y", "err_distance"]


def run(job):
    """One survey and its analysis. All inputs are in `job`, so any start method works."""
    case, variant, dist, ti, seed = job
    street_kw, channel_kw = MISMATCH[variant]
    rng = np.random.default_rng(1000 * seed + 7)
    channel = TrapezoidalChannel(**channel_kw)
    y_rock = rng.uniform(-0.8, 0.8)
    rock = HemisphereRock.from_protrusion(channel, 0.0, y_rock)
    street = KarmanStreet(channel, rock, phase=rng.uniform(0, 20), **street_kw)
    turb = RiverTurbulence(street.U, intensity=ti, scale=channel.depth, seed=seed) if ti > 0 else None
    survey = lateral_survey(street, dist - FISH_LENGTH / 3, STATIONS, turbulence=turb, seed=seed,
                            fish_length=FISH_LENGTH)
    fix = locate_rock(survey)
    return [case, variant, dist, ti, seed, round(y_rock, 4), round(street.f, 5), round(street.Uc, 4),
            round(street.D, 4), int(fix.detected), round(fix.f, 5), round(fix.Uc, 4),
            round(fix.diameter, 4), round(fix.distance, 3), round(fix.y, 3), round(fix.snr, 2),
            round(fix.y - y_rock, 3), round(fix.distance - dist, 3)]


def jobs(n_seeds):
    out = [("distance", "nominal", d, ti, s) for d in DISTANCES for ti in TURBULENCE
           for s in range(n_seeds)]
    out += [("mismatch", v, 6.0, 0.10, s) for v in MISMATCH for s in range(n_seeds)]
    return out


def summary(rows):
    def stats(sel):
        det = [r for r in sel if r["detected"] == 1]
        if not det:
            return f"{0:3d}/{len(sel):<3d} detected"
        ey = np.array([r["err_y"] for r in det])
        ed = np.array([r["err_distance"] for r in det])
        eD = np.array([r["D_est"] / r["D_true"] - 1 for r in det]) * 100
        return (f"{len(det):3d}/{len(sel):<3d} lateral {np.median(np.abs(ey)):5.2f} m (90% {np.percentile(np.abs(ey), 90):4.2f}), "
                f"distance bias {np.mean(ed):+5.2f} SD {np.std(ed):4.2f} median|e| {np.median(np.abs(ed)):4.2f} m, "
                f"D {np.mean(eD):+5.1f} % SD {np.std(eD):4.1f} %")

    print("\nDistance and turbulence (fish's model = true model)")
    for d in DISTANCES:
        for ti in TURBULENCE:
            sel = [r for r in rows if r["case"] == "distance" and r["distance_m"] == d and r["turbulence"] == ti]
            if sel:
                print(f"  {d:4.0f} m  Tu {ti:4.2f}: {stats(sel)}")
    print("\nModel mismatch (6 m, Tu 0.10)")
    for v in MISMATCH:
        sel = [r for r in rows if r["case"] == "mismatch" and r["variant"] == v]
        if sel:
            print(f"  {v:16s}: {stats(sel)}")


def main():
    todo = jobs(N_SEEDS)
    with Pool(int(os.environ.get("N_PROC", os.cpu_count() or 1))) as pool:
        rows = pool.map(run, todo, chunksize=1)
    out = OUT if N_SEEDS == N_SEEDS_FULL else OUT_PARTIAL
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(FIELDS)
        w.writerows(rows)
    print(f"wrote {len(rows)} runs to {out}")
    summary([dict(zip(FIELDS, r)) for r in rows])


if __name__ == "__main__":
    main()
