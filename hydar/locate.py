"""HYDAR: locate a wave source by time of flight, like radar or sonar.

Each sensor array works like a small sonar receiver:

1. Arrival picking. For every sensor the arrival time of the wave front is
   the first moment the dynamic pressure exceeds a fraction of the array's
   peak value (with linear interpolation between samples).

2. Direction (bearing) from time differences of arrival (TDOA). Across a
   small array the front is close to a plane wave, so the arrival times obey

       t_i = tau + s . (r_i - r_c)

   where r_c is the array centre, tau the arrival time at the centre and s
   the slowness vector (|s| = 1/c). A linear least-squares fit gives s; the
   source lies in the direction -s / |s| from the array centre.

3. Range from time of flight. When the emission (drop) time t0 is known,
   range = c * (tau - t0), exactly as in radar ranging, with
   c = sqrt(g h) the shallow-water wave speed. With the bearing this gives a
   position estimate from a single array.

4. Fusion. With several arrays, all arrival times are fitted together with a
   circular wave front, t_i = t0 + |r_i - x_s| / c, solving for the source
   position x_s and, if unknown, the emission time t0. Because t0 is fitted,
   any constant picking delay cancels out.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import least_squares


def pick_arrivals(t: np.ndarray, p: np.ndarray, fraction: float = 0.1,
                  t_start: float = 0.0, t_end: float | None = None) -> np.ndarray:
    """First-arrival time at each sensor [s].

    t : (n_samples,) sample times; p : (n_sensors, n_samples) dynamic pressure.
    A common threshold, `fraction` of the largest |p| on the array inside the
    window [t_start, t_end], is used for all sensors so their picks match.
    Returns NaN for sensors that never cross the threshold.
    """
    win = (t >= t_start) & (t <= (t[-1] if t_end is None else t_end))
    tw, pw = t[win], np.abs(p[:, win])
    thr = fraction * pw.max()
    out = np.full(p.shape[0], np.nan)
    if thr <= 0:
        return out
    for k in range(p.shape[0]):
        above = np.nonzero(pw[k] >= thr)[0]
        if above.size == 0:
            continue
        n = above[0]
        if n == 0:
            out[k] = tw[0]
        else:
            a, b = pw[k, n - 1], pw[k, n]
            out[k] = tw[n - 1] + (thr - a) / (b - a) * (tw[n] - tw[n - 1])
    return out


@dataclass
class ArrayFix:
    """Result for one sensor array."""

    name: str
    center: np.ndarray
    arrivals: np.ndarray        # per-sensor arrival times [s]
    tau: float                  # arrival at array centre [s]
    bearing_deg: float          # direction to source, degrees counter-clockwise from +x
    direction: np.ndarray       # unit vector towards the source
    apparent_speed: float       # 1/|s| from the TDOA fit [m/s]
    range: float | None         # c * (tau - t0) if t0 known [m]
    position: np.ndarray | None # centre + range * direction [m]

    def __str__(self) -> str:
        s = f"{self.name}: bearing {self.bearing_deg:6.1f} deg"
        if self.range is not None:
            s += f", range {self.range:6.2f} m, position ({self.position[0]:.2f}, {self.position[1]:.2f})"
        return s + f"  [apparent speed {self.apparent_speed:.2f} m/s]"


def locate_array(array, c: float, t0: float | None = None, fraction: float = 0.1,
                 t_start: float = 0.0, t_end: float | None = None) -> ArrayFix:
    """Bearing (TDOA) and, if t0 is given, range (time of flight) for one array."""
    t, p = array.time, array.dynamic_pressure
    arr = pick_arrivals(t, p, fraction, t_start, t_end)
    ok = np.isfinite(arr)
    if ok.sum() < 3:
        raise RuntimeError(f"{array.name}: fewer than 3 arrivals picked")
    A = np.column_stack([np.ones(ok.sum()), array.offsets[ok]])
    tau, sx, sy = np.linalg.lstsq(A, arr[ok], rcond=None)[0]
    s = np.array([sx, sy])
    d = -s / np.linalg.norm(s)
    bearing = float(np.degrees(np.arctan2(d[1], d[0])))
    rng = pos = None
    if t0 is not None:
        rng = float(c * (tau - t0))
        pos = array.center + rng * d
    return ArrayFix(array.name, array.center, arr, float(tau), bearing, d,
                    float(1.0 / np.linalg.norm(s)), rng, pos)


def intersect_bearings(fixes) -> np.ndarray:
    """Least-squares intersection of the bearing lines of several arrays."""
    A, b = np.zeros((2, 2)), np.zeros(2)
    for f in fixes:
        P = np.eye(2) - np.outer(f.direction, f.direction)  # projector normal to line
        A += P
        b += P @ f.center
    return np.linalg.solve(A, b)


def locate(arrays, c: float, t0: float | None = None, fraction: float = 0.1,
           t_start: float = 0.0, t_end: float | None = None) -> dict:
    """Full HYDAR estimate from one or more arrays.

    Returns a dict with
      'fixes'     : list of ArrayFix (bearing and TOF range per array)
      'bearings'  : intersection of the bearing lines (needs >= 2 arrays)
      'position'  : fused circular-front fit over all sensors
      't0'        : emission time used or estimated
    """
    fixes = [locate_array(a, c, t0, fraction, t_start, t_end) for a in arrays]
    out = {"fixes": fixes, "bearings": None}
    if len(fixes) >= 2:
        out["bearings"] = intersect_bearings(fixes)

    pos = np.vstack([a.positions for a in arrays])
    arr = np.concatenate([f.arrivals for f in fixes])
    ok = np.isfinite(arr)
    pos, arr = pos[ok], arr[ok]

    # initial guess: bearing intersection, else first array's TOF position
    if out["bearings"] is not None:
        x_init = out["bearings"]
    elif fixes[0].position is not None:
        x_init = fixes[0].position
    else:
        x_init = fixes[0].center + 5.0 * fixes[0].direction

    if t0 is None:
        t_init = arr.min() - np.linalg.norm(pos - x_init, axis=1).min() / c
        res = least_squares(lambda q: q[2] + np.linalg.norm(pos - q[:2], axis=1) / c - arr,
                            np.r_[x_init, t_init])
        out["position"], out["t0"] = res.x[:2], float(res.x[2])
    else:
        res = least_squares(lambda q: t0 + np.linalg.norm(pos - q, axis=1) / c - arr, x_init)
        out["position"], out["t0"] = res.x, t0
    out["residual_rms"] = float(np.sqrt(np.mean(res.fun ** 2)))
    return out
