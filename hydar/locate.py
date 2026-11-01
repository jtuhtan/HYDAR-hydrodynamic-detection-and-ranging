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


# ---------------------------------------------------------------------------
# Moving source: track a supercritical boat from its V-shaped wake
# ---------------------------------------------------------------------------
@dataclass
class BoatTrack:
    """Constant-velocity boat track estimated from wake arrivals."""

    p0: np.ndarray          # boat position at t = 0 [m]
    heading_deg: float      # direction of travel, degrees counter-clockwise from +x
    speed: float            # [m/s]
    mach: float             # speed / wave speed (depth Froude number)
    residual_rms: float     # [s]
    arrivals: list          # per-array arrival times [s]

    @property
    def direction(self) -> np.ndarray:
        th = np.radians(self.heading_deg)
        return np.array([np.cos(th), np.sin(th)])

    def position(self, t) -> np.ndarray:
        """Estimated boat position at time(s) t."""
        return self.p0 + np.multiply.outer(np.asarray(t, float) * self.speed, self.direction)

    @property
    def wake_half_angle_deg(self) -> float:
        return float(np.degrees(np.arcsin(1.0 / self.mach)))

    def __str__(self) -> str:
        return (f"speed {self.speed:.2f} m/s (Fr = {self.mach:.2f}), heading {self.heading_deg:.1f} deg, "
                f"position at t=0 ({self.p0[0]:.2f}, {self.p0[1]:.2f}) m, "
                f"wake half-angle {self.wake_half_angle_deg:.1f} deg, rms residual {1e3 * self.residual_rms:.0f} ms")


def wake_arrival(r: np.ndarray, p0, heading_deg: float, speed: float, c: float) -> np.ndarray:
    """Time the Mach-cone wake of a boat on a straight track first reaches points r.

    A boat at p(t) = p0 + U t e with U > c trails a cone of half-angle
    mu = arcsin(c/U). A point a distance s along the track and d to the side
    is reached at t = (s + d * sqrt(M**2 - 1)) / U, with M = U/c.
    """
    th = np.radians(heading_deg)
    e = np.array([np.cos(th), np.sin(th)])
    rel = np.atleast_2d(r) - np.asarray(p0)
    s = rel @ e
    d = np.abs(rel[:, 0] * e[1] - rel[:, 1] * e[0])
    M = speed / c
    return (s + d * np.sqrt(max(M * M - 1.0, 0.0))) / speed


def track_boat(arrays, c: float, fraction: float = 0.2, t_start: float = 0.0,
               t_end: float | None = None) -> BoatTrack:
    """HYDAR for a moving boat: fit a straight, constant-speed track.

    Picks the first wake arrival at every sensor of every array, then fits
    the boat's position at t = 0, heading and speed so that the predicted
    Mach-cone arrival times (wake_arrival) match the picks in a least-squares
    sense. Needs a supercritical boat (speed > c) seen by at least two arrays.
    """
    picks = [pick_arrivals(a.time, a.dynamic_pressure, fraction, t_start, t_end) for a in arrays]
    pos = np.vstack([a.positions for a in arrays])
    arr = np.concatenate(picks)
    ok = np.isfinite(arr)
    pos, arr = pos[ok], arr[ok]
    if arr.size < 6:
        raise RuntimeError("too few wake arrivals picked")

    def resid(q):
        x0, y0, th, U = q
        return wake_arrival(pos, (x0, y0), th, U, c) - arr

    centre, t_mean = pos.mean(axis=0), arr.mean()
    best = None
    for th in np.arange(0, 360, 30.0):
        for U in (1.2 * c, 1.6 * c, 2.2 * c):
            e = np.array([np.cos(np.radians(th)), np.sin(np.radians(th))])
            p_init = centre - e * U * t_mean
            res = least_squares(resid, [*p_init, th, U],
                                bounds=([-np.inf, -np.inf, -np.inf, 1.0001 * c],
                                        [np.inf, np.inf, np.inf, 10 * c]))
            if best is None or res.cost < best.cost:
                best = res
    x0, y0, th, U = best.x
    return BoatTrack(np.array([x0, y0]), float(th % 360), float(U), float(U / c),
                     float(np.sqrt(np.mean(best.fun ** 2))), picks)
