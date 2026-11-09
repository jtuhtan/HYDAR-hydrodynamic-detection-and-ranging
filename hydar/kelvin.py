"""Case study 3: localizing a passing ship from its Kelvin wake.

This follows the field method of Rätsep et al. (2020, 2021), who detected and
located ferries passing about 3 km from a 3 x 3 array of seabed pressure
sensors (2.5 m spacing, 100 Hz) in Tallinn Bay:

* the spectrogram of the wake gives the ship's speed and its distance;
* time delays of the same wave between sensors give the direction the wake
  travels;
* the Kelvin cusp angle (35 deg 16') turns that direction into the ship's
  course, and with the distance gives the ship's position.

Ships of this kind are subcritical and make a *dispersive* Kelvin wake, which
the shallow-water Solver cannot represent. The wake is therefore synthesized
here from linear deep-water theory, and added to a random wind sea and sensor
noise. The analysis functions only see the recorded pressures.

Theory used (deep water, ship speed U along a straight track)
-------------------------------------------------------------
In the ship's frame the wake is a sum of steady waves travelling at angle
theta to the course, with wavenumber k = g / (U cos theta)^2 (Havelock):

    eta(x, y) = Re  integral  A(theta) exp(i k(theta) (x cos theta + y sin theta)) d theta

A fixed observer a distance d from the track sees, once the ship has passed
by a distance x, the two stationary-phase branches with tan(theta) = T,

    2 T^2 - (x/d) T + 1 = 0          (real for x/d >= 2 sqrt 2: inside the wedge)

at the frequency f = g / (2 pi U cos theta). The smaller root is the
transverse branch (f falls towards f_t = g / (2 pi U)), the larger root the
divergent branch (f rises). Both start at the cusp, theta_c = arctan(1/sqrt2)
= 35.26 deg, with f_cusp = sqrt(3/2) f_t, when the ship is x = 2 sqrt(2) d past
the observer. Fitting these two branches to the spectrogram ridge gives U, d
and the time the ship passed abeam. At the cusp the waves travel at theta_c to
the course and with phase speed U cos theta_c.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares
from scipy.signal import butter, sosfiltfilt, spectrogram

G = 9.81
THETA_CUSP = np.arctan(1 / np.sqrt(2))   # 35.26 deg


# ----------------------------------------------------------------- the ship
@dataclass
class KelvinShip:
    """Ship on a straight track making a deep-water Kelvin wake.

    speed       : U [m/s]
    course_deg  : direction of travel, degrees counter-clockwise from +x (east)
    x0, y0      : ship position at t = 0 [m]
    amplitude   : wake amplitude scale [m]; the wave height at a few km is of
                  the same order
    length_scale: hull length scale [m] that sets how fast the short divergent
                  waves die out (larger = smoother hull, fewer short waves)
    """

    speed: float = 11.0
    course_deg: float = -69.0
    x0: float = 0.0
    y0: float = 0.0
    amplitude: float = 0.25
    length_scale: float = 4.0

    @property
    def e(self) -> np.ndarray:
        th = np.radians(self.course_deg)
        return np.array([np.cos(th), np.sin(th)])

    @property
    def n(self) -> np.ndarray:
        return np.array([-self.e[1], self.e[0]])

    def position(self, t) -> np.ndarray:
        return np.array([self.x0, self.y0]) + np.multiply.outer(np.asarray(t, float) * self.speed, self.e)

    def wake(self, points: np.ndarray, t: np.ndarray, n_theta: int = 12000,
             chunk: int = 200) -> np.ndarray:
        """Surface elevation of the wake at points (n_pts, 2) and times t (n_t,).

        Returns (n_pts, n_t). Evaluated by direct numerical integration of the
        Havelock integral, so the cusp and both branches are exact in linear
        deep-water theory.
        """
        U = self.speed
        k0 = G / U ** 2
        th = np.linspace(-np.radians(82), np.radians(82), n_theta)
        k = k0 / np.cos(th) ** 2
        A = (1 / np.cos(th) ** 3) * np.exp(-0.5 * (k * self.length_scale) ** 2)
        dth = th[1] - th[0]
        # normalise so that `amplitude` is the wave height scale at 1 km abeam
        norm = self.amplitude * np.sqrt(k0 * 1000.0) / (A.sum() * dth)
        out = np.zeros((len(points), len(t)))
        P = self.position(t)                                    # (n_t, 2)
        for i, r in enumerate(np.atleast_2d(points)):
            rel = r - P                                         # observer in ship frame
            xs = rel @ self.e                                   # >0 ahead, <0 behind
            ys = rel @ self.n
            behind = xs < 0
            idx = np.nonzero(behind)[0]
            for j0 in range(0, idx.size, chunk):
                jj = idx[j0:j0 + chunk]
                ph = k[None, :] * (xs[jj, None] * np.cos(th)[None, :]
                                   + np.abs(ys[jj, None]) * np.sin(th)[None, :])
                out[i, jj] = norm * np.real(np.exp(1j * ph) @ A) * dth
        return out


def wind_sea(points: np.ndarray, t: np.ndarray, hs: float = 0.15, f_peak: float = 0.4,
             mean_dir_deg: float = 112.5, spread_deg: float = 30.0, n: int = 80,
             seed: int = 0) -> np.ndarray:
    """Random short-crested deep-water wind sea (sum of n wave components)."""
    rng = np.random.default_rng(seed)
    f = rng.uniform(0.25, 0.7, n)
    # JONSWAP-like weighting around f_peak
    S = f ** -5 * np.exp(-1.25 * (f_peak / f) ** 4)
    a = np.sqrt(S / S.sum()) * hs / 4 * np.sqrt(2)
    d = np.radians(mean_dir_deg + spread_deg * rng.standard_normal(n))
    k = (2 * np.pi * f) ** 2 / G
    ph0 = rng.uniform(0, 2 * np.pi, n)
    pts = np.atleast_2d(points)
    proj = pts[:, 0, None] * np.cos(d) + pts[:, 1, None] * np.sin(d)      # (n_pts, n)
    arg = (k * proj)[:, :, None] - (2 * np.pi * f)[None, :, None] * t[None, None, :] + ph0[None, :, None]
    return np.einsum("j,ijk->ik", a, np.cos(arg))


def record_passage(arrays, ship: KelvinShip, t_end: float, depth: float, fs_model: float = 4.0,
                   wind: dict | None = None, noise_pa: float = 2.0, seed: int = 1) -> None:
    """Synthesize the wake (+ wind sea + sensor noise) and load it into arrays.

    The wave field is computed at fs_model and spline-interpolated to each
    array's sampling rate; noise_pa is white sensor noise [Pa] added afterwards
    as an equivalent surface elevation.
    """
    rng = np.random.default_rng(seed)
    tm = np.arange(0.0, t_end, 1.0 / fs_model)
    for a in arrays:
        eta_m = ship.wake(a.positions, tm)
        if wind is not None:
            eta_m += wind_sea(a.positions, tm, **wind)
        tf = np.arange(0.0, tm[-1], 1.0 / a.fs)
        eta = CubicSpline(tm, eta_m, axis=1)(tf)
        eta += rng.normal(0, noise_pa / (1000 * G), eta.shape)
        a.load(tf, eta, depth=depth)


# ----------------------------------------------------------------- theory
def branch_frequencies(t, U: float, d: float, t_abeam: float):
    """Transverse and divergent wake frequencies [Hz] seen at distance d.

    Returns (f_transverse, f_divergent); NaN outside the Kelvin wedge.
    """
    x = U * (np.asarray(t, float) - t_abeam)
    q = x / d
    disc = q ** 2 - 8
    ok = (q > 0) & (disc >= 0)
    sq = np.sqrt(np.where(ok, disc, 0))
    Tt = np.where(ok, (q - sq) / 4, np.nan)
    Td = np.where(ok, (q + sq) / 4, np.nan)
    f = lambda T: G * np.sqrt(1 + T ** 2) / (2 * np.pi * U)
    return f(Tt), f(Td)


def cusp_time(U: float, d: float, t_abeam: float) -> float:
    return t_abeam + 2 * np.sqrt(2) * d / U


# ----------------------------------------------------------------- analysis
def stacked_spectrogram(array, sensors=(0, 2, 4, 6, 8), fs_out: float = 4.0,
                        window_s: float = 64.0, overlap: float = 0.9, band=(0.05, 0.6)):
    """Spectrogram of the dynamic pressure, averaged over several sensors.

    Following Rätsep et al. (2020) the spectrograms of closely spaced sensors
    are combined; by default the four corners and the centre of a 3 x 3 array.
    Returns t [s], f [Hz], S (power, f x t), normalised to its maximum.
    """
    p = array.dynamic_pressure
    step = int(round(array.fs / fs_out))
    sos = butter(4, 0.4 * fs_out, fs=array.fs, output="sos")
    S = None
    for k in sensors:
        x = sosfiltfilt(sos, p[k])[::step]
        nper = int(window_s * fs_out)
        f, t, Sk = spectrogram(x - x.mean(), fs=fs_out, nperseg=nper,
                               noverlap=int(overlap * nper), nfft=4 * nper, window="hann")
        S = Sk if S is None else S + Sk
    keep = (f >= band[0]) & (f <= band[1])
    S = S[keep] / len(sensors)
    return t, f[keep], S / S.max()


def detect_wake(t, f, S, band=(0.08, 0.3), k_noise: float = 6.0):
    """Start and end time of the wake: band power above k_noise x its median."""
    sel = (f >= band[0]) & (f <= band[1])
    pw = S[sel].sum(axis=0)
    thr = k_noise * np.median(pw)
    on = np.nonzero(pw > thr)[0]
    if on.size == 0:
        return None
    return float(t[on[0]]), float(t[on[-1]]), pw, thr


def ridge_points(t, f, S, t_window, band=(0.08, 0.27), rel: float = 0.03, floor: float = 30.0,
                 max_peaks: int = 2):
    """Spectral peaks (time, frequency) of the wake inside t_window.

    Up to max_peaks peaks per time column, each at least `rel` of that column's
    strongest peak and `floor` times the band's median background. The upper
    band edge sits below the wind sea, which masks the upper divergent branch
    (as noted by Rätsep et al.).
    """
    from scipy.signal import find_peaks
    pts = []
    sel = (f >= band[0]) & (f <= band[1])
    fs_, Ss = f[sel], S[sel]
    bg = np.median(Ss)
    for j, tj in enumerate(t):
        if not (t_window[0] <= tj <= t_window[1]):
            continue
        col = Ss[:, j]
        pk, prop = find_peaks(col, height=max(rel * col.max(), floor * bg))
        if pk.size == 0:
            continue
        pk = pk[np.argsort(prop["peak_heights"])[::-1][:max_peaks]]
        pts += [(tj, fs_[i]) for i in pk]
    return np.array(pts)


@dataclass
class ChirpFit:
    speed: float        # U [m/s]
    distance: float     # perpendicular distance from the track d [m]
    t_abeam: float      # time the ship passed abeam [s]
    t_cusp: float       # time the cusp reached the array [s]
    rms_hz: float


def fit_chirp(points: np.ndarray, t_detect: float) -> ChirpFit:
    """Fit U, d and t_abeam so the Kelvin branches match the ridge points."""
    tp, fp = points[:, 0], points[:, 1]

    def resid(q):
        U, d, ta = q
        ft, fd = branch_frequencies(tp, U, d, ta)
        # the spectrogram window smears the cusp slightly earlier in time:
        # treat points just before the cusp as sitting at the cusp frequency
        f_c = np.sqrt(1.5) * G / (2 * np.pi * U)
        ft, fd = np.where(np.isnan(ft), f_c, ft), np.where(np.isnan(fd), f_c, fd)
        r = np.fmin(np.abs(np.nan_to_num(ft, nan=9) - fp), np.abs(np.nan_to_num(fd, nan=9) - fp))
        return np.minimum(r, 0.05)  # robust cap for stray peaks

    best = None
    U0 = G / (2 * np.pi * np.percentile(fp, 5))            # transverse asymptote
    for U in (0.9 * U0, U0, 1.1 * U0):
        for d in (500.0, 1000.0, 2000.0, 3000.0, 4500.0):
            ta = t_detect - 2 * np.sqrt(2) * d / U
            r = least_squares(resid, [U, d, ta], x_scale=[1.0, 500.0, 100.0],
                              bounds=([1.0, 50.0, -1e5], [40.0, 2e4, 1e5]))
            if best is None or r.cost < best.cost:
                best = r
    U, d, ta = best.x
    return ChirpFit(float(U), float(d), float(ta), float(cusp_time(U, d, ta)),
                    float(np.sqrt(np.mean(best.fun ** 2))))


def wave_direction(array, t0: float, t1: float, band=(0.1, 0.25)):
    """Propagation direction of the waves crossing the array in [t0, t1].

    Band-passes every sensor, cross-correlates all 36 sensor pairs to get
    their time delays (sub-sample peak interpolation), and fits a plane wave.
    Returns (unit propagation vector, apparent phase speed [m/s]).
    """
    p = array.dynamic_pressure
    t = array.time
    sos = butter(4, band, btype="band", fs=array.fs, output="sos")
    w = (t >= t0) & (t <= t1)
    x = np.array([sosfiltfilt(sos, pk)[w] for pk in p])
    n = len(x)
    span = max(np.linalg.norm(a - b) for a in array.positions for b in array.positions)
    max_lag = int(array.fs * span / 2.0)   # > largest inter-sensor delay at 2 m/s
    rows, lags = [], []
    for i in range(n):
        for j in range(i + 1, n):
            cc = np.correlate(x[j], x[i], mode="full")
            mid = len(x[i]) - 1
            seg = cc[mid - max_lag: mid + max_lag + 1]
            m = int(np.argmax(seg))
            if 0 < m < len(seg) - 1:   # parabolic sub-sample refinement
                y0, y1, y2 = seg[m - 1], seg[m], seg[m + 1]
                m = m + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2)
            lags.append((m - max_lag) / array.fs)
            rows.append(array.positions[j] - array.positions[i])
    s = np.linalg.lstsq(np.array(rows), np.array(lags), rcond=None)[0]   # slowness vector
    return s / np.linalg.norm(s), float(1 / np.linalg.norm(s))


@dataclass
class VesselFix:
    speed: float
    distance: float
    course_deg: float
    t_abeam: float
    t_cusp: float
    position_abeam: np.ndarray     # ship position when abeam of the array [m]
    cusp_direction: np.ndarray     # propagation direction of the cusp waves
    candidates_deg: tuple          # the two course candidates from the cusp angle
    speed_from_phase: float        # U = c_phase / cos(theta_c), independent check
    chirp: ChirpFit

    def position(self, t) -> np.ndarray:
        th = np.radians(self.course_deg)
        e = np.array([np.cos(th), np.sin(th)])
        return self.position_abeam + np.multiply.outer(np.asarray(t, float) - self.t_abeam, self.speed * e)

    def __str__(self) -> str:
        return (f"speed {self.speed:.2f} m/s, distance from track {self.distance:.0f} m, "
                f"course {self.course_deg:.1f} deg, abeam at t = {self.t_abeam:.0f} s, "
                f"cusp at t = {self.t_cusp:.0f} s (phase-speed check {self.speed_from_phase:.2f} m/s)")


def locate_vessel(array, spectro_kw: dict | None = None) -> VesselFix:
    """Full case-study-3 pipeline for one 3 x 3 array."""
    t, f, S = stacked_spectrogram(array, **(spectro_kw or {}))
    det = detect_wake(t, f, S)
    if det is None:
        raise RuntimeError("no wake detected")
    t_on, t_off = det[0], det[1]
    pts = ridge_points(t, f, S, (t_on, t_off))
    chirp = fit_chirp(pts, t_on)
    U, d = chirp.speed, chirp.distance

    # direction of the cusp waves (strongest, theta = theta_c) and of the late
    # transverse waves (theta -> 0, i.e. along the course)
    n_c, c_app = wave_direction(array, chirp.t_cusp - 10, chirp.t_cusp + 30)
    f_t = G / (2 * np.pi * U)
    t_late = chirp.t_abeam + 6 * d / U
    n_t, _ = wave_direction(array, t_late - 40, t_late + 40, band=(0.8 * f_t, 1.05 * f_t))
    rot = lambda v, a: np.array([v[0] * np.cos(a) - v[1] * np.sin(a), v[0] * np.sin(a) + v[1] * np.cos(a)])
    cands = [rot(n_c, -THETA_CUSP), rot(n_c, THETA_CUSP)]
    e = max(cands, key=lambda c: c @ n_t)
    course = float(np.degrees(np.arctan2(e[1], e[0])))
    # the array lies on the side of the track the cusp waves travel towards
    nperp = np.array([-e[1], e[0]])
    side = np.sign(n_c @ nperp)
    p_abeam = array.center - side * d * nperp
    return VesselFix(U, d, course, chirp.t_abeam, chirp.t_cusp, p_abeam, n_c,
                     tuple(float(np.degrees(np.arctan2(c[1], c[0]))) for c in cands),
                     float(c_app / np.cos(THETA_CUSP)), chirp)
