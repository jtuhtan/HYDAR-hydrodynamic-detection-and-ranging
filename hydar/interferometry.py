"""Case study 4: array design and ship tracking by correlative interferometry.

Case studies 2 and 3 locate a vessel with Kelvin wake theory: the 19.47 deg
cusp line (case 2) and the cusp angle and spectrogram chirp of the two Kelvin
branches (case 3). This module tracks a vessel without any of that. It uses
only two facts of linear wave physics that hold for any steady wave pattern
made by any source moving at constant velocity V, in any depth:

1. Stationarity. A wave that belongs to the steady pattern does not move
   relative to the source. A wave component with wavenumber vector
   k = k n (n a unit vector) and frequency omega therefore satisfies
   omega = k . V, i.e.

       n . V = c(omega)                                            (1)

   with c = omega / k the phase speed from the dispersion relation
   omega^2 = g k tanh(k h).

2. Rays. The energy of that component left the source when it was at
   P(t_e) = P0 + V t_e and travelled to the sensor at the group velocity
   c_g n. Observed at position r and time t:

       P0 + V t_e + c_g n (t - t_e) = r                             (2)

Equation (1) is linear in V, and with V known equation (2) is linear in P0
once t_e is eliminated. A 3 x 3 array measures n for every frequency bin of
every time window by correlative interferometry: the measured cross-spectral
phases between all 36 sensor pairs are correlated with the phases a plane
wave from each candidate direction would produce (the "reference table"
of correlative-interferometer direction finders in radio), and the best
match gives n and a correlation coefficient rho in [0, 1] that says how
plane-wave-like the field is. Each (t, f, n) observation from any array is
one equation (1) and one equation (2). A robust fit (RANSAC, then Huber
least squares) gives V and P0, so the track P(t) = P0 + V t.

The Kelvin angle, the cusp and the chirp are consequences of (1) and (2) in
deep water. They are not assumed here, so the method also applies in finite
depth (any Froude number below 1) and to hulls whose wake does not show a
clean cusp. Wind waves are free waves; they do not satisfy (1) for one V, so
the robust fit rejects them.

Array design
------------
How well an array measures n depends on its geometry. For a plane wave with
wavenumber k and direction theta, the Cramer-Rao bound of the bearing is

    var(theta) >= 1 / (2 SNR k^2 sum_i ((r_i - r_mean) . t_theta)^2)

with t_theta the unit vector normal to n: a wide aperture measures the
bearing better. But a wide, regular grid aliases: its correlation pattern
gets grating lobes once the spacing exceeds half a wavelength, and in noise
the wrong lobe wins. design_array() searches 3 x 3 grids with free column
and row positions (non-uniform spacing) for the geometry with the smallest
mean bearing bound over the frequency band and all directions, subject to
an upper limit on the highest sidelobe of the correlation pattern, a
maximum footprint and a minimum sensor spacing.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import differential_evolution, least_squares
from scipy.signal import butter, sosfiltfilt

from .sensors import wavenumber

G = 9.81


# ---------------------------------------------------------------- dispersion
def dispersion(f, depth: float = np.inf):
    """Wavenumber k [1/m], phase speed c and group speed c_g [m/s] at frequency f [Hz].

    depth = np.inf gives the deep-water relation omega^2 = g k.
    """
    w = 2 * np.pi * np.atleast_1d(np.asarray(f, float))
    if np.isinf(depth):
        k = w ** 2 / G
        return k, w / k, 0.5 * w / k
    k = wavenumber(w, depth)
    kh = k * depth
    c = w / k
    n = 0.5 * (1 + np.where(kh < 20, 2 * kh / np.sinh(np.minimum(2 * kh, 40)), 0.0))
    return k, c, n * c


def _unit(theta):
    return np.stack([np.cos(theta), np.sin(theta)], axis=-1)


# ---------------------------------------------------------------- array design
def ci_pattern(offsets: np.ndarray, k: float, theta0: float = 0.0, n_dir: int = 720) -> tuple:
    """Correlation pattern of correlative interferometry for one array.

    For a noise-free plane wave of wavenumber k from direction theta0, the
    correlation coefficient the direction search sees at every candidate
    direction (n_dir of them). Returns (theta, rho).
    """
    th = np.linspace(0, 2 * np.pi, n_dir, endpoint=False)
    i, j = np.triu_indices(len(offsets), 1)
    d = offsets[i] - offsets[j]                         # (n_pairs, 2)
    ph0 = k * d @ _unit(theta0)
    ph = k * d @ _unit(th).T                            # (n_pairs, n_dir)
    return th, np.cos(ph - ph0[:, None]).mean(axis=0)


def peak_sidelobe(offsets: np.ndarray, k: float, n_dir: int = 360, n_src: int = 24) -> float:
    """Highest secondary maximum of the correlation pattern, over n_src source directions.

    A sidelobe is any local maximum of rho(theta) other than the main peak.
    Close to 1 means an ambiguous array: in noise the wrong lobe can win
    (grating lobes of a sparse regular grid are the typical case).
    """
    th = np.linspace(0, 2 * np.pi, n_dir, endpoint=False)
    i, j = np.triu_indices(len(offsets), 1)
    d = offsets[i] - offsets[j]
    ph = k * d @ _unit(th).T                             # (n_pairs, n_dir)
    worst = -1.0
    for s in range(n_src):
        k0 = s * n_dir // n_src
        rho = np.cos(ph - ph[:, [k0]]).mean(axis=0)
        loc = (rho >= np.roll(rho, 1)) & (rho > np.roll(rho, -1))
        loc[k0] = False
        # drop the main peak itself if the grid put it next to k0
        near = np.abs((th - th[k0] + np.pi) % (2 * np.pi) - np.pi) < 2 * (th[1] - th[0])
        loc &= ~near
        if loc.any():
            worst = max(worst, float(rho[loc].max()))
    return max(worst, 0.0)


def bearing_crlb(offsets: np.ndarray, k: float, snr: float = 1.0, n_dir: int = 72) -> np.ndarray:
    """Cramer-Rao bound of the bearing [rad] for each of n_dir directions."""
    th = np.linspace(0, 2 * np.pi, n_dir, endpoint=False)
    r = offsets - offsets.mean(axis=0)
    t = _unit(th + np.pi / 2)
    moment = ((r @ t.T) ** 2).sum(axis=0)
    return 1.0 / np.sqrt(2 * snr * k ** 2 * moment)


def grid_layout(cols, rows) -> np.ndarray:
    """3 x 3 (or n x m) grid offsets from column x positions and row y positions."""
    X, Y = np.meshgrid(np.asarray(cols, float), np.asarray(rows, float), indexing="ij")
    out = np.column_stack([X.ravel(), Y.ravel()])
    return out - out.mean(axis=0)


@dataclass
class ArrayDesign:
    """Result of design_array."""

    offsets: np.ndarray          # (9, 2) sensor offsets from the centre [m]
    cols: np.ndarray             # column x positions [m]
    rows: np.ndarray             # row y positions [m]
    mean_crlb_deg: float         # mean bearing bound over band and directions (SNR = 1)
    peak_sidelobe: float         # worst sidelobe over the band
    aperture: float              # largest sensor separation [m]
    band: tuple
    depth: float

    def __str__(self) -> str:
        gx, gy = np.diff(self.cols), np.diff(self.rows)
        return (f"columns at {np.round(self.cols - self.cols.mean(), 2)} m (gaps {np.round(gx, 2)}), "
                f"rows at {np.round(self.rows - self.rows.mean(), 2)} m (gaps {np.round(gy, 2)}); "
                f"aperture {self.aperture:.1f} m, peak sidelobe {self.peak_sidelobe:.2f}, "
                f"mean bearing bound {self.mean_crlb_deg:.2f} deg at SNR 1")


def evaluate_layout(offsets: np.ndarray, band: tuple, depth: float = np.inf, n_freq: int = 6):
    """(mean bearing bound [deg], peak sidelobe) of a layout over a frequency band."""
    f = np.linspace(band[0], band[1], n_freq)
    k, _, _ = dispersion(f, depth)
    crlb = np.mean([bearing_crlb(offsets, kk).mean() for kk in k])
    psl = max(peak_sidelobe(offsets, kk) for kk in k)
    return float(np.degrees(crlb)), float(psl)


def design_array(band: tuple, depth: float = np.inf, max_aperture: float = 20.0,
                 min_spacing: float = 0.5, psl_max: float = 0.6, n_freq: int = 6,
                 seed: int = 0, maxiter: int = 60) -> ArrayDesign:
    """Optimal 3 x 3 grid for correlative interferometry over a frequency band.

    The design variables are the width and height of a 3 x 3 grid and where
    the middle column and the middle row sit between the outer ones, so the
    grid may be rectangular and non-uniform. The objective is the mean bearing
    Cramer-Rao bound over the band and over all source directions.
    Constraints, applied as penalties: highest sidelobe of the correlation
    pattern at most psl_max at every frequency of the band (ambiguity and
    aliasing), largest sensor separation (the diagonal) at most max_aperture
    (footprint), and every gap at least min_spacing. Solved with differential
    evolution; the returned design is re-evaluated at full resolution.
    """
    f = np.linspace(band[0], band[1], n_freq)
    k, _, _ = dispersion(f, depth)

    def layout(q):
        Lx, Ly, ux, uy = q
        cx = min_spacing + ux * (Lx - 2 * min_spacing)
        cy = min_spacing + uy * (Ly - 2 * min_spacing)
        return grid_layout([0, cx, Lx], [0, cy, Ly])

    def cost(q):
        off = layout(q)
        ap = np.hypot(q[0], q[1])
        crlb = np.mean([bearing_crlb(off, kk, n_dir=36).mean() for kk in k])
        psl = max(peak_sidelobe(off, kk, n_dir=360, n_src=12) for kk in k)
        pen = 200 * max(0.0, psl - psl_max) + 200 * max(0.0, ap / max_aperture - 1)
        return np.log(crlb) + pen

    bounds = [(2 * min_spacing, max_aperture), (2 * min_spacing, max_aperture), (0, 1), (0, 1)]
    res = differential_evolution(cost, bounds, seed=seed, maxiter=maxiter, tol=1e-8,
                                 polish=False, popsize=15, init="sobol")
    off = layout(res.x)
    crlb, psl = evaluate_layout(off, band, depth)
    ap = float(np.sqrt(((off[:, None] - off[None]) ** 2).sum(-1)).max())
    xs, ys = np.unique(np.round(off[:, 0], 9)), np.unique(np.round(off[:, 1], 9))
    return ArrayDesign(off, xs, ys, crlb, psl, ap, tuple(band), depth)


# ---------------------------------------------------------------- measurement
@dataclass
class Observations:
    """Plane-wave observations from correlative interferometry.

    One row per (array, time window, frequency bin) that passed the
    coherence test.
    """

    t: np.ndarray        # window centre [s]
    f: np.ndarray        # frequency [Hz]
    theta: np.ndarray    # propagation direction [rad, counter-clockwise from +x]
    rho: np.ndarray      # correlation coefficient of the best direction
    power: np.ndarray    # mean auto-spectral power at that bin [Pa^2/Hz]
    c: np.ndarray        # phase speed from the dispersion relation [m/s]
    cg: np.ndarray       # group speed [m/s]
    r: np.ndarray        # (n, 2) array centre [m]
    array: np.ndarray    # index of the array

    @property
    def n(self) -> np.ndarray:
        return _unit(self.theta)

    def __len__(self) -> int:
        return len(self.t)

    def subset(self, m) -> "Observations":
        return Observations(*(getattr(self, k)[m] for k in
                              ("t", "f", "theta", "rho", "power", "c", "cg", "r", "array")))

    @staticmethod
    def concat(obs) -> "Observations":
        keys = ("t", "f", "theta", "rho", "power", "c", "cg", "r", "array")
        return Observations(*(np.concatenate([getattr(o, k) for o in obs]) for k in keys))


def ci_scan(array, band: tuple, depth: float = np.inf, window: float = 120.0, hop: float = 30.0,
            segment: float | None = None, fs_work: float | None = None, n_dir: int = 720,
            pad: int = 4, time_stamp: str = "centre", t_start: float = 0.0,
            t_end: float | None = None):
    """Correlative-interferometry scan of one array.

    Sliding windows of `window` s (step `hop`); inside each window a Welch
    cross-spectral matrix (segments of `segment` s, half overlap, Hann,
    zero-padded `pad` times) is normalized to the complex coherency gamma_ij. For each frequency bin in
    band, every candidate direction theta gets the correlation

        rho(theta) = mean_{i<j} Re[ gamma_ij exp(i k n(theta) . (r_i - r_j)) ]

    with k from the dispersion relation, and the best theta is refined by a
    parabola.

    time_stamp sets the time given to each cell: "centre" (the window centre,
    right for a slowly varying wave field such as a far ship's wake) or
    "centroid" (the power-weighted mean time of the Welch segments at that
    frequency, for a near, fast-changing field, where a wave group crosses
    the array early or late in the window).
    Returns (t_cell, f, theta_best, rho_best, power), each (n_windows, n_freq)
    except f.
    """
    p = array.dynamic_pressure
    t = array.time
    fs = array.fs
    if fs_work is None:
        fs_work = min(fs, 8 * band[1])
    step = max(1, int(round(fs / fs_work)))
    fs_w = fs / step
    sos = butter(4, min(0.45 * fs_w, 0.99 * fs / 2), fs=fs, output="sos")
    x = sosfiltfilt(sos, p, axis=1)[:, ::step]
    tw = t[::step]
    if segment is None:
        segment = window / 2
    nseg = int(round(segment * fs_w))
    nwin = int(round(window * fs_w))
    nhop = max(1, int(round(hop * fs_w)))
    nfft = pad * nseg
    freqs = np.fft.rfftfreq(nfft, 1 / fs_w)
    fsel = (freqs >= band[0]) & (freqs <= band[1])
    freqs = freqs[fsel]
    k, _, _ = dispersion(freqs, depth)
    i, j = np.triu_indices(array.n_sensors, 1)
    d = array.positions[i] - array.positions[j]
    th = np.linspace(0, 2 * np.pi, n_dir, endpoint=False)
    # reference table: expected pair phases for each (frequency, direction)
    ref = np.exp(1j * k[:, None, None] * (d @ _unit(th).T)[None])    # (nf, npairs, ndir)
    hann = np.hanning(nseg)
    t_hi = tw[-1] if t_end is None else t_end
    starts = [s for s in range(0, len(tw) - nwin + 1, nhop)
              if tw[s] >= t_start and tw[s + nwin - 1] <= t_hi]
    T, TH, RHO, PW = [], [], [], []
    for s in starts:
        seg_starts = range(s, s + nwin - nseg + 1, max(1, nseg // 2))
        C = None
        tc_num, tc_den = 0.0, 0.0
        for q in seg_starts:
            xq = x[:, q:q + nseg]
            X = np.fft.rfft((xq - xq.mean(1, keepdims=True)) * hann, n=nfft, axis=1)[:, fsel]
            Cq = X[:, None, :] * np.conj(X[None, :, :])                 # (ns, ns, nf)
            C = Cq if C is None else C + Cq
            pq = (np.abs(X) ** 2).sum(axis=0)
            tc_num = tc_num + pq * (tw[q] + 0.5 * (nseg - 1) / fs_w)
            tc_den = tc_den + pq
        auto = np.real(np.einsum("iif->if", C))
        gam = C[i, j, :] / np.sqrt(auto[i] * auto[j] + 1e-300)          # (npairs, nf)
        R = np.real(np.einsum("pf,fpd->fd", gam, ref)) / len(i)          # (nf, ndir)
        m = np.argmax(R, axis=1)
        y0, y1, y2 = R[np.arange(len(m)), m - 1], R[np.arange(len(m)), m], R[np.arange(len(m)), (m + 1) % n_dir]
        den = y0 - 2 * y1 + y2
        dm = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
        if time_stamp == "centroid":
            T.append(tc_num / (tc_den + 1e-300))
        else:
            T.append(np.full(len(freqs), tw[s] + 0.5 * (nwin - 1) / fs_w))
        TH.append((th[m] + dm * (th[1] - th[0])) % (2 * np.pi))
        RHO.append(y1)
        PW.append(auto.mean(axis=0) / (fs_w * (hann ** 2).sum()) / len(seg_starts))
    return np.array(T), freqs, np.array(TH), np.array(RHO), np.array(PW)


def observe(array, band: tuple, depth: float = np.inf, rho_min: float = 0.9,
            excess: float = 3.0, peak_rel: float = 0.0, index: int = 0, **scan_kw) -> Observations:
    """Plane-wave observations of one array: ci_scan + two detection tests.

    A (window, frequency) cell is kept when
      * its correlation coefficient is at least rho_min (a plane wave crosses
        the array), and
      * its power is at least `excess` times the median power of the same
        frequency bin over the whole record (something transient arrived;
        a stationary wind sea sits at the median and is rejected), and
      * it is a local maximum of the spectrum, at least `peak_rel` times the
        strongest cell of its window (peak_rel = 0.01 drops the sidelobes
        that a strong, narrow peak leaks into neighbouring bins).
    """
    T, f, TH, RHO, PW = ci_scan(array, band, depth, **scan_kw)
    if len(T) == 0:
        raise RuntimeError("record shorter than one window")
    keep = (RHO >= rho_min) & (PW >= excess * np.median(PW, axis=0, keepdims=True))
    # one observation per spectral peak, its frequency refined by a parabola
    # through the log power (the wake chirps; bin centres would bias c, c_g)
    peak = np.zeros_like(keep)
    peak[:, 1:-1] = (PW[:, 1:-1] >= PW[:, :-2]) & (PW[:, 1:-1] > PW[:, 2:])
    peak &= PW >= peak_rel * PW.max(axis=1, keepdims=True)
    wi, fi = np.nonzero(keep & peak)
    L = np.log(PW + 1e-300)
    y0, y1, y2 = L[wi, fi - 1], L[wi, fi], L[wi, fi + 1]
    den = y0 - 2 * y1 + y2
    df = np.where(np.abs(den) > 1e-12, 0.5 * (y0 - y2) / den, 0.0)
    fr = f[fi] + np.clip(df, -0.5, 0.5) * (f[1] - f[0])
    _, c, cg = dispersion(fr, depth)
    return Observations(T[wi, fi], fr, TH[wi, fi], RHO[wi, fi], PW[wi, fi], c, cg,
                        np.tile(array.center, (len(wi), 1)), np.full(len(wi), index))


def detection_statistic(array, band: tuple, depth: float = np.inf, rho_min: float = 0.9,
                        **scan_kw):
    """Coherent excess power per window.

    Each cell's power is divided by the median power of its frequency bin
    over the record; the statistic sums that ratio over the bins of the
    window whose correlation is at least rho_min, divided by the number of
    bins. Stationary noise gives about rho_min-weighted 1 or less; a
    transient, coherent wake gives much more. Returns (t_window, D).
    """
    T, f, TH, RHO, PW = ci_scan(array, band, depth, **scan_kw)
    T = np.median(T, axis=1)
    ratio = PW / (np.median(PW, axis=0, keepdims=True) + 1e-300)
    D = np.where(RHO >= rho_min, ratio, 0.0).mean(axis=1)
    return T, D


# ---------------------------------------------------------------- tracking
@dataclass
class CITrack:
    """Constant-velocity track from correlative interferometry."""

    p0: np.ndarray            # position at t = 0 [m]
    velocity: np.ndarray      # [m/s]
    n_obs: int                # observations used (inliers)
    n_total: int              # observations offered
    rms_speed: float          # rms misfit of equation (1) [m/s]
    rms_ray: float            # rms misfit of equation (2) [m]
    observations: Observations = field(repr=False, default=None)
    inliers: np.ndarray = field(repr=False, default=None)

    @property
    def speed(self) -> float:
        return float(np.linalg.norm(self.velocity))

    @property
    def heading_deg(self) -> float:
        return float(np.degrees(np.arctan2(self.velocity[1], self.velocity[0])))

    def position(self, t) -> np.ndarray:
        return self.p0 + np.multiply.outer(np.asarray(t, float), self.velocity)

    def __str__(self) -> str:
        return (f"speed {self.speed:.2f} m/s, heading {self.heading_deg:.1f} deg, "
                f"position at t=0 ({self.p0[0]:.1f}, {self.p0[1]:.1f}) m, "
                f"{self.n_obs}/{self.n_total} observations, rms {self.rms_speed:.3f} m/s, "
                f"{self.rms_ray:.1f} m")


def _ray_rows(obs: Observations, V: np.ndarray):
    """Rows of the linear system m . P0 = m . (r - c_g n t) (equation 2)."""
    n = obs.n
    a = V[None, :] - obs.cg[:, None] * n
    m = np.column_stack([-a[:, 1], a[:, 0]])
    m /= np.linalg.norm(m, axis=1, keepdims=True) + 1e-300
    rhs = np.einsum("ij,ij->i", m, obs.r - obs.cg[:, None] * n * obs.t[:, None])
    return m, rhs


def _velocity_ls(obs, w=None):
    A, b = obs.n, obs.c
    if w is not None:
        A, b = A * w[:, None], b * w
    return np.linalg.lstsq(A, b, rcond=None)[0]


def fit_track(obs: Observations, tol: float = 0.04, tol_t: float = 60.0, n_iter: int = 3000,
              seed: int = 0, v_max: float = 40.0, weighted: bool = False) -> CITrack:
    """Robust constant-velocity track from plane-wave observations.

    1. RANSAC: two observations fix V (equation 1, 2 x 2) and P0
       (equation 2, 2 x 2). An observation is an inlier when its relative
       misfit of (1) is below `tol` and its misfit of (2), expressed as a
       travel time (metres / c_g), is below `tol_t` seconds. The candidate
       with the largest rho-weighted inlier count wins. Equation (2) is what
       rejects wind waves: the same wind-wave direction seen in many windows
       cannot come from one point of a moving source.
    2. Joint robust least squares (Cauchy loss, scales from the median
       absolute deviation of each misfit) of P0 and V on the inliers, each
       observation weighted by rho^2 / (1 - rho^2) if `weighted`.
    """
    if len(obs) < 4:
        raise RuntimeError("too few coherent observations to fit a track")
    rng = np.random.default_rng(seed)
    n, c, w = obs.n, obs.c, obs.rho
    N = len(obs)

    def inliers(P, V):
        r1 = np.abs(n @ V - c) / c
        m, rhs = _ray_rows(obs, V)
        r2 = np.abs(m @ P - rhs) / obs.cg
        return (r1 < tol) & (r2 < tol_t)

    best, best_score = None, -1.0
    for _ in range(n_iter):
        a, b = rng.choice(N, 2, replace=False)
        A = np.array([n[a], n[b]])
        if abs(np.linalg.det(A)) < 0.02:
            continue
        V = np.linalg.solve(A, [c[a], c[b]])
        if not 0.1 < np.linalg.norm(V) <= v_max:
            continue
        m, rhs = _ray_rows(obs.subset(np.array([a, b])), V)
        if abs(np.linalg.det(m)) < 1e-3:
            continue
        P = np.linalg.solve(m, rhs)
        score = w[inliers(P, V)].sum()
        if score > best_score:
            best, best_score = (P, V), score
    if best is None:
        raise RuntimeError("no consistent track found")
    P0, V = best

    def resid(q, sub):
        P, v = q[:2], q[2:]
        mm, rr = _ray_rows(sub, v)
        return (sub.n @ v - sub.c) / sub.c, (mm @ P - rr) / sub.cg

    def mad(x):
        return 1.4826 * np.median(np.abs(x - np.median(x)))

    inl = inliers(P0, V)
    sub = obs.subset(inl)
    if len(sub) < 4:
        raise RuntimeError("too few inliers")
    # weights: the bearing variance of a plane-wave fit grows like
    # (1 - rho^2) / rho^2, so weight each observation by its inverse
    if weighted:
        wt = np.clip(sub.rho ** 2 / np.maximum(1 - sub.rho ** 2, 1e-3), 0, 1e3)
        sw = np.sqrt(wt / wt.mean())
    else:
        sw = np.ones(len(sub))
    # robust scales of both misfits from the RANSAC solution; the Cauchy loss
    # then gives little weight to observations far outside them (e.g. near
    # the cusp, a caustic where ray theory and the window timing are poor)
    for _ in range(3):
        r1, r2 = resid(np.r_[P0, V], sub)
        s1, s2 = max(mad(r1), 1e-3), max(mad(r2), 0.02 * tol_t)
        q = least_squares(lambda q: np.concatenate([sw * a / b for a, b in zip(resid(q, sub), (s1, s2))]),
                          np.r_[P0, V], loss="cauchy", f_scale=2.0).x
        P0, V = q[:2], q[2:]
    inl = inliers(P0, V)
    sub = obs.subset(inl)
    mm, rr = _ray_rows(sub, V)
    return CITrack(P0, V, int(inl.sum()), N,
                   float(np.sqrt(np.mean((sub.n @ V - sub.c) ** 2))),
                   float(np.sqrt(np.mean((mm @ P0 - rr) ** 2))), obs, inl)


def track_vessel(arrays, band: tuple, depth: float = np.inf, rho_min: float = 0.9,
                 excess: float = 3.0, tol: float = 0.04, tol_t: float | None = None,
                 weighted: bool = False, **scan_kw) -> CITrack:
    """Full case-study-4 pipeline: observe every array, fuse, fit the track.

    tol_t defaults to half the analysis window.
    """
    obs = Observations.concat([observe(a, band, depth, rho_min, excess, index=i, **scan_kw)
                               for i, a in enumerate(arrays)])
    if tol_t is None:
        tol_t = 0.5 * scan_kw.get("window", 120.0)
    return fit_track(obs, tol=tol, tol_t=tol_t, weighted=weighted)


# ---------------------------------------------------------------- wind sea in finite depth
def wind_sea_head(points: np.ndarray, t: np.ndarray, depth: float, z: float, hs: float,
                  f_peak: float, mean_dir_deg: float = 112.5, spread_deg: float = 30.0,
                  f_range: tuple | None = None, n: int = 80, seed: int = 0) -> np.ndarray:
    """Random short-crested wind sea as pressure head p_dyn/(rho g) [m] at depth z.

    Like hydar.kelvin.wind_sea (JONSWAP-like weighting, n random components
    with a Gaussian directional spread), but with the finite-depth dispersion
    relation and the Airy pressure attenuation cosh(k(h+z))/cosh(kh) of each
    component, so it can be added to SpectralSolver records (pressure_model
    "solver"). hs is the significant wave height at the surface.
    """
    rng = np.random.default_rng(seed)
    lo, hi = f_range if f_range is not None else (0.6 * f_peak, 1.8 * f_peak)
    f = rng.uniform(lo, hi, n)
    S = f ** -5 * np.exp(-1.25 * (f_peak / f) ** 4)
    a = np.sqrt(S / S.sum()) * hs / 4 * np.sqrt(2)
    d = np.radians(mean_dir_deg + spread_deg * rng.standard_normal(n))
    k = dispersion(f, depth)[0]
    Kp = np.exp(k * z) * (1 + np.exp(-2 * k * (depth + z))) / (1 + np.exp(-2 * k * depth))
    ph0 = rng.uniform(0, 2 * np.pi, n)
    pts = np.atleast_2d(points)
    proj = pts[:, 0, None] * np.cos(d) + pts[:, 1, None] * np.sin(d)
    arg = (k * proj)[:, :, None] - (2 * np.pi * f)[None, :, None] * t[None, None, :] + ph0[None, :, None]
    return np.einsum("j,ijk->ik", a * Kp, np.cos(arg))
