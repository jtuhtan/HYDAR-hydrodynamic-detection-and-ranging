"""Case study 5: a fish locates a rock in a river from its Karman vortex street.

A hemispherical rock sits on the bed of a river with a trapezoidal cross
section and sticks out of the water. The current sheds a Karman vortex street
from it. A 25 cm fish holds station somewhere downstream and carries a lateral
line of 10 pressure sensors: 5 along the head third of its body, 5 along the
middle third and none on the tail third. From the pressures alone, the fish
estimates where the rock is.

The vortex street is not a wave. Its vortices are carried downstream by the
current, so HYDAR's time of flight becomes a time of convection: the delay of
the street between the fish's head and its mid-body gives the convection
speed, and the sign of that delay says the source is upstream.

Model (all of it linear in the geometry, nothing is solved on a grid)
--------------------------------------------------------------------
River. Straight channel along +x (downstream), y across, y = 0 on the
centreline. Flat bed of width b, banks of slope 1 : m (horizontal : vertical),
flow depth h on the bed. The mean velocity follows Manning,
U = R^(2/3) S^(1/2) / n; across the section the depth-averaged velocity of
each vertical strip follows Manning with the local depth, scaled to the same
discharge. The bed friction coefficient is c_f = g n^2 / h^(1/3).

Rock. A hemisphere of radius R on the bed. "Protruding 50 %" is read as half
its height above the surface, so R = 2 h. The rock's sectional diameter falls
from 2R at the bed to 2 sqrt(R^2 - h^2) at the surface; the street is assumed
to be shed at one frequency set by the depth-averaged diameter D.

Street. Vortices are shed alternately from the two sides at
f = St U / D (Strouhal number St = 0.2) and appear one formation length
(1.5 D) behind the rock centre. They drift downstream at U_c = kappa U
(kappa = 0.88) in two rows b_s = 0.281 lambda apart (von Karman's stable
ratio), lambda = U_c / f being the streamwise spacing of like vortices. The
upper row (y > y_rock) turns clockwise. The circulation follows from von
Karman's relation for the self-induced speed of the street,

    U - U_c = Gamma / (2 lambda) * tanh(pi b_s / lambda),

so Gamma is not a free parameter. Each vortex is a Lamb-Oseen vortex whose
core grows by turbulent diffusion, r_c^2 = (0.1 D)^2 + 4 nu_t tau with
nu_t = 0.01 U D and tau the vortex age, and whose circulation decays by bed
friction as exp(-c_f U tau / h) (shallow-wake damping). A vortex grows to full
strength over half a shedding period after it forms. The banks are vertical
walls at the waterline (method of images).

Pressure. Unsteady Bernoulli for the street moving at U_c in a stream U,

    p = -rho [ (U - U_c) u'_x + |u'|^2 / 2 ],

with u' the velocity induced by all vortices. Bernoulli does not hold inside
a vortex core, so each real vortex contributes its exact Lamb-Oseen pressure
instead of -rho |u_k|^2 / 2; only the cross terms between vortices come from
Bernoulli. The fish is a passive probe: its own flow field is left out.

Noise. Every sensor has white noise (2 Pa by default, as in case study 4) and
the river carries turbulence: a frozen random pressure field advected at U,
with an RMS of rho (0.1 U)^2 (turbulence intensity 10 %) and a spectrum with
a corner set by an integral scale equal to the depth and an f^(-7/3)
inertial range. It covers the shedding frequency, so it is the main noise.

The fish's method (locate_rock)
-------------------------------
1. Detection and frequency. The spectrum of the whole survey shows a line at
   the shedding frequency f and its harmonic 2 f. Near the street's
   centreline the two rows cancel at f and only 2 f is left, so a peak is
   taken as the harmonic when there is a line at half its frequency.
2. Convection speed. At each station the phase of f and 2 f along the 10
   sensors falls linearly with position; the slope gives U_c and its sign
   says the street comes from upstream.
3. Size. U = U_c / kappa and D = St U / f, with the same St and kappa as the
   model (the fish has to assume them).
4. Lateral position. The fish steps across the river and holds each station
   for a while. The amplitude of f peaks under each vortex row and vanishes
   on the centreline, where only 2 f is left. The rock lies on the centreline
   of this mirror-symmetric pattern.
5. Range. The vortex cores grow and the circulation decays as the street
   ages, so the amplitude and width of the lateral profile change with
   distance. Fitting the model's profiles (with D, U and U_c from steps 2
   and 3) to the measured ones gives the distance to the rock and refines the
   lateral position. This step depends on the wake-ageing model, unlike steps
   1 to 4; see the sensitivity results in the README.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.special import exp1

from .sensors import SensorArray

G = 9.81
RHO = 1000.0
_trapz = getattr(np, "trapezoid", None) or np.trapz   # numpy < 2 has only trapz


# ----------------------------------------------------------------- the river
@dataclass
class TrapezoidalChannel:
    """Straight river with a trapezoidal cross section (see module docstring).

    bed_width  : width of the flat bed [m]
    side_slope : bank slope m, horizontal per vertical (2 means 1 : 2)
    depth      : flow depth over the flat bed [m]
    slope      : bed (energy) slope S [-]
    manning_n  : Manning roughness [s m^-1/3] (0.035: gravel-bed stream)
    """

    bed_width: float = 4.0
    side_slope: float = 2.0
    depth: float = 0.25
    slope: float = 0.002
    manning_n: float = 0.035
    g: float = G

    @property
    def top_width(self) -> float:
        return self.bed_width + 2 * self.side_slope * self.depth

    @property
    def area(self) -> float:
        return (self.bed_width + self.side_slope * self.depth) * self.depth

    @property
    def wetted_perimeter(self) -> float:
        return self.bed_width + 2 * self.depth * np.sqrt(1 + self.side_slope ** 2)

    @property
    def hydraulic_radius(self) -> float:
        return self.area / self.wetted_perimeter

    @property
    def mean_velocity(self) -> float:
        """Section-mean velocity from Manning's equation [m/s]."""
        return self.hydraulic_radius ** (2 / 3) * np.sqrt(self.slope) / self.manning_n

    @property
    def discharge(self) -> float:
        return self.mean_velocity * self.area

    @property
    def froude(self) -> float:
        return self.mean_velocity / np.sqrt(self.g * self.area / self.top_width)

    def local_depth(self, y) -> np.ndarray:
        """Water depth at lateral position y [m] (0 outside the wetted section)."""
        y = np.abs(np.asarray(y, dtype=float))
        bank = self.depth - (y - 0.5 * self.bed_width) / self.side_slope
        return np.clip(np.where(y <= 0.5 * self.bed_width, self.depth, bank), 0.0, None)

    def local_velocity(self, y) -> np.ndarray:
        """Depth-averaged velocity at y: Manning per strip, scaled to the discharge."""
        yy = np.linspace(-0.5 * self.top_width, 0.5 * self.top_width, 2001)
        q_shape = _trapz(self.local_depth(yy) ** (5 / 3), yy)
        return self.discharge * self.local_depth(y) ** (2 / 3) / q_shape

    def friction_coefficient(self, y=0.0) -> float:
        """Bed friction coefficient c_f = g n^2 / h^(1/3) at y."""
        return float(self.g * self.manning_n ** 2 / self.local_depth(y) ** (1 / 3))


# ----------------------------------------------------------------- the rock
@dataclass
class HemisphereRock:
    """Hemispherical rock of radius R resting on the bed at (x, y)."""

    x: float
    y: float
    radius: float

    @classmethod
    def from_protrusion(cls, channel: TrapezoidalChannel, x: float, y: float,
                        protrusion: float = 0.5) -> "HemisphereRock":
        """Rock whose height sticks out of the water by the given fraction (0.5: half)."""
        if not 0 < protrusion < 1:
            raise ValueError("protrusion must lie between 0 and 1")
        h = float(channel.local_depth(y))
        if h <= 0:
            raise ValueError("the rock must stand in water")
        return cls(x, y, h / (1 - protrusion))

    def section_diameter(self, z) -> np.ndarray:
        """Diameter of the rock at height z above the bed [m]."""
        z = np.clip(np.asarray(z, dtype=float), 0.0, self.radius)
        return 2 * np.sqrt(self.radius ** 2 - z ** 2)

    def waterline_diameter(self, depth: float) -> float:
        return float(self.section_diameter(depth))

    def shedding_diameter(self, depth: float) -> float:
        """Sectional diameter averaged over the water column, the length scale of the street."""
        z = np.linspace(0.0, min(depth, self.radius), 401)
        return float(_trapz(self.section_diameter(z), z) / depth)


# ----------------------------------------------------------------- the street
def lamb_oseen_pressure(r, gamma, rc, rho: float = RHO) -> np.ndarray:
    """Pressure of a Lamb-Oseen vortex relative to far away [Pa] (cyclostrophic balance).

    p(r) = -rho Gamma^2 / (8 pi^2 rc^2) [ (1 - e^-s)^2 / s + 2 (E1(s) - E1(2 s)) ],  s = r^2/rc^2;
    at the centre it is -rho Gamma^2 ln 2 / (4 pi^2 rc^2).
    """
    s = np.maximum(np.asarray(r, dtype=float) ** 2 / rc ** 2, 1e-12)
    F = np.expm1(-s) ** 2 / s + 2 * (exp1(s) - exp1(2 * s))     # (1 - e^-s)^2 / s + ...
    return -rho * gamma ** 2 / (8 * np.pi ** 2 * rc ** 2) * F


@dataclass
class KarmanStreet:
    """Karman vortex street shed by a rock in a river (see module docstring).

    strouhal       : St = f D / U
    convection     : kappa = U_c / U
    spacing_ratio  : b_s / lambda, lateral over streamwise vortex spacing
    formation      : distance from the rock centre to where vortices form, in D
    core           : initial core radius, in D
    eddy_viscosity : nu_t / (U D), turbulent diffusion of the cores
    friction_decay : damp the circulation by bed friction
    phase          : time [s] at which a vortex forms on the upper row
    n_images       : bank reflections on each side (0: no banks)
    velocity       : stream velocity U at the rock [m/s]; None: from the channel
    friction       : bed friction coefficient c_f; None: from the channel
    """

    channel: TrapezoidalChannel
    rock: HemisphereRock
    strouhal: float = 0.2
    convection: float = 0.88
    spacing_ratio: float = 0.281
    formation: float = 1.5
    core: float = 0.1
    eddy_viscosity: float = 0.01
    friction_decay: bool = True
    phase: float = 0.0
    n_images: int = 2
    rho: float = RHO
    reach: float = 40.0        # vortices are followed this far downstream [m]
    velocity: float | None = None
    friction: float | None = None

    def __post_init__(self) -> None:
        ch, rk = self.channel, self.rock
        self.depth = float(ch.local_depth(rk.y))
        self.U = float(ch.local_velocity(rk.y)) if self.velocity is None else float(self.velocity)
        self.D = rk.shedding_diameter(self.depth)
        self.f = self.strouhal * self.U / self.D
        self.Uc = self.convection * self.U
        self.wavelength = self.Uc / self.f
        self.row_spacing = self.spacing_ratio * self.wavelength
        self.gamma0 = (2 * self.wavelength * (self.U - self.Uc)
                       / np.tanh(np.pi * self.row_spacing / self.wavelength))
        self.cf = ch.friction_coefficient(rk.y) if self.friction is None else float(self.friction)
        self.decay_rate = self.cf * self.U / self.depth if self.friction_decay else 0.0
        self.nu_t = self.eddy_viscosity * self.U * self.D

    # -------------------------------------------------- derived numbers
    @property
    def reynolds(self) -> float:
        return self.U * self.D / 1.0e-6

    @property
    def stability_number(self) -> float:
        """Shallow-wake stability number S = c_f D / (2 h); a vortex street needs S below ~0.2."""
        return self.cf * self.D / (2 * self.depth)

    @property
    def x_formation(self) -> float:
        return self.rock.x + self.formation * self.D

    # -------------------------------------------------- vortices
    def vortices(self, t: float, x_max: float | None = None):
        """Real vortices alive at time t: arrays x, y, Gamma, r_c, and the row sign."""
        x_max = self.rock.x + self.reach if x_max is None else x_max
        half = 0.5 / self.f
        tau_max = (x_max - self.x_formation) / self.Uc
        k_hi = np.floor((t - self.phase) / half)
        k_lo = np.ceil((t - tau_max - self.phase) / half)
        k = np.arange(k_lo, k_hi + 1)
        tau = t - (self.phase + k * half)
        upper = (k % 2) == 0
        x = self.x_formation + self.Uc * tau
        y = self.rock.y + np.where(upper, 0.5, -0.5) * self.row_spacing
        ramp = -np.expm1(-(tau / half) ** 2)
        gamma = np.where(upper, -1.0, 1.0) * self.gamma0 * ramp * np.exp(-self.decay_rate * tau)
        rc = np.sqrt((self.core * self.D) ** 2 + 4 * self.nu_t * tau)
        return x, y, gamma, rc

    def _with_images(self, y, gamma):
        W = self.channel.top_width
        ys, gs = [y], [gamma]
        for n in range(1, self.n_images + 1):
            # n odd: mirror in one bank (sign flips); n even: shift by whole periods 2W
            for side in (1, -1):
                if n % 2:
                    ys.append(side * n * W - y)
                    gs.append(-gamma)
                else:
                    ys.append(y + side * n * W)
                    gs.append(gamma)
        return np.concatenate(ys), np.concatenate(gs)

    def pressure(self, x, y, t) -> np.ndarray:
        """Dynamic pressure [Pa] of the street at points (x, y) and times t.

        x, y : arrays of the same shape (n_points,), or scalars
        t    : (n_times,) times [s]
        returns (n_points, n_times)
        """
        x = np.atleast_1d(np.asarray(x, dtype=float))
        y = np.atleast_1d(np.asarray(y, dtype=float))
        t = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.empty((x.size, t.size))
        x_max = max(x.max(), self.x_formation) + 3 * self.wavelength
        for j, tj in enumerate(t):
            out[:, j] = self._pressure_at(x, y, tj, x_max)
        return out

    def _pressure_at(self, x, y, t, x_max):
        xv, yv, gv, rv = self.vortices(t, x_max)
        if xv.size == 0:
            return np.zeros_like(x)
        n_real = xv.size
        yi, gi = self._with_images(yv, gv)
        xi = np.tile(xv, len(yi) // n_real)
        ri = np.tile(rv, len(yi) // n_real)
        dx = x[:, None] - xi[None, :]
        dy = y[:, None] - yi[None, :]
        r2 = np.maximum(dx ** 2 + dy ** 2, 1e-12)
        fac = gi / (2 * np.pi * r2) * (-np.expm1(-r2 / ri ** 2))    # u_theta / r
        ux, uy = -fac * dy, fac * dx
        Ux, Uy = ux.sum(1), uy.sum(1)
        self_sq = (ux[:, :n_real] ** 2 + uy[:, :n_real] ** 2).sum(1)
        p_lo = lamb_oseen_pressure(np.sqrt(r2[:, :n_real]), gv, rv, self.rho).sum(1)
        return (-self.rho * (self.U - self.Uc) * Ux
                - 0.5 * self.rho * (Ux ** 2 + Uy ** 2 - self_sq) + p_lo)

    def harmonic_profile(self, dx: float, y_rel, n_phase: int = 32) -> tuple:
        """Amplitude [Pa] of the f and 2 f components at a distance dx behind the rock.

        y_rel : lateral offsets from the rock [m]. Uses one shedding period, during
        which the pattern repeats exactly. Returns (A1, A2), each shaped like y_rel.
        """
        y_rel = np.asarray(y_rel, dtype=float)
        t = self.phase + np.arange(n_phase) / (n_phase * self.f)
        p = self.pressure(np.full(y_rel.size, self.rock.x + dx), self.rock.y + y_rel.ravel(), t)
        c = np.fft.rfft(p, axis=1) * 2 / n_phase
        return np.abs(c[:, 1]).reshape(y_rel.shape), np.abs(c[:, 2]).reshape(y_rel.shape)


# ----------------------------------------------------------------- the fish
def fish_sensor_positions(length: float = 0.25, n_head: int = 5, n_middle: int = 5) -> np.ndarray:
    """Distances of the lateral-line sensors from the snout [m].

    n_head sensors are spread evenly over the head third of the body and n_middle
    over the middle third (each at the centre of an equal sub-segment); the tail
    third carries none.
    """
    third = length / 3
    head = (np.arange(n_head) + 0.5) * third / n_head
    middle = third + (np.arange(n_middle) + 0.5) * third / n_middle
    return np.concatenate([head, middle])


@dataclass
class Fish:
    """A fish holding station with its head upstream (facing -x).

    x, y   : position of the snout [m]
    length : body length [m]
    """

    x: float
    y: float
    length: float = 0.25
    n_head: int = 5
    n_middle: int = 5

    @property
    def sensor_offsets(self) -> np.ndarray:
        """Sensor distances behind the snout [m]; the body lies along +x from the snout."""
        return fish_sensor_positions(self.length, self.n_head, self.n_middle)

    @property
    def sensor_positions(self) -> np.ndarray:
        s = self.sensor_offsets
        return np.column_stack([self.x + s, np.full(s.size, self.y)])

    @property
    def center(self) -> np.ndarray:
        """Centre of the sensor line (it lies between the head and middle thirds)."""
        return self.sensor_positions.mean(0)

    def array(self, name: str = "fish", fs: float = 100.0, noise_pa: float = 0.0,
              noise_seed: int = 0) -> SensorArray:
        """The lateral line as a SensorArray (a line of sensors along the body)."""
        c = self.center
        layout = self.sensor_positions - c
        return SensorArray(name, float(c[0]), float(c[1]), z=-0.01, layout=layout, fs=fs,
                           noise_pa=noise_pa, noise_seed=noise_seed)


# ----------------------------------------------------------------- turbulence
@dataclass
class RiverTurbulence:
    """Background turbulence: a frozen random pressure field advected at U.

    The pressure fluctuations of turbulence scale with rho sigma_u^2 (sigma_u =
    intensity * U); their spectrum is flat at low frequency and falls as
    f^(-7/3) in the inertial range, with the integral length scale `scale`
    (about the depth) setting the corner, as in a von Karman form. The field
    is one random series g(tau) (Gaussian, built by inverse FFT), seen at x and t
    as g(t - x / U).

    intensity : sigma_u / U
    scale     : integral length scale [m]
    duration  : length of the series [s]; t - x / U must stay inside [-60 s, duration - 60 s]
    fs        : sampling rate of the series [Hz]
    """

    U: float
    intensity: float = 0.1
    scale: float = 0.25
    seed: int = 0
    rho: float = RHO
    duration: float = 3600.0
    fs: float = 25.0

    def __post_init__(self) -> None:
        rng = np.random.default_rng(self.seed)
        n = int(round(self.duration * self.fs))
        f = np.fft.rfftfreq(n, 1 / self.fs)
        L, U = self.scale, self.U
        S = 1.0 / (1 + 70.8 * (f * L / U) ** 2) ** (7 / 6)
        S[0] = 0.0
        spec = np.sqrt(S) * (rng.standard_normal(f.size) + 1j * rng.standard_normal(f.size))
        g = np.fft.irfft(spec, n=n)
        self._g = g * self.rho * (self.intensity * U) ** 2 / g.std()
        self._tau0 = -60.0

    @property
    def rms(self) -> float:
        return float(self._g.std())

    def pressure(self, x, t) -> np.ndarray:
        """(n_points, n_times) pressure [Pa] at streamwise positions x and times t."""
        x = np.atleast_1d(np.asarray(x, dtype=float))
        t = np.atleast_1d(np.asarray(t, dtype=float))
        tau = t[None, :] - x[:, None] / self.U
        i = (tau - self._tau0) * self.fs
        if i.min() < 0 or i.max() > self._g.size - 2:
            raise ValueError("times outside the turbulence series; raise duration")
        i0 = np.floor(i).astype(int)
        w = i - i0
        return (1 - w) * self._g[i0] + w * self._g[i0 + 1]


# ----------------------------------------------------------------- the survey
@dataclass
class Survey:
    """Lateral-line records of a fish that steps across the river.

    x_snout : streamwise position of the snout [m]
    y       : (n_stations,) lateral positions
    t       : (n_stations, n_samples) sample times
    p       : (n_stations, n_sensors, n_samples) dynamic pressure [Pa]
    offsets : (n_sensors,) sensor positions along the body, relative to its centre [m]
    """

    x_snout: float
    y: np.ndarray
    t: np.ndarray
    p: np.ndarray
    offsets: np.ndarray
    fs: float
    length: float

    @property
    def x_center(self) -> float:
        return self.x_snout + self.length / 3


def lateral_survey(street: KarmanStreet, x_snout: float, ys, hold: float = 60.0,
                   fs: float = 100.0, t0: float = 0.0, noise_pa: float = 2.0,
                   turbulence: RiverTurbulence | None = None, seed: int = 0,
                   fish_length: float = 0.25, fs_model: float = 10.0) -> Survey:
    """Record the street with a fish that holds each lateral station for `hold` seconds.

    The stations follow one another in time (the fish moves between them
    instantly). The street is evaluated at fs_model and interpolated to the
    sampling rate fs (its content lies below 1 Hz); sensor noise and turbulence
    are added at fs.
    """
    ys = np.asarray(ys, dtype=float)
    n = int(round(hold * fs))
    T, P = [], []
    for i, y in enumerate(ys):
        fish = Fish(x_snout, float(y), fish_length)
        t = t0 + i * hold + np.arange(n) / fs
        pos = fish.sensor_positions
        tm = np.arange(t[0] - 1.0, t[-1] + 1.0 + 0.5 / fs_model, 1.0 / fs_model)
        p = CubicSpline(tm, street.pressure(pos[:, 0], pos[:, 1], tm), axis=1)(t)
        if turbulence is not None:
            p += turbulence.pressure(pos[:, 0], t)
        arr = fish.array(fs=fs, noise_pa=noise_pa, noise_seed=seed * 1000 + i)
        arr.load(t - t[0], p / (street.rho * G), depth=street.depth, rho=street.rho, g=G)
        T.append(t)
        P.append(arr.dynamic_pressure)
    fish = Fish(x_snout, 0.0, fish_length)
    return Survey(x_snout, ys, np.array(T), np.array(P), fish.sensor_offsets - fish.sensor_offsets.mean(),
                  fs, fish_length)


# ----------------------------------------------------------------- the fish's method
@dataclass
class RockFix:
    """What the fish concludes about the rock."""

    x: float                 # streamwise position of the rock [m]
    y: float                 # lateral position [m]
    distance: float          # from the rock centre to the centre of the sensor line [m]
    diameter: float          # shedding diameter D [m]
    f: float                 # shedding frequency [Hz]
    Uc: float                # convection speed [m/s]
    U: float                 # stream velocity [m/s]
    upstream: bool           # the street travels from head to tail
    cost: float              # sum of squared, noise-normalized profile residuals
    snr: float               # largest amplitude of f over its noise amplitude
    A1: np.ndarray = field(repr=False)   # measured amplitude of f per station [Pa]
    A2: np.ndarray = field(repr=False)   # measured amplitude of 2 f per station [Pa]
    noise1: float = field(repr=False)    # noise amplitude near f [Pa]
    noise2: float = field(repr=False)    # noise amplitude near 2 f [Pa]
    misfit: np.ndarray = field(repr=False)   # (n_dx, n_y) fit cost over the search grid
    dx_grid: np.ndarray = field(repr=False)
    yc_grid: np.ndarray = field(repr=False)
    model: "KarmanStreet" = field(repr=False)   # the fish's wake model, rock at (0, 0)
    detected: bool = True

    @property
    def wavelength(self) -> float:
        return self.Uc / self.f


def _spectral_peak(survey: Survey, band=(0.02, 2.0), pad: int = 16) -> float:
    """Strongest line of the spectrum averaged over stations and sensors [Hz].

    Each station is windowed on its own: the phase of the street jumps by half a
    period from one vortex row to the other, so the stations cannot be joined
    into one long record. The line is the shedding frequency f or, if the fish
    spent most of its time near the centreline, its harmonic 2 f.
    """
    x = survey.p - survey.p.mean(-1, keepdims=True)
    n = x.shape[-1]
    nf = pad * n
    P = (np.abs(np.fft.rfft(x * np.hanning(n), n=nf, axis=-1)) ** 2).sum((0, 1))
    fr = np.fft.rfftfreq(nf, 1 / survey.fs)
    sel = (fr >= band[0]) & (fr <= band[1])
    i = np.flatnonzero(sel)[np.argmax(P[sel])]
    a, b, c = np.log(P[i - 1:i + 2])
    return float(fr[i] + 0.5 * (a - c) / (a - 2 * b + c) * (fr[1] - fr[0]))


def _project(p, t, f):
    """Complex amplitude of frequency f in each record (Hann window, amplitude-normalized)."""
    w = np.hanning(p.shape[-1])
    e = np.exp(-2j * np.pi * f * (t - t[..., :1]))
    return 2 * (p * w * e).sum(-1) / w.sum()


def _convection_speed(coeffs, offsets) -> float:
    """Slowness [s/m] of the street along the body from the phase of each harmonic.

    coeffs : list of (c, f, noise) with c the (n_stations, n_sensors) complex
             amplitudes and noise the noise amplitude near f; each station is
             weighted by its power above the noise.
    The phase falls as -2 pi f s x along the body; a weighted fit over all stations
    and harmonics gives s. Positive s: the street travels from head to tail.
    """
    x = offsets
    num = den = 0.0
    for c, fh, noise in coeffs:
        ref = c[:, np.argmin(np.abs(x))][:, None]
        ph = np.angle(c * np.conj(ref))                   # phase relative to mid-line sensor
        w = np.maximum(np.mean(np.abs(c) ** 2, 1, keepdims=True) - noise ** 2, 0.0) + 1e-30
        w = np.broadcast_to(w, c.shape)
        xm = (w * x).sum(1, keepdims=True) / w.sum(1, keepdims=True)
        pm = (w * ph).sum(1, keepdims=True) / w.sum(1, keepdims=True)
        omega = 2 * np.pi * fh
        num -= omega * (w * (x - xm) * (ph - pm)).sum()
        den += omega ** 2 * (w * (x - xm) ** 2).sum()
    return num / den


def _parabolic(v, k, grid):
    if 0 < k < len(v) - 1 and np.isfinite(v[k - 1:k + 2]).all():
        a, b, c = v[k - 1:k + 2]
        den = a - 2 * b + c
        if den > 0:
            return grid[k] + 0.5 * (a - c) / den * (grid[1] - grid[0])
    return grid[k]


def locate_rock(survey: Survey, strouhal: float = 0.2, convection: float = 0.88,
                spacing_ratio: float = 0.281, formation: float = 1.5, core: float = 0.1,
                eddy_viscosity: float = 0.01, protrusion: float = 0.5,
                manning_n: float = 0.035, friction: float | None = None,
                dx_range=(1.0, 20.0), n_dx: int = 96, min_snr: float = 3.0,
                rho: float = RHO) -> RockFix:
    """Locate the rock from a lateral survey (steps 1 to 5 of the module docstring).

    strouhal ... eddy_viscosity : the wake model the fish assumes
    protrusion                  : fraction of the rock's height above water; with D it
                                  gives the depth at the rock (D = 3.83 h for one half)
    manning_n                   : roughness the fish assumes; gives the bed friction
                                  c_f = g n^2 / h^(1/3) that damps the street
    friction                    : c_f given directly instead (0: no damping)
    dx_range                    : distances searched, from the rock to the sensor line

    The strongest spectral line is either f or 2 f. Both readings are fitted; a
    reading counts only if its line at f stands out of the noise (peak amplitude
    at least `min_snr` times the noise near f, which a reading of noise as f fails),
    and of those the one whose lateral profiles match the model better is kept.
    If neither does, the street is not detected (RockFix.detected is False).
    """
    fp = _spectral_peak(survey)
    p = survey.p - survey.p.mean(-1, keepdims=True)
    t = survey.t[:, None, :]
    fits = [_fit(survey, p, t, f, strouhal, convection, spacing_ratio, formation, core,
                 eddy_viscosity, protrusion, manning_n, friction, dx_range, n_dx, rho)
            for f in (fp, 0.5 * fp)]
    ok = [fx for fx in fits if fx.snr >= min_snr]
    if not ok:
        return replace(fits[0], detected=False)
    return min(ok, key=lambda fx: fx.cost)


def _fit(survey, p, t, f, strouhal, convection, spacing_ratio, formation, core,
         eddy_viscosity, protrusion, manning_n, friction, dx_range, n_dx, rho) -> RockFix:
    c1, c2 = _project(p, t, f), _project(p, t, 2 * f)                  # (n_st, n_sens)
    A1 = np.sqrt(np.mean(np.abs(c1) ** 2, 1))
    A2 = np.sqrt(np.mean(np.abs(c2) ** 2, 1))
    # noise amplitude of each harmonic, from frequencies between the lines
    s1 = np.median(np.sqrt(np.mean(np.abs(_project(p, t, 1.5 * f)) ** 2, 1)))
    s2 = np.median(np.sqrt(np.mean(np.abs(_project(p, t, 2.5 * f)) ** 2, 1)))

    # 2 and 3: convection speed, stream velocity and size
    slowness = _convection_speed([(c1, f, s1), (c2, 2 * f, s2)], survey.offsets)
    Uc = 1.0 / abs(slowness)
    U = Uc / convection
    D = strouhal * U / f

    # 4 and 5: fit the model's lateral profiles for lateral position and distance.
    # The fish's model river: depth from D and the protrusion, no banks.
    unit = HemisphereRock(0.0, 0.0, 1 / (1 - protrusion))
    h = D / unit.shedding_diameter(1.0)
    cf = G * manning_n ** 2 / h ** (1 / 3) if friction is None else friction
    model = KarmanStreet(TrapezoidalChannel(bed_width=1e3, depth=h),
                         HemisphereRock(0.0, 0.0, h / (1 - protrusion)),
                         strouhal=strouhal, convection=convection, spacing_ratio=spacing_ratio,
                         formation=formation, core=core, eddy_viscosity=eddy_viscosity,
                         n_images=0, rho=rho, velocity=U, friction=cf)

    ys = survey.y
    span = ys.max() - ys.min()
    y_fine = np.linspace(-span - 1.0, span + 1.0, int(round((2 * span + 2) / 0.02)) + 1)
    dx_grid = np.linspace(*dx_range, n_dx)
    yc_grid = np.linspace(ys.min(), ys.max(), int(round(span / 0.01)) + 1)
    rel = ys[None, :] - yc_grid[:, None]                                # (n_yc, n_st)
    cost = np.empty((dx_grid.size, yc_grid.size))
    for i, dx in enumerate(dx_grid):
        m1, m2 = model.harmonic_profile(dx, y_fine)
        r1 = (A1 - np.sqrt(np.interp(rel, y_fine, m1) ** 2 + s1 ** 2)) / s1
        r2 = (A2 - np.sqrt(np.interp(rel, y_fine, m2) ** 2 + s2 ** 2)) / s2
        cost[i] = (r1 ** 2).sum(1) + (r2 ** 2).sum(1)
    i, j = np.unravel_index(np.argmin(cost), cost.shape)
    dx = _parabolic(cost[:, j], i, dx_grid)
    yc = _parabolic(cost[i, :], j, yc_grid)
    return RockFix(x=survey.x_center - dx, y=yc, distance=dx, diameter=D, f=f, Uc=Uc, U=U,
                   upstream=bool(slowness > 0), cost=float(cost[i, j]),
                   snr=float(A1.max() / s1), A1=A1, A2=A2,
                   noise1=s1, noise2=s2, misfit=cost, dx_grid=dx_grid, yc_grid=yc_grid,
                   model=model)
