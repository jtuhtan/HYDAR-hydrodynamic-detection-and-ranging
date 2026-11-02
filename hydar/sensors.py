"""Pressure sensor arrays.

A SensorArray is a regular nx-by-ny grid of pressure sensors with a known
spacing, centred at (x, y) and mounted at depth z below the still-water
surface (z is negative, -depth < z <= 0). Each sensor samples at `fs` Hz.

Pressure at depth (linear wave theory)
--------------------------------------
The solver gives the surface elevation eta(x, y, t). Gauge pressure at a
sensor at depth z is

    p(t) = rho * g * (-z) + p_dyn(t)

where the static part is the hydrostatic head and the dynamic part follows
linear (Airy) wave theory: every frequency component of eta is scaled by
the pressure response factor

    K_p(omega) = cosh(k (h + z)) / cosh(k h),   omega**2 = g k tanh(k h)

so short waves are attenuated with depth and long waves are not. In the
shallow-water limit (k h -> 0) K_p -> 1 and the pressure is hydrostatic.

Two pressure models are available (SensorArray.pressure_model):

"hydrostatic" (default)  p_dyn = rho g eta. Consistent with the shallow-water
                         solver, and the right choice for time-of-flight
                         picking.
"airy"                   p_dyn from the K_p filter above, applied with an FFT.
                         Shows how strongly each array's depth attenuates
                         short waves. The filter is zero-phase, so it smears
                         a sharp front slightly earlier in time; do not use
                         it for arrival picking.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def wavenumber(omega: np.ndarray, h: float, g: float = 9.81) -> np.ndarray:
    """Solve the linear dispersion relation omega^2 = g k tanh(k h) for k."""
    omega = np.abs(np.asarray(omega, dtype=float))
    k = np.maximum(omega / np.sqrt(g * h), omega ** 2 / g)  # good first guess
    k = np.where(k == 0, 1e-12, k)
    for _ in range(50):  # Newton iterations
        th = np.tanh(k * h)
        f = g * k * th - omega ** 2
        df = g * th + g * k * h * (1 - th ** 2)
        k = np.maximum(k - f / df, 1e-12)
    return np.where(omega == 0, 0.0, k)


def depth_response(eta: np.ndarray, fs: float, z: float, h: float,
                   g: float = 9.81, rho: float = 1000.0) -> np.ndarray:
    """Dynamic pressure [Pa] at depth z from a surface-elevation record."""
    eta = np.asarray(eta, dtype=float)
    n = eta.shape[-1]
    omega = 2 * np.pi * np.fft.rfftfreq(n, d=1.0 / fs)
    k = wavenumber(omega, h, g)
    # cosh(k(h+z))/cosh(kh) written to avoid overflow for large k h
    Kp = np.exp(k * z) * (1 + np.exp(-2 * k * (h + z))) / (1 + np.exp(-2 * k * h))
    return rho * g * np.fft.irfft(np.fft.rfft(eta, axis=-1) * Kp, n=n, axis=-1)


@dataclass
class SensorArray:
    """Regular grid of pressure sensors.

    name     : label used in plots and results
    x, y     : horizontal position of the array centre [m]
    z        : sensor depth [m], negative below the still-water surface
    nx, ny   : number of sensors along x and y (3 x 3 by default)
    spacing  : distance between neighbouring sensors [m]
    fs       : sampling rate [Hz]
    pressure_model : "hydrostatic" or "airy" (see module docstring)
    rotation_deg   : rotation of the grid, degrees counter-clockwise
    """

    name: str
    x: float
    y: float
    z: float = -0.5
    nx: int = 3
    ny: int = 3
    spacing: float = 1.0
    fs: float = 100.0
    pressure_model: str = "hydrostatic"
    rotation_deg: float = 0.0
    _t: list = field(init=False, default_factory=list, repr=False)
    _eta: list = field(init=False, default_factory=list, repr=False)

    def __post_init__(self) -> None:
        ox = (np.arange(self.nx) - (self.nx - 1) / 2) * self.spacing
        oy = (np.arange(self.ny) - (self.ny - 1) / 2) * self.spacing
        OX, OY = np.meshgrid(ox, oy, indexing="ij")
        r = np.radians(self.rotation_deg)
        R = np.array([[np.cos(r), -np.sin(r)], [np.sin(r), np.cos(r)]])
        #: sensor offsets from the centre, shape (n_sensors, 2)
        self.offsets = np.column_stack([OX.ravel(), OY.ravel()]) @ R.T
        #: absolute sensor positions, shape (n_sensors, 2)
        self.positions = self.offsets + [self.x, self.y]
        self._w = None

    @property
    def center(self) -> np.ndarray:
        return np.array([self.x, self.y])

    @property
    def n_sensors(self) -> int:
        return len(self.positions)

    # ------------------------------------------------ solver coupling
    def _attach(self, solver) -> None:
        gr = solver.grid
        if not -solver.depth < self.z <= 0:
            raise ValueError(f"{self.name}: z must lie in (-{solver.depth}, 0]")
        if solver.dt > 1.0 / self.fs:
            raise ValueError(f"{self.name}: solver dt {solver.dt:.4g}s is longer than "
                             f"the sampling interval 1/fs; refine dx")
        # bilinear interpolation weights on cell centres
        fx = (self.positions[:, 0] - gr.x[0]) / gr.dx
        fy = (self.positions[:, 1] - gr.y[0]) / gr.dx
        i0, j0 = np.floor(fx).astype(int), np.floor(fy).astype(int)
        if (i0 < gr.n_pml).any() or (i0 >= gr.n_pml + gr.nx_in - 1).any() or \
           (j0 < gr.n_pml).any() or (j0 >= gr.n_pml + gr.ny_in - 1).any():
            raise ValueError(f"{self.name}: a sensor lies outside the interior domain")
        ax, ay = fx - i0, fy - j0
        self._idx = (i0, j0)
        self._w = ((1 - ax) * (1 - ay), ax * (1 - ay), (1 - ax) * ay, ax * ay)
        self._depth, self._g, self._rho = solver.depth, solver.g, solver.rho

    def source(self, solver):  # sensors are passive
        if self._w is None:
            self._attach(solver)
        return None

    def advance(self, solver) -> None:
        eta = solver.eta
        i, j = self._idx
        w00, w10, w01, w11 = self._w
        val = w00 * eta[i, j] + w10 * eta[i + 1, j] + w01 * eta[i, j + 1] + w11 * eta[i + 1, j + 1]
        self._t.append(solver.t)
        self._eta.append(val)

    # ------------------------------------------------ recorded data
    @property
    def time(self) -> np.ndarray:
        """Sample times at fs [s]."""
        if not self._t:
            return np.empty(0)
        return np.arange(0.0, self._t[-1] + 1e-12, 1.0 / self.fs)

    @property
    def eta(self) -> np.ndarray:
        """Surface elevation above each sensor, shape (n_sensors, n_samples)."""
        t_raw = np.concatenate([[0.0], self._t])
        e_raw = np.vstack([np.zeros(self.n_sensors), np.array(self._eta)])
        return np.array([np.interp(self.time, t_raw, e_raw[:, k])
                         for k in range(self.n_sensors)])

    @property
    def dynamic_pressure(self) -> np.ndarray:
        """Wave-induced pressure at the sensor depth [Pa], (n_sensors, n_samples)."""
        if self.pressure_model == "hydrostatic":
            return self._rho * self._g * self.eta
        if self.pressure_model == "airy":
            return depth_response(self.eta, self.fs, self.z, self._depth, self._g, self._rho)
        raise ValueError("pressure_model must be 'hydrostatic' or 'airy'")

    @property
    def pressure(self) -> np.ndarray:
        """Total gauge pressure (hydrostatic + dynamic) [Pa]."""
        return self._rho * self._g * (-self.z) + self.dynamic_pressure

    def load(self, t: np.ndarray, eta: np.ndarray, depth: float, g: float = 9.81,
             rho: float = 1000.0) -> None:
        """Fill the array with an externally computed surface elevation record.

        t   : (n,) times [s], starting at or after 0
        eta : (n_sensors, n) surface elevation above each sensor [m]
        Used for wave fields that do not come from the Solver, e.g. hydar.kelvin.
        """
        if not -depth < self.z <= 0:
            raise ValueError(f"{self.name}: z must lie in (-{depth}, 0]")
        eta = np.asarray(eta, dtype=float)
        if eta.shape[0] != self.n_sensors:
            raise ValueError("eta must have one row per sensor")
        self._t = list(np.asarray(t, dtype=float))
        self._eta = list(eta.T)
        self._depth, self._g, self._rho = depth, g, rho
        self._w = ()  # mark as filled (no solver coupling)

    def to_dict(self) -> dict:
        """Recorded data as plain arrays (e.g. for np.savez)."""
        return {"time": self.time, "pressure": self.pressure,
                "positions": self.positions, "z": self.z, "fs": self.fs}
