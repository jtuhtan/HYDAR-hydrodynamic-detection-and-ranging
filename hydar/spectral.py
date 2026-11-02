"""Dispersive 2D wave solver: linear potential flow with exact dispersion.

The shallow-water Solver (hydar.solver) is non-dispersive: every wave moves
at sqrt(g h). That is right for long waves, but a boat in water deeper than
a fraction of its wavelength makes a Kelvin wake, whose shape comes entirely
from dispersion (short waves are slower than long ones). This solver keeps
the full linear free-surface physics instead.

Model
-----
For inviscid, irrotational flow of uniform depth h, linearised about still
water, the surface elevation eta(x, y, t) and the velocity potential at the
surface phi(x, y, t) obey, in Fourier space (wavenumber k = |k|),

    d(eta_k)/dt = k tanh(kh) * phi_k
    d(phi_k)/dt = -g * eta_k - P_k / rho

where P(x, y, t) is a pressure applied on the surface. Every Fourier mode is
then a harmonic oscillator with the exact dispersion relation

    omega^2 = g k tanh(kh),

valid from shallow to deep water. The oscillator is integrated exactly over a
step (an exponential integrator), with P held at its mid-step value, so the
linear part has no stability limit and no numerical dispersion; dt is set
only by the forcing and by the sensors' sampling rate.

A moving boat is represented by the classical pressure patch: the hull
pushes on the water with P = rho * g * D(x, y, t), D being the hull draft, so
that a stationary hull would sink the surface by exactly D. This is the
standard linear model of ship waves (Havelock, Wehausen and Laitone).

Pressure below the surface follows from the same potential:

    p_dyn(z)/(rho g) = F^-1[ (eta_k + P_k/(rho g)) cosh k(h+z) / cosh kh ],

which sensors with pressure_model="solver" record directly.

Open boundaries
---------------
The FFT makes the domain periodic. Waves are absorbed in a sponge layer of
n_pml cells along every wall: eta and phi are damped there by
exp(-sigma(r) dt) with sigma growing quadratically into the layer, so that
waves leaving the interior die out before they can wrap around. The sponge
should be at least about one wavelength of the longest important wave.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import fft as sfft

from .solver import G, Grid

_WORKERS = -1


@dataclass
class SpectralSolver:
    """Linear potential-flow wave solver on a Grid (see module docstring).

    grid      : Grid; n_pml is the sponge-layer width in cells
    depth     : water depth h [m]
    dt        : time step [s] (default 0.01 s, i.e. 100 Hz)
    sigma_max : damping rate at the outer edge of the sponge [1/s]
    smooth    : Gaussian smoothing length of the applied pressure, in cells
    """

    grid: Grid
    depth: float
    dt: float = 0.01
    sigma_max: float = 3.0
    smooth: float = 1.0
    g: float = G
    rho: float = 1000.0
    t: float = field(init=False, default=0.0)
    objects: list = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        gr = self.grid
        kx = 2 * np.pi * sfft.fftfreq(gr.nx, gr.dx)
        ky = 2 * np.pi * sfft.rfftfreq(gr.ny, gr.dx)
        KX, KY = np.meshgrid(kx, ky, indexing="ij")
        self.k = np.hypot(KX, KY)
        self.K = self.k * np.tanh(self.k * self.depth)           # k tanh kh
        self.omega = np.sqrt(self.g * self.K)
        w = self.omega * self.dt
        self._cos, self._sin = np.cos(w), np.sin(w)
        with np.errstate(divide="ignore", invalid="ignore"):
            self._K_over_w = np.where(self.omega > 0, self.K / self.omega, 0.0)
            self._w_over_K = np.where(self.omega > 0, self.g / self.omega, 0.0)
        self._filter = np.exp(-0.5 * (self.k * self.smooth * gr.dx) ** 2)
        self._damp = np.exp(-self._sponge() * self.dt)
        self.c = float(np.sqrt(self.g * self.depth))              # longest-wave speed
        self._eta = np.zeros((gr.nx, gr.ny))
        self._phi = np.zeros((gr.nx, gr.ny))
        self._P = np.zeros((gr.nx, gr.ny))                         # P/(rho g) [m]
        self._a_hat = None                                         # last eta_k + P_k/(rho g)
        self._zcache: dict = {}

    def _sponge(self) -> np.ndarray:
        gr, n = self.grid, self.grid.n_pml
        if n == 0:
            return np.zeros((gr.nx, gr.ny))

        def ramp(m):
            i = np.arange(m)
            d = np.maximum(np.maximum(n - i - 0.5, i + 0.5 - (m - n)), 0) / n
            return d

        dx, dy = np.meshgrid(ramp(gr.nx), ramp(gr.ny), indexing="ij")
        return self.sigma_max * np.maximum(dx, dy) ** 2

    # ---------------------------------------------------------- state
    @property
    def eta(self) -> np.ndarray:
        """Surface elevation at cell centres, shape (nx, ny), incl. sponge."""
        return self._eta

    def set_eta(self, eta: np.ndarray) -> None:
        """Set an initial surface elevation (full grid), fluid at rest."""
        self._eta = np.array(eta, dtype=float)
        self._phi[...] = 0.0

    @property
    def surface_pressure_head(self) -> np.ndarray:
        """Applied surface pressure P/(rho g) [m] at the current time."""
        return self._P

    def add(self, obj) -> None:
        """Attach an object.

        Forcing objects (e.g. hydar.Boat) provide pressure_head(solver, t)
        returning (slices, field) of P/(rho g) on a sub-box of the grid, or None.
        Passive objects (sensors) provide source(solver) and advance(solver).
        """
        self.objects.append(obj)

    # ---------------------------------------------------------- forcing
    def _pressure(self, t: float) -> np.ndarray:
        P = np.zeros((self.grid.nx, self.grid.ny))
        for obj in self.objects:
            if hasattr(obj, "pressure_head"):
                out = obj.pressure_head(self, t)
                if out is not None:
                    sl, val = out
                    P[sl] += val
        return P

    # ---------------------------------------------------------- stepping
    def step(self) -> None:
        """Advance one time step dt."""
        for obj in self.objects:
            if hasattr(obj, "source"):
                obj.source(self)          # lets sensors attach; no volume sources here
        dt = self.dt
        Pm = self._pressure(self.t + 0.5 * dt)
        Pm_hat = sfft.rfft2(Pm, workers=_WORKERS) * self._filter
        a = sfft.rfft2(self._eta, workers=_WORKERS) + Pm_hat
        f = sfft.rfft2(self._phi, workers=_WORKERS)
        a_new = a * self._cos + f * self._K_over_w * self._sin
        f_new = f * self._cos - a * self._w_over_K * self._sin
        shape = (self.grid.nx, self.grid.ny)
        self._eta = sfft.irfft2(a_new - Pm_hat, s=shape, workers=_WORKERS) * self._damp
        self._phi = sfft.irfft2(f_new, s=shape, workers=_WORKERS) * self._damp
        P_end = self._pressure(self.t + dt)
        self._P = sfft.irfft2(sfft.rfft2(P_end, workers=_WORKERS) * self._filter, s=shape,
                              workers=_WORKERS)
        self._a_hat = None
        self._zcache = {}
        self.t += dt
        for obj in self.objects:
            if hasattr(obj, "advance"):
                obj.advance(self)

    def pressure_head(self, z: float) -> np.ndarray:
        """Dynamic pressure head p_dyn/(rho g) [m] at depth z (< 0) on the grid."""
        if z in self._zcache:
            return self._zcache[z]
        if self._a_hat is None:
            self._a_hat = sfft.rfft2(self._eta + self._P, workers=_WORKERS)
        kh = self.k * self.depth
        # cosh k(h+z)/cosh kh, overflow-safe
        Kz = np.exp(self.k * z) * (1 + np.exp(-2 * (kh + self.k * z))) / (1 + np.exp(-2 * kh))
        out = sfft.irfft2(self._a_hat * Kz, s=(self.grid.nx, self.grid.ny), workers=_WORKERS)
        self._zcache[z] = out
        return out

    def run(self, t_end: float, snapshot_every: int | None = None, callback=None) -> dict:
        """Integrate to time t_end; optionally keep interior eta every N steps."""
        n = max(int(np.ceil((t_end - self.t) / self.dt - 1e-9)), 0)
        snaps, snap_t = [], []
        si, sj = self.grid.interior
        for k in range(n):
            self.step()
            if snapshot_every and (k + 1) % snapshot_every == 0:
                snaps.append(self._eta[si, sj].copy())
                snap_t.append(self.t)
            if callback:
                callback(self)
        return {"snapshots": snaps, "snapshot_t": np.array(snap_t)}

    def energy(self, interior_only: bool = True) -> float:
        """Wave energy [J]: potential 0.5 rho g eta^2 plus kinetic 0.5 rho phi K phi."""
        f = sfft.rfft2(self._phi, workers=_WORKERS)
        Kphi = sfft.irfft2(self.K * f, s=self._phi.shape, workers=_WORKERS)
        e = 0.5 * self.rho * (self.g * self._eta ** 2 + self._phi * Kphi)
        if interior_only:
            e = e[self.grid.interior]
        return float(e.sum() * self.grid.dx ** 2)

    def volume(self) -> float:
        """Integral of eta over the interior [m^3]."""
        return float(self._eta[self.grid.interior].sum() * self.grid.dx ** 2)

