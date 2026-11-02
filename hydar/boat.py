"""A small boat travelling through the domain.

Hull geometry: Series 60
------------------------
The default hull follows the Series 60 parent form (Todd, 1963; block
coefficient C_B = 0.60), the standard systematic hull series used in ship
hydrodynamics:

    L/B = 7.5,  B/T = 2.5,  C_B = 0.60,  C_P = 0.614,  C_M = 0.977,
    C_WP = 0.706,  LCB = 1.5 % L aft of midships

The published offset table is not reproduced here. Instead the hull is built
from smooth curves fitted to exactly these coefficients:

* sectional-area curve  A(xi) = A_M (1 - |xi|^n), with separate exponents
  for the fore and aft body so that C_P = 0.614 and the LCB lies 1.5 % L aft;
* design waterline      b(xi) = B/2 (1 - |xi|^m), with m set by C_WP = 0.706;
* each section is given the depth profile D = T (1 - |y/b|^k), with k chosen
  so the section has the area A(xi). Midships this gives C_M = 0.977.

xi runs from -1 (stern) to +1 (bow). The result has the Series 60 principal
dimensions, displacement, prismatic distribution and LCB; the finer details
of the section shapes are approximate. Hull.coefficients() reports the values
actually achieved on the grid.

Coupling to the solvers
-----------------------
Dispersive hydar.spectral.SpectralSolver (recommended, used in case study 2):
the hull pushes on the surface as a pressure patch, P/(rho g) = D(x, y, t),
the hull draft below the still surface (Boat.pressure_head). The wake is the
full linear Kelvin pattern: transverse and diverging waves inside a wedge of
half-angle 19.47 deg. HYDAR tracks it with hydar.locate.track_boat.

Shallow-water hydar.solver.Solver: the rate at which the hull displaces water
is fed in as the volume source q = dD/dt (Boat.source). This model is
non-dispersive, so with depth Froude number Fr = U / sqrt(g h):

* Fr > 1 (supercritical): the boat drags a Mach cone of half-angle
  arcsin(1/Fr); track it with track_boat(..., wake="mach").
* Fr < 1 (subcritical): the boat only carries a local depression with it and
  radiates no wake. Use the SpectralSolver instead.
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np
from scipy.optimize import fsolve

SERIES60 = {"L/B": 7.5, "B/T": 2.5, "CB": 0.60, "CP": 0.614, "CM": 0.977,
            "CWP": 0.706, "LCB": -0.015}  # LCB as a fraction of L, negative = aft


@lru_cache(maxsize=None)
def _exponents(CP: float, CWP: float, LCB: float) -> tuple[float, float, float]:
    """Fore/aft sectional-area exponents and waterline exponent."""
    lcb = 2 * LCB  # in half-lengths

    def eqs(v):
        nf, na = v
        cp = 0.5 * (nf / (nf + 1) + na / (na + 1))
        mom = (0.5 - 1 / (nf + 2)) - (0.5 - 1 / (na + 2))
        return [cp - CP, mom / (2 * cp) - lcb]

    nf, na = fsolve(eqs, [1.5, 1.5])
    m = CWP / (1 - CWP)
    return float(nf), float(na), float(m)


@dataclass
class Series60Hull:
    """Series 60 (C_B = 0.60) hull scaled to length L [m].

    beam and draft default to L/7.5 and beam/2.5.
    """

    length: float = 4.5
    beam: float | None = None
    draft: float | None = None
    CP: float = SERIES60["CP"]
    CM: float = SERIES60["CM"]
    CWP: float = SERIES60["CWP"]
    LCB: float = SERIES60["LCB"]

    def __post_init__(self) -> None:
        self.beam = self.beam or self.length / SERIES60["L/B"]
        self.draft = self.draft or self.beam / SERIES60["B/T"]
        self.nf, self.na, self.m = _exponents(self.CP, self.CWP, self.LCB)

    def depth(self, xi, zeta):
        """Hull depth below the waterline [m].

        xi   : longitudinal coordinate, -1 at the stern, +1 at the bow
        zeta : transverse coordinate in metres from the centreline
        """
        xi, zeta = np.broadcast_arrays(np.asarray(xi, float), np.asarray(zeta, float))
        a = np.abs(xi)
        inside = a < 1
        n = np.where(xi >= 0, self.nf, self.na)
        area = np.where(inside, 1 - np.minimum(a, 1) ** n, 0.0)          # A / A_M
        half_b = 0.5 * self.beam * np.where(inside, 1 - np.minimum(a, 1) ** self.m, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            s = np.clip(self.CM * area / np.where(half_b > 0, 2 * half_b / self.beam, 1), 0.05, 0.995)
            k = s / (1 - s)
            u = np.where(half_b > 0, np.abs(zeta) / half_b, 2.0)
        return np.where(inside & (u < 1), self.draft * (1 - np.minimum(u, 1) ** k), 0.0)

    def coefficients(self, n: int = 801) -> dict:
        """Form coefficients of the hull as built (numerical integration)."""
        xi = np.linspace(-1, 1, n)
        z = np.linspace(-self.beam / 2, self.beam / 2, n)
        XI, Z = np.meshgrid(xi, z, indexing="ij")
        D = self.depth(XI, Z)
        dA = (self.length / 2 * (xi[1] - xi[0])) * (z[1] - z[0])
        vol = D.sum() * dA
        sec = D.sum(axis=1) * (z[1] - z[0])
        x = xi * self.length / 2
        L, B, T = self.length, self.beam, self.draft
        return {"volume": vol, "CB": vol / (L * B * T), "CM": sec.max() / (B * T),
                "CP": vol / (sec.max() * L), "CWP": (D > 0).sum() * dA / (L * B),
                "LCB": float((sec * x).sum() / sec.sum() / L)}


@dataclass
class EllipseHull:
    """Simple smooth elliptical hull (cos^2 depth profile), for comparison."""

    length: float = 3.0
    beam: float = 1.2
    draft: float = 0.2

    def depth(self, xi, zeta):
        r = np.sqrt(np.asarray(xi) ** 2 + (np.asarray(zeta) / (0.5 * self.beam)) ** 2)
        return np.where(r < 1, self.draft * np.cos(0.5 * np.pi * r) ** 2, 0.0)

    def coefficients(self) -> dict:
        vol = self.draft * self.length * self.beam / 4 * (np.pi / 2 - 2 / np.pi)
        return {"volume": vol, "CB": vol / (self.length * self.beam * self.draft)}


@dataclass
class Boat:
    """Hull moving in a straight line at constant speed.

    hull        : Series60Hull (default, 4.5 m) or EllipseHull
    speed       : boat speed U [m/s]
    x, y        : position of the hull centre at time t_start [m]
    heading_deg : direction of travel, degrees counter-clockwise from +x
    t_start     : time [s] the boat starts moving
    ramp        : time [s] over which the draft grows from zero at the start
    """

    hull: object = field(default_factory=Series60Hull)
    speed: float = 5.0
    x: float = 0.0
    y: float = 0.0
    heading_deg: float = 0.0
    t_start: float = 0.0
    ramp: float = 1.0
    history: dict = field(init=False, default_factory=lambda: {"t": [], "x": [], "y": []})

    def __post_init__(self) -> None:
        th = np.radians(self.heading_deg)
        self.e = np.array([np.cos(th), np.sin(th)])          # along track (bow)
        self.n = np.array([-np.sin(th), np.cos(th)])         # to port
        self._ready = False

    @property
    def length(self) -> float:
        return self.hull.length

    @property
    def beam(self) -> float:
        return self.hull.beam

    @property
    def volume(self) -> float:
        """Displaced volume [m^3]."""
        return self.hull.coefficients()["volume"]

    def froude(self, depth: float, g: float = 9.81) -> float:
        """Depth Froude number U / sqrt(g h)."""
        return self.speed / np.sqrt(g * depth)

    def position(self, t) -> np.ndarray:
        """Hull centre at time(s) t [s]; shape (2,) or (len(t), 2)."""
        t = np.asarray(t, dtype=float)
        s = self.speed * np.clip(t - self.t_start, 0, None)
        return np.array([self.x, self.y]) + np.multiply.outer(s, self.e)

    def draft_field(self, X, Y, t):
        """Hull depth below the still surface at points (X, Y) and time t."""
        p = self.position(t)
        dx, dy = X - p[0], Y - p[1]
        xi = (dx * self.e[0] + dy * self.e[1]) / (0.5 * self.length)
        zeta = dx * self.n[0] + dy * self.n[1]
        amp = np.clip((t - self.t_start) / self.ramp, 0, 1) if self.ramp > 0 else 1.0
        return amp * self.hull.depth(xi, zeta)

    # ------------------------------------------------------------ coupling
    def _attach(self, solver) -> None:
        gr = solver.grid
        if self.beam / gr.dx < 8:
            warnings.warn(f"boat beam spans only {self.beam / gr.dx:.1f} cells; "
                          f"use dx <= beam/8 = {self.beam / 8:.3g} m")
        self._q = np.zeros_like(solver.eta)
        self._ready = True

    def source(self, solver):
        if not hasattr(solver, "u"):      # spectral solver: forced by pressure_head
            return None
        if not self._ready:
            self._attach(solver)
        t0, t1 = solver.t, solver.t + solver.dt
        if t1 <= self.t_start:
            return None
        gr = solver.grid
        # evaluate only on the box swept by the hull during this step
        p0, p1 = self.position(t0), self.position(t1)
        R = 0.5 * max(self.length, self.beam) + 2 * gr.dx
        lo, hi = np.minimum(p0, p1) - R, np.maximum(p0, p1) + R
        i = np.nonzero((gr.x >= lo[0]) & (gr.x <= hi[0]))[0]
        j = np.nonzero((gr.y >= lo[1]) & (gr.y <= hi[1]))[0]
        self._q[...] = 0.0
        if i.size == 0 or j.size == 0:
            return self._q
        sl = (slice(i[0], i[-1] + 1), slice(j[0], j[-1] + 1))
        X, Y = np.meshgrid(gr.x[sl[0]], gr.y[sl[1]], indexing="ij")
        self._q[sl] = (self.draft_field(X, Y, t1) - self.draft_field(X, Y, t0)) / solver.dt
        return self._q

    def pressure_head(self, solver, t):
        """Surface pressure head P/(rho g) = D for hydar.spectral.SpectralSolver.

        Returns (slices, values) on the box around the hull, or None.
        """
        if t <= self.t_start:
            return None
        gr = solver.grid
        p = self.position(t)
        R = 0.5 * max(self.length, self.beam) + 2 * gr.dx
        i = np.nonzero(np.abs(gr.x - p[0]) <= R)[0]
        j = np.nonzero(np.abs(gr.y - p[1]) <= R)[0]
        if i.size == 0 or j.size == 0:
            return None
        sl = (slice(i[0], i[-1] + 1), slice(j[0], j[-1] + 1))
        X, Y = np.meshgrid(gr.x[sl[0]], gr.y[sl[1]], indexing="ij")
        return sl, self.draft_field(X, Y, t)

    def advance(self, solver) -> None:
        p = self.position(solver.t)
        self.history["t"].append(solver.t)
        self.history["x"].append(p[0])
        self.history["y"].append(p[1])
