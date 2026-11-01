"""Drop a sphere into the water and let it generate surface waves.

The sphere moves vertically only. Its centre height z (z = 0 is the still
water surface, positive up) obeys

    (m + m_a) dw/dt = rho_w * g * V_sub - m * g - 0.5 * rho_w * Cd * A * w|w|
    dz/dt = w

with m = rho_s * V the sphere mass, V_sub the volume below the local water
surface, m_a = Ca * rho_w * V_sub the added mass (Ca = 0.5 for a sphere) and A
the frontal area of the submerged part. A sphere lighter than water settles
at a floating draft after a damped heave oscillation; a heavier one sinks
until it rests on the bottom (z = -h + R).

Coupling to the wave solver: while the sphere crosses the surface its
submerged volume changes, and exactly that volume of water is pushed
outward. The rate dV_sub/dt is injected into the continuity equation as a
source q(x, y) spread over the waterline disc, so the total volume of eta
equals V_sub. The local mean eta under the sphere is fed back into V_sub,
so the floating sphere feels the waves it radiates (two-way coupling).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np


def cap_volume(s: float, R: float) -> float:
    """Volume of a sphere of radius R submerged to depth s (0 <= s <= 2R)."""
    s = min(max(s, 0.0), 2 * R)
    return np.pi * s * s * (3 * R - s) / 3.0


def waterline_radius(s: float, R: float) -> float:
    """Radius of the circle where the water surface cuts the sphere."""
    s = min(max(s, 0.0), 2 * R)
    return np.sqrt(max(s * (2 * R - s), 0.0))


@dataclass
class DropObject:
    """A sphere released from rest above (or at) the water surface.

    diameter : sphere diameter D [m]
    density  : sphere density [kg/m^3] (water is 1000: lighter floats)
    x, y     : horizontal drop location [m] in interior coordinates
    z0       : initial height of the sphere's lowest point above the
               surface [m]; 0 means it starts just touching the water
    t_release: time [s] at which the sphere is let go (allows several drops
               at different times in one simulation)
    Cd, Ca   : drag and added-mass coefficients
    """

    diameter: float
    density: float
    x: float
    y: float
    z0: float = 0.0
    t_release: float = 0.0
    Cd: float = 0.5
    Ca: float = 0.5
    rho_w: float = 1000.0
    history: dict = field(init=False, default_factory=lambda: {"t": [], "z": [], "w": [], "V_sub": []})

    def __post_init__(self) -> None:
        self.R = 0.5 * self.diameter
        self.V = 4.0 / 3.0 * np.pi * self.R ** 3
        self.m = self.density * self.V
        self.z = self.z0 + self.R  # centre height
        self.w = 0.0
        self.V_sub = 0.0
        self.on_bottom = False
        self._s = 0.0
        self._mask = None

    # ---------------------------------------------------------- set-up
    def attach(self, solver) -> None:
        """Pre-compute geometry on the solver grid and check resolution."""
        gr = solver.grid
        X, Y = np.meshgrid(gr.x, gr.y, indexing="ij")
        self._r2 = (X - self.x) ** 2 + (Y - self.y) ** 2
        near = self._r2 <= (self.R + 2 * gr.dx) ** 2
        self._win = np.nonzero(near.any(axis=1))[0], np.nonzero(near.any(axis=0))[0]
        i, j = self._win
        self._sl = (slice(i[0], i[-1] + 1), slice(j[0], j[-1] + 1))
        self._r2 = self._r2[self._sl]
        self._q = np.zeros_like(solver.eta)

        cells = self.diameter / gr.dx
        if cells < 8:
            warnings.warn(f"sphere spans only {cells:.1f} cells; use dx <= D/8 "
                          f"(dx <= {self.diameter / 8:.4g} m) to resolve the waves")
        if self.density > self.rho_w and solver.depth < self.diameter:
            warnings.warn("water is shallower than the sphere; it will rest on the "
                          "bottom partly emerged")
        lx, ly = gr.nx_in * gr.dx, gr.ny_in * gr.dx
        if not (self.R <= self.x <= lx - self.R and self.R <= self.y <= ly - self.R):
            warnings.warn("drop location is too close to (or outside) the domain walls")

    def _footprint(self, a: float, dx: float) -> np.ndarray:
        """Smooth unit-integral weight over the waterline disc of radius a."""
        a_eff = max(a, 1.5 * dx)  # keep at least a few cells wide
        wgt = 0.5 * (1.0 - np.tanh((np.sqrt(self._r2) - a_eff) / (0.75 * dx)))
        return wgt / (wgt.sum() * dx * dx)

    def _local_level(self, solver) -> float:
        """Mean water elevation in a ring just outside the sphere."""
        eta = solver.eta[self._sl]
        dx = solver.grid.dx
        ring = (self._r2 > self.R ** 2) & (self._r2 <= (self.R + 2 * dx) ** 2)
        return float(eta[ring].mean()) if ring.any() else 0.0

    # ---------------------------------------------------------- coupling
    def source(self, solver):
        """Volume source q [m/s] from the change of submerged volume."""
        if self._mask is None:
            self.attach(solver)
            self._mask = True
        if solver.t < self.t_release:
            return None
        level = self._local_level(solver)
        s = level - (self.z - self.R)  # submergence depth of lowest point
        self._s = s
        V_new = cap_volume(s, self.R)
        dV = V_new - self.V_sub
        self.V_sub = V_new
        self._q[...] = 0.0
        if dV != 0.0:
            a = waterline_radius(s, self.R)
            self._q[self._sl] = (dV / solver.dt) * self._footprint(a, solver.grid.dx)
        return self._q

    def advance(self, solver) -> None:
        """Integrate the sphere's heave motion over one solver step."""
        dt, g = solver.dt, solver.g
        if solver.t < self.t_release:
            return
        if not self.on_bottom:
            # semi-implicit Euler with quadratic drag treated implicitly
            m_tot = self.m + self.Ca * self.rho_w * self.V_sub
            s = self._s
            A = np.pi * (self.R ** 2 if s >= self.R else waterline_radius(s, self.R) ** 2)
            F = self.rho_w * g * self.V_sub - self.m * g
            k = 0.5 * self.rho_w * self.Cd * A / m_tot
            self.w = (self.w + dt * F / m_tot) / (1.0 + dt * k * abs(self.w))
            self.z += dt * self.w
            z_bed = -solver.depth + self.R
            if self.z <= z_bed:
                self.z, self.w, self.on_bottom = z_bed, 0.0, True
        h = self.history
        h["t"].append(solver.t)
        h["z"].append(self.z)
        h["w"].append(self.w)
        h["V_sub"].append(self.V_sub)

    # ---------------------------------------------------------- helpers
    @property
    def equilibrium_draft(self) -> float | None:
        """Submerged depth of the lowest point when floating (None if it sinks)."""
        if self.density >= self.rho_w:
            return None
        target = self.density / self.rho_w * self.V
        lo, hi = 0.0, 2 * self.R
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if cap_volume(mid, self.R) < target else (lo, mid)
        return 0.5 * (lo + hi)


def recommend_grid(diameter: float, depth: float, g: float = 9.81,
                   cells_per_diameter: int = 10) -> dict:
    """Suggest dx and a domain size that resolve the waves from a drop.

    The wave pulse length is set by the sphere size and by the entry time
    T ~ sqrt(D/g); the domain should hold the wave front for a few T.
    """
    c = np.sqrt(g * depth)
    T = np.sqrt(diameter / g)
    dx = diameter / cells_per_diameter
    L = max(20 * diameter, 2 * c * 10 * T)
    return {"dx": float(dx), "domain": float(L), "wave_speed": float(c),
            "entry_time": float(T)}
