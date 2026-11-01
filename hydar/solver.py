"""2D linear shallow-water solver with non-reflecting (PML) walls.

The domain is a horizontal rectangle (x, y) of the free surface. Unknowns are
the surface elevation eta(x, y, t) [m] and the depth-averaged horizontal
velocity (u, v) [m/s]. For still water of uniform depth h:

    d(eta)/dt = -h * (du/dx + dv/dy) + q(x, y, t)
    du/dt     = -g * d(eta)/dx
    dv/dt     = -g * d(eta)/dy

q [m/s] is a volume source per unit area (used by objects that displace
water, see hydar.drop). Waves travel at c = sqrt(g*h). The model is valid for
waves much longer than the depth (non-dispersive, linear: amplitude << h).

The interior is wrapped in a split-field perfectly matched layer
(Berenger 1994): outgoing waves are absorbed with negligible reflection, so
the walls behave like open water. Discretisation is a staggered (Arakawa C /
Yee) grid, second order in space and time.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

G = 9.81


@dataclass
class Grid:
    """Rectangular domain.

    lx, ly : size of the physical (interior) region [m]; origin at its corner
    dx     : grid spacing [m], equal in x and y
    n_pml  : absorbing-layer thickness in cells, added outside every wall
    """

    lx: float
    ly: float
    dx: float
    n_pml: int = 40

    def __post_init__(self) -> None:
        self.nx_in = int(round(self.lx / self.dx))
        self.ny_in = int(round(self.ly / self.dx))
        self.nx = self.nx_in + 2 * self.n_pml
        self.ny = self.ny_in + 2 * self.n_pml
        # Cell-centre (eta) coordinates.
        self.x = (np.arange(self.nx) - self.n_pml + 0.5) * self.dx
        self.y = (np.arange(self.ny) - self.n_pml + 0.5) * self.dx

    def index(self, x: float, y: float) -> tuple[int, int]:
        """Nearest cell index (i, j) of physical position (x, y)."""
        return (int(np.argmin(np.abs(self.x - x))),
                int(np.argmin(np.abs(self.y - y))))

    @property
    def interior(self) -> tuple[slice, slice]:
        """Slices selecting the physical (non-PML) part of an eta-array."""
        return (slice(self.n_pml, self.n_pml + self.nx_in),
                slice(self.n_pml, self.n_pml + self.ny_in))

    @property
    def extent(self) -> list[float]:
        """[xmin, xmax, ymin, ymax] of the interior, for plt.imshow."""
        return [0.0, self.nx_in * self.dx, 0.0, self.ny_in * self.dx]


def pml_profile(n: int, n_pml: int, dx: float, c: float, order: int = 2,
                reflection: float = 1e-6, staggered: bool = False) -> np.ndarray:
    """Damping coefficient sigma [1/s] along one axis.

    sigma(d) = sigma_max * (d/L)**order, with sigma_max set for a theoretical
    normal-incidence reflection coefficient `reflection`.
    """
    if n_pml == 0:
        return np.zeros(n + (1 if staggered else 0))
    L = n_pml * dx
    sigma_max = -(order + 1) * c * np.log(reflection) / (2.0 * L)
    pos = np.arange(n + (1 if staggered else 0), dtype=float)
    if not staggered:
        pos += 0.5  # cell centres
    d = np.maximum(np.clip(n_pml - pos, 0, None), np.clip(pos - (n - n_pml), 0, None))
    return sigma_max * (d * dx / L) ** order


@dataclass
class Solver:
    """Time-domain shallow-water solver.

    grid  : Grid
    depth : still-water depth h [m]
    g     : gravity [m/s^2]
    cfl   : Courant number c*dt/dx, must be < 1/sqrt(2)
    """

    grid: Grid
    depth: float = 1.0
    g: float = G
    cfl: float = 0.5
    reflection: float = 1e-6
    rho: float = 1000.0
    t: float = field(init=False, default=0.0)

    def __post_init__(self) -> None:
        gr = self.grid
        if self.cfl >= 1.0 / np.sqrt(2.0):
            raise ValueError("cfl must be below 1/sqrt(2) on a 2D staggered grid")
        self.c = np.sqrt(self.g * self.depth)
        self.dt = self.cfl * gr.dx / self.c

        # eta is split into ex + ey inside the PML; velocities on cell faces.
        self.ex = np.zeros((gr.nx, gr.ny))
        self.ey = np.zeros((gr.nx, gr.ny))
        self.u = np.zeros((gr.nx + 1, gr.ny))
        self.v = np.zeros((gr.nx, gr.ny + 1))

        def ab(sigma):  # Crank-Nicolson damping factors
            den = 1.0 + 0.5 * self.dt * sigma
            return (1.0 - 0.5 * self.dt * sigma) / den, self.dt / den

        prof = lambda n, s: pml_profile(n, gr.n_pml, gr.dx, self.c,
                                        reflection=self.reflection, staggered=s)
        ax, bx = ab(prof(gr.nx, False))
        ay, by = ab(prof(gr.ny, False))
        aux, bux = ab(prof(gr.nx, True))
        avy, bvy = ab(prof(gr.ny, True))
        self._ax, self._bx = ax[:, None], bx[:, None]
        self._ay, self._by = ay[None, :], by[None, :]
        self._aux, self._bux = aux[1:-1, None], bux[1:-1, None]
        self._avy, self._bvy = avy[None, 1:-1], bvy[None, 1:-1]

        self.objects: list = []

    # ------------------------------------------------------------------
    @property
    def eta(self) -> np.ndarray:
        """Surface elevation at cell centres, shape (nx, ny), incl. PML."""
        return self.ex + self.ey

    def set_eta(self, eta: np.ndarray) -> None:
        """Set an initial surface elevation (full grid, incl. PML)."""
        self.ex = 0.5 * np.asarray(eta, dtype=float)
        self.ey = self.ex.copy()

    def add(self, obj) -> None:
        """Attach an object with methods `source(solver) -> q-array` and
        optionally `advance(solver)`; see hydar.drop.DropObject."""
        self.objects.append(obj)

    def step(self) -> None:
        """Advance one time step dt."""
        dx, dt = self.grid.dx, self.dt
        eta = self.ex + self.ey

        self.u[1:-1, :] = self._aux * self.u[1:-1, :] - self._bux * self.g * (eta[1:] - eta[:-1]) / dx
        self.v[:, 1:-1] = self._avy * self.v[:, 1:-1] - self._bvy * self.g * (eta[:, 1:] - eta[:, :-1]) / dx

        h = self.depth
        self.ex = self._ax * self.ex - self._bx * h * (self.u[1:] - self.u[:-1]) / dx
        self.ey = self._ay * self.ey - self._by * h * (self.v[:, 1:] - self.v[:, :-1]) / dx

        for obj in self.objects:
            q = obj.source(self)
            if q is not None:
                self.ex += 0.5 * dt * q
                self.ey += 0.5 * dt * q
        self.t += dt
        for obj in self.objects:
            if hasattr(obj, "advance"):
                obj.advance(self)

    def run(self, t_end: float, gauges=None, snapshot_every: int | None = None,
            callback=None) -> dict:
        """Integrate to time t_end.

        gauges         : list of (x, y) wave-gauge positions; eta recorded each step
        snapshot_every : keep a copy of the interior eta every N steps
        callback       : called as callback(solver) after every step

        Returns dict: 't', 'gauges' (n_steps, n_gauges), 'snapshots', 'snapshot_t'.
        """
        gauges = gauges or []
        idx = [self.grid.index(x, y) for x, y in gauges]
        n = max(int(np.ceil((t_end - self.t) / self.dt)), 0)
        times, rec = np.empty(n), np.empty((n, len(idx)))
        snaps, snap_t = [], []
        si, sj = self.grid.interior
        for k in range(n):
            self.step()
            times[k] = self.t
            eta = self.ex + self.ey
            for m, (i, j) in enumerate(idx):
                rec[k, m] = eta[i, j]
            if snapshot_every and (k + 1) % snapshot_every == 0:
                snaps.append(eta[si, sj].copy())
                snap_t.append(self.t)
            if callback:
                callback(self)
        return {"t": times, "gauges": rec, "snapshots": snaps,
                "snapshot_t": np.array(snap_t)}

    def energy(self, interior_only: bool = True) -> float:
        """Wave energy [J]: 0.5*rho*(h*|u|^2 + g*eta^2) integrated over area."""
        uc = 0.5 * (self.u[1:] + self.u[:-1])
        vc = 0.5 * (self.v[:, 1:] + self.v[:, :-1])
        e = 0.5 * self.rho * (self.depth * (uc ** 2 + vc ** 2) + self.g * self.eta ** 2)
        if interior_only:
            e = e[self.grid.interior]
        return float(e.sum() * self.grid.dx ** 2)

    def volume(self) -> float:
        """Net displaced water volume in the interior, integral of eta [m^3]."""
        return float(self.eta[self.grid.interior].sum() * self.grid.dx ** 2)
