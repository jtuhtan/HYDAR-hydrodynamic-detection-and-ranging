"""The compiled kernels (Numba) and their NumPy fallbacks give bitwise identical results."""

import numpy as np
import pytest

from hydar import _kernels as K
from hydar.solver import Grid, Solver
from hydar.spectral import SpectralSolver


def _gauss_solver(**kw):
    sim = Solver(Grid(lx=2.0, ly=1.5, dx=0.02, n_pml=10), depth=0.5)
    X, Y = np.meshgrid(sim.grid.x, sim.grid.y, indexing="ij")
    sim.set_eta(0.01 * np.exp(-((X - 0.8) ** 2 + (Y - 0.6) ** 2) / 0.02))
    return sim


def test_numpy_fallback_matches_the_original_expressions():
    sim = _gauss_solver()
    ex, ey, u, v = (a.copy() for a in (sim.ex, sim.ey, sim.u, sim.v))
    dx, g, h = sim.grid.dx, sim.g, sim.depth
    kern = K.ShallowWaterKernel(sim.grid.nx, sim.grid.ny, sim._aux, sim._bux, sim._avy, sim._bvy,
                                sim._ax, sim._bx, sim._ay, sim._by, use_numba=False)
    for _ in range(20):
        # the update as it was written before the kernels (v1.3.0)
        eta = ex + ey
        u[1:-1, :] = sim._aux * u[1:-1, :] - sim._bux * g * (eta[1:] - eta[:-1]) / dx
        v[:, 1:-1] = sim._avy * v[:, 1:-1] - sim._bvy * g * (eta[:, 1:] - eta[:, :-1]) / dx
        ex = sim._ax * ex - sim._bx * h * (u[1:] - u[:-1]) / dx
        ey = sim._ay * ey - sim._by * h * (v[:, 1:] - v[:, :-1]) / dx
        kern(sim.ex, sim.ey, sim.u, sim.v, g, h, dx)
    for a, b in ((ex, sim.ex), (ey, sim.ey), (u, sim.u), (v, sim.v)):
        assert np.array_equal(a, b)


@pytest.mark.skipif(not K.HAVE_NUMBA, reason="numba not installed")
def test_numba_shallow_water_kernel_is_bitwise_identical():
    a, b = _gauss_solver(), _gauss_solver()
    b._kernel = K.ShallowWaterKernel(b.grid.nx, b.grid.ny, b._aux, b._bux, b._avy, b._bvy,
                                     b._ax, b._bx, b._ay, b._by, use_numba=False)
    assert a._kernel.numba and not b._kernel.numba
    for _ in range(50):
        a.step()
        b.step()
    for name in ("ex", "ey", "u", "v"):
        assert np.array_equal(getattr(a, name), getattr(b, name))


@pytest.mark.skipif(not K.HAVE_NUMBA, reason="numba not installed")
def test_numba_oscillator_is_bitwise_identical():
    rng = np.random.default_rng(3)
    sim = SpectralSolver(Grid(lx=8.0, ly=6.0, dx=0.25, n_pml=4), depth=2.0, dt=0.05)
    shape = sim.k.shape
    cplx = lambda: rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    e, f, p = cplx(), cplx(), cplx()
    args = (sim._cos, sim._sin, sim._K_over_w, sim._w_over_K)
    ref = K._oscillate_numpy(e, f, p, *args)
    out = K._oscillate_numba(e.copy(), f.copy(), p, *args)
    for r, o in zip(ref, out):
        assert np.array_equal(r, o)
