"""Compiled inner loops, with pure-NumPy fallbacks that give the same numbers.

Numba is optional (``pip install hydar[fast]``). When it is installed the
kernels below are compiled to machine code on first use and cached on disk;
without it, or with the environment variable ``HYDAR_NO_NUMBA=1``, the NumPy
versions run instead. Both evaluate every expression in the same order as the
original NumPy code, without fast-math, so the results are bitwise identical
(checked in tests/test_kernels.py); only the speed differs.

Measured on a 460 x 460 grid (4 m at 1 cm with a 30-cell PML), one core:
NumPy as written before 8.1 ms per step, in-place NumPy 3.8 ms, Numba 0.7 ms.
A Rust (PyO3) version of the same loop ran at 0.7 ms too: the loop is memory
bound, so any compiled language lands at the same speed, and Numba needs no
compiler on the user's machine.
"""

from __future__ import annotations

import os

import numpy as np

try:
    if os.environ.get("HYDAR_NO_NUMBA", "") not in ("", "0"):
        raise ImportError("disabled by HYDAR_NO_NUMBA")
    import numba as _nb
except ImportError:  # pragma: no cover - depends on the environment
    _nb = None

HAVE_NUMBA = _nb is not None


def _sw_step_numpy(ex, ey, u, v, aux, bux, avy, bvy, ax, bx, ay, by, g, h, dx, work):
    """One shallow-water step in place (NumPy, no full-size temporaries but `work`)."""
    eta, du, dv, t = work
    np.add(ex, ey, out=eta)
    np.subtract(eta[1:], eta[:-1], out=du)
    np.multiply(bux * g, du, out=du)
    np.divide(du, dx, out=du)
    ui = u[1:-1]
    np.multiply(aux, ui, out=ui)
    np.subtract(ui, du, out=ui)
    np.subtract(eta[:, 1:], eta[:, :-1], out=dv)
    np.multiply(bvy * g, dv, out=dv)
    np.divide(dv, dx, out=dv)
    vi = v[:, 1:-1]
    np.multiply(avy, vi, out=vi)
    np.subtract(vi, dv, out=vi)
    np.subtract(u[1:], u[:-1], out=t)
    np.multiply(bx * h, t, out=t)
    np.divide(t, dx, out=t)
    np.multiply(ax, ex, out=ex)
    np.subtract(ex, t, out=ex)
    np.subtract(v[:, 1:], v[:, :-1], out=t)
    np.multiply(by * h, t, out=t)
    np.divide(t, dx, out=t)
    np.multiply(ay, ey, out=ey)
    np.subtract(ey, t, out=ey)


if HAVE_NUMBA:
    @_nb.njit(cache=True)
    def _sw_step_numba(ex, ey, u, v, aux, bux, avy, bvy, ax, bx, ay, by, g, h, dx):
        nx, ny = ex.shape
        for i in range(1, nx):
            a, b = aux[i - 1], bux[i - 1]
            for j in range(ny):
                u[i, j] = a * u[i, j] - b * g * ((ex[i, j] + ey[i, j]) - (ex[i - 1, j] + ey[i - 1, j])) / dx
        for i in range(nx):
            for j in range(1, ny):
                v[i, j] = (avy[j - 1] * v[i, j]
                           - bvy[j - 1] * g * ((ex[i, j] + ey[i, j]) - (ex[i, j - 1] + ey[i, j - 1])) / dx)
        for i in range(nx):
            a, b = ax[i], bx[i]
            for j in range(ny):
                ex[i, j] = a * ex[i, j] - b * h * (u[i + 1, j] - u[i, j]) / dx
                ey[i, j] = ay[j] * ey[i, j] - by[j] * h * (v[i, j + 1] - v[i, j]) / dx


class ShallowWaterKernel:
    """Advances (ex, ey, u, v) of hydar.Solver by one step, in place.

    aux, bux : (nx-1,) PML factors of the u faces; avy, bvy : (ny-1,) of the v faces
    ax, bx   : (nx,) of ex; ay, by : (ny,) of ey
    """

    def __init__(self, nx, ny, aux, bux, avy, bvy, ax, bx, ay, by, use_numba=None):
        vec = lambda a: np.ascontiguousarray(np.ravel(a), dtype=float)
        self.coef = tuple(vec(a) for a in (aux, bux, avy, bvy, ax, bx, ay, by))
        self.numba = HAVE_NUMBA if use_numba is None else (use_numba and HAVE_NUMBA)
        if not self.numba:
            aux, bux, avy, bvy, ax, bx, ay, by = self.coef
            self.coef = (aux[:, None], bux[:, None], avy[None, :], bvy[None, :],
                         ax[:, None], bx[:, None], ay[None, :], by[None, :])
            self.work = (np.empty((nx, ny)), np.empty((nx - 1, ny)),
                         np.empty((nx, ny - 1)), np.empty((nx, ny)))

    def __call__(self, ex, ey, u, v, g, h, dx):
        if self.numba:
            _sw_step_numba(ex, ey, u, v, *self.coef, g, h, dx)
        else:
            _sw_step_numpy(ex, ey, u, v, *self.coef, g, h, dx, self.work)


# ------------------------------------------------------------- spectral solver
def _oscillate_numpy(e_hat, f, Pm_hat, cos, sin, K_over_w, w_over_K):
    """Exact step of every Fourier mode; returns (eta_k - Pm_k, phi_k) after the step."""
    a = e_hat + Pm_hat
    a_new = a * cos + f * K_over_w * sin
    f_new = f * cos - a * w_over_K * sin
    return a_new - Pm_hat, f_new


if HAVE_NUMBA:
    @_nb.njit(cache=True)
    def _oscillate_numba(e_hat, f, Pm_hat, cos, sin, K_over_w, w_over_K):
        nx, ny = e_hat.shape
        for i in range(nx):
            for j in range(ny):
                a = e_hat[i, j] + Pm_hat[i, j]
                fij = f[i, j]
                c, s = complex(cos[i, j]), complex(sin[i, j])
                a_new = a * c + fij * complex(K_over_w[i, j]) * s
                f[i, j] = fij * c - a * complex(w_over_K[i, j]) * s
                e_hat[i, j] = a_new - Pm_hat[i, j]
        return e_hat, f


def oscillate(e_hat, f, Pm_hat, cos, sin, K_over_w, w_over_K):
    """Advance the Fourier modes of SpectralSolver by one step (may reuse e_hat and f)."""
    if HAVE_NUMBA:
        return _oscillate_numba(e_hat, f, Pm_hat, cos, sin, K_over_w, w_over_K)
    return _oscillate_numpy(e_hat, f, Pm_hat, cos, sin, K_over_w, w_over_K)
