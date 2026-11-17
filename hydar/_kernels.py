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
    @_nb.njit(cache=True, parallel=True)
    def _sw_step_numba(ex, ey, u, v, aux, bux, avy, bvy, ax, bx, ay, by, g, h, dx):
        # rows are independent within each sweep, so the threads (prange) change no number
        nx, ny = ex.shape
        for i in _nb.prange(1, nx):
            a, b = aux[i - 1], bux[i - 1]
            for j in range(ny):
                u[i, j] = a * u[i, j] - b * g * ((ex[i, j] + ey[i, j]) - (ex[i - 1, j] + ey[i - 1, j])) / dx
        for i in _nb.prange(nx):
            for j in range(1, ny):
                v[i, j] = (avy[j - 1] * v[i, j]
                           - bvy[j - 1] * g * ((ex[i, j] + ey[i, j]) - (ex[i, j - 1] + ey[i, j - 1])) / dx)
        for i in _nb.prange(nx):
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


# ------------------------------------------------------------- vortex street
_EULER = 0.5772156649015329


def _street_numpy(x, y, xi, yi, gi, ri, n_real, rho, dU):
    """Pressure of a Kármán street with bank images at points (x, y) (see KarmanStreet)."""
    from .karman import lamb_oseen_pressure
    dx = x[:, None] - xi[None, :]
    dy = y[:, None] - yi[None, :]
    r2 = np.maximum(dx ** 2 + dy ** 2, 1e-12)
    fac = gi / (2 * np.pi * r2) * (-np.expm1(-r2 / ri ** 2))    # u_theta / r
    ux, uy = -fac * dy, fac * dx
    Ux, Uy = ux.sum(1), uy.sum(1)
    self_sq = (ux[:, :n_real] ** 2 + uy[:, :n_real] ** 2).sum(1)
    p_lo = lamb_oseen_pressure(np.sqrt(r2[:, :n_real]), gi[:n_real], ri[:n_real], rho).sum(1)
    return -rho * dU * Ux - 0.5 * rho * (Ux ** 2 + Uy ** 2 - self_sq) + p_lo


if HAVE_NUMBA:
    @_nb.njit(cache=True)
    def _e1(x):
        """Exponential integral E1(x), x > 0, to about 1e-15 (series below 1, continued fraction above)."""
        if x <= 1.0:
            term, total, k = 1.0, 0.0, 1
            while True:
                term *= -x / k
                d = -term / k
                total += d
                if abs(d) < 1e-17 * abs(total) or k > 60:
                    break
                k += 1
            return -np.log(x) - _EULER + total
        # modified Lentz for E1(x) = e^-x / (x + 1 / (1 + 1 / (x + 2 / (1 + 2 / (x + ...)))))
        b = x + 1.0
        c = 1.0 / 1e-300
        d = 1.0 / b
        h = d
        for i in range(1, 200):
            an = -float(i) * i
            b += 2.0
            d = 1.0 / (an * d + b)
            c = b + an / c
            de = c * d
            h *= de
            if abs(de - 1.0) < 1e-16:
                break
        return h * np.exp(-x)

    @_nb.njit(cache=True)
    def _street_numba(x, y, xi, yi, gi, ri, n_real, rho, dU, s_far):
        n, m = x.size, xi.size
        out = np.empty(n)
        c0 = 1.0 / (8 * np.pi ** 2)
        for p in range(n):
            Ux = Uy = self_sq = p_lo = 0.0
            for q in range(m):
                dx = x[p] - xi[q]
                dy = y[p] - yi[q]
                r2 = max(dx * dx + dy * dy, 1e-12)
                fac = gi[q] / (2 * np.pi * r2) * (-np.expm1(-r2 / ri[q] ** 2))
                ux, uy = -fac * dy, fac * dx
                Ux += ux
                Uy += uy
                if q < n_real:
                    self_sq += ux * ux + uy * uy
                    s = max(r2 / ri[q] ** 2, 1e-12)
                    F = np.expm1(-s) ** 2 / s
                    if s <= s_far:
                        F += 2 * (_e1(s) - _e1(2 * s))
                    p_lo += -rho * gi[q] ** 2 * c0 / ri[q] ** 2 * F
            out[p] = -rho * dU * Ux - 0.5 * rho * (Ux * Ux + Uy * Uy - self_sq) + p_lo
        return out


def street_pressure(x, y, xi, yi, gi, ri, n_real, rho, dU, s_far=50.0):
    """Street pressure [Pa] at points (x, y); vortices (xi, yi, gi, ri), the first n_real real."""
    if HAVE_NUMBA:
        return _street_numba(x, y, xi, yi, gi, ri, n_real, rho, dU, s_far)
    return _street_numpy(x, y, xi, yi, gi, ri, n_real, rho, dU)
