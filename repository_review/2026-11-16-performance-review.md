# HYDAR: performance review of the solvers and calculations, 2026-11-16

Builds on the comprehensive review of 2026-11-14 (`2026-11-14-comprehensive-repository-review.md`),
which covered correctness, reproducibility and CI; this one covers only speed.

| Item | Value |
|---|---|
| Reviewed commit | `be4e828` (main, v1.3.0) |
| Environment | Python 3.13, NumPy 2.5.3, SciPy 1.18.1, Numba 0.68.1, Rust 1.x with PyO3 0.23, Linux, 4 cores |
| Question (Jeff) | Can the calculations and the numerical solver be made more efficient, perhaps with Rust and/or Fortran? |
| Rule kept | Every change must leave the results unchanged: tests pass and the benchmarks reproduce |

## 1. Answer in short

Rust or Fortran would not pay for themselves here. The code is already vectorised NumPy, and the time
goes to four places, of which only one is a loop that a compiled language speeds up:

1. **The shallow-water `Solver`** (case study 1, drops) is a memory-bound stencil. A compiled loop is
   11 to 12 times faster than NumPy. Numba and Rust (PyO3) were both measured, and both took 0.7 ms
   per step against 8.1 ms for NumPy, with bitwise identical results. Fortran was not built: the
   container has no Fortran compiler. The same loop compiled by gfortran would almost certainly run at
   the same speed, since the step is limited by memory bandwidth, not by the language.
2. **The `SpectralSolver`** (case studies 2 and 4, boats) spends about 75 % of each step in FFTs, which
   already run in SciPy's multi-threaded C++ (pocketfft/ducc). No language beats that. Removing repeated
   work gives 1.35 times.
3. **The vortex street** (case study 5) spent 70 % of its time in the special function E1 (`scipy.special.exp1`).
   It is now evaluated only where it can change the result, which makes it 1.7 times faster.
4. **Array design** (case study 4) spent its time in one cosine per (pair, direction, source). One matrix
   product now replaces them, which makes it 4 times faster and returns the same designs.

Numba gives the compiled speed of item 1 without a compiler on the user's machine and without
platform wheels. It is an optional extra (`pip install hydar[fast]`); without it the same kernels run in
NumPy and give the same numbers. Rust or Fortran would mean a build toolchain, one compiled wheel per
platform and Python version (cibuildwheel), and on Windows a Fortran compiler next to MSVC. That buys
nothing over Numba for this code. Rust becomes worth it if HYDAR's processing (not the simulator) has to
run on a sensor node or another embedded device, where Python is not available.

## 2. Where the time goes (profiled on main)

| Workload | Wall time on main | Hot spot | Share |
|---|---|---|---|
| Fast tests (46) | 118 s | `peak_sidelobe` (array design), shallow-water `step`, FFTs, `lamb_oseen_pressure` | 35 s, 30 s, 15 s, 12 s |
| Boat solver run of benchmark 07 (4400 steps, 1120 x 560) | about 4.3 min | FFTs 46 %, element-wise complex arithmetic 22 %, `depth_attenuation` recomputed every step 14 % | |
| `design_array` (3 calls in benchmark 07) | 52 s each | `peak_sidelobe` | 80 % |
| Ferry wake synthesis (3 geometries in benchmark 07) | 28 s each | complex exponential of (times x 12000 angles) in `KelvinShip.wake` | 100 % |
| One run of benchmark 09 (625 runs) | 10 s | `lamb_oseen_pressure`, mostly `exp1` | 70 % |
| One ferry Monte Carlo run of benchmark 07 | 3 s | spread over FFTs, correlations and `dynamic_pressure`; no single hot spot | |

## 3. What was changed (all on branch `claude/project-thread-oebx14`)

| Change | Where | Result identical? | Speed-up (measured) |
|---|---|---|---|
| Shallow-water step as one compiled loop (Numba), or in-place NumPy without Numba | `hydar/_kernels.py`, `Solver.step` | bitwise | 8.1 → 0.7 ms per step (11x); 4.0 ms with NumPy only (2x) |
| `Solver.run` reads the gauges from `ex + ey` at their cells instead of summing the whole grid | `Solver.run` | bitwise | part of the above |
| `DropObject` sums `ex + ey` on its own box only | `drop.py` | bitwise | drop test 21 → 6 s |
| `depth_attenuation(k, z)` cached per sensor depth (it never changes) | `SpectralSolver.pressure_head` | bitwise | 14 % of a step with sensors |
| Mode update (harmonic oscillator) fused in one Numba loop; damping in place | `_kernels.oscillate`, `SpectralSolver.step` | bitwise | with the cache: 300 boat steps with 4 arrays 14.5 → 10.7 s (1.35x) |
| No FFTs of an all-zero pressure field when nothing presses on the water | `SpectralSolver._pressure` | bitwise | sponge test 4.2 → 2.3 s |
| E1 terms of the Lamb-Oseen pressure only where r^2/rc^2 <= 50 (beyond, they are below half an ulp of the result, checked from 40 to 1000) | `karman.lamb_oseen_pressure` | bitwise | one run of benchmark 09: 10.0 → 5.8 s (1.7x) |
| All source directions of the sidelobe search from one matrix product, `cos(a - b) = cos a cos b + sin a sin b` | `interferometry.peak_sidelobe` | to 1e-15 in rho; `design_array` returns the identical layout | 1.9 → 0.4 ms per call; `design_array` 52 → 13 s (4x) |

Identity was checked by running each hot path on main and on the branch and comparing the arrays:
fish benchmark run, street pressure, `design_array`, 300 boat steps with sensors, the shallow-water
solver with and without a dropped sphere, and the spectral solver without forcing are bitwise
identical, with and without Numba. A partial run of benchmark 07 (ferry, 2 seeds, all geometries and
methods, 60 rows) writes a CSV that is byte-identical to the one main writes.
`tests/test_kernels.py` keeps the Numba kernels and the NumPy fallback bitwise equal to the
original expressions. All 49 fast and 5 slow tests pass, with and without Numba.

**Tried and dropped:** `KelvinShip.wake` can use `cos(ph) @ A` instead of `Re(exp(i ph) @ A)` (A is real),
which is 2.3 times faster. Because a real matrix product rounds differently from a complex one, the wake then differs by 2e-16 m,
and that moved the ferry benchmark results in the 9th digit (track error 401.464704
against 401.464707 m). The CSV would no longer reproduce bit for bit, for a saving of about 50 s in a
30 minute benchmark, so the change was not kept.

The fast tests take 45 to 55 s instead of 118 s. `examples/10_speed.py` prints the per-call times above,
so the numbers can be repeated on Jeff's PC (with `HYDAR_NO_NUMBA=1` for NumPy only).

## 4. Rust and Fortran, measured and judged

The shallow-water step was written three ways and timed on the same 460 x 460 grid (one core):

| Implementation | ms per step | Identical to NumPy |
|---|---|---|
| NumPy as on main | 8.1 | (reference) |
| NumPy in place (no temporaries) | 3.8 to 4.0 | bitwise |
| Numba `@njit` loop | 0.69 | bitwise |
| Rust, PyO3 + `numpy` crate, `opt-level=3` | 0.70 | bitwise |
| Fortran (f2py or meson) | not built: no compiler in the container | expected to be the same as Rust and Numba |

Numba and Rust both go through LLVM and are limited by memory bandwidth, so they tie. What separates
them is the cost of shipping:

| | Numba (chosen) | Rust (PyO3, maturin) | Fortran (f2py, meson) |
|---|---|---|---|
| User installs | `pip install numba` (wheels for Windows, Linux, macOS) | a compiled wheel of HYDAR per platform and Python version | the same, plus a Fortran runtime |
| Developer needs | nothing | Rust toolchain | gfortran or flang; on Windows also MSVC, the hardest of the three |
| CI | one more matrix entry | cibuildwheel for every platform and Python version | the same, with Fortran on Windows |
| Without the accelerator | falls back to NumPy, same numbers | the package does not build from source | the same |
| Zenodo source archive | still pure Python | needs a compiler to run | needs a compiler to run |

## 5. What was not done, and what would come next

These keep the results only to rounding (about 1e-15) or change them, so they need Jeff's decision:

- **Fewer FFTs in the spectral solver.** By linearity, rfft2(eta + P) = rfft2(eta) + rfft2(P), so the
  spectra of the step can be reused for the sensors: two of the eleven FFTs per step with four sensor
  depths can go, about 15 %. The results then change in the last bits.
- **Single precision FFTs** (`float32`) halve the FFT time of the spectral solver. The results change
  by about 1e-7 relative, which is far below the sensor noise but would change the benchmark CSV.
- **Threads in the shallow-water loop** (`numba.prange` over rows) keep the results bitwise and should
  give another 2 to 3 times on 4 cores, but compete with the process pools of the benchmarks. Only
  case study 1 would gain.
- **Vortex street over all times at once** instead of a Python loop over time steps; about 20 % of a fish
  run is that loop's overhead.
- **The Monte Carlo runs of benchmark 07** (ferry and boat, 2700 runs) were not touched. They have no
  single hot spot: FFTs, correlations, `dynamic_pressure` and a 50-iteration Newton solve of the
  dispersion relation in `wavenumber` (within 1e-15 of its 50-iteration answer after 4 iterations for 0.01 to 5 Hz in 3 m water; stopping on convergence would change the last bits) share the time. This is the largest
  remaining part of the 30 minutes of that benchmark, so it is the next place to look.

## 6. Effect on the benchmarks

| Benchmark | Before | After (estimate from the parts above) |
|---|---|---|
| 07 interferometry (100 seeds) | about 30 min on 4 cores | about 27 min: the boat solver run (4.3 → 3.2 min) and array design (2.6 → 0.7 min) are faster; the wake synthesis (1.4 min) and the Monte Carlo runs are not |
| 09 fish and rock (25 seeds) | about 25 min on 4 cores | about 15 min |

## 7. Second round: speed-ups that change the numbers (Jeff chose "More speed, numbers shift", 2026-11-16)

| Change | Speed (measured) | Effect on results |
|---|---|---|
| `numba.prange` over rows in the shallow-water step | 0.69 → 0.17 ms per step (48x v1.3.0) | none: bitwise identical (rows are independent) |
| `SpectralSolver`: sensors reuse rfft2(eta) and the pressure spectrum (linearity); P on the grid only on request | 2 of 11 FFTs per step saved; boat without sensors 25 → 16 ms per step | to 1e-18 m (rounding of the sum) |
| `SpectralSolver(precision="single")`, used by benchmark 07 | boat with 4 arrays 41 → 16.5 ms per step (2.5x) | sensor records 5e-5 relative, at most 2.6e-6 m = 0.025 Pa (sensor noise 2 Pa) over the 44 s run |
| Vortex street as one Numba loop with its own E1 (series below 1, continued fraction above) | a run of benchmark 09 10.0 → 3.4 s (2.9x v1.3.0) | 6e-14 Pa on 26 Pa; the benchmark's rounded CSV is byte-identical |
| `KelvinShip.wake`: `cos(ph) @ A` | 2.3x | 2e-16 m |
| `wavenumber` stops Newton at round-off (4 to 6 instead of 50 iterations) | small | to round-off |

Both benchmarks were rerun in full on 4 cores:

- **Benchmark 07:** 21.5 min (it was quoted as about 30). The CSV changes in the last digits. Medians move by about 1e-5 m. A few boat runs in wind (H_s 0.1 m and more), where `track_boat` or the interferometry fit sit between two solutions, jump to the other one: by up to 143° in course and 608 m in track error. No run gains or loses a fix. Every ferry cell of the summary table is unchanged. Four cells of the README's boat table moved: calm 1 m interferometry median 0.73 → 0.72 m; H_s 0.05 m optimal interferometry course bias +0.86 → +0.85°; H_s 0.1 m Kelvin course bias +26.1 ± 61.1 → +24.7 ± 61.7° and speed SD 2.06 → 2.05 m/s. In the seed-by-seed comparison of boat geometries, the H_s 0.05 m median difference moved from −0.13 to −0.12 m (p 0.131 → 0.132). No conclusion changes. Figure 7 was not regenerated; its changes are far below a line width.
- **Benchmark 09:** 7.5 min (it was about 25). The CSV and the results tables are byte-identical.

The library default stays `precision="double"`. Only benchmark 07 switches to single precision, where the change is 80 times below the sensor noise.
