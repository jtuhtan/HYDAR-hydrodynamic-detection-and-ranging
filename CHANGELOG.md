# Changelog

## v1.1.0 (2026-11-09)

Adds case study 4: a method that designs the 3 × 3 array geometry for correlative interferometry and detects and tracks a vessel **without Kelvin wake theory**, with a Monte Carlo benchmark against the Kelvin methods of case studies 2 and 3.

### Added

- **Case study 4: array design and vessel tracking by correlative interferometry** (`hydar.interferometry`). Tracks a vessel without Kelvin wake theory, from two general facts of linear waves: a wave of the steady pattern does not move relative to the source (n · V = c), and its energy travels from the source along a ray at the group velocity. Each 3 × 3 array measures wave directions by correlative interferometry (all 36 sensor pairs against a table of plane-wave phases); a robust fit (RANSAC, then a Cauchy-loss least squares; optional rho weighting) gives the track and rejects the wind sea.
- `design_array`: optimal 3 × 3 grid (free column and row positions) for a frequency band, minimizing the bearing Cramér-Rao bound subject to a sidelobe (ambiguity) limit, a footprint limit and a minimum spacing.
- `wind_sea_head`: finite-depth wind sea as pressure head at the sensor depth, to add wind noise to `SpectralSolver` records.
- `SensorArray(layout=...)` for arrays with any sensor positions.
- `examples/06_array_design.py`: designs the optimal grids for the ferry band (same 7.1 m footprint as the Rätsep et al. grid, and a 20 m footprint), records the case study 3 ferry passage with each in a 0.3 m wind sea, and compares the Kelvin method with interferometry (about 2 minutes).
- `examples/07_benchmark_interferometry.py`: 540 Monte Carlo runs over wind-sea heights for both vessel scenarios, geometries and methods; results in `docs/data/ci_benchmark.csv`.
- README case study 4 with Figure 7, five tests in `tests/test_interferometry.py`, and three references (Fischler & Bolles 1981, Lighthill 1978, Van Trees 2002).

### Case study 4 results

| Scenario | Result |
|---|---|
| Ferry, one array 2.5 km away, wind sea 0 to 0.6 m | Interferometry better on every measure: course bias below 0.1° (Kelvin method about 2°), track error 3 to 9 times smaller (Hs 0.3 m: median 103 m with the optimal 20 m array vs 537 m), fixes in 300/300 runs (Kelvin 293/300). The optimal 20 m array lowers the error a further 9 to 22 % (median) |
| Boat, four arrays 8 to 18 m away, wind sea 0 to 0.2 m | Mixed. Calm and light wind: interferometry roughly halves track error and speed bias. Hs 0.05 m: equal median track error, but the cusp-line method is more precise in course (0.23° vs 0.8 to 1.3° spread). Hs ≥ 0.1 m: both fail |

### Changed

- `kelvin.wave_direction` sizes its correlation lag window from the largest sensor separation instead of `spacing`, so it also works for non-uniform layouts.
- **Kelvin method of case study 3 retuned** after the repository review of 2026-11-14 (`repository_review/`), so that the benchmark compares against a fair baseline:
  - `locate_vessel` measures the cusp-wave direction in the 30 s before the cusp (`cusp_window`), not from 10 s before to 30 s after it. Behind the cusp line the two wake branches overlap and bias the course.
  - `ridge_points` interpolates each spectrogram peak between frequency bins.
  - `wave_direction` uses the unbiased cross-correlation. The raw correlation pulled the delays of a narrow-band signal towards zero and overstated the phase speed by 5 to 45 %.
  - Case study 3 now reads: distance 2403 m (was 2390), course 0.7° off (was 1.9°), position 365 m off (was 463 m).
- Sensor noise is a property of the array (`SensorArray(noise_pa=..., noise_seed=...)`) and is added after the pressure model. Before, `record_passage` added it to the surface elevation, so the depth filter reduced the stated 2 Pa to 0.14 Pa.
- `wind_sea` takes `f_range`, `depth` and `z`; `wind_sea_head` now calls it. `stacked_spectrogram`, `fit_chirp` and `wave_direction` expose constants that were fixed for the ferry scenario (`sensors`, `d_guesses`, `min_speed`). `track_vessel` passes `peak_rel` on to `observe`.
- The benchmark (`examples/07_benchmark_interferometry.py`) runs with any multiprocessing start method (it failed on macOS, Windows and Python 3.14), uses 100 seeds instead of 20, records why a run gave no fix (`reason` column), prints paired seed-by-seed comparisons, and writes a partial run to a separate file.

### Fixed

- `SpectralSolver` raises a `TypeError` for a volume source such as `DropObject` instead of ignoring it.
- Clear errors for inputs that used to fail obscurely: `locate_vessel` with an analysis window outside the record or no spectrogram ridge (`RuntimeError`), a drop outside the grid, parallel bearings in `intersect_bearings` (`ValueError`). `pick_wake` accepts a single record, and `pick_arrivals` returns no picks for an empty time window.

### Tests and CI

- 13 new tests (40 in total), among them the boat path of case study 4, `pick_wake`, `wave_direction`, the wind sea height, the sensor noise level and the benchmark under `spawn` and `forkserver`.
- GitHub Actions runs the fast tests on every pull request, on Python 3.9 with the oldest supported numpy, scipy and matplotlib and on Python 3.13 with the newest; all tests run weekly.

## v1.0.0 (2026-11-02)

First release of HYDAR (Hydrodynamic Detection and Ranging), a Python simulator for locating objects in water from the surface waves they make, using arrays of pressure sensors.

### Highlights

- **Two wave solvers.**
  - `Solver`: linear shallow-water equations on a staggered grid, with split-field PML walls that do not reflect.
  - `SpectralSolver`: linear potential flow with the exact dispersion relation, an exponential time integrator and absorbing sponge walls. It produces true Kelvin wakes.
- **Drop objects.** `DropObject` drops a sphere of any diameter and density. It sinks or floats, and is two-way coupled to the water surface.
- **Boats.** `Boat` with a Series 60 hull (`Series60Hull`, fitted to C_B 0.60, C_P 0.614, C_M 0.977, C_WP 0.706 and LCB 1.5 % L aft) or an elliptical hull. It acts as a volume source in the shallow-water solver and as a pressure patch in the dispersive solver.
- **Pressure sensor arrays.** `SensorArray` places n × n grids at any x, y and depth z, with any rotation, sampling at 100 Hz. Pressure models: `hydrostatic`, `airy` (depth filter) and `solver` (pressure read directly from the dispersive solver at the sensor depth).
- **HYDAR localization.**
  - Drops: arrival picking, bearing from time differences of arrival, range from time of flight, and fusion across arrays, with the drop time known or unknown.
  - Boats: Kelvin cusp-line picking (`pick_wake`) and fitting of a constant-velocity track (`track_boat`). Mach-cone tracking for supercritical boats in the shallow-water solver.
  - Ships far offshore (`hydar.kelvin`): synthetic Havelock Kelvin wake, spectrogram chirp fit, wave direction from cross-correlation, and vessel fix, after Rätsep et al. (2020, 2021).

### Case studies and results

| Case | Result |
|---|---|
| 1. 0.5 m sphere dropped among four 3 × 3 arrays | Bearings within 0.4°, ranges within 0.2 m, fused fix within 0.05 m |
| 2. 4.5 m Series 60 boat at 3 m/s, Kelvin wake from the dispersive solver | Speed 3.06 vs 3.00 m/s, heading 3.7° vs 4.0°, track within 0.7 to 2.0 m of the hull centre |
| 3. Ferry at 11 m/s, 2.5 km from one seabed array | Speed 10.99 vs 11.00 m/s, distance 2390 vs 2500 m, course 1.9° off, position about 460 m off |

### Also included

- Five runnable examples, plus scripts that render the README hero animation and the boat animation.
- README figures in publication style (PNG at 600 dpi, PDF, SVG), regenerated by `docs/make_figures.py`.
- A test suite of 22 tests, covering absorption at the walls, dispersion, Archimedes equilibrium, the Series 60 coefficients, and localization accuracy for all three case studies.

### Known limitations

- The drop case uses the non-dispersive shallow-water solver.
- All models are linear and assume uniform depth. Case study 3 ignores refraction and shoaling.
- No licence has been chosen yet.
