# Changelog

## Unreleased

### Added

- **Case study 5: a fish finds a rock from its Kármán vortex street** (`hydar.karman`). A hemispherical rock, half of its height above the water, in a river with a trapezoidal cross section sheds a vortex street; a 25 cm fish with 10 lateral-line pressure sensors (5 on the head third, 5 on the middle third, none on the tail third) surveys the river across the flow and estimates the rock's lateral position, size and distance. `TrapezoidalChannel` (Manning flow), `HemisphereRock`, `KarmanStreet` (kinematic street: Strouhal shedding, von Kármán spacing and circulation, Lamb-Oseen cores with eddy-viscosity growth, bed-friction decay, banks as images, unsteady Bernoulli pressure), `Fish`, `RiverTurbulence`, `lateral_survey` and `locate_rock`.
- `examples/08_fish_rock.py` (about 10 s) and `examples/09_benchmark_fish_rock.py` (625 runs, about 30 minutes on 4 cores; results in `docs/data/fish_rock_benchmark.csv`).
- README case study 5 with Figure 8 (`docs/make_figures.py 7`), 10 tests in `tests/test_karman.py`, and five references (Chen & Jirka 1995, Chow 1959, Liao et al. 2003, von Kármán 1911, Williamson 1996).
- Case study 5 animation `docs/img/hydar_fish.gif` / `.mp4` (`examples/make_fish_animation.py`), rendered from a real survey in the style of the other animations.
- Figure 9 (`fig10_fish_rock_results`): the benchmark in seven panels (lateral, distance and diameter errors and detection rate against distance and turbulence; errors under a wrong wake model). The benchmark panels of Figure 8 moved there.
- Results tables `docs/data/fish_rock_summary.md` (all 25 settings), written by the benchmark and rebuilt from the CSV with `examples/09_benchmark_fish_rock.py --summary`.

### Case study 5 results

| Quantity | Result |
|---|---|
| Lateral position | Median error 1 to 4 cm from 3 to 9 m behind the rock in turbulence up to 10 %, also when the wake model is wrong |
| Rock size D | Within 1 to 3 % at 6 to 9 m; 23 to 29 % too large at 3 m, where the street is still forming; carries any error of the assumed Strouhal number |
| Distance | Median error 0.13 to 0.29 m at 6 to 9 m with the right wake model; biased by 1.3 to 2.8 m when the eddy viscosity or bed roughness of the true wake differ by a factor of 1.5 to 2; detection fails in half the runs at 12 m in 15 % turbulence |

## v1.2.0 (2026-11-14)

Applies the repository review of 2026-11-14 (`repository_review/`). The Kelvin method of case study 3 is retuned, so that the benchmark of case study 4 compares against a fair baseline, and the benchmark is rerun with 100 seeds. The results of case studies 3 and 4 change; the conclusions for the ferry are weaker than in v1.1.0 and are stated below.

### Case study 4 results

| Scenario | Result |
|---|---|
| Ferry, one array 2.5 km away, wind sea 0 to 0.6 m | Interferometry more accurate in every wind: course bias below 0.1° (Kelvin method +0.5 to +0.7°), median track error 5.9 times smaller in calm water, 3.0 times at Hs 0.3 m (131 m vs 397 m) and 1.9 times at 0.6 m, fixes in 1500/1500 runs (Kelvin 1448/1500). The optimal 20 m array lowers the median error a further 13 to 28 % (better in 242 of 400 paired runs with wind, p < 0.001). The optimal array with the same footprint helps only in calm water and at Hs 0.6 m |
| Boat, four arrays 8 to 18 m away, wind sea 0 to 0.2 m | Mixed. Calm and light wind: interferometry halves track error and speed bias. Hs 0.05 m: similar median track error, but the cusp-line method is more precise in course (0.21° vs 1.07° spread). Hs ≥ 0.1 m: both degrade; interferometry returns fewer but mostly better tracks at 0.1 m |

100 seeds per wind height, 2 Pa sensor noise, Kelvin baseline retuned (see Changed). The benchmark of v1.1.0 (20 seeds, before the repository review) reported a Kelvin course bias of 2° and a factor of 3 to 9; most of that bias came from the baseline's analysis window.

### Changed

- **Kelvin method of case study 3 retuned** after the repository review of 2026-11-14 (`repository_review/`), so that the benchmark compares against a fair baseline:
  - `locate_vessel` measures the cusp-wave direction in the 30 s before the cusp (`cusp_window`), not from 10 s before to 30 s after it. Behind the cusp line the two wake branches overlap and bias the course.
  - `ridge_points` interpolates each spectrogram peak between frequency bins.
  - `wave_direction` uses the unbiased cross-correlation. The raw correlation pulled the delays of a narrow-band signal towards zero and overstated the phase speed by 5 to 46 %. It now fits in two passes (short baselines first, then all pairs within half a period of the predicted delay), so that long baselines of a wide array do not lock one period off.
  - Case study 3 now reads: distance 2403 m (was 2390), course 0.7° off (was 1.9°), position 365 m off (was 463 m).
- Sensor noise is a property of the array (`SensorArray(noise_pa=..., noise_seed=...)`) and is added after the pressure model. Before, `record_passage` added it to the surface elevation, so the depth filter reduced the stated 2 Pa to 0.14 Pa.
- `wind_sea` takes `f_range`, `depth` and `z`; `wind_sea_head` now calls it. `stacked_spectrogram`, `fit_chirp` and `wave_direction` expose constants that were fixed for the ferry scenario (`sensors`, `d_guesses`, `min_speed`). `track_vessel` passes `peak_rel` on to `observe`.
- The benchmark (`examples/07_benchmark_interferometry.py`) runs with any multiprocessing start method (it failed on macOS, Windows and Python 3.14), uses 100 seeds instead of 20, records why a run gave no fix (`reason` column), prints paired seed-by-seed comparisons, and writes a partial run to a separate file.

### Fixed

- `SpectralSolver` raises a `TypeError` for a volume source such as `DropObject` instead of ignoring it.
- Clear errors for inputs that used to fail obscurely: `locate_vessel` with an analysis window outside the record or no spectrogram ridge (`RuntimeError`), a drop outside the grid, parallel bearings in `intersect_bearings` (`ValueError`). `pick_wake` accepts a single record, and `pick_arrivals` returns no picks for an empty time window.

### Tests and CI

- 14 new tests (41 in total), among them the boat path of case study 4, `pick_wake`, `wave_direction`, the wind sea height, the sensor noise level and the benchmark under `spawn` and `forkserver`.
- GitHub Actions runs the fast tests on every pull request, on Python 3.9 with the oldest supported numpy, scipy and matplotlib and on Python 3.13 with the newest; all tests run weekly.

## v1.1.0 (2026-11-09)

Prepared but not released on its own; v1.2.0 is the first release that contains it.

Adds case study 4: a method that designs the 3 × 3 array geometry for correlative interferometry and detects and tracks a vessel **without Kelvin wake theory**, with a Monte Carlo benchmark against the Kelvin methods of case studies 2 and 3.

### Added

- **Case study 4: array design and vessel tracking by correlative interferometry** (`hydar.interferometry`). Tracks a vessel without Kelvin wake theory, from two general facts of linear waves: a wave of the steady pattern does not move relative to the source (n · V = c), and its energy travels from the source along a ray at the group velocity. Each 3 × 3 array measures wave directions by correlative interferometry (all 36 sensor pairs against a table of plane-wave phases); a robust fit (RANSAC, then a Cauchy-loss least squares; optional rho weighting) gives the track and rejects the wind sea.
- `design_array`: optimal 3 × 3 grid (free column and row positions) for a frequency band, minimizing the bearing Cramér-Rao bound subject to a sidelobe (ambiguity) limit, a footprint limit and a minimum spacing.
- `wind_sea_head`: finite-depth wind sea as pressure head at the sensor depth, to add wind noise to `SpectralSolver` records.
- `SensorArray(layout=...)` for arrays with any sensor positions.
- `examples/06_array_design.py`: designs the optimal grids for the ferry band (same 7.1 m footprint as the Rätsep et al. grid, and a 20 m footprint), records the case study 3 ferry passage with each in a 0.3 m wind sea, and compares the Kelvin method with interferometry (about 2 minutes).
- `examples/07_benchmark_interferometry.py`: 540 Monte Carlo runs (20 seeds; 2700 records and 100 seeds from v1.2.0) over wind-sea heights for both vessel scenarios, geometries and methods; results in `docs/data/ci_benchmark.csv`.
- README case study 4 with Figure 7, five tests in `tests/test_interferometry.py`, and three references (Fischler & Bolles 1981, Lighthill 1978, Van Trees 2002).

### Changed

- `kelvin.wave_direction` sizes its correlation lag window from the largest sensor separation instead of `spacing`, so it also works for non-uniform layouts.

### Case study 4 results

Superseded by v1.2.0. The 20-seed benchmark of this version reported a Kelvin course bias of about 2° on the ferry and an interferometry track error 3 to 9 times smaller.

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
