# HYDAR: comprehensive repository review, 2026-11-14

## 1. Review metadata

| Item | Value |
|---|---|
| Repository | `jtuhtan/hydar-hydrodynamic-detection-and-ranging` |
| Reviewed commit | `8d45565` (`8d455650c60661640ae570423823837de91fd94a`, 2026-11-09), `main`, clean working tree |
| Package version | 1.1.0 (prepared, not yet released) |
| Review date | 2026-11-14 |
| Repository type | Code-first research simulator (Python package, tests, examples, one benchmark data file) |
| Environment used | Python 3.13.18, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, Linux, 4 cores |
| Finding prefix | `HYDAR-###` (permanent; a re-review keeps the IDs and updates the status) |
| Source changes made | None. This review adds only this file and one JOURNAL.md entry. |

## 2. Executive assessment

**Overall risk: low to moderate.** The code is correct where it was tested, the published numbers are real, and the benchmark is reproducible to the last digit. No Critical finding. The weaknesses are in how far two conclusions of case study 4 can be pushed, and in a benchmark script that only runs on some platforms.

What held up:

- All 27 tests pass (23 fast in 104 s, 4 slow in 182 s).
- Every number in the two README benchmark tables was recomputed from `docs/data/ci_benchmark.csv` and matches.
- A fresh run of 60 benchmark rows (ferry, seeds 0 and 1, all geometries, both methods) is bit-identical to the committed CSV.
- The three array designs in the README reproduce exactly, and examples 03 and 05 print the numbers the README quotes.

The headline finding is **HYDAR-001**: the Kelvin baseline in the ferry benchmark carries a fixed 2.0° course bias that comes from its analysis window, not from Kelvin theory or from noise. Measuring the wave direction in the 30 s before the cusp instead of from 10 s before to 30 s after it cuts the bias to 0.3° in calm water. The README sentence "interferometry wins on every measure" is therefore a win against an untuned baseline. The second High finding, **HYDAR-002**, is that `examples/07_benchmark_interferometry.py` fails on macOS, Windows and Python 3.14 and later, because it relies on the `fork` start method.

| Severity | Count | IDs |
|---|---|---|
| Critical | 0 | |
| High | 2 | HYDAR-001, HYDAR-002 |
| Medium | 7 | HYDAR-003 to HYDAR-009 |
| Low | 9 | HYDAR-010 to HYDAR-018 |

## 3. Scope and method

The review covers the whole tracked tree at `8d45565`: 65 files (28 Python, 24 figure files, 4 animations, 4 Markdown, 1 CSV, packaging and citation metadata). The method follows the `repository-review` skill: inventory, the project's own gates, independent probes, then findings with stable IDs. The skill's reference files (template, severity rubric, probe script) were not present in this session, so the structure follows the skill's outline and the probes were written by hand.

Severity meanings used here, for a code-first research repository:

- **Critical**: a published result is wrong, or data can be lost.
- **High**: a published conclusion is not supported as worded, or a main workflow is broken for a typical user.
- **Medium**: a correctness or validity risk that passes the existing tests silently, or a missing safeguard.
- **Low**: drift, hygiene and maintainability.

All nine modules of `hydar/` were read in full, together with the tests, the examples for case studies 3 and 4 and the benchmark script. Every finding below cites a file, a command or a measured number.

## 4. Repository inventory

| Area | Files | Lines | Notes |
|---|---|---|---|
| `hydar/` | 10 | 2407 | two solvers, drop, boat, sensors, three localization modules, plotting |
| `tests/` | 7 | 382 | 27 tests, 4 marked `slow` |
| `examples/` | 9 | 1137 | 01 to 07 and two animation scripts |
| `docs/` | 2 scripts, 28 media, 1 CSV | 1428 | figures in PNG, PDF and SVG; `ci_benchmark.csv` (1080 rows) |
| Root | 8 | | README, CHANGELOG, JOURNAL, CLAUDE, CITATION.cff, LICENSE, pyproject, .gitignore |

There is no CI workflow, no lint or format configuration and no lockfile. `.git` is 49 MB, of which `docs/img` accounts for 32 MB in the working tree.

## 5. Verification record: the project's own gates

| Gate | Command | Result |
|---|---|---|
| Fast tests | `python -m pytest -m "not slow" -q` | 23 passed, 4 deselected, 104 s |
| Slow tests | `python -m pytest -m slow -q` | 4 passed, 23 deselected, 182 s |
| Install | `pip install -e . --no-deps` | succeeds, `hydar.__version__` is 1.1.0 |
| Example 03 | `python examples/03_hydar_localization.py` | fused error 0.04 m (README: 4 cm), 47 s |
| Example 05 | `python examples/05_vessel_localization.py` | 10.99 m/s, 2390 m, course −67.1°, 463 m (README: identical), 38 s |

The project defines no other gates. Examples 04, 06 and the animation and figure scripts were not executed (see section 12).

## 6. Independent probes

| Probe | Result |
|---|---|
| README benchmark tables against the CSV | All 21 table rows match (fixes, bias, SD, median, 90th percentile). Totals match: Kelvin 293/300, interferometry 300/300 for the ferry. |
| Benchmark reproduction | `N_SEEDS=2`, ferry only, written to a scratch file: 60 rows, 180 values, all bit-identical to the committed CSV. |
| Array designs | `design_array` for the three README cases returns the same gaps, sidelobes and bound ratios (0.907, 0.337, 0.954). |
| Paired statistics on the CSV | See HYDAR-003. |
| Benchmark under other start methods | `fork` runs; `spawn` and `forkserver` raise `KeyError: 'ferry'`. See HYDAR-002. |
| Cause of the 7 Kelvin "no fix" rows | Unhandled `ValueError` from `np.correlate` on an empty window. See HYDAR-004. |
| Kelvin course bias against the direction window | See HYDAR-001. |
| `design_array` against `min_spacing` | See HYDAR-005. |
| Sensor noise through the Airy filter | See HYDAR-011. |
| Edge inputs (drop in `SpectralSolver`, empty ridge, drop outside the grid, parallel bearings, 1-D record) | See HYDAR-007 and HYDAR-012. |
| Wind-sea realism (80 against 1000 components) | No effect found. See section 10. |
| `wind_sea` significant height | 4σ = 0.300 to 0.305 m for `hs=0.3`: correct. |
| Local links and images in README | All resolve. |
| Secret patterns, TODO and FIXME markers | None. |

## 7. Findings summary

| ID | Severity | Status | Finding |
|---|---|---|---|
| HYDAR-001 | High | Open | The ferry benchmark's Kelvin baseline has a 2.0° course bias set by its analysis window; a different window gives 0.3° |
| HYDAR-002 | High | Open | The benchmark script only works with the `fork` start method |
| HYDAR-003 | Medium | Open | The array-geometry claims rest on 20 seeds and are not significant seed by seed |
| HYDAR-004 | Medium | Open | The benchmark counts any exception as "no fix"; the 7 Kelvin misses are an unguarded crash |
| HYDAR-005 | Medium | Open | Same-footprint "optimal" designs sit on the `min_spacing` bound |
| HYDAR-006 | Medium | Open (tracked in JOURNAL) | The ferry case combines deep-water wave kinematics with a 3 m deep site |
| HYDAR-007 | Medium | Open | `SpectralSolver` silently ignores volume sources such as `DropObject` |
| HYDAR-008 | Medium | Open | No CI; the declared minimum versions are untested |
| HYDAR-009 | Medium | Open | The boat half of case study 4 and several analysis functions have no tests |
| HYDAR-010 | Low | Open (tracked in JOURNAL) | CITATION.cff pairs version v1.1.0 with the v1.0.0 DOI |
| HYDAR-011 | Low | Open | Sensor noise in `airy` records is 0.14 Pa, not the stated 2 Pa |
| HYDAR-012 | Low | Open | Edge inputs raise unhelpful errors or return the wrong shape |
| HYDAR-013 | Low | Open | Documentation drift (test count, Huber against Cauchy, package description) |
| HYDAR-014 | Low | Open | A partial benchmark run overwrites the full CSV; examples write PNGs into the working directory |
| HYDAR-015 | Low | Open | Scenario-specific constants are fixed inside library functions |
| HYDAR-016 | Low | Open | Duplicated constants and formulas; `peak_rel` is unreachable from `track_vessel` |
| HYDAR-017 | Low | Open | Repository weight and three figure-numbering schemes |
| HYDAR-018 | Low | Accepted risk | CC BY 4.0 is an unusual licence for software |

## 8. Findings in detail

### HYDAR-001 (High): the Kelvin baseline's course bias is a window artefact

**Where**: `hydar/kelvin.py:342`, `locate_vessel` calls `wave_direction(array, chirp.t_cusp - 10, chirp.t_cusp + 30)`.

**Evidence**: the README reports a Kelvin course bias of +2.03 ± 0.00° in calm water. A spread of zero over 20 runs means the bias is deterministic. On a noise-free, wind-free record of the benchmark scenario, the cusp-wave direction was measured in different windows and converted to a course:

| Window relative to the cusp | Course error, estimated cusp time | Course error, true cusp time |
|---|---|---|
| −10 to +30 s (as shipped) | +2.03° | +1.40° |
| −10 to +10 s | +1.51° | +0.94° |
| −5 to +5 s | +1.25° | +0.61° |
| −30 to 0 s | +0.32° | −0.31° |
| 0 to +60 s | +3.68° | +2.88° |
| +20 to +60 s | +5.83° | +5.39° |

After the cusp, the transverse and divergent branches coexist and their mean direction drifts away from the 35.26° cusp direction, so a window that extends past the cusp is biased. The same run also shows a distance bias (2375 m against 2500 m) and a cusp time 8 s late, which feed the 540 m track error.

**Why it matters**: the README and CHANGELOG say interferometry is "better on every measure" for the ferry and that the Kelvin method "carries a fixed course bias of about 2°". The numbers are correct, but the comparison is between a carefully tuned new method (Huber rejected, Cauchy chosen, frequency refinement added) and a baseline whose main error could be removed by changing one window. The benchmark has not shown that the bias is intrinsic to the Kelvin method.

**Recommendation**: retune the baseline (direction window ending at the cusp, and a check of the chirp fit's distance bias), rerun the ferry benchmark, and reword the conclusion to match. If interferometry still wins, the claim is stronger. The effect of the retuned window in wind is untested here.

### HYDAR-002 (High): the benchmark script needs the `fork` start method

**Where**: `examples/07_benchmark_interferometry.py:60`, `:81`, `:107`, `:192`. The wakes and solver records are stored in the module-level dict `STATE` in the parent process, then read inside `Pool` workers.

**Evidence**: with one dummy job, `multiprocessing.get_context("fork")` runs, while `"spawn"` and `"forkserver"` both raise `KeyError: 'ferry'`. `spawn` is the default on macOS and Windows; `forkserver` is the default on Linux from Python 3.14.

**Why it matters**: this script is the only way to regenerate `docs/data/ci_benchmark.csv`, the data behind case study 4. It works on Linux with Python 3.13 and earlier, and nowhere else.

**Recommendation**: pass the shared state through the `Pool(initializer=...)` argument, or build it in each worker behind a cache. Add a `N_SEEDS=1` smoke test.

### HYDAR-003 (Medium): geometry claims rest on 20 seeds without uncertainty

**Where**: README, "What the benchmark shows", and CHANGELOG v1.1.0.

**Evidence**: the benchmark uses the same seeds for every geometry, so the runs can be compared in pairs. For interferometry on the ferry, the optimal 20 m array against the Rätsep grid:

| H_s | Seeds where the 20 m array is better | Sign test p | Median ratio |
|---|---|---|---|
| 0 | 20 of 20 | < 0.001 | 0.91 |
| 0.15 m | 14 of 20 | 0.12 | 0.79 |
| 0.3 m | 12 of 20 | 0.50 | 0.78 |
| 0.45 m | 14 of 20 | 0.12 | 0.91 |
| 0.6 m | 9 of 20 | 0.82 | 0.85 |

Pooled over the four windy levels the 20 m array wins 49 of 80 runs (p = 0.057). The README says it "cuts the interferometry error further at every wind height, by 9 to 22 % in the median". The medians are right, but in wind the benefit is not significant at any single level. In calm water the 20 "runs" differ only by sensor noise and are effectively one run. Conversely, the README says the same-footprint optimal array helps "in calm water but not in wind", yet at H_s = 0.6 m it wins 15 of 20 runs (Wilcoxon p = 0.002). For the boat at H_s = 0.1 m, the quoted 4.4 m against 22 m rests on 11 of 15 paired runs (p = 0.12).

**Recommendation**: raise `N_SEEDS` to 100 or more for the geometry comparison, report paired differences with confidence intervals, and soften the wording until then. The method comparison (interferometry against Kelvin) is not affected: those differences are large and consistent.

### HYDAR-004 (Medium): the benchmark hides the reason for a failed fix

**Where**: `examples/07_benchmark_interferometry.py:100` and `:127` (`except Exception`), `hydar/kelvin.py:342` to `:345`.

**Evidence**: the CSV has 7 ferry rows where the Kelvin method has `ok = 0`. Two were rerun (Rätsep grid, H_s 0.3 m seed 11, and H_s 0.45 m seed 17). Both fail with `ValueError: first array argument cannot be empty` from `np.correlate` inside `wave_direction`: the chirp fit put the cusp or the late transverse window outside the 1800 s record, and `locate_vessel` does not check. A related path, `fit_chirp` with no ridge points, raises `IndexError`.

**Why it matters**: these are real failures of the Kelvin pipeline, so the count "293 of 300" stands. But a blanket `except` would count a programming error in either method the same way, and the CSV does not record why a run failed.

**Recommendation**: make `locate_vessel` raise a `RuntimeError` with a clear message when a window falls outside the record or no ridge is found; catch only `RuntimeError` in the benchmark; add a `reason` column.

### HYDAR-005 (Medium): same-footprint designs are set by `min_spacing`

**Where**: `hydar/interferometry.py:180`, `design_array`; README array-design table.

**Evidence**: for the ferry band with a 7.08 m footprint:

| `min_spacing` | Column gaps (m) | Row gaps (m) | Bearing bound relative to the Rätsep grid |
|---|---|---|---|
| 0.25 | 0.25, 4.69 | 4.82, 0.25 | 0.887 |
| 0.5 (README) | 0.50, 4.45 | 4.55, 0.50 | 0.907 |
| 1.0 | 1.00, 3.96 | 4.05, 1.00 | 0.944 |

The optimizer always pushes the middle column and row to the minimum spacing. The "optimal" geometry is therefore the constraint, not an interior optimum, and the 5 to 9 % gain quoted in the README scales with a default that has no stated physical basis. The ferry designs use 0.5 m and the boat design 0.25 m; the README does not say so. The bound also assumes noise that is independent between sensors, which a wind sea is not, and this agrees with the benchmark finding that the same-footprint design does not help in wind.

**Recommendation**: state in the README that the design is constraint-bound, give the minimum spacing used in each case and a reason for it (sensor housing, or the coherence length of the noise), and consider an objective that includes coherent noise.

### HYDAR-006 (Medium, tracked): deep-water waves at a 3 m deep site

**Where**: `hydar/kelvin.py` (`KelvinShip.wake`, `wind_sea`), `record_passage(..., depth=3.0)`; JOURNAL Open threads, "Ferry wake synthesis vs depth"; README Assumptions for case studies 3 and 4.

**Evidence**: at 0.1 Hz the synthesized wave has the deep-water wavenumber k = ω²/g = 0.040 m⁻¹, a wavelength of 156 m and a phase speed of 15.6 m/s. In 3 m of water a 0.1 Hz wave travels at about √(gh) = 5.4 m/s with a wavelength near 54 m. The array therefore sees a spatial phase pattern that cannot exist at that depth, while the pressure attenuation uses the 3 m dispersion relation. Both trackers succeed because they analyse the record with the same deep-water relation that made it.

**Why it matters**: this is the largest gap between the ferry results and field data. On a real record at this site, the delays across the array would correspond to about 5 m/s, and the waves would have refracted on the way in.

**Recommendation**: as the journal notes, a finite-depth wake synthesis. Until then, keep the limitation next to the ferry results, as the README does now.

### HYDAR-007 (Medium): `SpectralSolver` ignores volume sources without warning

**Where**: `hydar/spectral.py:155` to `:158` (`step` calls `obj.source(self)` and discards the result).

**Evidence**: a `DropObject` (0.5 m, 2500 kg/m³) added to a `SpectralSolver` and run for 1 s sinks (`V_sub` = 0.0654 m³) while `max|eta|` stays exactly 0. No error or warning is raised.

**Why it matters**: the README says "the dispersive `SpectralSolver` could also model the drop case". A user who tries it gets a flat surface and no explanation.

**Recommendation**: raise a `TypeError` when an object returns a non-`None` source in `SpectralSolver`, or implement the volume source (it enters the η equation directly).

### HYDAR-008 (Medium): no continuous integration; minimum versions untested

**Where**: no `.github/`; `pyproject.toml` declares `python >= 3.9`, `numpy >= 1.22`, `scipy >= 1.8`, `matplotlib >= 3.5`.

**Evidence**: the tests have only been run by hand, on numpy 2.5.3 and scipy 1.18.1 (JOURNAL, and this review). Nothing checks the lower bounds, and nothing runs when a pull request is opened.

**Recommendation**: a GitHub Actions workflow that runs the fast tests on the oldest and newest supported versions for every pull request, and the slow tests weekly.

### HYDAR-009 (Medium): published results without a regression test

**Where**: `tests/`.

**Evidence**: no test exercises the boat path of case study 4 (finite depth, `time_stamp="centroid"`, four arrays fused), although it produces half of the README benchmark. `pick_wake`, `wave_direction`, `detection_statistic`, `wind_sea`, `wind_sea_head`, `SensorArray(layout=...)` attached to a solver, and `hydar.plotting` have no direct tests.

**Recommendation**: add a slow test for the boat interferometry track in calm water, and unit tests for `pick_wake` and `wave_direction` on synthetic plane waves.

### HYDAR-010 (Low, tracked): citation metadata between releases

`CITATION.cff` has `version: v1.1.0` and `date-released: 2026-11-09` with the v1.0.0 DOI, whose description reads "Zenodo archive of HYDAR v1.0.0". GitHub's "Cite this repository" button currently produces a v1.1.0 citation with the v1.0.0 DOI. This is the open thread waiting on the release. When the new DOI arrives, consider using the Zenodo concept DOI (all versions) for the README badge and `pyproject.toml`, so that they do not go stale again.

### HYDAR-011 (Low): sensor noise is filtered by the depth response

`record_passage` (`hydar/kelvin.py:151`) and `run_ferry` in the benchmark add white noise as an equivalent surface elevation before the `airy` pressure model applies the depth filter. A record of pure noise with the stated 2 Pa comes out at 0.139 Pa RMS; within the analysis band the amplitude ratio is 0.97 at 0.07 Hz, 0.87 at 0.15 Hz and 0.51 at 0.3 Hz. A real sensor's noise is not attenuated by the water above it. The effect on the results is negligible, because the wind sea gives about 100 Pa RMS at the sensor (H_s = 0.15 m) and the wake about 500 Pa, but the documented noise level is not the applied one. The boat benchmark uses the hydrostatic model and applies the noise correctly.

### HYDAR-012 (Low): edge inputs

| Call | Result |
|---|---|
| `fit_chirp(np.array([]), ...)` | `IndexError: too many indices` |
| `DropObject` placed outside the grid | `IndexError` in `attach`, before the intended warning |
| `intersect_bearings` with parallel bearings | `LinAlgError: Singular matrix` |
| `pick_wake(t, p)` with a 1-D record of 2000 samples | returns 2000 NaNs instead of one pick |
| `pick_arrivals` with an empty time window | `ValueError` from `max` of an empty array |

### HYDAR-013 (Low): documentation drift

- README quick start: "22 tests, about 4 minutes". There are 27 tests and the full suite took 4.8 minutes here.
- `hydar/interferometry.py:33`: the module docstring says "RANSAC, then Huber least squares"; the code uses a Cauchy loss, and the journal records Huber as rejected.
- `pyproject.toml` description and the `hydar/__init__.py` docstring describe only the drop case.
- The README does not state that the boat design uses `min_spacing=0.25` (see HYDAR-005).
- The 20 m design reports a peak sidelobe of 0.61 against a limit of 0.6, because the optimizer evaluates the constraint at a coarser resolution than the final report. Reported honestly, but worth a sentence.

### HYDAR-014 (Low): outputs written over tracked or into untracked places

`examples/07_benchmark_interferometry.py` accepts `ferry` or `boat` and `N_SEEDS`, but always writes `docs/data/ci_benchmark.csv`, so a partial run replaces the full committed file. Example 03 run from the repository root leaves `hydar_map.png` and `hydar_timeseries.png` there; `.gitignore` covers only `examples/*.png`.

### HYDAR-015 (Low): scenario constants inside library functions

`fit_chirp` starts from distances of 500 to 4500 m; `wind_sea` draws frequencies from 0.25 to 0.7 Hz whatever `f_peak` is; `stacked_spectrogram` defaults to sensors (0, 2, 4, 6, 8) of a 3 × 3 grid; `locate_vessel` fixes its direction windows; `wave_direction` sizes its lag window for 2 m/s. These fit the ferry scenario and would fail quietly on others. Expose them as arguments or document the range of validity.

### HYDAR-016 (Low): duplication

`G = 9.81` is defined in `solver.py`, `kelvin.py` and `interferometry.py`, and again as default arguments. The overflow-safe cosh ratio is written three times (`sensors.depth_response`, `spectral.pressure_head`, `interferometry.wind_sea_head`). `wind_sea` and `wind_sea_head` are near copies. `kelvin.record_passage` hard-codes ρ = 1000. `observe(peak_rel=...)` cannot be set through `track_vessel`, whose extra keywords go to `ci_scan`.

### HYDAR-017 (Low): repository weight and figure numbering

`docs/img` holds 32 MB, including a 12.6 MB GIF, and each figure in PNG, PDF and SVG; every regeneration adds the full size to the history. Git LFS or release assets would keep clones small. Separately, figures carry three numbers: the `make_figures.py` argument (1 to 6), the file name (`fig01` to `fig08`) and the README caption (Figure 1 to 7). For example `python docs/make_figures.py 4` writes `fig05` and `fig06`, which the README calls Figures 4 and 5.

### HYDAR-018 (Low, accepted risk): licence

CC BY 4.0 is the author's choice and matches the Zenodo record. Creative Commons advises against its licences for software, because they do not address source and object code or patents. Recorded here for completeness; no action unless the author wants to revisit it.

## 9. Positive observations

- **The benchmark is exactly reproducible.** All seeds are explicit, and 60 rerun rows match the committed CSV in every digit.
- **The README reports results honestly**, including where the new method loses (the boat in moderate wind) and where both fail. Every table value traces to the CSV.
- **The numerics are sound.** The exponential integrator in `SpectralSolver` is exact for each mode, the pressure at depth follows from the same potential, and the tests check dispersion, absorption and volume conservation against analytic values. A control test shows that rigid walls keep the energy, which validates the PML test itself.
- **The boat case of the benchmark is a fair test**: its data come from the dispersive solver, a model independent of the tracker's assumptions.
- **The journal records rejected alternatives** (Huber loss, bin-centre frequencies, RANSAC on stationarity alone), which made this review faster and several suspicions easy to rule out.
- Local links resolve, no secrets, no stale TODO markers, clean packaging.

## 10. Hypotheses tested and not confirmed

- **A wind sea of 80 constant-amplitude sinusoids might flatter the detectors**, because each spectral bin then holds a nearly pure plane wave of constant power. Test: the ferry scenario on the Rätsep grid, 6 seeds, with 80 and with 1000 components. Interferometry track error, median: 171 m against 132 m at H_s = 0.3 m, and 206 m against 225 m at 0.6 m; the same number of observations and inliers in both. No effect.
- **The pair-phase sign convention in `ci_scan`** was checked by derivation against numpy's FFT convention. It is correct.
- **`SensorArray.load` with a record starting at t = 0** duplicates the first time stamp internally; the interpolated output is nevertheless correct.

## 11. Recommended remediation order

1. HYDAR-002 (benchmark start method), so that the benchmark can be rerun anywhere.
2. HYDAR-001 and HYDAR-004 together: retune and guard the Kelvin baseline, rerun the ferry benchmark, update the README conclusion.
3. HYDAR-003: more seeds and paired statistics in the same rerun.
4. HYDAR-005 and HYDAR-013: README wording.
5. HYDAR-007, HYDAR-008, HYDAR-009: guard, CI and tests.
6. The Low findings as convenient; HYDAR-010 when the Zenodo DOI exists.

## 12. Not checked

- `examples/04_boat_tracking.py`, `examples/06_array_design.py`, the two animation scripts and `docs/make_figures.py` were not run (about 30 minutes in total). The slow tests cover the same code paths as examples 04 and 06 at smaller size.
- The boat half of the benchmark was not rerun (about 6 minutes of solver time plus the Monte Carlo runs); its CSV rows were checked against the README only.
- The references were not verified against the publishers' records (the journal lists three as unverified by the author).
- The Zenodo record and the GitHub release state were not inspected beyond confirming that only v1.0.0 is published.
- No field data exist in the repository, so nothing here speaks to performance on measurements.
- The retuned Kelvin window of HYDAR-001 was tested in calm water only.

## 13. Review history

| Date | Commit | Reviewer | Critical | High | Medium | Low | Notes |
|---|---|---|---|---|---|---|---|
| 2026-11-14 | `8d45565` | Claude (repository-review skill) | 0 | 2 | 7 | 9 | First review |

## 14. Status conventions

`Open`, `Partially resolved`, `Resolved` (with the fixing commit), `Accepted risk`, `Not applicable`. IDs are never reused or renumbered.

## 15. Source changes made

None. A remediation series is a separate step, on request.
