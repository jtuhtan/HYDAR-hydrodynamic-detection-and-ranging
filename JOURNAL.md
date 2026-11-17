# HYDAR journal

The development thought tree of this repository: decisions, rejected alternatives, and what is still open.
Rules: entries are never rewritten; the newest entry sits first, directly under Open threads; only the
Open threads block (and the Current focus line in CLAUDE.md) is edited in place. A corrected entry gets a
new entry that names it, plus one line "*Superseded by the YYYY-MM-DD entry*" under the old heading.

## Open threads

- **GitHub repository description is out of date (Jeff to edit).** It still names only drops and boats ("shallow-water and dispersive wave solvers, drop tests and Kelvin wake boat tracking"); it should also mention ship localization, correlative interferometry (case study 4) and the fish finding a rock from its Kármán vortex street (case study 5). Only Jeff can change it, on the repository page.
- **Speed: what is left** (`repository_review/2026-11-16-performance-review.md`, sections 5 and 7). Done on 2026-11-16: fewer FFTs, single precision (benchmark 07 only), threads, the vortex street loop, `wavenumber`. Left: the Monte Carlo runs of benchmark 07 (now most of its 21 min), not profiled in depth. The boat fits in wind H_s >= 0.1 m flip between two solutions under tiny perturbations (2026-11-16 entry).
- **Zenodo for minor releases.** Jeff (2026-11-15): v1.3.0 gets a GitHub release but no new Zenodo
  version, as it is a minor release; README, CITATION.cff and the badge keep the v1.0.0 DOI
  10.5281/zenodo.23881055. If the Zenodo GitHub integration is switched on for the repository, publishing
  the release archives it automatically; check that before the release. A major release gets a new DOI;
  then update README "How to cite", CITATION.cff (`doi`, `identifiers`, `url`) and the bibtex block.
- **Boat case: interferometry loses course precision in moderate wind** (case study 4, Hs 0.05 m: course
  SD 0.9 to 1.1 deg vs 0.21 deg for `track_boat`, 100 seeds). Only ~60 coherent observations from 6 s windows.
  Untried ideas: multi-taper cross-spectra, fusing interferometry bearings with cusp-line picks, longer
  records. Start at `hydar/interferometry.py` (`ci_scan`, `fit_track`) and `examples/07_benchmark_interferometry.py`.
- **Ferry wake synthesis vs depth.** `hydar/kelvin.py` synthesizes a deep-water wake but sensors sit in
  3 m water (Airy filter with h = 3 m). Case study 4 therefore uses deep-water dispersion for the ferry.
  A finite-depth wake would need its own synthesis.
- **Chirp fit of the Kelvin method is late and short (HYDAR-020).** `fit_chirp` puts the ferry abeam about
  31 s late and 4 % too close; this leaves the Kelvin baseline with a +0.6 deg course bias and a 400 m track
  error, so the ferry factor of 2 to 6 in case study 4 is against a baseline that can still improve. A 128 s
  spectrogram window halves the error but loses the late wake. Start at `hydar/kelvin.py` (`ridge_points`,
  `fit_chirp`).
- **Repository review 2026-11-14: remaining items** (`repository_review/2026-11-14-comprehensive-repository-review.md`).
  Open: HYDAR-006 (ferry depth, above), 010 (Zenodo DOI, above), 020 (chirp fit, above); partially resolved:
  016 (`KelvinShip.wake` hard-codes rho = 1000), 017 (32 MB of figures in git history; LFS or release assets
  is Jeff's decision; case study 5 added two versions of `hydar_fish.gif`, 15 and 25 MB, on 2026-11-15). Review and remediation were merged with PR #4 (cf4f657, merge commit, so the
  commit hashes cited in the review stay valid).
- **Case study 5: distance depends on the wake-ageing model.** The fish's range comes from core growth and bed-friction decay, so a wrong eddy viscosity or roughness biases it by 1.3 to 2.8 m at 6 m; D is 23 to 29 % too large at 3 m (near wake). Untried: a fish that swims upstream along the found line, two surveys at different distances (ratio method), a body-bound flow model for the fish (it is a passive probe now), a near-wake model from the rock to the formation point at 1.5 D (the 3 m rows and the animation's straight dye both show it is missing). Start at `hydar/karman.py` (`locate_rock`, `_fit`) and `examples/09_benchmark_fish_rock.py`.
- **New references unverified by Jeff (case study 5):** Chen & Jirka (1995, doi 10.1016/0169-5983(95)00053-G), Chow (1959), Liao et al. (2003, doi 10.1126/science.1088295), von Kármán (1911), Williamson (1996, doi 10.1146/annurev.fl.28.010196.002401).
- **New references unverified by Jeff:** Fischler & Bolles (1981, doi 10.1145/358669.358692), Lighthill
  (1978), Van Trees (2002, doi 10.1002/0471221104) in README References. Jeff said "fine for now".

## 2026-11-16 — A second round of speed-ups changes results slightly; both benchmarks rerun, no conclusion changes

**Decided** (Jeff chose "More speed, numbers shift" after PR #12 merged): threads in the shallow-water loop (`numba.prange`, bitwise identical, 0.17 ms per step); two FFTs fewer per spectral step by linearity; `SpectralSolver(precision="single")`, used only by benchmark 07 (2.5x, boat records within 5e-5 relative, 0.025 Pa at most); the vortex street as one Numba loop with its own E1 (to 1e-14); `cos(ph) @ A` in `KelvinShip.wake` (the change rejected in the entry below for bit-for-bit reproducibility, now accepted); `wavenumber` stops at round-off. Benchmark 07 rerun in full (21.5 min): ferry summary cells unchanged; four boat cells in the README moved (calm interferometry median 0.73 to 0.72 m, H_s 0.05 m optimal course bias +0.86 to +0.85°, H_s 0.1 m Kelvin course bias +26.1 ± 61.1 to +24.7 ± 61.7°, speed SD 2.06 to 2.05 m/s). A few boat runs in wind H_s >= 0.1 m jump between two fit solutions (up to 143° and 608 m); no run gains or loses a fix. Benchmark 09 rerun in full (7.5 min, was about 25): CSV and tables byte-identical.
**Rejected**: single precision as the library default: keep "double" so that library users get the old numbers unless they ask; regenerating Figure 7 (changes far below a line width).
**Open**: the boat runs in wind that flip between two solutions under 1e-5 perturbations say those fits are ill-conditioned; the medians and fix counts are stable, the means and 90th percentiles of the Kelvin method at H_s >= 0.1 m are not.

## 2026-11-16 — The solvers are faster with unchanged results; Numba, not Rust or Fortran, is the accelerator

**Decided** (Jeff: review the current version, make the calculations and the numerical solver more efficient, maybe with Rust and/or Fortran): profile first, then speed up only what keeps the results bitwise identical, so the committed benchmark CSVs still reproduce. Review in `repository_review/2026-11-16-performance-review.md`. Changes: the shallow-water step is one compiled loop with Numba (0.7 instead of 8.1 ms per step on 460 x 460) and in-place NumPy without it (4.0 ms); `SpectralSolver` caches the depth factor per sensor depth, fuses the mode update and skips FFTs of an all-zero pressure (1.35x for the boat); `lamb_oseen_pressure` evaluates E1 only where r^2/rc^2 <= 50 (a run of benchmark 09 1.7x); `peak_sidelobe` uses one matrix product for all source directions (`design_array` 4x, same designs). Numba is an optional extra (`hydar[fast]`), off with `HYDAR_NO_NUMBA=1`; `hydar/_kernels.py` keeps a NumPy fallback for every kernel. CI gains a Windows job and runs with and without Numba. Fast tests: 118 s → about 50 s.
**Rejected**:
- Rust (PyO3, maturin): built and timed for the shallow-water step, 0.70 ms per step against Numba's 0.69, both bitwise identical to NumPy. The step is memory bound, so the language does not matter; Rust would add a toolchain and a compiled wheel per platform and Python version. Worth reconsidering only for processing on a sensor node, not for the simulator.
- Fortran (f2py, meson): not built (no compiler in the container); expected to tie with the others, and the hardest of the three to build on Windows.
- Compiled code for the spectral solver: about 75 % of a step is FFTs that already run in SciPy's multi-threaded C++.
- `cos(ph) @ A` in `KelvinShip.wake` (2.3x): it changes the wake by 2e-16 m and the ferry benchmark CSV in the 9th digit; not worth losing bit-for-bit reproducibility for 50 s.
**Confirmed by**: hot paths compared array by array with main (fish benchmark run, street pressure, `design_array`, 300 boat steps with sensors, shallow-water solver with and without a sphere, unforced spectral solver): bitwise identical, with and without Numba. A partial benchmark 07 run (ferry, 2 seeds) writes a byte-identical CSV. 49 fast and 5 slow tests pass.

## 2026-11-15 — Case study 5 is released as v1.3.0; v1.2.0 is folded into it like v1.1.0

**Decided** (Jeff: "Go ahead and update to 1.3"): version 1.3.0 in `pyproject.toml`, `hydar/__init__.py` and CITATION.cff (date 2026-11-15); the CHANGELOG "Unreleased" section becomes v1.3.0, and v1.2.0 is marked "prepared but not released on its own", as v1.1.0 was, since it was never tagged. README "How to cite" names v1.3.0 (this replaces the "case study 5 is unreleased" sentence of the entry below). CITATION, README and the badge keep the v1.0.0 DOI: Jeff will not archive v1.3.0 on Zenodo, as it is a minor release.
**Rejected**: keeping v1.2.0 as a separate release before v1.3.0: it was never tagged, so a second release would only add a Zenodo version nobody cites. No tag or GitHub release is made by Claude sessions.

## 2026-11-15 — The README describes all five case studies; no numbers changed

**Decided** (Jeff: go through the README and make sure it is up to date): the introduction said "four parts" and described only waves from drops and boats; it now lists five parts (models, sources, sensors, localization, benchmarks) and a table that links the five case studies. The quick-start test count is 46 fast and 51 in all (was 36 and 41), the "Idealized water" limitation covers the wind sea of case study 4 and the river of case study 5, the Theory section mentions `SensorArray(layout=...)`, and "How to cite" says case study 5 is on main but unreleased.
**Checked, unchanged**: example 8 rerun gives the README's values (f 0.0999 Hz, D 0.970 m, distance 5.99 m, y 0.39 m); the case study 5 tables agree with `docs/data/fish_rock_summary.md`; figure numbers 1 to 9 match `docs/make_figures.py`; version 1.2.0 and the v1.0.0 DOI stay until Jeff releases.
**Rejected**: changing the version or the citation to anticipate the next release (Jeff's call), and rerunning benchmarks (no number was edited).

## 2026-11-15 — Case study 5 and the README limitations wording are on main; PR #6 was carried over, not merged

**Decided** (Jeff: wrap up, commit and merge any remaining changes to main): case study 5 merged in three PRs (#7 code and benchmark, #8 animation, Figure 9 and results tables, #9 dye and translucent river in the animation). PR #6 (README limitations wording from the review thread) conflicted with the vortex-street bullet; its commit was cherry-picked onto `claude/project-thread-qmz92f` with both bullets kept and merged from there, so PR #6 itself can be closed unmerged.
**Rejected**: merging main into PR #6's branch, the usual fix: pushing to a branch another session opened was blocked without Jeff's explicit approval for that branch, and carrying the commit over reaches the same result.

## 2026-11-15 — The fish animation shows the vortices being shed from the rock with dye advected by the street

**Decided** (Jeff: show the vortex street, make the river a bit transparent so it looks like flowing water, and make it clear that the vortices are shed from the rock): `examples/make_fish_animation.py` now releases dye at four ports on the rock's flanks and advects it (RK2, 0.1 s steps, a little random diffusion, 45 s of warm-up so the street is full when the video starts) with the river velocity U plus the velocity the street's vortices and their bank images induce, the same kinematic street that gives the pressure. The dye leaves the rock as two shear layers and rolls up into the vortices, which a pressure glow alone did not show. The water surface is translucent (alpha about 0.35, opaque where a vortex core glows) over a gravel bed and the submerged half of the rock; ripples drift at U and sparse white streaks move with the flow.
**Rejected**:
- Drawing the vortices as glyphs or spirals: they would show the model's vortex list, not what a flume with dye shows; advected dye also shows the induced velocity field.
- Bending the flow around the rock (potential flow) or adding shear-layer flapping near it: the street model starts at the formation point (1.5 D), so the dye runs nearly straight for about 3 m before it rolls up; adding near-wake motion would show something the model does not have.
**Open**: the near wake (rock to about 3 D) is not modelled, which also limits the 3 m benchmark rows (see the 2026-11-14 entry).

## 2026-11-15 — Case study 5 gets an animation, a results figure and full results tables; the benchmark moves out of Figure 8

**Decided** (Jeff: an animation in the style of the others, plus a multi-plot figure and a table of the results, made with the paper-figures skill): `examples/make_fish_animation.py` renders `docs/img/hydar_fish.gif` / `.mp4` from a real survey (seed 1, Tu 10 %, rock 0.4 m off the centreline, fish 6 m behind it), in the dark style of the hero and boat animations (cyan glow of the street's pressure, purple sensors, red laser to the fix, 15 fps, 960 × 540). It opens on a close-up of the fish at the station with the strongest street, pulls back to the whole river, plays the 19 stations at speed while an inset builds the amplitude profile across the river, then fits it and fires the laser. Figure 9 (`fig10_fish_rock_results`, 164 mm, seven panels) holds the benchmark: lateral, distance and D errors and the detection rate against distance for four turbulence levels, and the three errors under each wrong wake model. Figure 8 keeps the single-survey panels a to f (132 mm high). `examples/09_benchmark_fish_rock.py --summary` writes `docs/data/fish_rock_summary.md` (Table A, 16 settings; Table B, 9 wake models) from the CSV, so the tables can be rebuilt without rerunning 625 surveys; the README keeps its shorter tables.
**Rejected**:
- Keeping the benchmark as panels g and h of Figure 8: two panels could show only the distance error, not lateral, D and detection; one figure per question reads better.
- A synthetic, idealised street for the animation: the animation uses the same survey and fix as `locate_rock`, so what it shows is what the method does.
- Constrained layout for Figure 9: it collapsed the panel grid; manual gridspec margins with the layout engine off.
**Confirmed by**: fast tests (46 passed); figures and animation frames checked by eye. `paperfig` prints no `[layout]` warning for Figure 9; Figure 8 still prints six, all legends within 0.3 em of data in panels c and e, which read cleanly at print size and were left as they are.

## 2026-11-14 — A fish finds the line to a rock from its vortex street to a few centimetres; the distance is only as good as the wake model

**Decided** (Jeff: river with a trapezoidal section, hemisphere rock 50 % above the surface, 25 cm fish, 10 sensors, 5 on the head third and 5 along the middle third, none on the tail third, rock found from its Kármán vortex): case study 5, `hydar/karman.py`, examples 08 and 09, Figure 8 (`fig09_fish_rock`), PR on branch `claude/project-thread-qmz92f`. Assumptions Claude chose (not Jeff's): 4 m bed, 1 : 2 banks, 0.25 m depth, slope 0.002, Manning n 0.035 (U = 0.48 m/s at the rock); "50 % above the surface" read as half the height, so R = 2 h = 0.5 m; shedding at the depth-averaged diameter (D = 3.83 h); St 0.2, U_c = 0.88 U, b/λ = 0.281, circulation from von Kármán's self-induced speed, formation at 1.5 D, Lamb-Oseen cores r0 = 0.1 D growing with ν_t = 0.01 U D, decay c_f U / h; banks as images; sensors evenly spaced 1.67 cm apart; turbulence RMS ρ (0.1 U)², integral scale h. The fish holds 19 stations across the river, 0.2 m apart, 60 s each.
Method: frequency from the averaged spectrum; U_c and the upstream direction from the phase slope along the body; D = St U / f; lateral position and distance from a fit of the model's lateral profiles of f and 2f (both readings f and f/2 of the strongest line are fitted; a reading needs its line at f at least 3 times the noise).
**Rejected**:
- A stationary fish: one line of sensors along the body cannot tell which side of the street it is on (the street is mirror symmetric with a half-period shift), so the fish steps across the river.
- Range from bed-friction decay alone: its length scale is about 20 m here, too slow to range over a few metres; core growth by turbulent diffusion carries most of the distance information.
- Log-amplitude residuals in the profile fit: noise-only stations dominated and the fit picked the near-field minimum (3.4 m for a true 6 m); residuals in amplitude, normalized by the noise, fixed it.
- Taking the strongest line as the fundamental after a subharmonic test on the spectrum: river turbulence near f/2 made the test unreliable; fitting both readings and comparing the misfit works.
- Turbulence pressure ρ U u′ (23 Pa RMS at 10 %): far above the ρ u′² scaling of turbulent pressure; replaced by an RMS of ρ (Tu U)² (2.3 Pa).
**Found during the work**: the Lamb-Oseen pressure first had a sign error in its (1 − e^−s)²/s term (unary minus applied after the square), which made 2f peak on the centreline; `test_lamb_oseen_pressure_is_cyclostrophic` caught it, and the benchmark was rerun after the fix. With the correct pressure, 2f is weak on the centreline and strongest near the rows.
**Confirmed by**: `tests/test_karman.py` (10 tests); benchmark of 625 runs, 25 seeds per setting (`docs/data/fish_rock_benchmark.csv`). Lateral median error 1 to 4 cm at 3 to 9 m and Tu up to 10 %, and 4 cm or less under every wrong-model variant except U_c/U = 0.95 (weak street, 0.26 m). D within 1 to 3 % at 6 to 9 m, but +23 to +29 % at 3 m (near wake). Distance: median error 0.13 to 0.29 m at 6 to 9 m with the right model; biased by −1.6/+2.2 m for ν_t × 0.5/× 2, −1.3/+2.8 m for n 0.025/0.05, −1.9 m for U_c/U 0.80; 13 of 25 detections at 12 m in 15 % turbulence.
**Open**: see Open threads (ranging model dependence, fish body flow, new references).

## 2026-11-14 — Version is 1.2.0; v1.1.0 was never released on its own

**Decided** (Jeff: "update the version as well to reflect the changes"): 1.2.0 in `pyproject.toml`,
`hydar/__init__.py` and CITATION.cff (date 2026-11-14). CHANGELOG has a v1.2.0 section with the review
changes and the new results; the v1.1.0 section keeps what was added on 2026-11-09 and is marked as not
released on its own. CITATION and README keep the v1.0.0 DOI until Zenodo archives v1.2.0.
**Rejected**: keeping 1.1.0 with a new date (the changed numbers and behaviour would hide under a version
that was already announced in the journal and in PR #3); 1.1.1 (new arguments and changed defaults are more
than a patch).

## 2026-11-14 — Review remediated: with a fair Kelvin baseline the ferry advantage is a factor of 2 to 6, not 3 to 9

**Decided** (Jeff: "Go ahead and update the repo based on the review"; PR #4):
- `locate_vessel` measures the cusp-wave direction in the 30 s before the cusp; `ridge_points` interpolates
  peaks; `wave_direction` uses the unbiased correlation and a two-pass fit (short baselines, then all pairs
  within half a period of the predicted delay).
- Benchmark: 100 seeds, any start method, `reason` column, paired seed-by-seed statistics, 2 Pa sensor noise
  added after the pressure model (`SensorArray(noise_pa=...)`).
- CI on GitHub Actions: Python 3.9 with the oldest supported numpy, scipy, matplotlib, and Python 3.13.
- 16 of the 18 findings resolved or partially resolved; two findings added (HYDAR-019 resolved, HYDAR-020 open).
**Confirmed by**: 41 tests pass (36 fast, 5 slow); full benchmark, 2700 records, no run ended in an
unexpected error. Ferry, Rätsep grid, median track error Kelvin against interferometry: 402 against 68 m
(calm), 397 against 131 m (Hs 0.3 m), 457 against 235 m (Hs 0.6 m); Kelvin course bias +0.68 deg (was +2.03).
Optimal 20 m array better than the Rätsep grid in 242 of 400 paired runs with wind (p < 0.001). Optimal
same-footprint array: 211 of 400 (p = 0.29), so no benefit in wind. Boat: optimal grid worse in calm water,
better at Hs 0.025 m (67 of 100), not significant at 0.05 and 0.1 m.
**Rejected**:
- Keeping the tapered ("biased") cross-correlation in `wave_direction`: it hid cycle skipping on wide arrays
  but overstated the phase speed by 5 to 46 % on plane waves.
- A 128 s spectrogram window for the chirp fit: halves the abeam-time error, loses the late wake.
- Rejecting every record whose analysis window touches the record edge: broke `test_locate_passing_ferry`;
  the guard now needs half of the window inside the record.
- `actions/setup-python` for Python 3.9: not available on the current runner image; CI uses uv instead.
**Found during the work**: the first 100-seed run gave kilometre-scale Kelvin errors on the 20 m array
(cycle skipping after the taper was removed); fixed by the two-pass fit and rerun. The full benchmark takes
about 30 minutes on 4 cores.
**Open**: see Open threads (chirp fit, release version).

## 2026-11-14 — Whole-repository review: results reproduce; the ferry comparison uses an untuned Kelvin baseline

**Decided**: review recorded at commit 8d45565 with no source changes; remediation is a separate step on
Jeff's request. 0 Critical, 2 High, 7 Medium, 9 Low.
**Confirmed by**: all 27 tests pass; every README benchmark table value recomputed from
`docs/data/ci_benchmark.csv`; 60 rerun ferry rows (seeds 0 and 1) bit-identical to the CSV; the three array
designs reproduce.
**Found**:
- `locate_vessel` measures the cusp-wave direction from 10 s before to 30 s after the cusp; in calm water
  that window gives the +2.03 deg course bias, and a window of 30 s before the cusp gives +0.32 deg. The
  claim "interferometry wins on every measure" for the ferry is against this untuned baseline (HYDAR-001).
- `examples/07_benchmark_interferometry.py` shares state with `Pool` workers through a module global, so it
  raises KeyError under `spawn` and `forkserver` (macOS, Windows, Python 3.14 and later) (HYDAR-002).
- The 7 Kelvin "no fix" rows are a ValueError from an empty window in `wave_direction`, caught by a blanket
  `except` (HYDAR-004). Same-footprint designs sit on the `min_spacing` bound (HYDAR-005). The 20 m array's
  benefit in wind is not significant seed by seed with 20 seeds (HYDAR-003).
**Rejected**: the suspicion that a wind sea of 80 constant-amplitude components flatters the detectors;
1000 components gave the same track errors and inlier counts (6 seeds, ferry, Rätsep grid).
**Open**: all 18 findings (Open threads).

## 2026-11-09 — Release v1.1.0 prepared; version bumped and changelog notes written

**Decided**: version 1.1.0 in `pyproject.toml`, `hydar/__init__.py`, CITATION.cff (date 2026-11-09);
CHANGELOG v1.1.0 with the case study 4 results table. CITATION and README keep the v1.0.0 DOI until
Zenodo archives v1.1.0 (see Open threads). Jeff asked for the bump and notes in the case study 4 thread.

## 2026-11-08 — Interferometry tracks a ship without Kelvin theory; it wins on the ferry, not in general on the boat

**Decided**: case study 4 (`hydar/interferometry.py`, PR #2, merged f6c938a). Tracking uses only wave
stationarity (n · V = c) and group-velocity rays; directions come from correlative interferometry over all
36 pairs; RANSAC on both equations then a Cauchy-loss fit. `design_array` optimizes a 3 x 3 grid's
width, height and middle column/row position (bearing CRLB objective, peak sidelobe <= 0.6, footprint and
min-spacing limits, differential evolution). Ferry designs at the Rätsep 7.1 m footprint and at 20 m (the
20 m limit is an assumption, not Jeff's spec); boat design at the 2.8 m footprint of the 1 m grid.
Results reported honestly in README: ferry, interferometry better on bias, precision and wind robustness;
boat, better only in calm/light wind.
**Rejected**:
- RANSAC on the stationarity equation alone: discrete wind-sea components repeat across windows and
  fit a false velocity. Adding the ray equation to the inlier test rejects them.
- Bin-centre frequencies: 1/60 Hz bins biased c and c_g by up to 5 % on the chirping transverse branch;
  now zero-padding x4 plus parabolic peak refinement.
- Huber refinement: near-cusp observations (a caustic, ray theory poor) pulled the ferry track ~250 m
  along-track; Cauchy loss with MAD scales fixed it.
- rho^2/(1 - rho^2) weighting as default: hurt calm cases (boat 0.72 -> 2.28 m); kept as `weighted=False` option.
- Power-centroid time stamps as default: helped the boat (source moves 18 m per window) but doubled the
  ferry error; default is window centre, boat uses `time_stamp="centroid"`.
- Boat apertures of 4 to 6 m: no feasible grid meets the sidelobe limit up to 1.2 Hz.
**Confirmed by**: `tests/test_interferometry.py`; 540-run benchmark `examples/07_benchmark_interferometry.py`
-> `docs/data/ci_benchmark.csv` (20 seeds per wind level, numpy 2.5.3, scipy 1.18.1).
**Open**: boat course precision in wind; ferry depth mismatch (Open threads).

## 2026-11-02 — v1.0.0 released with three case studies

Drop localization (shallow-water `Solver`, PML walls), boat tracking from the Kelvin cusp line in the
dispersive `SpectralSolver` (case study 2 redone from the Mach-cone model at Jeff's request for the
correct wave physics), and ferry localization after Rätsep et al. 2020/2021 from a synthesized Havelock
wake (case study 3). Zenodo DOI 10.5281/zenodo.23881055, licence CC BY 4.0 (Jeff's choice).
(a910bd4, c5b34af, 4dc3733, 7634d58)
