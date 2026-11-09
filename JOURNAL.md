# HYDAR journal

The development thought tree of this repository: decisions, rejected alternatives, and what is still open.
Rules: entries are never rewritten; the newest entry sits first, directly under Open threads; only the
Open threads block (and the Current focus line in CLAUDE.md) is edited in place. A corrected entry gets a
new entry that names it, plus one line "*Superseded by the YYYY-MM-DD entry*" under the old heading.

## Open threads

- **Zenodo DOI for v1.1.0.** Jeff publishes the GitHub release v1.1.0 (Claude sessions cannot create
  releases or tags); Zenodo then mints a new version DOI. Add it to README "How to cite", CITATION.cff
  (`doi`, `identifiers`, `url`) and the bibtex block. v1.0.0 DOI is 10.5281/zenodo.23881055.
- **Boat case: interferometry loses course precision in moderate wind** (case study 4, Hs 0.05 m: course
  SD 0.8 to 1.3 deg vs 0.23 deg for `track_boat`). Only ~60 coherent observations from 6 s windows.
  Untried ideas: multi-taper cross-spectra, fusing interferometry bearings with cusp-line picks, longer
  records. Start at `hydar/interferometry.py` (`ci_scan`, `fit_track`) and `examples/07_benchmark_interferometry.py`.
- **Ferry wake synthesis vs depth.** `hydar/kelvin.py` synthesizes a deep-water wake but sensors sit in
  3 m water (Airy filter with h = 3 m). Case study 4 therefore uses deep-water dispersion for the ferry.
  A finite-depth wake would need its own synthesis.
- **New references unverified by Jeff:** Fischler & Bolles (1981, doi 10.1145/358669.358692), Lighthill
  (1978), Van Trees (2002, doi 10.1002/0471221104) in README References. Jeff said "fine for now".

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
