# HYDAR — orientation

HYDAR (Hydrodynamic Detection and Ranging, Jeff Tuhtan) is a Python simulator that locates objects in
water from the surface waves they make, using pressure-sensor arrays: drops (case study 1), boats and
ships from their wakes (case studies 2 to 4), a rock in a river from its vortex street (case study 5). See README.md for the science and examples.

**Before working, read the top of `JOURNAL.md`** (Open threads and the newest one or two entries).
Project state and decisions live in JOURNAL.md; do not duplicate them into agent-side memory or Claude
Project memory. When compacting, preserve decisions made and alternatives rejected — they feed JOURNAL.md.

**Current focus**: v1.3.0 (case study 4, correlative interferometry; the fixes from the repository review
of 2026-11-14; case study 5, a fish finds a rock from its Kármán vortex street, `hydar/karman.py`) is the
version on main (Jeff chose 1.3, 2026-11-15); v1.1.0 and v1.2.0 were never released on their own. Jeff
publishes the GitHub release v1.3.0 without a new Zenodo version (minor release), so the v1.0.0 DOI stays.

## Rules

- Releases, tags and Zenodo records are made by Jeff, not by Claude sessions.
- Licence is CC BY 4.0 (Jeff's choice); keep CITATION.cff and README "How to cite" in sync with releases.
- Report benchmark results honestly, including where a new method loses.
- Figures: `docs/make_figures.py` with `docs/paperfig.py` (170 mm text width); run `python docs/make_figures.py N`.
- Tests: `python -m pytest -m "not slow"` before every commit; slow tests run full simulations.

## Commits

Narrative style: the subject states the conclusion, the body gives the reasoning and rejected alternatives.
