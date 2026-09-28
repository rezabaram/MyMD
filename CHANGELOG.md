# Changelog

Notable changes, newest first.  The repository had no changelog before this
point; the entries below cover the port and modernisation work.

## Unreleased

### Fixed

- **The translational velocity corrector in the Beeman integrator was wrong.**
  The `a_{n+1}` term had the wrong sign and the `a_{n−1}` term was missing, so
  translation was first-order accurate where rotation was second-order.  On a
  perfectly elastic bounce this lost 5% of the energy per bounce at `dt=1e-4`,
  converging only linearly with `dt`.  Now exact to four decimals.  Guarded by
  the new `elastic_bounce` regression case.
- `method restart` was broken: `celllist.build()` ran before `celllist.setup()`,
  so `which()` divided by an uninitialised `dx` and the run died with
  `(i,j,k): 2147483647 2147483647 2147483647`.
- Deposition aborted with `Point out of grid` once the box filled up: a new
  layer was seeded at a height that could exceed the lid.  It now reports that
  the box is full and places what fits.
- `zeta` and `zetaWidth` were ignored for `particleType general` — `zeta` was
  computed from `eta`, so the two shape parameters were not independent.
- The first line of `log_energy` was uninitialised garbage (`~1e-314`); the
  energies had never been computed when it was written.
- The solver printed `Relaxation criterion reached ... KE < 1e-8` on every exit
  including a normal `maxTime` exit, with KE orders of magnitude above 1e-8.
- `outDensity`, `e`, `friction_threshold` and `read_radii` were registered
  parameters that nothing read.  Removed.
- `cohesion` and `static_friction` are accepted, plumbed into the material and
  used by no force law.  Setting either now warns instead of silently doing
  nothing.
- `config_base` and `config_relax` set `particleShape` and `read_radii`, and
  `config_base_periodic` set `initialization` — none of them registered
  parameters, all of which warned on every run.
- `config_relax` pointed at an absolute `/home/reza/...` path for its radii file.

### Removed

- The raster3d and POV-Ray pipeline (`bin/coord2pov`, `bin/coord2pr3d`,
  `bin/genFrames.sh`, `bin/encodejpg.sh` and the `animate`/`movie`/`pov`
  Makefile targets), superseded by the tooling in `docs/VISUALIZATION.md`.
- The Condor/PBS job scripts, which had the author's absolute paths and a
  cluster-specific queue baked in.
- `tools/correlation.cc`, `tools/coord2xdr.cc`, `tools/asphericity.cc` and
  `tools/sphere_map.cc`, none of which could build.
- `CComposite` (composite particles) and `CCylinder`, both unfinished and never
  instantiated, together with their entries in the `GType` enum and the
  collision dispatch.
- The compiled-out Verlet neighbour-list subsystem (`CVerletManager` and every
  `#ifdef WITH_VERLET` block) and the orphaned `verletfactor` parameter.
  `CVerletList` is kept — despite the name it is the per-particle contact cache,
  and it is load-bearing.

### Added

- `bench/check_physics.py` and `make check`: three physics regression cases with
  calibrated tolerances, plus invariants and a case-specific energy check.
- `bench/run_bench.py`, `bench/configs/` and `bench/results.json`, so the
  performance baseline can be re-measured and diffed.
- `tools/web_viewer.py`, `tools/ovito_reader.py`, `tools/snapshot_to_dump.py`,
  `tools/ellipmd_io.py` and `tools/orientation_sample.py` — interactive HTML,
  OVITO and LAMMPS-dump visualisation, replacing the POV-Ray route.
- `docs/`: `ARCHITECTURE.md`, `PHYSICS.md`, `FORMAT.md`, `PORTING-NOTES.md`,
  `BENCHMARK.md`, `VISUALIZATION.md`.
- `README.md`, replacing the plain-text `README`; `ROADMAP.md`, the improvement
  plan; this file.

### Changed

- Build ported to macOS/GCC 14.  `Makefile.inc` now detects GSL (bundled
  `.deps/gsl` first, then Homebrew), selects `g++-14`, builds as `-std=gnu++98`
  and separates compile flags from link libraries.  `make gsl` fetches and
  builds the dependency.
- `include/size_dist.h`: nine functions were declared to return a value and fell
  off the end, which modern GCC turns into a trap — any configuration with a
  `SizeDistribution` line died with `SIGTRAP`.
- `include/packing.h`, `include/celllist.h`: qualified dependent-base calls with
  `this->` for modern two-phase name lookup.

---

## Pre-history

Everything before the above is the original research code, as committed by the
author between 2011 and 2013.  `git log dff95dd` is the last commit of that era.
