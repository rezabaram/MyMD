# MyMD — molecular dynamics of ellipsoidal packings

A soft-particle molecular-dynamics code for simulating packings of ellipsoids
under gravity.  Written for postdoc research on the structure of granular and
colloidal packings of non-spherical particles; released in the hope it is
useful to somebody else.

Particles interact through a viscoelastic (Hertzian) contact law with dynamic
friction, are integrated with a Beeman scheme including rotation via
quaternions, and can be deposited under gravity or compressed from a lattice
into a dense packing.  Boundaries may be solid walls or periodic in any
combination of axes.

```
   deposition                              Stillinger
   ──────────                              ──────────
   particles rain onto a floor             lattice, then compressed
   and settle into a packing               until the target density
```

---

## Quick start

```sh
make gsl          # one-off: fetch and build the GSL dependency (~2 min)
make ellipmd      # build the solver
make check        # physics regression tests (fast, run this after any change)

make run CONFIG=config_quick    # a ~1 s smoke test: 100 spheroids settling
open viz/trajectory_large.html  # look at it  (see docs/VISUALIZATION.md)
```

`make run` without `CONFIG=` looks for a file called `config` and falls back to
the compiled-in defaults, which are not a working configuration — always pass
`CONFIG=`, or copy one of the `config_*` files to `config`.

## Requirements

| | |
|---|---|
| compiler | **GCC 14** (Homebrew `g++-14` on macOS).  The code is C++98-era and uses `<tr1/random>` and dynamic exception specifications, which libc++/C++17 removed.  See [`docs/PORTING-NOTES.md`](docs/PORTING-NOTES.md). |
| GSL | for the 4×4 eigensolver used in ellipsoid contact detection.  `make gsl` builds it into `.deps/`; `brew install gsl` works too and is picked up automatically. |
| Python 3 | for the analysis, viewer and test tooling.  Standard library only, no numpy needed. |

Nothing else: no CMake yet, no exotic dependencies, no raster3d or POV-Ray.

## Configuration

A run is controlled by a flat `key value...` text file, passed as the second
argument.  `#` starts a comment.  The first argument is a seed:

```sh
./ellipmd <seed> <config-file>
```

`config_quick` is the smallest working example; `config_deposition`,
`config_stillinger`, `config_base` and `config_base_periodic` are the author's
originals.

### Box and boundaries

| key | default | meaning |
|---|---|---|
| `boxsize` | `1 1 2` | box edge lengths `Lx Ly Lz` |
| `boxcorner` | `0 0 0` | corner of the box |
| `boundary` | `solid` | `solid`, `periodic_x`, `periodic_xy` or `periodic_xyz`.  A periodic axis is implemented by marking that pair of walls non-solid and letting the cell list use image shifts. |

### Particles

| key | default | meaning |
|---|---|---|
| `nParticle` | `5` | target particle count |
| `particleSize` | `1` | nominal radius `r0`.  The shape parameters below redistribute the semi-axes while holding the **volume** at `(4/3)πr0³` constant |
| `particleSizeWidth` | `0` | relative width of the volume distribution |
| `particleType` | `general` | see below |
| `eta` | `1.0` | shape parameter, `c/a` |
| `etaWidth` | `0` | spread of `eta` (relative, truncated Gaussian) |
| `zeta` | `1.0` | shape parameter, `a/b` |
| `zetaWidth` | `0` | spread of `zeta` |
| `asphericity` | `-0.5` | used by `prolate` / `oblate` only |
| `asphericityWidth` | `0.1` | ditto |
| `rmin`, `rmax` | `0.05` | used by `sandstone` only |
| `radii` | `radii.dat` | semi-axis triples, one per line; used by `gen1`…`gen4` only |
| `SizeDistribution` | `mono 1.0` | `mono <r>`, `uniform <min> <max>` or `read <file>`.  **Only the `Stillinger` method reads this** |
| `scaling` | `1.0` | `Stillinger`: linear scale applied to every particle at each output step, compressing the lattice into a packing |

`particleType` selects how shapes are drawn:

| value | behaviour |
|---|---|
| `general` | polydisperse spheroids from `eta`/`zeta` and their widths — the normal choice |
| `prolate`, `oblate` | monodisperse spheroids whose aspect ratio is mapped from `asphericity` through a lookup table (`include/map_asph_aspect.h`) |
| `sandstone` | semi-axes drawn independently from `uniform(rmin, rmax)` |
| `gen1`…`gen4` | shapes read from the `radii` file |

Note the shape-parameter naming is the reverse of what one might expect:
algebraically `c/a = eta` and `a/b = zeta`.  (The top-level `README` in older
checkouts claimed `eta = a/b`; that was wrong.  See
[`docs/PHYSICS.md`](docs/PHYSICS.md).)

### Material

| key | default | meaning |
|---|---|---|
| `density` | `1.0` | particle density; mass = density × volume |
| `stiffness` | `5e2` | Hertzian contact stiffness `k` |
| `damping` | `5` | viscoelastic damping `γ` |
| `friction` | `0.2` | Coulomb **dynamic** friction coefficient |
| `fluiddampping` | `0.05` | bulk drag, proportional to `\|g\|·m·v` |
| `softwalls` | `false` | walls exert no torque on particles |
| `spherize_on` | `false` | gradually round the particles during the run |
| `cohesion` | `0` | **not implemented** — accepted, warns, has no effect |
| `static_friction` | `0` | **not implemented** — accepted, warns, has no effect |

### Solver and run

| key | default | meaning |
|---|---|---|
| `method` | `deposition` | `deposition`, `Stillinger` or `restart` |
| `timeStep` | `1e-5` | integration step `dt` |
| `maxTime` | `10.0` | stop at this simulated time |
| `Gravity` | `0 0 -10` | gravitational acceleration vector |

`method`:

* **`deposition`** — particles are rained onto the floor a layer at a time until
  `nParticle` is reached or the box fills up.  Needs `Gravity`.
* **`Stillinger`** — particles are placed on a lattice and then compressed by
  `scaling` at every output step.  Intended for periodic boundaries.
* **`restart`** — read a previous snapshot from `input` and continue.  Positions,
  orientations and shapes are restored; **velocities are not** (the snapshot
  format does not carry them), so a restarted run begins from rest.

### Output

| key | default | meaning |
|---|---|---|
| `output` | `out` | filename prefix |
| `outDt` | `0.02` | simulated time between snapshots; the number of steps between writes is `outDt / timeStep` |
| `outStart`, `outEnd` | `0.0`, `1000.0` | window of simulated time during which to write |

`outDt` and `timeStep` are independent: `timeStep` sets how finely the physics
is integrated, `outDt` how often it is written out.  Reducing `timeStep` by 10
gives ten times the work and **not** ten times the frames — if you want a smooth
animation, reduce `outDt`.

Snapshots are written as `out00000`, `out00001`, … plus a final `outend`.  The
format is specified in [`docs/FORMAT.md`](docs/FORMAT.md); `log_energy` holds
`t  E_total  E_kin  E_pot  E_rot`.

## Looking at the results

The old POV-Ray and raster3d pipelines have been removed.  In their place:

```sh
make viewer                     # self-contained interactive HTML, no installs
make dump                       # LAMMPS dump, for dragging into the OVITO GUI
make ovito-render FILE=out00010 PNG=frame.png
```

[`docs/VISUALIZATION.md`](docs/VISUALIZATION.md) covers the options, including
how to open snapshots in OVITO without OVITO Pro.

## Repository layout

```
main.cc              entry point; everything else is header-only and included
                     from include/main.h — one translation unit by design
include/             the solver: mdsys.h holds CSys, the rest is components
                     (see include/README for a map)
tools/               analysis utilities (C++) and the Python tooling
bench/               physics regression cases and the timing benchmark
data/                aspiration-ratio lookup tables and sample radii files
docs/                this documentation
config_*             example configurations
```

[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md) explains how the pieces fit
together and how a timestep flows.  [`ROADMAP.md`](ROADMAP.md) is the
improvement plan and known-issues list.

## Development

```sh
make check          # physics regression tests — run before and after any change
make bench          # full timing benchmark (~16 min, sequential)
make bench BENCH_ARGS="--only B2"     # ~80 s
make tools          # build the analysis utilities into bin/
```

`bench/check_physics.py` compares against committed reference snapshots with
tight tolerances, plus a set of physics invariants.  The tolerances are
calibrated: recompiling unchanged source with `-O3 -march=native` perturbs
positions by ~1e-12, while changing the contact force by 0.01% moves them by
~9e-6.  See [`docs/BENCHMARK.md`](docs/BENCHMARK.md) for the performance
baseline and [`ROADMAP.md`](ROADMAP.md) for what is known to be wrong.

## Known limitations

These are real and worth knowing before you trust a number:

* `cohesion` and `static_friction` are not implemented (the code warns).
* A restarted run loses velocities.
* The integrator's stability limit is not enforced — reducing `particleSize`
  without reducing `timeStep` will blow up.
* The deposition layer gate uses a hardcoded box height; see `ROADMAP.md`.
* The code is single-threaded and the hot path is dominated by a vendored
  heap-allocating matrix class; `ROADMAP.md` Phase 6 has a profile.

## Citing

<!-- TODO: add the paper(s) this code produced. -->

If you use this code in academic work, please ask the author for the
appropriate citation.

## License

There is no `LICENSE` file yet.  Every source file carries the author's notice:

> Do whatever you want with this code.  You can even replace my name with
> yours.  But you may not change the copyright itself.

Adding an explicit `LICENSE` file is on the roadmap; until then, treat the
above as the terms.

## Contact

Reza Baram — reza.baram@gmail.com
