# Changelog

Notable changes, newest first.  The repository had no changelog before this
point; the entries below cover the port and modernisation work.

## Unreleased

### Fixed

- **`CEllipsoid::gradient()` allocated a matrix on every call.**  It returned
  `2.0*ellip_mat*(X-Xc)`, which forms the scalar product as a matrix temporary
  before multiplying by the vector; `setcontact()` evaluates two per contact.
  `2.0*(ellip_mat*(X-Xc))` is the same arithmetic and allocates nothing.
- **`characteristic_polynomial()` called `matrixT::Det()`**, which copies the
  matrix and clones the copy to pivot on -- five heap allocations, once per
  candidate contact pair, roughly 420,000 times in one deposition run.  `det4()`
  expands along the first row instead and agrees with `Det()` to 5.6e-16
  relative over 20,000 random 4x4 matrices.

- **A missing config file was a warning, not an error.**  The run then
  continued with the compiled-in defaults, which are not a working
  configuration (`particleSize=1` in a `1x1x2` box), and died later with a
  confusing `Point out of grid`.  It now fails immediately and says what to
  pass.  A working `config` is also committed, so `make run` works on a fresh
  checkout.

- **`CParticle::`** carried an `avgforces` member that `addforce` overwrote
  with the last contact force (rather than accumulating), using a `static vec
  prev` shared by every particle; `avgtorque` was declared and never used.  The
  only reader was an `if(0)` block.  Removed.
- **`TNode::normal_fabric_tensor()`'s cache guard was dead and wrong**: it
  tested a local `bool calculated=false` and returned `branch_fabric_M` rather
  than `normal_fabric_M`.  There is no cache; the dead guard is gone.
- **`CException` was caught by value** at all eight call sites, slicing the
  type.  Now caught by reference.
- **The deposition layer gate used a literal `1.` as the box height**, left
  over from a 1x1x1 box.  It is now the box top from the geometry.  Same
  behaviour for a 1x1x1.2 box, correct for anything else.
- **Function-local `static` scratch variables** in `polynom.h` (five in the
  quadratic/cubic/quartic solvers) and `particle.h` (the rotational integrator)
  are now locals.  Nothing about them needed to outlive the call, and sharing
  them is what stood in the way of ever parallelising the force loop.

- **Heap buffer overflow in `DisBetaDistribution`.**  `bins` was allocated with
  `nbins` slots but every loop in the constructor indexed `0..nbins` inclusive,
  so the last write and read were one `double` past the end of the array.  This
  ran on every initialisation of the general particle type; it is the first
  thing AddressSanitizer reports.
- **Invalid downcast in `CBaseConfig::get_param<T>`.**  An unchecked
  `static_cast<CParam<T>*>` reinterpreted whatever was stored, so asking for
  `int` where the parameter was registered as `unsigned int` relied on the two
  layouts agreeing.  It is now a checked `dynamic_cast` that fails loudly.
  UBSan flags the original; the Stillinger path was hitting it via
  `get_param<int>("nParticle")`.
- **`CSizeDistribution` had inverted ownership.**  The destructor was
  `if(!p_dist) delete p_dist` (so it leaked), the compiler-generated copy shared
  the pointer (so the obvious fix would have double-freed), and
  `CBaseDistribution` had no virtual destructor (so deleting through the base
  pointer was undefined behaviour).  Now: virtual destructor, a `clone()`
  protocol, and real rule-of-three copy semantics.
- **`CUniformDist` was constructed from uninitialised members.**  The
  member-initialiser list built `unif(min, max)`, but members initialise in
  *declaration* order and `unif` was declared first, so it was built from
  garbage.  Harmless in practice only because `parse()` immediately reassigns it.
- **Destructors could terminate the process.**  `CSys::~CSys` used the
  TRY/CATCH macros, whose handler rethrows; a destructor is implicitly
  `noexcept` since C++11, so any exception there calls `std::terminate` instead
  of propagating.

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
- The regression reference snapshots were themselves gitignored (the repo's
  `*out*` rule matched `expected/outend`), so `make check` passed locally while
  having no reference data in a fresh clone.  They are now stored as
  `expected/final.snapshot` and `bench/reference/**` is exempted from the
  output rules.

### Removed

- Three dead globals: `gout` (an `ofstream*` nothing ever wrote through), a
  global `double friction` (nothing read it — the coefficient that matters is
  `CMaterial::friction`), and a global `vec G` that `CSys::G` shadowed inside
  the solver, making it invisible as well as unused.
- `include/eigen.h` renamed to `include/gsl_eigen.h`.  It is a wrapper around
  GSL's eigensolver, and the old name collided with the Eigen library.

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

- **Snapshots carry velocities and a header.**  The `id 14` record gained six
  fields (`vx vy vz wx wy wz`) and each file starts with
  `# ellipmd <version>  t=<time>`.  Appended and prepended respectively, so a
  reader that takes the first ten fields and skips non-`6`/`14` lines still
  works on both old and new files.
- **`method restart` now works.**  It restores the velocities, resumes the clock
  from the header (it used to be handed the whole of `maxTime` again, so a
  restart from `t=0.25` with `maxTime=0.3` ran six times too long) and primes the
  accelerations.  Not bit-exact: Beeman needs `a_{n-1}`, which nothing stores.
  Measured 8.3e-8 of divergence after one step, 3.9e-4 after 500 -- see
  `docs/FORMAT.md`.

- **A command line.**  `ellipmd --help`, `--version`, `-c/--config`,
  `-s/--seed`, `-D/--set KEY=VALUE` (repeatable), `-o/--output`,
  `--print-config`, `--save-config` and `--no-save-config`.  The original
  positional form still works, since `bin/run.sh` and the Makefile use it.
  `--set` is verified bit-identical to editing the config file.
- **Provenance.**  Every run writes `config.used` next to its output: the
  effective parameters (defaults included) under a header recording the
  version, seed, source config file and timestamp.  A result can now be traced
  back to the inputs that produced it.

- `bench/check_physics.py` and `make check`: three physics regression cases with
  calibrated tolerances, plus invariants and a case-specific energy check.  The
  `elastic_bounce` case exists specifically to catch integrator errors, which a
  snapshot comparison cannot see.
- A clean-checkout verification: `git archive HEAD` into an empty directory,
  `make ellipmd`, `make tools`, `make check` and the Python tooling all succeed
  with nothing but a compiler and GSL.
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

- **C++17.**  `<tr1/random>` is gone (TR1's `ranlux64_base_01` becomes
  `std::ranlux48_base`) and `matrix.h`'s dynamic exception specification is
  removed, which together mean the code no longer needs GCC specifically --
  **clang builds it now**, which is what made the sanitizer builds possible.
  The engine swap changes the random *realisation*, not the distribution: over
  2e6 draws the two engines give mean 0.499999/0.500084 and sd
  0.288646/0.288765 against a theoretical 0.5/0.2887.  The `stillinger`
  reference was regenerated accordingly and says so in its config.
- **Zero compiler warnings.**  The remaining 24 (`-Wreorder`, `-Wunused*`,
  `-Wsign-compare`, and 22 `register` keywords in the vendored Mersenne Twister)
  are fixed.  Two of them were real bugs, listed above.
- Added `make asan` and `make asan-check`.
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
