# Roadmap: from research code to a shareable project

Goal: a state that is correct, documented, buildable by someone else in one
command, and fast enough to be useful — without silently changing the physics.

Everything below is grounded in the code as it stands. Numbers come from
`docs/BENCHMARK.md` / `bench/results.json` and from a `sample` profile of run B2
(2500 particles).

---

## Where we are

**Working:** builds and runs on macOS with GCC 14 + GSL; deposition and
Stillinger initialisation; solid and periodic boundaries; snapshots, energy log,
fabric tensors; modern interactive visualisation (HTML viewer, OVITO, LAMMPS
dump).

**Measured baseline** (Apple M1 Max, single-threaded, `bench/results.json`):

| run | configuration | wall | ms/step |
|---|---|---|---|
| B1 | 250 p, dt=1e-4 | 5.43 s | 1.975 |
| B2 | 2500 p, dt=1e-4 | 78.58 s | 28.573 |
| B3 | 250 p, dt=1e-5 | 58.52 s | 2.128 |
| B4 | 2500 p, dt=1e-5 | 782.82 s | 28.466 |

**Profile of B2** (self time, leaf frames):

| % | what |
|---|---|
| **49.8%** | `math::matrix<double>` — the vendored linear algebra class |
| **18.0%** | GSL, essentially all of it the 4×4 non-symmetric eigensolver |
| 10.8% | ellipsoid contact geometry (`intersect`, `toBody`, `gradient`, `findMin`) |
| 4.3% | allocator (`operator new[]` / `delete[]`) |
| 16.8% | everything else |

Drilling into the matrix block: destructor 16.4%, `operator*=` 8.6%, `clone()`
5.1%, constructor 4.7%, `Inv()` 4.3%. That is not a mysterious hotspot — it is
heap traffic. `matrix<double>` allocates a row-pointer array plus one array per
row (five allocations for a 4×4), and `CEllipsoid` carries **six** 4×4 matrices
per particle, three of which are rebuilt from scratch every timestep in
`update_tranlation_mat()`.

**The thing that makes all of this affordable to touch:** the code is currently
*unverified*. There is no test that says "the physics is still the same". So
Phase 0 comes first — without it, every later change is a gamble.

---

## Principles

1. **Physics before speed.** A regression harness lands before any refactor.
2. **Every change measured.** `bench/results.json` is the record; a change that
   is not measured is not an optimisation.
3. **Delete before refactor.** Dead code has no tests and no users; removing it
   is the cheapest possible improvement.
4. **Upstream first.** A modern toolchain solves problems we should not be
   solving by hand (the matrix class, the eigensolver, TR1).
5. **Presentable means a stranger can build it.** `git clone && make && make
   check` on a clean machine, with a README that says why the project exists.

---

## Phase 0 — Safety net  *(must land first)*  — **DONE**

- [x] `bench/check_physics.py`: run a fixed config with a fixed seed, compare
      every snapshot against a committed reference within a tolerance
      (max |Δposition|, |Δquaternion|, |Δsemi-axis|).
- [x] Two reference cases, to cover both initialisation paths:
      - deposition + solid walls (`config_quick`)
      - Stillinger + `periodic_xyz`
- [x] `make check` target; non-zero exit on failure.
- [x] Record why a tolerance rather than a hash: recompilation with different
      flags may perturb the last bits, and a bit-exact test would block
      legitimate optimisation.

**Acceptance:** `make check` passes on the current code, and fails if the
physics is perturbed (verify by nudging a constant and watching it fail).

---

## Phase 1 — Delete the dead weight  — **mostly DONE**

Nothing here changes behaviour.

**Vestigial translation units** (not even compiled — the Makefile builds only
`main.cc`):
- [x] `CConfig.cc` — includes only.
- [x] `grid.cc` — a single `#include`.

**Tools that cannot build** (already excluded from `tools/Makefile`):
- [x] `tools/asphericity.cc` — includes a nonexistent `<CStat.h>` and a
      pre-`include/` path.
- [x] `tools/sphere_map.cc` — includes `include/define_params.h`, renamed to
      `config.h` years ago.
- [x] `tools/correlation.cc`, `tools/coord2xdr.cc` — need HDF5 and Sun RPC/XDR.
- [x] `bin/coord2xdr.sh` — driver for the above.

**Superseded pipeline:**
- [x] `bin/coord2pov`, `make pov` — the POV-Ray path, replaced by
      `docs/VISUALIZATION.md`.
- [x] `bin/coord2pr3d`, `bin/genFrames.sh`, `bin/encodejpg.sh`,
      `make animate` / `movie` — the raster3d path.

**Cluster scripts with the author's absolute paths baked in** (`/home/reza/...`),
useless to anyone else:
- [x] `bin/jobs`, `bin/jobs_abc`, `bin/jobs_gen_asp`, `bin/jobs_relax`,
      `bin/lastarg.sh`.
- [x] `bin/avgdensity.sh` — shells out to an `avg` binary that is not in the
      repo and assumes a directory layout that no longer exists.

**Half-finished features** that are referenced but never instantiated:
- [x] `include/cylinder.h` (`CCylinder`) — self-described as "not complete".
- [x] `include/composite.h` (`CComposite`) — "maybe not fully working"; needs
      the `interaction.h` overloads and `shapes.h` include removed with it.
- [x] `include/verlet.h` — the Verlet list is `#undef`'d out and has therefore
      never run in any build we have. Decide: repair and enable, or delete.
      Deleting is defensible until Phase 6 makes the cell list fast enough that
      it would only be needed for very large N.

**To keep, but reclassify:**
- [ ] `tools/generate_aspects.nb`, `tools/map_asph_aspect.nb` — the Mathematica
      provenance for `include/map_asph_aspect.h`. Move to `tools/notebooks/`
      and say so in the README.
- [ ] `bin/density.sh`, `bin/multidensity.sh` — still meaningful drivers for
      `packing_density`, but have hardcoded paths; rewrite or fold into an
      `analysis/` tool.
- [ ] `include/grid.h`, `include/slice.h` — only used by the analysis tools;
      move under a clearly tool-facing umbrella.

**Acceptance:** `make`, `make tools`, `make check` all still pass; the diff
contains no functional change; the file count drops materially.

---

## Phase 2 — Documentation a stranger can use  — **DONE**

- [x] **`README.md`** (replace the plain-text `README`). Sections: what the
      problem is, one-paragraph physics summary, build in one command,
      dependencies, quick start, where things live, how to cite, license.
- [x] **`docs/ARCHITECTURE.md`** — the object graph (`CSys` → `BoxContainer` /
      `CPacking<CParticle>` / `CCellList`), one translation unit and why,
      where each physics concept lives, the data flow of a timestep.
- [x] **`docs/FORMAT.md`** — the snapshot format as a specification: `id 6`
      planes, `id 14` ellipsoids, quaternion convention (scalar-first — this
      genuinely bites downstream), what `log_energy` columns mean, and the
      "first line of `log_energy` is uninitialised" caveat until it is fixed.
- [x] **`docs/PHYSICS.md`** — the force law, the contact algorithm, the
      integrator, and the assumptions (no static friction, no cohesion
      actually wired up, `zetaWidth` ignored for `particleType general`).
      Reference the original paper.
- [x] Fold `PORTING-NOTES.md` / `BENCHMARK.md` / `VISUALIZATION.md` into
      `docs/` (keeping the filenames) so the root stays clean.
- [x] A short `CHANGELOG.md`, since the roadmap is about tracking improvements.

**Acceptance:** someone who has never seen the repo can build it and produce a
picture from the README alone.

---

## Phase 3 — Modernise the build and toolchain  — **mostly DONE**

- [ ] **CMake** alongside (or replacing) the hand-written Makefiles:
      `cmake -B build && cmake --build build`. This is what makes the project
      installable and IDE-friendly, and it makes the GSL/Eigen dependency
      explicit via `find_package`.
- [x] **C++17** (drop `-std=gnu++98`).  The blockers were the vendored
      `matrix.h` exception specification and `<tr1/random>`; both are gone and
      **clang now builds the code**, which is what unlocked the sanitizers.
- [x] Replace `<tr1/random>` with `<random>`; `ranlux64_base_01` →
      `ranlux48_base` (note: changes the RNG stream, so it must be done with
      Phase 0 in place, and the reference outputs regenerated deliberately).
- [x] Compiler flags: `-O3 -march=native` for release, `-Wall -Wextra
      -Wpedantic` for development, and fix the 24 existing warnings
      (`-Wreorder`, `-Wunused*`, `-Wuninitialized`, `-Wdelete-non-virtual-dtor`).
- [x] AddressSanitizer / UBSan build target (`make asan`, `make asan-check`).
      It found a heap-buffer-overflow and an invalid downcast on its first run
      -- see CHANGELOG.  Replacing the remaining `static` scratch buffers is
      still Phase 5, and that is what stands between this and a thread-safe
      build.
- [ ] **CI** (GitHub Actions): build + `make check` on Linux and macOS. This is
      most of what "presentable" means in practice.

**Acceptance:** a clean checkout builds on a machine with no prior setup beyond
a compiler, CMake and GSL/Eigen.

*Verified so far, without CMake:* `git archive HEAD` into an empty directory,
then `make ellipmd`, `make tools`, `make check` and the Python tooling all
succeed with only a compiler and GSL.  Worth wiring into CI, which is the point
of this phase.

---

## Phase 4 — Bug fixes (correctness)  — **in progress**

From `docs/PORTING-NOTES.md`, plus what the port surfaced. Ordered by how much they
can silently mislead a result.

- [x] **The translational velocity corrector was wrong.**  The `a_{n+1}` term
      had the wrong sign and the `a_{n-1}` term was missing, so translation was
      first-order where rotation was second-order: 5% energy loss per elastic
      bounce at dt=1e-4, converging only linearly.  Found by building an
      energy-conservation test, and now guarded by
      `bench/reference/elastic_bounce`.
- [x] **`zeta` / `zetaWidth` are ignored for `particleType general`.** Reads
      `zeta0`/`zetaW`, then computes `zeta = eta0 * TruncGaussRand(1, etaW)`.
      The two shape parameters are therefore not independent and `zetaWidth`
      does nothing. Any published result using unequal `eta`/`zeta` is suspect.
- [x] **`CBaseConfig::get_param<T>` did an unchecked `static_cast`.** Asking
      for `int` where the parameter was registered as `unsigned int` was UB
      that happened to work.  Now a checked `dynamic_cast`; UBSan flagged it.
- [x] **`DisBetaDistribution` wrote one `double` past the end of its `bins`
      array** on every construction, and leaked it.  Found by ASan.
- [x] **`CSizeDistribution` ownership** -- inverted destructor, no deep copy,
      and a base class with no virtual destructor.  Fixed as one problem.
- [x] **`CUniformDist` was built from uninitialised members.**
- [x] **`CSys::~CSys` used TRY/CATCH**, which rethrows; a destructor is
      implicitly `noexcept`, so that would have called `std::terminate`.
- [x] **First line of `log_energy` was uninitialised garbage** — energies are
      accumulated at the end of `forward()`, so on the first snapshot they had
      never been computed.  Added `CSys::computeEnergies()`.
- [x] **The "Relaxation criterion reached ... KE < 1e-8" message was printed on
      every exit**, including a normal `maxTime` exit.
- [x] **Deposition aborted with "Point out of grid" once the box filled up** —
      a new layer was seeded unconditionally at a z that could exceed the lid.
      Only went unnoticed because `nParticle` was usually exhausted first.
- [x] **`method restart` was broken** — `celllist.build()` ran before
      `celllist.setup()`, so `which()` divided by an uninitialised `dx`.
- [x] **`nRadii` was a dead guard** — a function-local static nothing ever
      incremented.
- [ ] **Default parameters are unusable** — `particleSize=1` in a `1×1×2` box
      puts every particle outside the grid. Ship a `config` or fail loudly.
- [ ] **The deposition gate uses a literal `1.` as the box height**
      (`maxh < 1 + 2*maxRadii`), a leftover from a 1×1×1 box.  It should be
      `walls.L(2)`, but changing it changes how many particles get placed, so it
      needs a deliberate decision rather than a drive-by fix.
- [ ] **`README` documents the shape parameters wrongly.** It says
      `eta = a/b, xi = b/c`; the code computes `a/b = zeta` and `c/a = eta`.
      The config-file names do not mean what the README says.
- [ ] **The snapshot format carries no velocities**, so `method restart` starts
      from rest.  Needs a format change (extra columns) or a separate state file.
- [ ] `CParticle::addforce` overwrites `avgforces` with the *last* contact
      force rather than accumulating, and uses a `static prev` that is shared
      across all particles — so "average force" is neither averaged nor
      per-particle.  Its only consumer is an `if(0)` block in `forward()`, so
      the honest fix is probably deletion.
- [ ] `TNode::normal_fabric_tensor()` returns `branch_fabric_M` from its dead
      early return, and `bool calculated=false` is a local so the cache never
      works.
- [ ] `CPolynom<order,T>::operator()` uses function-local `static` accumulators
      — not reentrant, and wrong if anything nests.
- [ ] `CException` is thrown and caught **by value** everywhere, slicing the
      type and making the `catch(...)` fallbacks load-bearing.

**Acceptance:** each fix is a self-contained commit with a note in the
CHANGELOG, and `make check` still passes (except where the fix deliberately
changes a reference, which must be called out).

---

## Phase 5 — Design and robustness

The goal is fewer concepts, and failures that are loud and early.

- [ ] **Split the single translation unit.** Everything being header-only made
      sense for a prototype; for sharing it means every edit recompiles the
      world and every symbol is implicitly inline. Split into
      `src/` + `include/` with real `.cpp` files.
- [ ] **Replace `include/matrix.h`** (see Phase 6) — it is 1153 lines of
      borrowed code with copy-on-write reference counting, used only for
      fixed-size 3×3 and 4×4 matrices.
- [ ] **Rename `include/eigen.h`** — it is a GSL eigensolver wrapper, and
      collides conceptually with the Eigen library.
- [ ] **Kill the global mutable state**: globals `config`, `rgen`, `eng`,
      `particle_material`, `G`, `friction`, and the `static` scratch buffers in
      `Test::interact`, `CSys::interact`, `eigens`, `polynom` and `CParticle`.
      None of it is thread-safe, some of it is shared across particles in ways
      that are outright wrong.
- [ ] **Configuration**: a typed schema with validation and a `--print-config`
      that echoes the effective parameters into the run directory, so a result
      can be traced to the exact inputs that produced it.
- [ ] **Errors**: throw by reference, one exception hierarchy, and stop using
      exceptions for control flow in the inner loop (the `TRY`/`CATCH` macros
      wrap almost every function).
- [ ] **CLI**: replace the positional `ellipmd <seed> <config>` with proper
      argument parsing, plus `--help` and `--version`.
- [ ] **Output**: configurable precision, and a self-describing header in each
      snapshot so downstream tools do not have to re-derive the format.
- [ ] **Restart that actually works** — currently `method restart` reloads
      positions and orientations but not velocities.

**Acceptance:** no globals in the physics path; a run is reproducible from
`config` + seed alone; ASan/UBSan clean on a short run.

---

## Phase 6 — Performance  — **6a DONE, 6b open**

Ordered by (expected win) / (risk). Each step re-runs `bench/run_bench.py` and
records the delta.  See `docs/BENCHMARK.md` Part 2 for the full log.

**Phase 6a, done:** `-O3` (-14%), deleting work whose result was discarded in
`doOverlap` (-19%), one pose update per step instead of two (-24%), and
allocation-free matrix arithmetic (-34%).  Full re-baseline: B1 -32.0%, B2
-29.0%, B3 -35.6%, B4 -28.8%, with `make check` unchanged throughout.

Two things measured and rejected: `-march=native` (no effect at all) and a
conservative inscribed-sphere rejection test (correct, but it fired on 0.09% of
candidate pairs).  Both are written up in `docs/BENCHMARK.md`.

`bench/run_bench.py` gained `--repeat N` after discovering that a single run
varies by up to ~20%.

### Low-hanging fruit — no algorithmic change

- [x] **`-O3`** instead of `-O2` (-14% on B2).  `-march=native` measured at
      zero effect; it stays opt-in behind `NATIVE=1`.
- [x] **Hoist work out of the timestep.** `CSys::forward` calls
      `config.get_param<string>("method")` and `add_particle_layer` on *every*
      step even after all particles are placed; `config.get_param` is a `map`
      lookup by string plus a `static_cast`.
- [ ] **Stop rebuilding the cell list from scratch every step.**
      `celllist.build()` clears and refills every one of `nx*ny*nz` cells
      (4050 of them for the fine run, against 2500 particles) and
      `interact()` walks all of them. Track occupied cells, or move particles
      incrementally (`CCellList::update` already exists and is unused).
- [ ] **`CPacking` is a `std::list<CParticle*>`** — pointer chasing in the
      hottest traversal. `std::vector<CParticle*>` (or `vector<CParticle>`)
      would be markedly more cache-friendly; the list is only used for stable
      iterators during removal.
- [ ] **Re-measure the cell-list rebuild.**  It clears and refills every cell
      each step, but the profile never showed it above noise, so it is not a
      priority despite looking like one.
- [ ] **Avoid `findMin`'s 500 iterations when it has already converged** — it
      checks convergence, but the loop bound is unconditional in the callers.

### The main event — the matrix class

This is ~50% of the runtime and it is a data-structure problem, not a maths
problem.

- [ ] **Replace `matrix<double>` with fixed-size stack-allocated types** for
      3×3 and 4×4 (or use **Eigen**, which is already installed on this
      machine, `-I/opt/homebrew/include/eigen3`). `CEllipsoid` holds six 4×4
      matrices; each one heap-allocates five blocks today.
- [ ] **Stop recomputing `ellip_mat` in `update_tranlation_mat()` every step.**
      It is `(~tempmat)*scale_mat*tempmat` — three 4×4 multiplies plus
      temporaries and allocations per particle per step, to produce a matrix
      whose only needed content is a symmetric 3×3 block and a centre vector.
      Cache it against the last `(q, Xc)`, or reformulate the hot predicates
      (`operator()`, `gradient`, `toBody`) to work from the quaternion and
      centre directly — the quadratic form is
      `(R(x-Xc))ᵀ S (R(x-Xc))` and never needs a 4×4.
- [ ] Make `CEllipsoid`'s matrices `const`-after-setup so the compiler can see
      what actually changes per step.

### The eigensolver — 18%

- [ ] **Add cheap rejection tests before `doOverlap`.** Today every candidate
      pair pays `Inv()` + two 4×4 multiplies + a full GSL non-symmetric
      eigensolve, before any bounding-sphere or separating-plane test. Reject
      on the inscribed-sphere bound (`|ΔXc| > R₁+R₂`) first, then on a
      circumscribed-sphere test; only then do the eigen work. Contact detection
      is exactly the place where a cheap conservative test pays for itself.
- [ ] For the 4×4 problem the code only reads the real parts of the
      eigenvalues and the real eigenvector. The characteristic quartic is
      already available (`characteristicPolynom`, `include/polynom.h`) and
      `CQuartic::solve` exists — using it would remove GSL from the hot path
      entirely, and possibly from the dependency list.
- [ ] Reuse workspace/scratch instead of allocating the GSL workspace and
      complex vectors per call.

### Then, if still needed

- [ ] Verlet list on top of the (now cheap) cell list, with a skin distance —
      this is what makes 10⁴–10⁵ particles tractable.
- [ ] OpenMP over particles for the force loop (needs Phase 5's static-state
      removal first — this is the reason that item is on the list).
- [ ] Optional single-precision contact geometry, if the stability analysis
      supports it.

**Acceptance:** each step is a commit with a `bench/results.json` delta in the
message, and `make check` passes throughout. Target: at least 3× on B2 without
changing the physics, which the profile suggests is achievable from the matrix
work alone.

---

## Suggested order

```
Phase 0  safety net                    <-- first, non-negotiable
Phase 1  delete dead code              <-- cheap, shrinks everything after
Phase 2  documentation                 <-- while the design is fresh
Phase 4  bug fixes                     <-- correctness before speed
Phase 3  build modernisation           <-- unlocks C++17 for Phase 5/6
Phase 6a low-hanging performance
Phase 5  design refactor
Phase 6b matrix + eigensolver          <-- the big win
Phase 6c Verlet / OpenMP               <-- only if needed
```

Phases 3 and 4 can swap if the matrix replacement is done first — that is the
single change that unlocks C++17 *and* removes half the runtime, so it is
tempting to pull forward. The reason it sits late is risk: it rewrites the
geometry every contact calculation depends on, and it should not be attempted
before `make check` exists and the cheap wins are banked.
