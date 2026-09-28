# macOS build notes

> **Historical document.**  These are the notes from porting the original
> research code to macOS in 2026, and they describe the state *at that time*:
> `-std=gnu++98`, `CC=g++-14`, no CMake, and GSL required for contact detection.
> All of that has since changed -- the code is C++17, builds with GCC or clang
> under both a Makefile and CMake, and the contact test is a quartic root solve
> that does not use GSL at all.  Kept because the *why* behind each porting
> decision is still the useful part, and because some of them (the `push_back`
> in a dependent base, the `size_dist.h` fall-off-the-end bugs) are the kind of
> thing worth recognising again.  See `CHANGELOG.md` for the current state.


Status: **builds and runs** on macOS 14 (arm64) with Homebrew GCC 14.

For replacing the old POV-Ray / raster3d visualisation pipeline, see
[`VISUALIZATION.md`](VISUALIZATION.md) (built-in interactive HTML viewer plus
OVITO integration).

Verified commands:

```sh
make gsl                       # one-off: build the GSL dependency into .deps/gsl
make ellipmd                   # build ./ellipmd
make run CONFIG=config_quick   # 1000-step smoke test, ~1 s
make tools                     # build the post-processing utilities into bin/
```

---

## 1. What the program is

`ellipmd` is a soft-particle molecular-dynamics code for packings of
ellipsoids under gravity (the author's postdoc work, ~2012).

* `main.cc` sets up `CConfig` (key/value parameter file) and a `CSys`
  (`include/mdsys.h`), which owns the box (`BoxContainer`), the particles
  (`CPacking<CParticle>`), the cell list and the integrator. Everything is
  header-only apart from `main.cc`, `CConfig.cc` and `grid.cc`, so the build
  is a single translation unit.
* **Integration** — Beeman predictor/corrector for translation and rotation,
  quaternion orientation update, contact torques, uniform-ellipsoid inertia
  tensor (`CParticle::calPos/calVel`, `include/particle.h`).
* **Force law** — viscoelastic Hertzian normal force
  `Fn = -(k·δ + γ·v_n)·sqrt(δ)·n` plus Coulomb *dynamic* friction
  `Ft = -μ|Fn|·v_t` (no static friction), `Test::contactForce` in
  `include/interaction_force.h`.
* **Ellipsoid–ellipsoid contact** — `doOverlap` in `include/ellips_contact.h`
  builds `M = -(E1⁻¹·E2)` (4×4, homogeneous form), diagonalises it with the
  GSL non-symmetric eigensolver (`include/eigen.h`); the real eigenvector
  gives the pole direction through which the two ellipsoids interpenetrate.
  A ray between poles is intersected with each ellipsoid, then
  `findMin()` refines the pair of surface points with up to 500 iterations of
  a λ fixed-point minimisation of the quadratic potential. That is the
  expensive part of the code, and why it needs GSL.
* **Neighbour search** — `CCellList` (`include/celllist.h`), half-neighbour
  scheme (14 of the 26 neighbours) over a uniform grid of size
  `2·maxRadius`, with image shifts for periodic faces. The optional Verlet
  list (`WITH_VERLET`, `include/verlet.h`) is compiled out.
* **Initialisation** — `method deposition` rains particles in layer by layer
  onto the floor (`add_particle_layer`); `method Stillinger` places them on a
  lattice and compresses by repeatedly scaling all particles
  (`scaling` parameter). `restart` reads a previous snapshot.
* **Boundaries** — `CBox` gives six planes; `periodic_x`, `periodic_xy` and
  `periodic_xyz` mark the corresponding faces non-solid and make the cell
  list use periodic images. `solid` is the plain-wall case.
* **Output** — `outNNNNN` snapshots. `id 6` lines are the wall planes,
  `id 14` lines are ellipsoids: `14 x y z a b c q0 q1 q2 q3`.
  `log_energy` holds `t  E_total  E_kin  E_pot  E_rot`; `out<name>end` is the
  final configuration.

---

## 2. Toolchain requirements

Two things make this a 2012 codebase that a stock macOS toolchain will not
compile:

1. `<tr1/random>` (`include/size_dist.h`, `include/mdsys.h`). TR1 was never
   part of libc++, so **Apple clang cannot build this at all**.
2. Dynamic exception specifications (`throw (matrix_error)` in the borrowed
   `include/matrix.h`) — removed from the language in C++17.

A third problem is toolchain-independent and needed a source fix: two derived
classes called `push_back` on their dependent base without `this->`, which
modern GCC rejects (see §4).

So the build uses **Homebrew GCC 14 (`g++-14`)**, which still ships
`<tr1/random>`, and compiles with **`-std=gnu++98`** (the era-appropriate
standard). `gnu++11` and `gnu++14` also work; `gnu++17` and newer do not,
because of the `throw (matrix_error)` in `include/matrix.h`.

`Makefile.inc` picks `g++-14` automatically when `CC` is not given on the
command line.

If you ever want to build with the system compiler, the required source
changes are: replace `<tr1/random>` with `<random>`, `tr1::` with `std::`,
`ranlux64_base_01` with `ranlux48_base`, and drop the `throw (...)` from
`matrix.h:153`. Note that this changes the random stream, so results will not
be bit-identical to a GCC/TR1 run.

## 3. GSL

GSL is required and was not installed. Rather than touching the system,
`make gsl` downloads GSL 2.8 and builds it statically into `.deps/gsl`
(git-ignored). `Makefile.inc` then prefers `.deps/gsl` and falls back to
Homebrew's `/opt/homebrew` (or `/usr/local`) if you would rather run
`brew install gsl`.

The original `LDFLAGS` pointed at `/usr/lib`, `/sw/lib` (Fink) and
`-I/sw/include`; those paths are gone.

---

## 4. Changes made

### Build system

| File | Change |
| --- | --- |
| `Makefile.inc` | rewritten: GSL autodetection, `-std=gnu++98`, `CC=g++-14`; compile flags and link libraries separated so `-I` and `-l` flags each end up in the conventional place |
| `Makefile` | added `$(LDFLAGS)` to the link line; added a `make gsl` target; added `.PHONY` (the `tools` target was being shadowed by the `tools/` directory, so `make tools` silently did nothing) |
| `tools/Makefile` | `-O2` moved before the source, `$(CC)` instead of hard-coded `g++`, GSL libs added, `correlation`/`coord2xdr` removed from `all`, `clean` extended |
| `.gitignore` | ignore `.deps/`, `ellipmd`, `log`, `log_energy` and the compiled tool binaries |

### Source — required to compile

* `include/packing.h`, `include/celllist.h`: `push_back(x)` →
  `this->push_back(x)`. A derived class calling a member of its dependent
  base is no longer accepted by modern GCC's two-phase lookup.

### Source — required to *run* correctly (crash bugs)

`include/size_dist.h` had nine functions declared to return a value that fell
off the end — `CSizeDistribution::get()` and the `>>`/`<<` operators for
`CSizeDistribution`, `CMonoDist`, `CUniformDist` and `CReadDist`. Modern GCC
emits a trap (`brk`/`SIGTRAP`, exit 133) on that path, so **any configuration
containing a `SizeDistribution` line died immediately** — that includes the
shipped `config_deposition` and `config_stillinger`.

* `CSizeDistribution::get()` now returns `p_dist->get()` (it was discarding
  the value; the result was used as a particle radius).
* `operator>>`/`operator<<` for `CSizeDistribution`, `CMonoDist`,
  `CUniformDist`, `CReadDist` now return their stream.
* `CReadDist` stored the values read from file in a `vector<int>`; changed to
  `vector<double>` (radii were being truncated to 0/1).

---

## 5. Bugs found but deliberately *not* fixed

* **`CSizeDistribution` has no deep copy and an inverted destructor**:
  `~CSizeDistribution(){ if(!p_dist) delete p_dist; }`. It leaks, but simply
  flipping the condition causes a double free, because `CBaseConfig::add_param`
  passes the distribution by value at least twice and the default copy
  constructor shares the `p_dist` pointer. A correct fix needs a `clone()`
  on `CBaseDistribution`.
* **The defaults are unusable.** `main.cc` falls back to `config` when no
  file is given; the README says it will "run with default parameters", but
  `particleSize=1` inside the default `boxsize 1 1 2` puts every particle
  outside the grid and the run aborts with
  `Point out of grid: (x,y,z): ... 3.22114`. Either ship a `config` file or
  pass `CONFIG=`.
* **`config_quick` is new** — a small config I added so `make run` has
  something valid to run.
* **Misleading message**: `CSys::solve()` prints
  `"Relaxation criterion reached at time=... KE= <x> < 1e-8"` whenever it
  stops, including when it simply reached `maxTime` with `KE` three orders of
  magnitude larger. The message should be printed only for the real
  convergence branch.
* **First line of `log_energy` is garbage** (`~1e-314`). `output()` is called
  before the energies are computed on the first pass, so uninitialised
  doubles get written.
* **`CSys::initialize` reads `nParticle` as `int`** while it is registered as
  `unsigned int`; `CBaseConfig::get_param<T>` does an unchecked
  `static_cast<CParam<T>*>`, so a mismatched `T` is undefined behaviour that
  happens to work here.
* **`zeta` and `zetaWidth` are ignored for `particleType general`.** In
  `CSys::initialize` (include/mdsys.h) the code reads `zeta0`/`zetaW` and then
  never uses them:

  ```cpp
  double zeta0=config.get_param<double>("zeta");
  double zetaW=config.get_param<double>("zetaWidth");
  double eta0 =config.get_param<double>("eta");
  double etaW =config.get_param<double>("etaWidth");
  ...
  double zeta=eta0*TruncGaussRand(1, etaW);   // should be zeta0, zetaW
  double eta =eta0*TruncGaussRand(1, etaW);
  ```

  So the two shape parameters are not independent: `zeta` is forced to follow
  `eta`, and `zetaWidth` does nothing. `config_quick` and `config_demo` set
  `zeta == eta` anyway, so they are unaffected.
* **`TNode::normal_fabric_tensor()` is missing its cache check** and returns
  `branch_fabric_M` from its (dead) early return — a copy/paste slip. Harmless
  as written, wrong if the cache is ever made real.
* **`static` scratch buffers throughout** (`Test::interact`,
  `CSys::interact`, `eigens`, the `static double data[16]` in `include/eigen.h`)
  make the code non-reentrant and unsafe to thread.
* **`config_base` has `read_radii 1`**, but `particleType` defaults to
  `general`, which ignores it; and `particleShape` / `initialization` in
  several shipped configs are not registered parameters and only produce
  `Warning: ... is not a valid parameter or keyword`.
* **Performance**: with the Verlet list compiled out, the cost is ~1 ms/step
  for a 100–200 particle system while the particles are still falling, rising
  to ~11 ms/step once the packing has settled and contacts are numerous.
  A production run of 10⁴ particles is not something you will want to do on a
  laptop.

---

## 6. Smoke test results

`config_quick` (100 prolate ellipsoids, `η=ζ=1.4`, solid box 1×1×1.2,
gravity −z, `dt=1e-4`, `t_max=0.1`, seed 51):

```
Constructing the grid ... done: 7 X 7 X 8
Number of Particles: 0
Relaxation criterion reached at time=0.1: KE= 0.0162507 < 1e-8
```

`log_energy` (`t  E_total  E_kin  E_pot  E_rot`) shows the particles falling
and settling; note the garbage first line described above:

```
0                    7.3e-314  2.2e-314  3.0e-314  2.2e-314   <- uninitialised
0.02   0.085131  0.0038815  0.080683  0.00056612
0.04   0.084692  0.0067382  0.077339  0.00061507
0.06   0.083518  0.0102595  0.072354  0.00090450
0.08   0.081115  0.0129265  0.066514  0.00167482
```

The `Stillinger` + `periodic_xyz` path was checked separately: it placed 149
of 200 requested particles (the lattice loop stops at the top face, which is
the intended behaviour) and ran to completion.

## 7. Tools

`make tools` builds `coord_convert`, `packing_density`, `fabric`, `compact`
and `periodic2full` into `bin/`. `correlation` and `coord2xdr` need HDF5 and
Sun RPC/XDR (`hdf5.h`, `rpc/xdr.h`), which a stock macOS toolchain does not
have; they are kept as targets but left out of `all`.

`bin/genFrames.sh` and `make pov` additionally need `raster3d` (`render`) and
`povray`, and `bin/coord2pr3d` is a Linux-era script; none of that is needed
to run the simulation itself.
