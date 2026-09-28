# Architecture

How the solver is put together, and what happens during one timestep.

## One translation unit

`main.cc` includes `include/main.h`, which pulls in `mdsys.h`, and everything
else is header-only.  The Makefile compiles exactly one file:

```sh
$(CC) main.cc $(FLAGS) -DDEBUGFLAGS -o ellipmd $(LDFLAGS)
```

This is a deliberate prototype structure: no build graph, no link step to get
wrong, and the compiler can inline across the whole program.  The costs are that
every edit recompiles everything, every function is implicitly `inline`, and
there is no natural place to hide an implementation.  Splitting it into
`src/` + `include/` is on [`ROADMAP.md`](../ROADMAP.md) Phase 5.

## Ownership

```
CSys                                    the simulation
├── BoxContainer walls                  the box; six CPlane faces
├── CPacking<CParticle> particles       the particle container
├── CCellList<Packing, CParticle>       uniform grid, half-neighbour scheme
├── CSizeDistribution size_dist         radii source (Stillinger only)
├── vec G                               gravity
└── double t, dt, DT, tMax, outDt, ...  time bookkeeping
```

`CPacking<T>` is a `std::list<T*>` plus a contact network and a packing-fraction
helper.  `CCellList` owns a flat array of `CCell` (each a `std::list<TParticle*>`)
and precomputes which of the 13 "forward" neighbour cells each cell should check,
so no pair is visited twice.

Shapes form their own small hierarchy, rooted at `GeomObjectBase`:

```
GeomObjectBase
├── CSphere      radius
├── CPlane       point + normal; used for the six box faces
├── CBox         six CPlane faces; knows which are solid
└── CEllipsoid   semi-axes, quaternion, and six 4x4 matrices
```

`CParticle` owns a `GeomObjectBase*` (so a particle is a shape plus state) and
holds the derivative arrays `x` and `w` — `CDFreedom<3>` instances, where index
0 is position, 1 velocity, 2 acceleration.  `CParticle::vlist` is the
**contact cache** (`CVerletList`, despite the name): a map from a neighbouring
particle to the contact geometry computed for that pair this step.

### Global mutable state

These are global and are on the cleanup list:

| symbol | where | what |
|---|---|---|
| `config` | `main.cc` | the parsed configuration, read from everywhere |
| `rgen` | `common.h` | Mersenne Twister, the main RNG |
| `eng`, `eng0` | `mdsys.h`, `size_dist.h` | TR1 ranlux engines for the size distributions |
| `particle_material` | `mdsys.h` | the material applied to every particle |
| `G` | `common.h` | a stray global gravity vector, shadowed by `CSys::G` |
| `friction` | `particle.h` | a stray global, unrelated to the material's friction |
| `gout` | `common.h` | an unused `ofstream*` |

Several functions also use function-local `static` scratch buffers, which makes
them non-reentrant and is why the force loop cannot be parallelised as it
stands.  See `ROADMAP.md` Phase 5.

## A timestep

`CSys::solve()` runs `forward(dt)` until `t` reaches `maxTime` (clamping the
last step so it lands exactly on `tMax`), then writes `outend`.

`CSys::forward(dt)` does, in order:

1. **Retire expired particles** — anything with `expired` set is deleted from
   the container.
2. **Deposit a layer** if `method == deposition` and there is still room; a new
   layer is seeded at `maxh + 1.02·maxRadii` (see the caveat in
   `docs/FORMAT.md` and `ROADMAP.md`).
3. **Write output** if this step is due: a snapshot `outNNNNN`, a line in
   `log_energy`, and any per-output housekeeping (shrinking gravity, applying
   `scaling`, `spherize_on`).
4. **`calPos(dt)`** for every particle: Beeman position update, rotational
   update, quaternion advance and renormalisation, and the shape's matrices
   refreshed via `rotateTo`/`moveto`.
5. **`calForces()`**:
   * reset each particle's force to `m·g − fluid_damping·|g|·m·v`, and its
     torque to zero;
   * wall contacts: `interact(p, &walls)`, which tests the particle against each
     non-solid-skipped face;
   * `celllist.build(particles)` then `celllist.interact()`, which for each cell
     tests the pairs inside it and the pairs against its 13 forward neighbours,
     shifting a particle's position by the image vector when the neighbour
     crosses a periodic boundary and shifting it back afterwards.
6. **`calVel(dt)`** for every particle: Beeman velocity correction and the
   angular-velocity update from the torque and the inertia tensor, accumulating
   `rEnergy`, `pEnergy` and `kEnergy` on the way.

Note that a particle only ever *adds* force; the equal-and-opposite force on the
neighbour is applied inside `Test::interact`, which is why pair contacts are
handled in the cell list and wall contacts separately in `CSys::interact`.

## Contact detection

`CInteraction::overlaps` is a double-dispatch typed on `GeomObjectBase::type`,
with overloads for sphere/sphere, sphere/box, ellipsoid/plane, ellipsoid/box and
ellipsoid/ellipsoid.  The last is the interesting one and is described in
[`PHYSICS.md`](PHYSICS.md).

## Where the time goes

A `sample` profile of the 2500-particle benchmark run attributes roughly:

| % | |
|---|---|
| 68 | `math::matrix<double>`, the vendored 4×4 matrix class — `operator()`, `~matrix()`, `clone()` and `invert_into`, i.e. the *type* rather than the arithmetic |
| 6 | `std::vector` reallocation, in the polynomial solvers |
| 4 | the polynomial solvers themselves |
| 3 | ellipsoid contact geometry |

The eigensolver is no longer in this table: the contact test is a quartic root
solve (see [`PHYSICS.md`](PHYSICS.md)), which took GSL off the hot path.  What
dominates now is that `matrix<double>::operator()` bounds-checks and
reference-count-checks *every element access*, and there are sixteen of them per
4×4 multiply.  Replacing the type with a fixed-size stack array is the largest
remaining optimisation and is deliberately left as future work --.  See
[`BENCHMARK.md`](BENCHMARK.md) for the log and `ROADMAP.md` Phase 5.

`CEllipsoid` holds six `matrix<double>` members and rebuilds three of them from
the quaternion and position on every step in `update_tranlation_mat()`.  Each
`matrix<double>` heap-allocates a row-pointer array plus one array per row — five
allocations for a 4×4.  Replacing this with fixed-size stack types (or Eigen) is
the single largest win available; see `ROADMAP.md` Phase 6.

## Things that look like they do something and do not

* `CMaterial::color` is written (for frozen particles) and never read.
* `CMaterial::cohesion`, `static_friction` and `friction_threshold` are set from
  the config and consulted by no force law.  The first two now warn.
* `CParticle::avgforces` is set to the *last* contact force from a `static`
  scratch variable shared across all particles, and its only reader is an
  `if(0)` block.
* `CSys::interactions()` is a debug helper that is never called.
