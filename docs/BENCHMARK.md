# Timing benchmark

Two things live here: the scaling study that shows how cost moves with particle
count and time step, and the log of what the optimisation work actually bought.
Numbers are from `bench/results.json`, which is committed.

**A single run varies by up to ~20% on a laptop**, which is more than most
optimisations are worth.  Use `--repeat N` (the median is recorded, along with
every sample) before believing any delta below that.

`--fast` skips the longest configuration (B4, about 9 of the ~16 minutes) for
use while iterating, and merges rather than overwrites, so the skipped entry
keeps its last full-sweep number.  Run a full sweep before quoting the file as
a whole — the `repeats` field on each entry says whether it was re-measured.

---

# Part 1 — scaling: 10x particles x 1/10 time step

Machine: Apple M1 Max, 10 cores, 32 GB, macOS 14.4.  Single-threaded run, no
other load (the four runs were executed strictly one after another — running
them concurrently would have made the numbers meaningless).

## What changed

| | baseline | new |
|---|---|---|
| particles | 250 | **2500** |
| `particleSize` | 0.05 | **0.0232079** |
| `timeStep` | 1e-4 | **1e-5** |
| `maxTime` | 0.275 | 0.275 (unchanged) |
| box | 1 x 1 x 1.2 | unchanged |
| everything else | `eta=zeta=1.4` with 25% polydispersity, `stiffness=5e2`, `damping=5`, `friction=0.2`, `gravity=-10 z` | unchanged |

`particleSize` is scaled by `10^(-1/3) = 0.46416`, not by 1/10, so the **volume
fraction is exactly preserved**:

```
250  particles at r = 0.050000  ->  total volume 0.130900,  phi = 0.10908
2500 particles at r = 0.023208  ->  total volume 0.130900,  phi = 0.10908
```

`maxTime` was deliberately left alone, so B4 is the same physical interval as B1
with 10x the particles and 10x the steps — a clean 100x-work comparison.

Configs: `config_large` (this is B4) and `bench/configs/{B1..B4}`.

## Reproducing / re-measuring

The driver and the exact configurations are in `bench/`, and the numbers below
are committed as `bench/results.json` so a later optimisation can be diffed
against this baseline:

```sh
make ellipmd
python3 bench/run_bench.py            # all four, sequential; ~16 min
python3 bench/run_bench.py --only B2  # ~80 s, the cheapest useful signal
git diff bench/results.json
```

Runs are sequential on purpose: overlapping them invalidates the timings.
Per-run `out*` directories and `stdout.log` land in `bench/runs/` (git-ignored).

## Results

| run | configuration | steps | wall time | ms/step | program CPU | particles in final snapshot |
|---|---|---|---|---|---|---|
| B1 | 250 particles, dt=1e-4 | 2 750 | **5.43 s** | 1.975 | 5.30 s | 250 |
| B2 | 2500 particles, dt=1e-4 | 2 750 | **78.58 s** | 28.573 | 76.95 s | 2500 |
| B3 | 250 particles, dt=1e-5 | 27 500 | **58.52 s** | 2.128 | 54.82 s | 250 |
| B4 | 2500 particles, dt=1e-5 | 27 500 | **782.82 s (13 min 3 s)** | 28.466 | 744.30 s | 2500 |

All four exited cleanly (`rc=0`), reached `t = 0.275`, and produced the requested
number of particles.

## Reading the numbers

**Particle count is slightly super-linear.**  B2/B1 = 78.58 / 5.43 = **14.5x**
for 10x the particles.  Per step, 28.573 / 1.975 = 14.5x.  The cell list makes
the neighbour search O(N), but the ellipsoid pair test (`doOverlap` plus a
500-iteration `findMin`) scales with the number of *contacts*, and a denser
packing has more of them, so the constant grows a little faster than N.

**Time step is exactly linear.**  B3/B1 = 58.52 / 5.43 = **10.8x** for 10x the
steps, and the per-step cost is essentially unchanged (2.128 vs 1.975 ms/step,
+7.7%).  Halving dt does not make a step cheaper or dearer — it just buys more
of them.  B4 confirms it at the larger size: 28.466 vs 28.573 ms/step.

**Combined, B4 costs 144x B1** (782.82 / 5.43), against the 100x you would get
if both scalings were linear.  The extra 1.44x is the super-linear particle
term.

**Where the time goes.**  The per-step cost is not constant over a run: it
roughly doubles as the packing densifies.  B4 was running at ~13 ms/step up to
`t = 0.13` and averaged 28.5 ms/step over the whole interval; the settled part of
B2 ran at ~42 ms/step.  So the numbers above are averages over a transient, and
a run that reaches a settled packing will be slower per step than these.

## What the run physically covers

`maxTime = 0.275` is short: with gravity -10, a particle falling from the top of
the box (`z = 1.2`) needs `sqrt(2*1.2/10) = 0.49` to reach the floor, so nothing
has had time to fully settle.  In the final state:

```
z 0.0-0.1 :  965   <- deposited layer
z 0.1-0.2 :  261
z 0.2-0.3 :  109   <- still falling
...
z 1.0-1.1 :  134
zmax = 1.08,  <z> = 0.356
```

965 of 2500 particles have reached the floor and the rest are still in flight.
Raising `maxTime` is where the cost really bites:

| target | steps at dt=1e-5 | projected wall time |
|---|---|---|
| t = 0.275 (this run) | 27 500 | 13 min |
| t = 1.0 | 100 000 | ~47 min |
| t = 1.5 (settled) | 150 000 | ~71 min |

At the original `dt = 1e-4` the same settled state is ~7 min, at 10x less
integration accuracy.  Note also that the code has no true restart: `method
restart` reloads positions and orientations but not velocities, so a long run has
to be done in one go.

## Visualisation of this run

Built from the B4 snapshots, all with the random per-particle colouring:

| file | what |
|---|---|
| `viz/trajectory_large.html` | 12 frames, interactive, 2500 particles (4.7 MB, self-contained) |
| `viz/ovito_movie.mp4` | OVITO, 1280x960, 12 frames, 2 fps, dark outlines (1.7 MB) |
| `viz/ovito_final_state.png` | OVITO still of the final state, 1600x1200 |
| `viz/large_run/` | the raw `out*` snapshots, for re-rendering |
| `viz/config_large` | the config that produced them |

Render commands (the whole thing takes a few seconds — OVITO's renderer is not
the bottleneck, the MD integrator is):

```sh
.deps/venv/bin/python tools/ovito_reader.py 'viz/large_run/out0*' viz/large_run/outend \
    --random-colors --anim --fps 2 --outlines --out viz/ovito_movie.mp4 --size 1280x960
python3 tools/web_viewer.py viz/large_run/out0* viz/large_run/outend \
    --color random -o viz/trajectory_large.html
```

One gotcha found while doing this: OVITO's `Viewport.camera_up` defaults to
**+Y**, but gravity here is along -Z, so without setting `camera_up = (0,0,1)`
the box renders lying on its side.  `tools/ovito_reader.py` now sets it.


---

# Part 2 — optimisation log

Starting point: the ported code before any performance work, built `-O2`.
Everything below is measured on B2 (2500 particles, 2750 steps) unless stated,
and `make check` reports the same 1e-12 delta throughout, so no contact decision
moved.

| step | B2 wall | ms/step | gain |
|---|---|---|---|
| baseline (`-O2`) | 78.58 s | 28.573 | — |
| `-O3` | 67.34 s | 24.486 | -14% |
| + drop dead work from `doOverlap` | 63.64 s | 23.142 | -19% |
| + one pose update per step instead of two | 59.98 s | 21.811 | -24% |
| + allocation-free pose matrices | 54.37 s | 19.771 | -31% |
| + allocation-free matrix inverse | 52.06 s | 18.929 | -34% |
| + GSL workspace reuse | ~52 s | ~18.9 | no measurable change |

**The C++17 migration is performance-neutral** — B2 55.76 s before, 55.22 s
after, both within the noise floor — which is what you would expect, since it
changed an RNG engine and an exception specification rather than anything on
the hot path.

Full re-baseline after the optimisation work (`bench/results.json`):

| run | before | after | delta |
|---|---|---|---|
| B1 | 5.43 s | 3.69 s | **-32.0%** |
| B2 | 78.58 s | 55.76 s | **-29.0%** |
| B3 | 58.52 s | 37.66 s | **-35.6%** |
| B4 | 782.82 s | 557.27 s | **-28.8%** |

## What worked

**Allocation-free matrix arithmetic.**  This was the whole game.  `matrixT`
owns a row-pointer array plus one array per row -- five heap allocations for a
4x4 -- and its operators each build a temporary.  `update_tranlation_mat()` ran
the obvious two-line expression once per particle per step, so ~40 malloc/free
pairs per particle-step, and a profile put about half the runtime in
`matrixT`'s destructor, `clone()` and constructor.

Two helpers, `matmul()` and `transpose_into()`, do the same arithmetic into
caller-owned storage with the accumulation order copied from
`matrixT::operator*=`, so results are bit-identical.  A third, `invert_into()`,
replaces `operator!`, which takes its argument by value and then inverts in
place -- so it cloned a 4x4 on every candidate contact pair.

**Calling `update_tranlation_mat()` once instead of twice per step.**  `calPos`
called `rotateTo()` then `moveto()`, and both end in that function, so the first
of the two full rebuilds was thrown away one line later.  `GeomObjectBase` gained
`setPose(q, x)` for this.

**Deleting work whose result was discarded.**  When `doOverlap` finds a
separating axis it used to carry on building the pole ray, intersecting it with
both ellipsoids and constructing a separating plane, then return -- and the
caller skips the contact block entirely in that case, so every write was dead.
That is the *common* path, because the bounding-sphere test that gates the call
is deliberately loose.

## What did not work

**`-march=native`.**  67.34 s with, 67.35 s without.  No effect on this
workload; it now lives behind `NATIVE=1` rather than being a default, which is
also the right call for a binary you might hand to somebody else.

**A conservative inscribed-sphere rejection test.**  An ellipsoid contains its
largest inscribed sphere, so if two particles' inscribed spheres overlap, the
ellipsoids overlap and the eigensolve cannot say otherwise.  It was implemented,
and then checked rather than assumed: a build with `-DVERIFY_INSCRIBED_SHORTCUT`
ran the eigensolve anyway and compared verdicts.

```
420 000 candidate pairs, shortcut fired 383 times (0.09%), 0 disagreements
```

Correct, and useless here: the inscribed radius is `min(a,b,c)` (0.036 for these
spheroids) against a circumscribed 0.07, so almost every overlapping pair sits
between the two bounds.  Removed again.  Worth recording because it is exactly
the kind of plausible optimisation that costs readability for nothing.

**Reusing the GSL eigensolver workspace.**  Strictly fewer allocations, no
measurable effect.  Kept, but not claimed as a win.

## Where the time is now

Re-profiled after the matrix work, B2:

| % | |
|---|---|
| 47 | `matrix<double>` -- still the destructor and `clone()` |
| 29 | the eigensolver (GSL 26.6 + the code's own Schur decomposition 2.6) |
| 9 | the allocator |
| rest | ellipsoid geometry, containers, everything else |

The next real win has to come from the eigensolver, which currently runs once
per candidate pair -- about 1500 times per step at 2500 particles.  The
characteristic polynomial is a quartic and `characteristicPolynom()` and
`CQuartic::solve()` already exist in the codebase (the call is commented out in
`doOverlap`), so replacing a 4x4 non-symmetric eigendecomposition with a quartic
root solve is the obvious next move.  It is a numerical change rather than a
mechanical one, so it needs the regression harness to demonstrate that the
contact decisions do not move near the complex/real boundary.

Further out: the remaining `matrix<double>` cost is in the type itself, and
replacing it with fixed-size stack types (or Eigen) would remove it.
