# Snapshot file format

A specification of what `ellipmd` writes, so that downstream tools do not have
to reverse-engineer it.  The reference reader is `tools/ellipmd_io.py`, which is
standard-library Python only.

## Files in a run directory

| file | contents |
|---|---|
| `out00000`, `out00001`, … | snapshots, one every `outDt` of simulated time |
| `outend` | the final configuration, written after `maxTime` |
| `log_energy` | one line per snapshot: `t  E_total  E_kin  E_pot  E_rot` |
| `stdout.log` | only if you redirected it; the solver writes progress to `cerr` |

The filename is `<output><NNNNN>`, where `output` is the `output` config key
(default `out`) and the number is zero-padded to five digits.

## Snapshot files

A snapshot is a flat text file, one record per line, whitespace-separated.
Records are identified by their first field.  A snapshot is written in this
order: one header line, five wall-plane lines, then one line per particle.
Lines whose first field is neither `6` nor `14` are ignored by the reader, which
is what keeps the header (and the `id 5` lines the raster3d writer used to emit)
from breaking anything.

### Header line

```
# ellipmd 0.1.0  t=0.25
```

The version that wrote the file and the simulated time it represents.  The time
is what lets `method restart` resume the clock rather than restarting it at
zero — without it a restarted run was given the whole of `maxTime` over again.
A file with no header is treated as being at `t=0`, which is what snapshots
written before the header existed are.

### `id 6` — wall planes

Five lines, each `6` followed by three 3-D points:

```
6  x1 y1 z1  x2 y2 z2  x3 y3 z3
```

Together they describe the simulation box.  Each line gives three points on one
face; the five lines cover the box's extent.

> **Caveat.** `CBox::print` builds the second and third points with hard-coded
> unit offsets rather than the box lengths, so for anything other than a
> `1×1×1` box the points do not lie on the actual faces:
>
> ```
> 6   0  0  0     1  0  0     1  1  0
> 6   0  0  0     0  0  1     0  1  1
> 6   1  1  1.2   1  1  2.2   0  1  1
> ```
>
> above, a box of `L = (1, 1, 1.2)`.  The **first** point of each line, however,
> is always exactly either the box corner (lines 1–2) or `corner + L` (lines
> 3–5), so the bounding box of those five points recovers the true box.  That is
> what `ellipmd_io.bounding_box()` does.  Fixing the writer is on
> [`ROADMAP.md`](../ROADMAP.md).

### `id 14` — ellipsoid

One line per particle: eleven fields for the geometry, then six for the
velocities.

```
14  x y z  a b c  q0 q1 q2 q3   vx vy vz   wx wy wz
```

| field | meaning |
|---|---|
| `x y z` | centroid, in world coordinates |
| `a b c` | semi-axis lengths along the ellipsoid's **body** axes, in that order.  `a·b·c` is constant for a given particle size, so changing the shape parameters redistributes the axes at constant volume |
| `q0 q1 q2 q3` | orientation quaternion, **scalar first**: `q0` is the real part `w`, and `(q1,q2,q3)` the vector part |
| `vx vy vz` | linear velocity |
| `wx wy wz` | angular velocity |

Written with `setprecision(12)`.

The velocities were appended rather than inserted, so a reader that only wants
the geometry can take the first ten fields and ignore the rest — which is what
`tools/ellipmd_io.py` does, and what makes it work on files written before they
existed.  The author had left `//out<<"  "<<p.x(1);` commented out in
`CParticle`'s output operator, so this was the original intent.

#### Orientation convention

This is the part that bites.  Two conventions are in play:

* **This format stores `(w, x, y, z)`** — scalar first.  This comes from
  `Quaternion::print`, which writes `q.u` and then `q.v`.
* **OVITO and three.js both want `(x, y, z, w)`** — scalar last.

Any consumer must reorder.  Getting it wrong rotates every particle and fails
silently, so it is worth verifying against
`tools/orientation_sample.py`, which emits ellipsoids whose long axes are known
by construction (see [`VISUALIZATION.md`](VISUALIZATION.md)).

The rotation itself is a standard **active** rotation.  The world-space
direction of body axis `eᵢ` is column `i` of

```
R = ⎡ 1-2(y²+z²)   2(xy-wz)    2(xz+wy) ⎤
    ⎢ 2(xy+wz)    1-2(x²+z²)   2(yz-wx) ⎥
    ⎣ 2(xz-wy)     2(yz+wx)   1-2(x²+y²)⎦
```

which is what `Quaternion::toWorld` computes.  Equivalently, the ellipsoid is the
set of points `X` with `(R·(X−Xc))ᵀ · diag(1/a², 1/b², 1/c²) · (R·(X−Xc)) = 1`.

### Shape parameter algebra

For `particleType general` the semi-axes are built from the config's `eta` and
`zeta` and a nominal radius `r`:

```
a = r · ζ^(1/3) / η^(1/3)
b = r / (ζ^(2/3) · η^(1/3))
c = r · η · (ζ/η)^(1/3)
```

whose ratios are simply

```
c/a = eta        a/b = zeta        a·b·c = r³
```

So `eta` controls the elongation and `zeta` the flattening of the cross-section,
and the volume is exactly `(4/3)πr³` whatever they are.  (The old top-level
`README` claimed `eta = a/b`; that was wrong.)

## `log_energy`

```
t  E_total  E_kin  E_pot  E_rot
```

one line per snapshot, written with `setprecision(14)`.

* `E_kin` = Σ ½m|v|²
* `E_pot` = Σ −m (g · x)
* `E_rot` = Σ ½ (I₁ω₁² + I₂ω₂² + I₃ω₃²), with ω in body coordinates
* `E_total` = the sum

The first line corresponds to the initial state, before any force has been
evaluated.  It used to be uninitialised garbage (`~1e-314`); it now holds real
values computed from the initial configuration.

## Reading a run

```python
import sys; sys.path.insert(0, "tools")
from ellipmd_io import read_snapshot, expand_paths, bounding_box, quat_to_matrix

snap = read_snapshot("out00010")
print(len(snap), "ellipsoids, box", bounding_box(snap))
print("particle 0:", snap.positions[0], snap.axes[0], snap.quats[0])
R = quat_to_matrix(snap.quats[0])          # body -> world rotation
```

Or from the shell, for a summary of a whole trajectory:

```sh
python3 tools/ellipmd_io.py --summary out0* outend
```

## Continuations

`method restart` reads a snapshot back in through `CPacking::parse`, which takes
any file in this format (it skips the header and the `id 6` lines).  Positions,
orientations, shapes and velocities are restored, the clock is resumed from the
header, and the accelerations are primed by evaluating the forces once.

It is **not bit-exact**, and the reason is worth knowing: Beeman advances the
position from `a_n` and `a_{n-1}`, and the format carries neither.  Priming gives
both, but `a_{n-1}` properly belongs to the previous step, so the first step
differs by the `(1/6)(a_n − a_{n-1})dt²` term.  Measured on the deposition case
with `dt=1e-4`:

| | |
|---|---|
| one step past the restart point | 8.3e-8 |
| 500 steps (0.05 of simulated time) | 3.9e-4, about 1% of a particle radius |

That is the amplification of a one-step difference by dense contact dynamics,
not a defect in the restart.  Storing the acceleration history too would make it
exact, at the cost of six more columns per particle.
