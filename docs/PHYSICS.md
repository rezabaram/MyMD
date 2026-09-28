# Physics

What the solver actually computes.  Every statement here was read off the code,
and the places where the code and the documentation disagreed are called out.

## Geometry

A particle is an ellipsoid defined in its own body frame by

```
x²/a² + y²/b² + z²/c² = 1
```

placed at a centre `Xc` and rotated by a unit quaternion `q`.  Writing `R(q)` for
that rotation, the implicit surface in world coordinates is

```
(R·(X − Xc))ᵀ · diag(1/a², 1/b², 1/c²) · (R·(X − Xc)) = 1
```

The code stores the quadratic form as a 4×4 homogeneous matrix, `ellip_mat`,
rebuilt from `q` and `Xc` whenever either changes
(`CEllipsoid::update_tranlation_mat`).  `CEllipsoid::operator()(X)` returns the
left-hand side minus one, so it is negative inside and positive outside; that is
the containment test used everywhere.  `gradient(X)` returns `2·ellip_mat·(X−Xc)`,
which is the outward normal up to a scale factor.

Inertia is that of a uniform solid ellipsoid, `I₁ = m(b²+c²)/5` and cyclic.

## Force law

Contacts are viscoelastic Hertzian, evaluated per contact as
(`Test::contactForce`, `include/interaction_force.h`):

```
δ     = penetration depth (Contact::dx_n)
v_n   = relative normal velocity of the contact points, v·n
n     = contact normal

F_n   = −(k·δ + γ·v_n)·√δ · n            (normal)
F_t   = −μ·|F_n|·(v − v_n·n)             (tangential, dynamic friction only)

F     = F_n + F_t
```

with `k = stiffness`, `γ = damping`, `μ = friction` from the config.  The
normal force is repulsive only: if the bracket would go negative, `δ` is clamped
so the force is zero rather than attractive.  There is **no static friction** and
**no cohesion** — see "Things that do not exist" below.

For an ellipsoid the force is applied at the contact point, giving the torque
`r × F` with `r` the contact point relative to the centre.  For a sphere the
contact normal passes through the centre and the torque vanishes, which is what
makes the elastic-bounce regression test a clean 1-D problem.

Gravity enters as a body force `m·g`, plus an optional bulk drag

```
F_drag = −fluiddampping · |g| · m · v
```

**Note that `fluiddampping` defaults to 0.05, not zero.**  Every run therefore has
a velocity-proportional damping term unless it is explicitly disabled, and it is
strong enough to dominate a small elastic bounce (it moves the effective
acceleration by several per cent at typical speeds).  Set `fluiddampping 0` when
you want to look at energy conservation.

## Contact detection

`CInteraction::overlaps` dispatches on the pair of shape types.  Three of the
four paths are analytic:

* **sphere–sphere**: overlap along the centre line.
* **sphere–box** and **ellipsoid–plane**: the extreme point of the particle
  toward the plane, giving `δ` directly.  (`CEllipsoid::point_to_plane` picks
  whichever of the two poles is nearer the plane, which selects the penetrating
  one.)
* **ellipsoid–box**: the six faces are tested as planes, skipping the non-solid
  ones.

Ellipsoid–ellipsoid is the hard one and is worth describing.  It follows the
"geometric potential" idea: two ellipsoids `E1`, `E2` overlap iff a certain
4×4 matrix pencil has the right eigenvalue structure, and the real eigenvector
gives the *pole* direction along which they interpenetrate.

```
doOverlap(ovs, E1, E2):
    M = −(!E1.ellip_mat) · E2.ellip_mat     // 4x4, one inverse + a multiply
    roots = solve_quartic(characteristic_polynomial(M))
    for r in roots:
        if |Im(r)| > eps: return true       // a complex root -> they interpenetrate
    return false                            // all real -> a separating axis exists
```

The eigenvalue test is a *quartic root solve*, not an eigensolver.  The two agree
exactly -- a build that ran both on the real matrices of a full deposition run
reported zero disagreements over 420,000 candidate pairs -- and the quartic is
much cheaper, which is what let GSL leave the contact path entirely (Phase 6b).
`characteristicPolynom()` and `CQuartic::solve()` were already in the repository;
the author had left the call commented out.

then, for a genuine contact, `intersect()` casts the ray between the two poles
and intersects it with each ellipsoid, and `findMin()` refines the contact points
by up to 500 iterations of a λ fixed point:

```
λ  = |(x−Xc₁)·E2·(x−Xc₂)|
x ← (Em2 + λ·Em1)⁻¹ · (Em2·Xc₂ + λ·Em1·Xc₁)
```

converged when the point moves less than 1e-13 and λ less than 1e-10.  Written in
caller-owned scratch, because as written above the loop built about five matrix
temporaries per iteration -- and it runs up to 500 times, twice per contact.  That
single change halved the runtime; see [`BENCHMARK.md`](BENCHMARK.md).

This is paid for *every candidate pair* before any cheap rejection test, which is
the cost of the method.  A conservative inscribed-sphere rejection test was
implemented and measured, and it fires on 0.09% of pairs -- the shape is too
elongated for its inscribed sphere to be useful; see `BENCHMARK.md`.

## Integrator

Beeman, with the position and rotation advanced in `CParticle::calPos` and the
velocities corrected in `calVel` after the forces are evaluated.

Translation and rotation use the same standard Beeman coefficients:

```
position:  r_{n+1} = r_n + v_n·dt + (2/3) a_n dt² − (1/6) a_{n−1} dt²
velocity:  v_{n+1} = v_n + (1/3) a_{n+1} dt + (5/6) a_n dt − (1/6) a_{n−1} dt
```

implemented as a predictor in `calPos` followed by a correction in `calVel`.
The rotation is advanced by integrating `q̇ = ½ ω ⊗ q` with one step and
renormalising, and the angular acceleration comes from Euler's equations in the
body frame:

```
I₁·ω̇₁ = τ₁ + ω₂ω₃(I₂−I₃)     and cyclic
```

For constant acceleration the scheme is exact, which is why free fall and the
ballistic part of a bounce have no energy drift.

Until recently the translational corrector was wrong — the `a_{n+1}` term had
the wrong sign and the `a_{n−1}` term was missing — making translation
first-order while rotation was second-order.  It cost 5% of the energy per
elastic bounce, and the error fell off only linearly as `dt` was reduced rather
than quadratically.  That is fixed, and
[`bench/reference/elastic_bounce`](../bench/reference/elastic_bounce) now guards
it.

## Initialisation

**`deposition`.**  Particles are rained onto the floor.  `celllist.setup()` sizes
a grid with cell size `2·maxRadii`; a layer adds one particle per grid cell in
the x–y plane at height `maxh + 1.02·maxRadii`, with a small random offset, and
another layer is added on a later step while the pile top is still below
`1 + 2·maxRadii`.  Each particle gets small random linear and angular velocity.

Note the `1.` in that gate is a hardcoded box height left over from a 1×1×1 box;
it should be `walls.L(2)`.  See `ROADMAP.md` Phase 4.

**`Stillinger`.**  Particles are placed on a simple cubic lattice of spacing
`dl = nParticle^(−1/3)`, with radius `0.2·dl·SizeDistribution.get()`, and random
orientations.  The packing is then compressed by multiplying every particle's
size by `scaling` at each write, so a run that writes more often compresses
faster.  Intended for periodic boundaries.

**`restart`.**  Reads a snapshot; positions, orientations and shapes are
restored, velocities are not.

## Shape parameters

For `particleType general` the semi-axes come from a nominal radius `r` and the
two shape parameters:

```
a = r · ζ^(1/3) / η^(1/3)
b = r / (ζ^(2/3) · η^(1/3))
c = r · η · (ζ/η)^(1/3)
```

giving the clean relations

```
c/a = eta        a/b = zeta        a·b·c = r³
```

so the volume `(4/3)πr³` is independent of the shape.  `eta` is the elongation
and `zeta` the flattening of the cross-section.  Both may be drawn from a
truncated Gaussian of relative width `etaWidth` / `zetaWidth`.

This is the reverse of the naming in the old top-level `README`, which claimed
`eta = a/b`.  The code is the authority.

For `prolate` and `oblate` the aspect ratio instead comes from `asphericity`
through the lookup tables in `include/map_asph_aspect.h`; for `sandstone` the
three semi-axes are drawn independently.

## Things that do not exist

Worth stating plainly, because the configuration surface suggests otherwise:

| | |
|---|---|
| `cohesion` | accepted, never used by any force law.  Now warns. |
| `static_friction` | accepted, never used.  Now warns. |
| `CMaterial::friction_threshold`, `CMaterial::color` | never read |
| `CParticle::avgforces` | set to the last contact force from a `static` shared across all particles; the only reader is an `if(0)` block |

So the model is: gravity, viscoelastic normal contact, Coulomb *dynamic*
friction, solid or periodic walls.  Nothing else.

## Validation

`make check` runs three reference cases:

| case | what it covers |
|---|---|
| `deposition` | 250 polydisperse spheroids settling into a dense layer — contact-dominated, so a force-law or contact-geometry change is caught immediately |
| `stillinger` | Stillinger lattice + `periodic_xyz` — the size-distribution parser, periodic image shifts, the compression path |
| `elastic_bounce` | a single sphere dropped onto the floor with no damping — **energy conservation**, which a snapshot comparison cannot see |

Tolerances are calibrated rather than guessed: recompiling the unchanged source
with `-O3 -march=native` perturbs positions by ~1e-12, while changing the contact
force by 0.01% moves them by ~9e-6, so the 1e-9 threshold separates noise from
signal with room to spare on both sides.

The tests were validated by deliberately breaking the physics and watching them
fail — a 1.0001 factor on the Hertzian force is caught by `deposition`, and the
old velocity corrector is caught by `elastic_bounce` on both the snapshot and
the energy check.
