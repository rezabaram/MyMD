#!/usr/bin/env python3
"""Write a single synthetic snapshot whose orientation is known by
construction, so the viewer (or OVITO) can be sanity-checked at a glance.

    python3 tools/orientation_sample.py data/orientation_sample
    python3 tools/web_viewer.py data/orientation_sample -o orientation_sample.html

The file contains seven ellipsoids in a 1 x 1 x 1 box.  The first four are a
compass: all have semi-axes ``a=0.18, b=0.06, c=0.03``, so the *body X* axis
is the long one, and each is rotated to point somewhere different:

  ============  =========================  ========================
  position      orientation               long axis should point
  ============  =========================  ========================
  (0.25,0.20)   identity                   +X  (right)
  (0.25,0.50)   90 deg about Z             +Y  (into the screen)
  (0.25,0.80)   90 deg about Y             -Z  (down, if Z is up)
  (0.72,0.50)   identity, a=0.20 b=0.12    +X long, +Y medium,
                c=0.05                     +Z thin
  ============  =========================  ========================

The last three are a ``c/a`` ladder (``1, 2, 3``) at ``y = 0.15, 0.50, 0.85``
with the body Z axis rotated onto world X, so they lie horizontally in a row
with increasing length; they exist so the viewer's *colour by aspect ratio*
mode has a spread to show.

If the long axes do not come out as listed, the quaternion is being applied
in the wrong order or with the wrong handedness.  The C++ side stores
``(q0,q1,q2,q3) = (w,x,y,z)`` (scalar first, see ``Quaternion::print``), while
OVITO and three.js both want ``(x,y,z,w)``.
"""

from __future__ import annotations

import math
import os
import sys

BOX = (1.0, 1.0, 1.0)

PLANE_LINES = [
    ((0, 0, 0), (1, 0, 0), (1, 1, 0)),
    ((0, 0, 0), (0, 0, 1), (0, 1, 1)),
    ((1, 1, 1), (1, 1, 2), (1, 2, 2)),
    ((1, 1, 1), (1, 1, 2), (2, 1, 2)),
    ((1, 1, 1), (2, 1, 1), (2, 2, 1)),
]


def q_axis_angle(axis, degrees):
    """Quaternion (w, x, y, z) for a rotation about ``axis``."""
    n = math.sqrt(sum(c * c for c in axis))
    ax = tuple(c / n for c in axis)
    h = math.radians(degrees) / 2.0
    s = math.sin(h)
    return (math.cos(h), ax[0] * s, ax[1] * s, ax[2] * s)


def q_mul(p, q):
    pw, px, py, pz = p
    qw, qx, qy, qz = q
    return (
        pw * qw - px * qx - py * qy - pz * qz,
        pw * qx + px * qw + py * qz - pz * qy,
        pw * qy - px * qz + py * qw + pz * qx,
        pw * qz + px * qy - py * qx + pz * qw,
    )


def build_lines():
    out = []
    for pts in PLANE_LINES:
        flat = "  ".join("%g" % v for p in pts for v in p)
        out.append("6   " + flat)
    return out


def build_ellipsoids():
    ident = (1.0, 0.0, 0.0, 0.0)
    rz90 = q_axis_angle((0, 0, 1), 90)
    ry90 = q_axis_angle((0, 1, 0), 90)

    rows = [
        # x,    y,    z,    a,    b,    c,    q
        # --- orientation compass (body X is the long axis) ---
        (0.25, 0.20, 0.30, 0.18, 0.06, 0.03, ident),   # long axis -> +X
        (0.25, 0.50, 0.30, 0.18, 0.06, 0.03, rz90),    # long axis -> +Y
        (0.25, 0.80, 0.30, 0.18, 0.06, 0.03, ry90),    # long axis -> -Z
        # triaxial: shows all three semi-axes at once, unrotated
        (0.72, 0.50, 0.30, 0.20, 0.12, 0.05, ident),
        # --- c/a ladder, body Z (the long axis) rotated onto world X ---
        (0.30, 0.15, 0.75, 0.06, 0.06, 0.06, ry90),    # c/a = 1
        (0.30, 0.50, 0.75, 0.06, 0.06, 0.12, ry90),    # c/a = 2
        (0.30, 0.85, 0.75, 0.06, 0.06, 0.18, ry90),    # c/a = 3
    ]

    lines = []
    for x, y, z, a, b, c, q in rows:
        lines.append("14   %s   %s   %s" % (
            "  ".join("%.6f" % v for v in (x, y, z)),
            "  ".join("%.6f" % v for v in (a, b, c)),
            "  ".join("%.9f" % v for v in q)))
    return lines


def main(argv):
    path = argv[1] if len(argv) > 1 else "orientation_sample"
    text = "\n".join(build_lines() + build_ellipsoids()) + "\n"
    with open(path, "w") as fh:
        fh.write(text)
    print("wrote %s (%d ellipsoids)" % (path, len(build_ellipsoids())))
    print("view with:  python3 tools/web_viewer.py %s -o orientation_sample.html"
          % path)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
