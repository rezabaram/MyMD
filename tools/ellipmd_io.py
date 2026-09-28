#!/usr/bin/env python3
"""Reader for the ``out*`` snapshot format written by ellipmd (MyMD).

Only the Python standard library is used, so this runs with the system
``python3`` -- no numpy, no virtualenv.

A snapshot file contains

  * ``id 6``  lines -- five lines describing the wall planes as three points
    each: ``6  x y z  x y z  x y z``
  * ``id 14`` lines -- one per ellipsoid:
    ``14  x y z  a b c  q0 q1 q2 q3``

where ``(x,y,z)`` is the centroid, ``(a,b,c)`` the three semi-axis lengths and
``(q0,q1,q2,q3)`` the orientation quaternion **in scalar-first order**
(``q0`` is the real part ``w``; this is ``Quaternion::print`` in
``include/quaternion.h``).

The code's ``Quaternion::toWorld`` is a standard active rotation
``v' = (w^2-|v|^2) v + 2 v (v.v) + 2 w (v x v)``, so the world-space
principal axes of an ellipsoid are the columns of ``R(q)`` scaled by
``(a, b, c)``.  ``quat_to_matrix`` below returns that ``R``.

Run this file directly for a quick look at a file or a whole trajectory::

    python3 tools/ellipmd_io.py out00000
    python3 tools/ellipmd_io.py --summary out0*
"""

from __future__ import annotations

import argparse
import glob
import math
import os
import re
import sys

PLANE_ID = "6"
ELLIPSOID_ID = "14"

__all__ = [
    "Snapshot",
    "read_snapshot",
    "read_times",
    "expand_paths",
    "quat_to_matrix",
    "bounding_box",
    "principal_axes",
    "particle_color",
    "hsl_to_rgb",
]


class Snapshot:
    """One configuration: wall planes plus a list of ellipsoids."""

    __slots__ = ("path", "label", "planes", "positions", "axes", "quats")

    def __init__(self, path="", label="", planes=None, positions=None,
                 axes=None, quats=None):
        self.path = path
        self.label = label
        # planes: list of 3 points, each a 3-tuple
        self.planes = planes if planes is not None else []
        # positions/axes: N x 3 lists of floats, quats: N x 4 (w, x, y, z)
        self.positions = positions if positions is not None else []
        self.axes = axes if axes is not None else []
        self.quats = quats if quats is not None else []

    def __len__(self):
        return len(self.positions)

    def __repr__(self):
        return "Snapshot(%s, %d ellipsoids, %d planes)" % (
            os.path.basename(self.path), len(self), len(self.planes))


def _to_floats(tokens):
    out = []
    for t in tokens:
        try:
            out.append(float(t))
        except ValueError:
            return None
    return out


# --------------------------------------------------------------- colours
# A deterministic pseudo-random colour per particle, stable across frames so
# you can follow one particle through the trajectory.  Implemented with
# splitmix32 so that Python, the JavaScript in tools/web_viewer.py (Math.imul
# with the same constants) and OVITO all produce the identical colour.

def _splitmix32(i):
    x = (i + 0x9E3779B9) & 0xFFFFFFFF
    x = ((x ^ (x >> 16)) * 0x21F0AAAD) & 0xFFFFFFFF
    x = ((x ^ (x >> 15)) * 0x735A2D97) & 0xFFFFFFFF
    return (x ^ (x >> 15)) & 0xFFFFFFFF


def particle_color(index):
    """``(h, s, l)`` with h, s, l in [0, 1] for particle ``index``."""
    x = _splitmix32(index)
    h = (x & 0xFFFF) / 65536.0
    s = 0.55 + ((x >> 16) & 0xFF) / 255.0 * 0.35
    l = 0.42 + ((x >> 24) & 0xFF) / 255.0 * 0.22
    return h, s, l


def _hue_to_rgb(p, q, t):
    if t < 0.0:
        t += 1.0
    if t > 1.0:
        t -= 1.0
    if t < 1.0 / 6.0:
        return p + (q - p) * 6.0 * t
    if t < 1.0 / 2.0:
        return q
    if t < 2.0 / 3.0:
        return p + (q - p) * 6.0 * (2.0 / 3.0 - t)
    return p


def hsl_to_rgb(h, s, l):
    """sRGB triple in [0, 1].  Same algorithm as three.js ``Color.setHSL``."""
    h = h % 1.0
    s = min(1.0, max(0.0, s))
    l = min(1.0, max(0.0, l))
    if s == 0.0:
        return (l, l, l)
    q = l * (1.0 + s) if l < 0.5 else l + s - l * s
    p = 2.0 * l - q
    return (_hue_to_rgb(p, q, h + 1.0 / 3.0),
            _hue_to_rgb(p, q, h),
            _hue_to_rgb(p, q, h - 1.0 / 3.0))


def read_snapshot(path):
    """Parse one ``out*`` file.  Unknown lines (including the ``id 5`` ray
    lines produced by the raster3d writer) are ignored."""
    planes = []
    positions = []
    axes = []
    quats = []

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            tokens = line.split()
            if not tokens:
                continue
            tag = tokens[0]
            if tag == PLANE_ID:
                vals = _to_floats(tokens[1:10])
                if vals and len(vals) == 9:
                    planes.append((tuple(vals[0:3]), tuple(vals[3:6]),
                                   tuple(vals[6:9])))
            elif tag == ELLIPSOID_ID:
                vals = _to_floats(tokens[1:11])
                if vals and len(vals) == 10:
                    positions.append(tuple(vals[0:3]))
                    axes.append(tuple(vals[3:6]))
                    quats.append(tuple(vals[6:10]))

    return Snapshot(path, _label_from_path(path), planes, positions, axes, quats)


def _label_from_path(path):
    base = os.path.basename(path)
    m = re.search(r"(\d+)$", base)
    if m:
        return str(int(m.group(1)))
    return base


def read_times(directory="."):
    """Times from ``log_energy`` (first column), one per written snapshot.

    Returns an empty list if the file is missing.  Note that the *first* line
    of ``log_energy`` has correct time but uninitialised energies (the code
    writes it before computing them)."""
    path = os.path.join(directory, "log_energy")
    times = []
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                tokens = line.split()
                if not tokens:
                    continue
                try:
                    times.append(float(tokens[0]))
                except ValueError:
                    pass
    except OSError:
        return []
    return times


def expand_paths(patterns):
    """Expand shell globs, drop ``log*``/``config*`` noise and sort the
    ``out00000``-style names numerically, with ``outend`` last."""
    files = []
    for pattern in patterns:
        hits = sorted(glob.glob(pattern))
        if not hits and os.path.isfile(pattern):
            hits = [pattern]
        for h in hits:
            if h not in files:
                files.append(h)

    def sort_key(p):
        base = os.path.basename(p)
        m = re.search(r"(\d+)$", base)
        if m:
            return (0, int(m.group(1)), base)
        return (1, 0, base)

    return sorted(files, key=sort_key)


def quat_to_matrix(q):
    """Rotation matrix (row-major 3x3 tuple of tuples) for ``q = (w,x,y,z)``.

    Matches ``Quaternion::toWorld`` in include/quaternion.h: the world-space
    direction of body axis ``e_i`` is column ``i`` of the returned matrix."""
    w, x, y, z = q
    n = math.sqrt(w * w + x * x + y * y + z * z)
    if n == 0.0:
        return ((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0))
    w, x, y, z = w / n, x / n, y / n, z / n
    return (
        (1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)),
        (2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)),
        (2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)),
    )


def principal_axes(quat, axes):
    """The three world-space semi-axis vectors of an ellipsoid."""
    a, b, c = axes
    r = quat_to_matrix(quat)
    out = []
    for col, length in enumerate((a, b, c)):
        out.append(tuple(length * r[row][col] for row in range(3)))
    return out


def bounding_box(snapshot, fallback=((0.0, 0.0, 0.0), (1.0, 1.0, 1.0))):
    """Axis-aligned box of the simulation cell.

    ``CBox::print`` (include/box.h) emits five lines and uses hard-coded unit
    offsets instead of ``L`` for the second and third points, so the AABB of
    *all* the points it writes overshoots whenever ``L > 1``.  The *first*
    point of each line, however, is exactly ``corner`` (lines 1-2) or
    ``corner + L`` (lines 3-5), so the AABB of just those anchors recovers the
    real box."""
    anchors = [tri[0] for tri in snapshot.planes]
    if len(anchors) >= 2:
        lo = tuple(min(p[i] for p in anchors) for i in range(3))
        hi = tuple(max(p[i] for p in anchors) for i in range(3))
        if all(hi[i] > lo[i] for i in range(3)):
            return (lo, hi)

    pts = [p for tri in snapshot.planes for p in tri] or list(snapshot.positions)
    if not pts:
        return fallback
    lo = tuple(min(p[i] for p in pts) for i in range(3))
    hi = tuple(max(p[i] for p in pts) for i in range(3))
    if any(hi[i] - lo[i] <= 0 for i in range(3)):
        return fallback
    return (lo, hi)


def _main(argv=None):
    ap = argparse.ArgumentParser(
        description="Inspect ellipmd 'out*' snapshot files.")
    ap.add_argument("files", nargs="+", help="snapshot files or globs")
    ap.add_argument("--summary", action="store_true",
                    help="one line per file instead of full detail")
    args = ap.parse_args(argv)

    paths = expand_paths(args.files)
    if not paths:
        ap.error("no matching files")

    times = read_times(os.path.dirname(os.path.abspath(paths[0])) or ".")

    for i, path in enumerate(paths):
        snap = read_snapshot(path)
        if args.summary:
            lo, hi = bounding_box(snap)
            print("%-16s N=%-6d planes=%d  box=%s..%s" % (
                os.path.basename(path), len(snap), len(snap.planes),
                tuple(round(v, 3) for v in lo), tuple(round(v, 3) for v in hi)))
            continue

        print("%s  ->  %r" % (path, snap))
        if i < len(times):
            print("  t = %g" % times[i])
        lo, hi = bounding_box(snap)
        print("  box: %s .. %s" % (
            tuple(round(v, 4) for v in lo), tuple(round(v, 4) for v in hi)))
        for n in range(min(3, len(snap))):
            a, b, c = snap.axes[n]
            pa = principal_axes(snap.quats[n], snap.axes[n])
            print("  #%d x=%s abc=(%.5f, %.5f, %.5f) q=%s" % (
                n, tuple(round(v, 5) for v in snap.positions[n]), a, b, c,
                tuple(round(v, 5) for v in snap.quats[n])))
            print("      axis0=%s" % (tuple(round(v, 5) for v in pa[0]),))
    return 0


if __name__ == "__main__":
    sys.exit(_main())
