#!/usr/bin/env python3
"""Convert ellipmd ``out*`` snapshots into a LAMMPS dump file.

Why: OVITO's GUI cannot pick up a *new* file format without OVITO Pro, but it
reads LAMMPS dump files natively, including the column names that carry
non-spherical particle data.  Writing the columns as
``Orientation.X/.Y/.Z/.W`` and ``AsphericalShape.X/.Y/.Z`` makes OVITO map them
straight onto the ``Orientation`` and ``Aspherical Shape`` particle properties,
so the snapshots come up as oriented ellipsoids with no scripting at all.

Two naming traps, both avoided here:

* The LAMMPS-native names ``shapex shapey shapez`` carry an **automatic division
  by 2** (LAMMPS stores ellipsoid *diameters*), which would silently halve every
  particle.  The explicit ``AsphericalShape.*`` names have no such scaling.
* ``quati quatj quatk quatw`` do map to Orientation X, Y, Z, W, but the explicit
  ``Orientation.*`` names say the same thing without relying on that ordering.

    python3 tools/snapshot_to_dump.py out0* outend -o viz/large_run.dump
    open -a Ovito viz/large_run.dump          # or just drag it onto the window

The whole trajectory goes into one file (concatenated TIMESTEP blocks), so OVITO
shows it as a single animation with a frame slider.
"""

from __future__ import annotations

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ellipmd_io import (  # noqa: E402
    bounding_box, expand_paths, hsl_to_rgb, particle_color, read_snapshot,
)

PBC_FLAGS = {
    "solid": "ff ff ff",
    "periodic_x": "pp ff ff",
    "periodic_xy": "pp pp ff",
    "periodic_xyz": "pp pp pp",
}

COLUMNS = ("id type x y z "
           "Orientation.X Orientation.Y Orientation.Z Orientation.W "
           "AsphericalShape.X AsphericalShape.Y AsphericalShape.Z "
           "Color.R Color.G Color.B")


def write_dump(paths, out, pbc="solid", random_colors=True):
    flags = PBC_FLAGS[pbc]
    frames = 0
    with open(out, "w") as fh:
        for index, path in enumerate(paths):
            snap = read_snapshot(path)
            lo, hi = bounding_box(snap)

            fh.write("ITEM: TIMESTEP\n%d\n" % index)
            fh.write("ITEM: NUMBER OF ATOMS\n%d\n" % len(snap))
            fh.write("ITEM: BOX BOUNDS %s\n" % flags)
            fh.write("%.10g %.10g\n" % (lo[0], hi[0]))
            fh.write("%.10g %.10g\n" % (lo[1], hi[1]))
            fh.write("%.10g %.10g\n" % (lo[2], hi[2]))
            fh.write("ITEM: ATOMS %s\n" % COLUMNS)

            for i in range(len(snap)):
                x, y, z = snap.positions[i]
                a, b, c = snap.axes[i]
                w, qx, qy, qz = snap.quats[i]        # file order: scalar first
                if random_colors:
                    cr, cg, cb = hsl_to_rgb(*particle_color(i))
                else:
                    cr, cg, cb = 0.5, 0.6, 0.8
                # the columns are named Orientation.X/.Y/.Z/.W, so the values
                # must be written in X,Y,Z,W order -- i.e. scalar LAST
                fh.write("%d 1 %.10g %.10g %.10g %.10g %.10g %.10g %.10g "
                         "%.10g %.10g %.10g %.6f %.6f %.6f\n"
                         % (i + 1, x, y, z, qx, qy, qz, w, a, b, c,
                            cr, cg, cb))
            frames += 1
    return frames


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Convert ellipmd 'out*' snapshots to a LAMMPS dump file.")
    ap.add_argument("files", nargs="+", help="snapshot files or globs")
    ap.add_argument("-o", "--out", default="trajectory.dump",
                    help="output dump file (default: trajectory.dump)")
    ap.add_argument("--pbc", default="solid", choices=sorted(PBC_FLAGS),
                    help="boundary condition to record in BOX BOUNDS")
    ap.add_argument("--no-random-colors", action="store_true",
                    help="use a single uniform colour instead of the "
                         "per-particle random colours")
    args = ap.parse_args(argv)

    paths = expand_paths(args.files)
    if not paths:
        ap.error("no matching snapshot files")

    frames = write_dump(paths, args.out, args.pbc,
                        random_colors=not args.no_random_colors)
    size = os.path.getsize(args.out)
    print("wrote %s" % args.out)
    print("  %d frames, %.1f MB" % (frames, size / 1e6))
    print("  open it with:  open -a Ovito %s" % args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
