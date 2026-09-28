#!/usr/bin/env python3
"""OVITO integration for ellipmd ``out*`` snapshots.

OVITO (https://ovito.org) is the modern replacement for the
``coord_convert | render`` (raster3d) pipeline: real-time OpenGL/Vulkan
viewing, a proper camera, ray-traced rendering, colour coding, slicing,
movies, and analysis.  OVITO Basic is free and MIT-licensed, and the
``ovito`` Python module on PyPI is MIT-licensed too.

Install::

    pip install ovito

Two ways to use this file
-------------------------

1. Standalone -- render an image (or a movie) from the command line::

       python3 tools/ovito_reader.py out0* --out frame.png
       python3 tools/ovito_reader.py out0* --out movie.mp4

2. As a reader inside your own OVITO Python scripts::

       from ovito_reader import EllipMDFileReader
       from ovito.io import import_file
       pipeline = import_file("/path/to/out0*", input_format=EllipMDFileReader)

   Dropping the reader into OVITO's GUI as an auto-detected extension needs
   OVITO Pro; from the free Python module the explicit ``input_format=``
   argument above is the supported route.

Conventions
-----------

* ``Aspherical Shape`` gets the file's ``(a, b, c)`` semi-axes.
* ``Orientation`` gets the quaternion reordered from the file's scalar-first
  ``(q0,q1,q2,q3) = (w,x,y,z)`` (``Quaternion::print``) to OVITO's Shoemake
  order ``(x,y,z,w)``.  ``Quaternion::toWorld`` is an active rotation, which
  matches how OVITO applies ``Orientation``.
* ``Radius`` is ``max(a,b,c)``, so an OVITO build that falls back to spheres
  still shows something sane.

Periodic boundaries are not recorded in the snapshot files, so the cell is
imported as non-periodic by default; set ``EllipMDFileReader.pbc`` (or pass
``--pbc`` on the command line) if you want the box drawn as periodic.
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ellipmd_io import (  # noqa: E402
    bounding_box, expand_paths, hsl_to_rgb, particle_color, read_snapshot,
)

try:
    from ovito.data import DataCollection
    from ovito.io import FileReaderInterface, import_file
    from ovito.vis import Viewport
    HAVE_OVITO = True
except ImportError:  # keep the module importable without ovito installed
    HAVE_OVITO = False
    FileReaderInterface = object


PBC_OPTIONS = {
    "solid": (False, False, False),
    "periodic_x": (True, False, False),
    "periodic_xy": (True, True, False),
    "periodic_xyz": (True, True, True),
}


class EllipMDFileReader(FileReaderInterface):
    """OVITO file reader for the ellipmd ``out*`` snapshot format."""

    #: set to "solid", "periodic_x", "periodic_xy" or "periodic_xyz"
    pbc = "solid"

    #: write a per-particle `Color` property (deterministic random colours)
    random_colors = False

    @staticmethod
    def detect(filename: str) -> bool:
        """Claim files whose first records contain an ``id 14`` ellipsoid."""
        try:
            with open(filename, "r", errors="replace") as fh:
                for _ in range(50):
                    line = fh.readline()
                    if not line:
                        break
                    tokens = line.split()
                    if tokens and tokens[0] == "14":
                        return True
        except OSError:
            return False
        return False

    def scan(self, filename, register_frame):
        """One frame per file.  Supports a glob (``out0*``) so a whole
        trajectory can be opened as a single OVITO pipeline."""
        paths = expand_paths([filename]) if _has_magic(filename) else [filename]
        for path in paths:
            register_frame(frame_info=path, label=os.path.basename(path))
        if not paths:
            register_frame(frame_info=filename, label=os.path.basename(filename))

    def parse(self, data: DataCollection, filename, frame_index=0,
              frame_info=None, **kwargs):
        path = frame_info if isinstance(frame_info, str) else filename
        snapshot = read_snapshot(path)
        n = len(snapshot)

        lo, hi = bounding_box(snapshot)
        length = [max(hi[k] - lo[k], 1e-9) for k in range(3)]
        cell = data.create_cell(
            matrix=[[length[0], 0.0, 0.0, lo[0]],
                    [0.0, length[1], 0.0, lo[1]],
                    [0.0, 0.0, length[2], lo[2]]],
            pbc=PBC_OPTIONS.get(self.pbc, (False, False, False)))

        particles = data.create_particles(count=n)
        positions = particles.create_property("Position")
        shapes = particles.create_property("Aspherical Shape")
        orientations = particles.create_property("Orientation")
        radius = particles.create_property("Radius")

        # Random colour per particle.  Keyed on the particle index, which the
        # solver assigns in insertion order, so a particle keeps its colour
        # across frames.  Same hash as tools/web_viewer.py.
        colors = None
        if self.random_colors:
            colors = particles.create_property("Color")

        for i in range(n):
            positions[i] = snapshot.positions[i]
            shapes[i] = snapshot.axes[i]
            w, x, y, z = snapshot.quats[i]          # file order: scalar first
            orientations[i] = (x, y, z, w)          # OVITO order: scalar last
            radius[i] = max(snapshot.axes[i])
            if colors is not None:
                colors[i] = hsl_to_rgb(*particle_color(i))

        # OVITO's default `Sphere` shape mode renders ellipsoids automatically
        # as soon as the `Aspherical Shape` property is present, so there is no
        # visual element to configure here.
        return None


def _has_magic(pattern):
    return glob.has_magic(pattern)


def _main(argv=None):
    ap = argparse.ArgumentParser(
        description="Import ellipmd snapshots into OVITO and render them.")
    ap.add_argument("files", nargs="+", help="snapshot files or a glob")
    ap.add_argument("--out", default=None,
                    help="output image/movie file (omit to just test the import)")
    ap.add_argument("--size", default="1200x900",
                    help="image size, e.g. 1200x900 (default)")
    ap.add_argument("--pbc", default="solid", choices=sorted(PBC_OPTIONS),
                    help="periodic boundary condition to draw (default solid)")
    ap.add_argument("--random-colors", action="store_true",
                    help="give every particle its own deterministic colour")
    ap.add_argument("--frame", type=int, default=0,
                    help="frame index to work on (default 0)")
    ap.add_argument("--anim", action="store_true",
                    help="render the whole trajectory to a movie (implies "
                         "that --out should be e.g. movie.mp4)")
    ap.add_argument("--fps", type=float, default=10.0,
                    help="movie playback rate (default 10)")
    ap.add_argument("--quality", default="high",
                    choices=["low", "medium", "high"],
                    help="movie encoding quality (default high)")
    ap.add_argument("--outlines", action="store_true",
                    help="draw dark outlines around particles, which makes "
                         "touching particles easier to tell apart")
    ap.add_argument("--list-frames", action="store_true",
                    help="print the frames the reader found and stop")
    args = ap.parse_args(argv)

    if not HAVE_OVITO:
        print("error: the 'ovito' module is not installed.\n"
              "       pip install ovito\n"
              "       (MIT-licensed; macOS arm64 wheels available)",
              file=sys.stderr)
        return 2

    EllipMDFileReader.pbc = args.pbc
    EllipMDFileReader.random_colors = args.random_colors

    location = args.files[0] if len(args.files) == 1 else args.files
    pipeline = import_file(location, input_format=EllipMDFileReader)
    pipeline.add_to_scene()

    if args.list_frames:
        paths = expand_paths(args.files)
        print("%d frame(s) registered by the reader" % pipeline.source.num_frames)
        for i, path in enumerate(paths):
            print("%4d  %s" % (i, os.path.basename(path)))
        return 0

    data = pipeline.compute(args.frame)
    count = data.particles.count
    print("frame %d: %d ellipsoids, cell %s"
          % (args.frame, count,
             tuple(round(v, 4) for v in data.cell[:, :3].diagonal())))

    if not count:
        print("error: no particles in this frame", file=sys.stderr)
        return 1

    if args.out is None:
        print("import OK (use --out to render)")
        return 0

    vp = Viewport(type=Viewport.Type.Perspective,
                  camera_dir=(1.15, -1.75, 1.0))
    # OVITO's default camera_up is +Y, but gravity here points along -Z, so
    # without this the box would be rendered lying on its side.
    try:
        vp.camera_up = (0.0, 0.0, 1.0)
    except Exception as exc:                        # pragma: no cover
        print("note: could not set camera_up to Z (%s)" % exc, file=sys.stderr)
    vp.zoom_all()
    if args.anim:
        # OVITO insists on all-or-none of the FFmpeg parameters
        vp.render_anim(filename=args.out, size=_parse_size(args.size),
                       fps=args.fps,
                       ffmpeg_executable="ffmpeg",
                       ffmpeg_codec="libx264",
                       ffmpeg_quality=args.quality,
                       outlines_enabled=args.outlines)
    else:
        pipeline.compute(args.frame)   # make sure this frame is loaded
        vp.render_image(filename=args.out, size=_parse_size(args.size),
                        frame=args.frame, outlines_enabled=args.outlines)
    print("wrote %s" % args.out)
    return 0


def _parse_size(text):
    m = re.match(r"^\s*(\d+)\s*[xX*]\s*(\d+)\s*$", text)
    if not m:
        raise argparse.ArgumentTypeError("size must look like 1200x900")
    return (int(m.group(1)), int(m.group(2)))


if __name__ == "__main__":
    sys.exit(_main())
