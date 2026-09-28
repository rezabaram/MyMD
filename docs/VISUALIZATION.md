# Visualising the output

## Open something right now

If a run has already produced `out*` snapshots, these are the three fastest
routes from "files on disk" to "picture on screen":

```sh
# 1. Interactive, in the browser.  No installs.  One self-contained file.
open viz/trajectory_large.html

# 2. Just watch a movie.
open viz/ovito_movie.mp4

# 3. Real interactive 3D in OVITO (orbit, slice, colour by property, ...).
brew install --cask ovito                                   # one-off
python3 tools/snapshot_to_dump.py 'viz/large_run/out0*' viz/large_run/outend \
        -o viz/large_run.dump                               # one-off per run
open -a Ovito viz/large_run.dump                            # or drag it in
```

Route 3 is the one that matters: `tools/snapshot_to_dump.py` writes a LAMMPS dump
whose column names OVITO recognises, so the snapshots open as properly oriented
ellipsoids **with no scripting on your side**.  See "Opening snapshots in the
OVITO GUI" below for why the column names are what they are.

---

POV-Ray still works (`make pov FILE=out00010`), but it is the slowest possible
way to look at a trajectory: no interactivity, no colour coding, one image per
invocation.  There are three better options here, in increasing order of how
much you get and how much you have to install.

A snapshot contains only geometry, so every tool below shows the same thing:
the wall planes (`id 6`) and the ellipsoids (`id 14`, i.e.
`x y z a b c q0 q1 q2 q3`).  There are no velocities, forces or per-particle
stress in the file — see "Showing something other than shape" at the end if you
want to colour by a physical quantity.

| | install | interact | web/shareable | ray-traced stills | analysis |
|---|---|---|---|---|---|
| `make viewer` (built in) | none | yes | **yes, single HTML file** | no | no |
| **OVITO** | `brew install --cask ovito` (or `pip install ovito`) | yes | no | yes (Tachyon/OSPRay) | yes |
| PyVista + trame | `pip install pyvista trame …` | yes | yes (served) | via VTK | some |
| Blender | big download | yes | no | best-looking | no |

---

## 0. The live dashboard (watch a run as it happens)

```sh
make live                       # http://127.0.0.1:8770/
make live LIVE_CONFIG=config_quick PORT=8800
python3 tools/live_viewer.py --config config_rain --rundir viz/live
```

`tools/live_viewer.py` is a local web app: it serves a page with a **parameter
form**, Start and Stop, a progress bar, live statistics and a 3D view that picks
up each new snapshot as the solver writes it.  Edit the config in the page, press
Start, and watch the box fill.

Nothing is installed and nothing leaves the machine -- the server is the Python
standard library's `http.server`, and three.js comes from the same CDN the other
viewers use.

| endpoint | |
|---|---|
| `GET /` | the dashboard |
| `GET /api/status` | running, t, particles, energies, frame list, log tail |
| `GET /api/frame/<name>` | one snapshot as a float32 blob |
| `POST /api/run` | `{"config": "...", "seed": 1}` |
| `POST /api/stop` | terminate the run |

The page shares its scene, colour maps and per-frame update with
`web_viewer.py` (they live in `tools/viewer_common.py`), so the two render
identically.

`follow` keeps the newest frame on screen; turn it off, or drag the slider, to
scrub back through what has been produced so far while the run carries on.

---

## 1. The built-in web viewer (no installs)

```sh
make viewer                       # reads out0*  -> trajectory.html
make viewer OUT='out0001*' HTML=part1.html
python3 tools/web_viewer.py out0* -o trajectory.html   # same thing
make viewer COLOR=random          # every particle gets its own colour
open trajectory.html
```

`trajectory.html` is a **single self-contained file**: the whole trajectory is
embedded in it, so it can be e-mailed, put on a share, or attached to a paper's
supplementary material.  It gives you orbit/pan/zoom, a frame slider with the
simulation time from `log_energy`, play/pause (space bar), step with the arrow
keys, and colour-by-uniform / random / `c` / `a` / `c/a` / volume.  Z is up,
because that is where gravity points.  A frame with 10^4 particles is one draw
call.

### Exporting a video

Both viewers have a **record** button next to the colour selector, with a
resolution choice (540p / 720p / 1080p).  Pressing it plays the whole
trajectory through once and downloads a video file.

* **Only the 3D view is recorded.**  The capture is taken from the WebGL
  canvas, and the panels, buttons and sliders are separate DOM elements, so
  they cannot appear in the file.
* **The simulation box and the camera spin are always included**, whatever the
  view toggles happen to be set to -- a video is a rendered artefact, not a
  live view, and both are part of the render.
* **The camera turns one fixed step per frame** while recording.  The on-screen
  spin advances with wall-clock time, which would give a different angle per
  frame depending on how long each frame took to render.

#### For GitHub: use `make live`

GitHub plays **H.264 in an MP4 container** inline and nothing else, and a
browser's `MediaRecorder` will not necessarily give that -- Chrome offers WebM,
which GitHub will not play in a README.  So when the page is served by
`tools/live_viewer.py` and `ffmpeg` is on `PATH`, the record button renders each
frame to a PNG in the browser and pipes it to ffmpeg on the server, which
produces a file that works:

```
make live          # then press "record"
```

The encoder is invoked as

```
ffmpeg -f image2pipe -framerate 30 -i - -c:v libx264 -preset slow -crf 23 \
       -pix_fmt yuv420p -movflags +faststart out.mp4
```

`yuv420p` because anything else will not play in most players, and
`+faststart` to put the index at the front so it starts playing before it has
finished downloading.  The button shows the resulting size when it is done.

**Sizing.**  GitHub's web uploader takes files up to 10 MB, warns above 50 MB
and refuses above 100 MB -- and a large binary in git history is permanent.  For
a README, aim under 10 MB.  The lever is the resolution and the length: 540p or
720p keeps a full filling animation in single-digit megabytes at CRF 23, and
1080p does not.  If it comes out too big, record at 540p rather than shortening
the run.

To embed it, drag the file into the README editor on GitHub; it writes the
`<video>` or link markup for you.  A raw link to a `.mp4` in the repository
also works, but GitHub only serves it as a download rather than playing it
inline.

#### Without a server

The standalone `web_viewer.py` page has no server, so its button falls back to
the browser's own recorder.  On Safari that produces MP4; on Chrome, WebM --
the button says which before you press it, rather than naming a file `.mp4`
when it is not.

For a ray-traced mp4 instead of a screen recording, OVITO can do it offline and
uses the same ffmpeg:

```sh
make ovito
.deps/venv/bin/python tools/ovito_reader.py 'out0*' outend --anim --out movie.mp4 --fps 30
```

### Random colours

`--color random` (or picking *random per particle* in the page) gives every
particle its own colour.  The colour is a **deterministic hash of the particle
index**, not a fresh draw, so a particle keeps the same colour in every frame and
you can follow it through the trajectory.  The solver assigns indices in
insertion order (`CSys::add`), so this is stable across snapshots.

The same hash is implemented in `tools/ellipmd_io.py` (`particle_color`), in the
JavaScript of the generated page, and in the OVITO reader — all three produce
the identical colour, verified bit-for-bit over 5000 particles.

The page pulls three.js from a CDN, so the *first* open needs network access;
after that the browser caches it.  If you need it fully offline, download
`three.module.js` and `examples/jsm/` for version 0.186.1 next to the HTML and
edit the import map at the top.

Big trajectories make big pages (~40 bytes per particle per frame).  Use
`--stride` or `--first/--last` to keep the file manageable:

```sh
python3 tools/web_viewer.py out0* --stride 5 --last 60 -o movie.html
```

### Is the viewer lying to me?

Orientation is the one thing that is easy to get wrong (the C++ code stores
quaternions scalar-first `(w,x,y,z)`, while OVITO and three.js want
`(x,y,z,w)`), so there is a known-answer page:

```sh
make viewer-check
open orientation_sample.html
```

It draws seven ellipsoids whose long axes are known by construction: `+X`,
`+Y`, `-Z`, one triaxial ellipsoid showing all three different semi-axes, and a
`c/a = 1, 2, 3` ladder.  If that page looks right, the geometry pipeline is
right.  (`tools/orientation_sample.py` generates it.)

---

## 2. OVITO — the real replacement for the raster3d pipeline

[OVITO](https://ovito.org) is the standard tool for particle simulations now.
OVITO Basic is free and MIT-licensed; the `ovito` Python module on PyPI is
MIT-licensed too.  It handles ellipsoids natively: writing the
`Aspherical Shape` and `Orientation` particle properties is enough, the default
`Sphere` render mode switches to ellipsoid geometry automatically.

### Opening snapshots in the OVITO GUI

OVITO cannot auto-detect a *new* file format without OVITO Pro (that needs a
registered Python extension).  The way around it is to convert to a format OVITO
already reads, and LAMMPS dump is the natural one because its column names can
carry orientation and shape:

```sh
python3 tools/snapshot_to_dump.py 'out0*' outend -o run.dump
open -a Ovito run.dump        # or drag run.dump onto the OVITO window
```

The converter writes exactly these columns, which OVITO maps by name onto its
standard particle properties with no further configuration:

```
id type x y z  Orientation.X Orientation.Y Orientation.Z Orientation.W
               AsphericalShape.X AsphericalShape.Y AsphericalShape.Z
               Color.R Color.G Color.B
```

Two traps this avoids, both of which fail *silently*:

* `shapex shapey shapez` are the LAMMPS-native names but OVITO applies an
  **automatic division by 2** to them (LAMMPS stores ellipsoid diameters).
  Every particle would come out half-size.  `AsphericalShape.*` has no scaling.
  (`c_shape[1..3]` is unscaled; `c_diameter[1..3]` is halved.)
* A column named `Orientation.X` must contain the X component, so the values
  have to be written in `X, Y, Z, W` order.  The solver stores them scalar-first
  `(w, x, y, z)`, so writing them straight through shifts every component by
  one and rotates every particle wrongly.

The round trip is verified exact: position, shape, orientation and colour all
agree with the source snapshots to float32 precision (≤ 4e-9 for the axes,
≤ 1.5e-7 for the rotation matrix).

The whole trajectory goes into the one `.dump` file (concatenated `TIMESTEP`
blocks), so OVITO shows it as an animation with a frame slider, and
`--pbc periodic_xy` records the boundary condition in `BOX BOUNDS`.

**Python module** (scriptable, and what `tools/ovito_reader.py` targets):

```sh
make ovito                                    # creates .deps/venv, installs ovito (~72 MB)
make ovito-render FILE=out00010 PNG=frame.png
```

or directly:

```sh
.deps/venv/bin/python tools/ovito_reader.py 'out0*' --list-frames
.deps/venv/bin/python tools/ovito_reader.py out00010 --out frame.png --size 1600x1200
.deps/venv/bin/python tools/ovito_reader.py out00010 --out frame.png --random-colors
.deps/venv/bin/python tools/ovito_reader.py 'out0*' --out movie.mp4 --anim --fps 2 --outlines
.deps/venv/bin/python tools/ovito_reader.py out00010 --pbc periodic_xy
```

`--random-colors` writes a per-particle `Color` property using the same hash as
the HTML viewer, so the two show the same colouring.  `--outlines` draws dark
edges around particles, which helps tell touching neighbours apart; `--fps`
controls playback speed.

Two OVITO-specific gotchas this wrapper already handles:

* `Viewport.camera_up` defaults to **+Y**.  Gravity here is along -Z, so it is
  forced to `(0, 0, 1)`; otherwise the box renders lying on its side.
* OVITO rejects a partial set of FFmpeg options, so the codec, quality and
  executable are always passed together.

And in your own scripts, which is where OVITO really pays off:

```python
import sys; sys.path.insert(0, "tools")
from ovito_reader import EllipMDFileReader
from ovito.io import import_file
from ovito.modifiers import ColorCodingModifier

pipeline = import_file("out0*", input_format=EllipMDFileReader)
pipeline.modifiers.append(ColorCodingModifier(property="Aspherical Shape.X"))
data = pipeline.compute(200)
print(data.particles.count, data.cell[:])
```

That gets you colour coding by any property, slicing, the coordination
analysis, RDF, Wigner-Seitz/Voronoi, displacement vectors between frames, and
ray-traced output — none of which the POV-Ray pipeline could do.

---

## 3. Browser-based pipeline (PyVista + trame)

If you want a *hosted* viewer rather than a single file, ParaView's web stack
is the usual answer: [trame](https://kitware.github.io/trame/) runs a VTK
render in the browser and can [monitor a running
simulation](https://kitware.github.io/trame/blogs/monitor-your-simulation-in-your-web-browser-with-trame-and-catalyst.html).
PyVista is the friendly wrapper around the same machinery.  Sketch:

```python
import numpy as np, pyvista as pv
from ellipmd_io import read_snapshot, quat_to_matrix

snap = read_snapshot("out00010")
cloud = pv.PolyData(np.array(snap.positions))
cloud["shape"] = np.array(snap.axes)          # semi-axes, for scaling
# glyph a unit sphere per particle; orientation via the rotation matrices
# from quat_to_matrix(), scaling via "shape", then
#   plotter.export_html("scene.html")        # or plotter.show(...)
```

I have not run this path here (it needs a VTK install), so treat it as a
starting point rather than tested code.  `tools/ellipmd_io.py` has the parsing
and the quaternion→matrix conversion you need.

## 4. Blender

If the goal is a figure for a paper or a talk, Blender with the
[Molecular Nodes](https://www.blender.org/) ecosystem gives the best-looking
result.  Export the ellipsoids as a mesh (a unit sphere instanced with
position/orientation/scale — the same data `trajectory.html` embeds) and render
with Cycles.  This is the most work of the four options and I would only reach
for it once you know which frame you want.

---

## Showing something other than shape

Right now `CSys::output` writes only `x y z a b c q0 q1 q2 q3`, so colouring by
speed, contact force or coordination is not possible without changing the C++
side.  The cheap way to add it is to append extra columns in
`CSys::output(ostream&)` in `include/mdsys.h`:

```cpp
out << **it
    << "  " << (*it)->x(1).abs()                  // speed
    << "  " << (*it)->avgforces.abs()             // last contact force
    << endl;
```

then extend the column list in `tools/ellipmd_io.py` (`_to_floats(tokens[1:11])`
and the `STRIDE`-10 packing in `tools/web_viewer.py`, or `create_property` in
`tools/ovito_reader.py`).  OVITO's colour-coding modifier will then pick the new
property up for free, which is the main reason to prefer it over the built-in
viewer for actual analysis.

Two things that *are* already dumped and worth knowing about: `method Stillinger`
writes an `orient` file of contact normals in spherical coordinates, and the
`fabric` tool (`make tools && bin/fabric <snapshot>`) prints fabric tensor
eigenvalues per particle.
