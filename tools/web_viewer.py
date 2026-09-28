#!/usr/bin/env python3
"""Turn a sequence of ellipmd ``out*`` snapshots into one self-contained
interactive HTML page (three.js, orbit/zoom, frame slider).

    python3 tools/web_viewer.py out0* --out trajectory.html
    open trajectory.html

No server, no Python packages, no OVITO/ParaView install.  The page pulls
three.js from a CDN, so it needs network access the first time it is opened;
everything else (the whole trajectory) is embedded in the file, so the HTML
can be e-mailed or dropped on a shared drive as-is.

The page draws each ellipsoid with an ``InstancedMesh``: position, quaternion
and the three semi-axes are uploaded per instance, so a frame with 10^4
particles is one draw call.
"""

from __future__ import annotations

import argparse
import base64
import os
import struct
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ellipmd_io import (  # noqa: E402
    bounding_box, expand_paths, read_snapshot, read_times,
)

THREE_VERSION = "0.186.1"

from viewer_common import (  # noqa: E402
    CSS, ERROR_BLOCK, THREE_VERSION, colour_ranges, importmap, note_ranges,
    pack_snapshot, scene_js,
)

PAGE_HEAD = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>__CSS__
</style>
</head>
<body>
<div id="view"></div>
<div id="hud" class="panel">
  <b>__TITLE__</b><br>
  <span id="count"></span> &middot; box __BOXSIZE__<br>
  <span class="dim">drag = orbit &middot; shift+drag = pan &middot; wheel = zoom &middot; z is up</span>
</div>
<div id="ui" class="panel">
  <button id="play" title="play/pause (space)">&#9654;</button>
  <input id="slider" type="range" min="0" max="__LAST__" value="0" step="1">
  <span id="label"></span>
  <label>colour
    <select id="color">
      <option value="uniform">uniform</option>
      <option value="random">random per particle</option>
      <option value="c">long semi-axis c</option>
      <option value="a">short semi-axis a</option>
      <option value="aspect">aspect ratio c/a</option>
      <option value="volume">volume</option>
    </select>
  </label>
  <label><input id="box" type="checkbox" checked> box</label>
  <label><input id="spin" type="checkbox"> spin</label>
  <select id="vres" title="video resolution">
      <option value="960x540">540p</option>
      <option value="1280x720" selected>720p</option>
      <option value="1920x1080">1080p</option>
    </select>
  <button id="export" title="record the 3D view to a video file">record</button>
  <button id="savepng" title="save the current frame as a PNG">save png</button>
</div>
""" + ERROR_BLOCK + r"""
<script id="trajectory" type="application/json">__DATA__</script>
"""

PAGE_TAIL = r"""
</script>
</body>
</html>
"""

PLAYBACK_JS = r"""
const DATA = JSON.parse(document.getElementById('trajectory').textContent);
const FRAMES = DATA.frames;
const cache = new Array(FRAMES.length).fill(null);

function frameData(i) {
  if (!cache[i]) cache[i] = decodeB64(FRAMES[i].b64);
  return cache[i];
}

let current = 0;
let pending = 0;

async function showFrame(i) {
  current = i;
  const token = ++pending;
  const f = FRAMES[i];
  applyFrame(frameData(i), f.n);
  if (token !== pending) return;                 // a newer request won
  slider.value = i;
  label_el.textContent = f.label;
  count.textContent = f.n + (f.n === 1 ? ' particle' : ' particles');
}

// ------------------------------------------------------------------- UI
const slider = document.getElementById('slider');
const label_el = document.getElementById('label');
const count = document.getElementById('count');
const colorSelect = document.getElementById('color');
colorSelect.value = colorMode;
const playBtn = document.getElementById('play');

let playing = FRAMES.length > 1;
playBtn.textContent = playing ? '⏸' : '▶';

slider.addEventListener('input', () => {
  playing = false;
  playBtn.textContent = '▶';
  showFrame(+slider.value);
});
playBtn.addEventListener('click', () => {
  playing = !playing;
  playBtn.textContent = playing ? '⏸' : '▶';
});
colorSelect.addEventListener('change', (e) => {
  setColorMode(e.target.value);
  showFrame(current);
});
document.getElementById('box').addEventListener('change', (e) => {
  cell.visible = e.target.checked;
});
addEventListener('keydown', (e) => {
  if (e.code === 'Space') { e.preventDefault(); playBtn.click(); }
  if (e.code === 'ArrowRight') { slider.value = +slider.value + 1; slider.dispatchEvent(new Event('input')); }
  if (e.code === 'ArrowLeft') { slider.value = +slider.value - 1; slider.dispatchEvent(new Event('input')); }
});

let last = performance.now();
let acc = 0;
function animate(now) {
  requestAnimationFrame(animate);
  const dt = (now - last) / 1000;
  last = now;
  if (document.getElementById('spin').checked) spinCamera(dt);
  if (playing && FRAMES.length > 1) {
    acc += dt;
    if (acc >= 1 / DATA.fps) { acc = 0; showFrame((current + 1) % FRAMES.length); }
  }
  drawScene();
}
initVideoExport('export', 'vres', { frames: () => FRAMES, show: showFrame });
initFrameExport('savepng', 'vres', () => (FRAMES[current] ? FRAMES[current].label : 'frame'));
showFrame(0);
requestAnimationFrame(animate);
"""
def build(paths, out_path, title, fps=4.0, color="uniform"):
    times = read_times(os.path.dirname(os.path.abspath(paths[0])) or ".")

    frames = []
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    # colour-mode ranges, so the scale is stable across the whole trajectory
    rng = colour_ranges()

    for i, path in enumerate(paths):
        snap = read_snapshot(path)
        blobb64 = base64.b64encode(pack_snapshot(snap)).decode("ascii")
        if i < len(times):
            text = "%s   t = %g" % (os.path.basename(path), times[i])
        else:
            text = os.path.basename(path)
        frames.append({"label": text, "n": len(snap), "b64": blobb64})
        for a, b, c in snap.axes:
            note_ranges(rng, a, b, c)
        blo, bhi = bounding_box(snap)
        for k in range(3):
            lo[k] = min(lo[k], blo[k])
            hi[k] = max(hi[k], bhi[k])

    import json
    payload = json.dumps({"box": {"min": lo, "max": hi},
                          "ranges": rng, "fps": fps, "frames": frames})
    # </script> inside the JSON would terminate the host <script> element
    payload = payload.replace("</", "<\\/")

    html = (PAGE_HEAD
            .replace("__TITLE__", title)
            .replace("__CSS__", CSS)
            .replace("__LAST__", str(max(0, len(frames) - 1)))
            .replace("__BOXSIZE__", "&times;".join(
                "%g" % (hi[k] - lo[k]) for k in range(3)))
            .replace("__DATA__", payload)
            + importmap(THREE_VERSION)
            + '<script type="module">\n'
            + scene_js({"min": lo, "max": hi}, rng, color)
            + PLAYBACK_JS
            + PAGE_TAIL)

    with open(out_path, "w") as fh:
        fh.write(html)
    return frames, os.path.getsize(out_path)


def _main(argv=None):
    ap = argparse.ArgumentParser(
        description="Build a self-contained interactive HTML viewer for "
                    "ellipmd 'out*' snapshots.")
    ap.add_argument("files", nargs="+", help="snapshot files or globs")
    ap.add_argument("-o", "--out", default="trajectory.html",
                    help="output HTML file (default: trajectory.html)")
    ap.add_argument("--first", type=int, default=0,
                    help="index of the first snapshot to include")
    ap.add_argument("--last", type=int, default=None,
                    help="index of the last snapshot to include")
    ap.add_argument("--stride", type=int, default=1,
                    help="keep every Nth snapshot")
    ap.add_argument("--fps", type=float, default=4.0,
                    help="playback rate in frames per second (default 4)")
    ap.add_argument("--title", default=None,
                    help="title shown in the page (default: output file name)")
    ap.add_argument("--color", default="uniform",
                    choices=["uniform", "random", "a", "c", "aspect", "volume"],
                    help="initial colour mode (default: uniform); 'random' "
                         "gives every particle its own deterministic colour")
    args = ap.parse_args(argv)

    paths = expand_paths(args.files)
    if not paths:
        ap.error("no matching snapshot files")
    paths = paths[args.first:args.last:args.stride]
    if not paths:
        ap.error("the --first/--last/--stride selection is empty")

    title = args.title or ("ellipmd trajectory \u2014 %d snapshots" % len(paths))
    frames, size = build(paths, args.out, title, args.fps, args.color)

    print("wrote %s" % args.out)
    print("  %d snapshots, %d particles in the largest frame"
          % (len(frames), max(f["n"] for f in frames)))
    print("  %.1f MB" % (size / 1e6))
    if size > 60e6:
        print("  note: that is a large page; consider --stride or --first/--last")
    return 0


if __name__ == "__main__":
    sys.exit(_main())
