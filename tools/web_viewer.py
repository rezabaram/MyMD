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

PAGE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
  :root { color-scheme: dark; }
  * { box-sizing: border-box; }
  html, body { margin: 0; height: 100%; overflow: hidden;
               background: #14161a; color: #e8e8ea;
               font: 13px/1.4 ui-sans-serif, -apple-system, "Segoe UI", sans-serif; }
  #view { position: fixed; inset: 0; }
  #ui { position: fixed; left: 12px; bottom: 12px; right: 12px;
        display: flex; flex-wrap: wrap; align-items: center; gap: 10px;
        padding: 9px 12px; border-radius: 10px;
        background: rgba(24,26,31,.88); backdrop-filter: blur(8px);
        border: 1px solid rgba(255,255,255,.10); }
  #ui button { width: 34px; height: 30px; border-radius: 7px; cursor: pointer;
               border: 1px solid rgba(255,255,255,.16);
               background: rgba(255,255,255,.07); color: inherit; font-size: 13px; }
  #ui button:hover { background: rgba(255,255,255,.14); }
  #slider { flex: 1 1 220px; min-width: 140px; accent-color: #6aa9ff; }
  #ui select, #ui label { background: rgba(255,255,255,.07); color: inherit;
        border: 1px solid rgba(255,255,255,.16); border-radius: 7px;
        padding: 4px 7px; font: inherit; }
  #ui label { display: inline-flex; align-items: center; gap: 5px; cursor: pointer; }
  #ui label input { accent-color: #6aa9ff; margin: 0; }
  #label { font-variant-numeric: tabular-nums; min-width: 150px; opacity: .9; }
  #hud { position: fixed; left: 12px; top: 12px; padding: 8px 11px;
         border-radius: 9px; background: rgba(24,26,31,.8);
         border: 1px solid rgba(255,255,255,.10); line-height: 1.5;
         max-width: 340px; }
  #hud b { font-weight: 600; }
  #hud .dim { opacity: .6; }
  #err { position: fixed; inset: 0; display: none; place-items: center;
         padding: 40px; text-align: center; }
  #err code { background: rgba(255,255,255,.1); padding: 2px 5px; border-radius: 4px; }
</style>
</head>
<body>
<div id="view"></div>
<div id="hud">
  <b>__TITLE__</b><br>
  <span id="count"></span> &middot; box __BOXSIZE__<br>
  <span class="dim">drag = orbit &middot; shift+drag = pan &middot; wheel = zoom &middot; z is up</span>
</div>
<div id="ui">
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
</div>
<div id="err">
  <div>
    <p><b>Could not load three.js from the CDN.</b></p>
    <p>This page needs network access the first time it is opened.
       If you are offline, download <code>three.module.js</code> and the
       <code>examples/jsm</code> directory for version __THREE__ next to this
       file and edit the import map below.</p>
  </div>
</div>

<script type="importmap">
{
  "imports": {
    "three": "https://unpkg.com/three@__THREE__/build/three.module.js",
    "three/addons/": "https://unpkg.com/three@__THREE__/examples/jsm/"
  }
}
</script>

<script id="trajectory" type="application/json">__DATA__</script>

<script type="module">
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

const DATA = JSON.parse(document.getElementById('trajectory').textContent);
const FRAMES = DATA.frames;
const STRIDE = 10;                     // px,py,pz, qx,qy,qz,qw, a,b,c

// ---------------------------------------------------------------- decoding
const cache = new Array(FRAMES.length).fill(null);

async function decode(b64) {
  try {
    const res = await fetch('data:application/octet-stream;base64,' + b64);
    if (res.ok) return new Float32Array(await res.arrayBuffer());
  } catch (e) { /* fall through to atob */ }
  const bin = atob(b64);
  const u8 = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  return new Float32Array(u8.buffer);
}

function frameData(i) {
  if (!cache[i]) cache[i] = decode(FRAMES[i].b64);
  return cache[i];
}

// ------------------------------------------------------------- colormaps
const VIRIDIS = [
  [0.267, 0.005, 0.329], [0.283, 0.141, 0.458], [0.254, 0.265, 0.530],
  [0.207, 0.372, 0.553], [0.164, 0.471, 0.558], [0.128, 0.567, 0.551],
  [0.135, 0.659, 0.518], [0.267, 0.749, 0.441], [0.478, 0.821, 0.318],
  [0.741, 0.873, 0.150], [0.993, 0.906, 0.144],
];
function viridis(out, t) {
  t = Math.min(1, Math.max(0, t));
  const x = t * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x));
  const f = x - i, a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return out.setRGB(a[0] + f * (b[0] - a[0]),
                    a[1] + f * (b[1] - a[1]),
                    a[2] + f * (b[2] - a[2]));
}

// ------------------------------------------------------------------ scene
const view = document.getElementById('view');
const renderer = new THREE.WebGLRenderer({ antialias: true });
renderer.setPixelRatio(Math.min(devicePixelRatio, 2));
view.appendChild(renderer.domElement);

const scene = new THREE.Scene();
scene.background = new THREE.Color(0x14161a);

const boxMin = new THREE.Vector3(...DATA.box.min);
const boxMax = new THREE.Vector3(...DATA.box.max);
const boxCentre = boxMin.clone().add(boxMax).multiplyScalar(0.5);
const boxSize = boxMax.clone().sub(boxMin);
const diagonal = boxSize.length();

const camera = new THREE.PerspectiveCamera(38, 1, diagonal / 500, diagonal * 40);
camera.up.set(0, 0, 1);                       // gravity is along -z here
camera.position.copy(boxCentre).add(
  new THREE.Vector3(1.15, -1.75, 1.0).normalize().multiplyScalar(diagonal * 0.95));

const controls = new OrbitControls(camera, renderer.domElement);
controls.target.copy(boxCentre);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

scene.add(new THREE.HemisphereLight(0xdfe8ff, 0x2a2c33, 1.6));
const key = new THREE.DirectionalLight(0xffffff, 2.6);
key.position.set(1, -1.3, 2).multiplyScalar(diagonal);
scene.add(key);
const fill = new THREE.DirectionalLight(0x9fb8ff, 1.1);
fill.position.set(-1.4, 1.0, -0.6).multiplyScalar(diagonal);
scene.add(fill);

// simulation cell
const cellBox = new THREE.Box3(boxMin, boxMax);
const cell = new THREE.LineSegments(
  new THREE.EdgesGeometry(new THREE.BoxGeometry(
    boxSize.x || 1, boxSize.y || 1, boxSize.z || 1)),
  new THREE.LineBasicMaterial({ color: 0x6f7683 }));
cell.position.copy(boxCentre);
scene.add(cell);

// -------------------------------------------------------------- particles
const maxN = FRAMES.reduce((m, f) => Math.max(m, f.n), 1);
const geometry = new THREE.SphereGeometry(1, 24, 16);
const material = new THREE.MeshStandardMaterial({
  roughness: 0.45, metalness: 0.12, flatShading: false,
});
const mesh = new THREE.InstancedMesh(geometry, material, maxN);
mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
mesh.frustumCulled = false;
mesh.count = 0;
scene.add(mesh);

const _m = new THREE.Matrix4();
const _p = new THREE.Vector3();
const _q = new THREE.Quaternion();
const _s = new THREE.Vector3();
const _c = new THREE.Color();
const plain = new THREE.Color(0x7fb0e8);

// --- per-frame scalar ranges, so the colour scale does not jump around
const ranges = {};
for (const k of ['a', 'c', 'aspect', 'volume']) ranges[k] = [Infinity, -Infinity];
for (const f of FRAMES) {
  for (const p of f.stats) {
    for (const k of Object.keys(ranges)) {
      ranges[k][0] = Math.min(ranges[k][0], p[k]);
      ranges[k][1] = Math.max(ranges[k][1], p[k]);
    }
  }
}

let colorMode = '__COLORMODE__';

// Deterministic pseudo-random colour per particle, matching particle_color()
// in tools/ellipmd_io.py (splitmix32) so the same particle gets the same
// colour here and in OVITO.  It is keyed on the particle index, which the
// solver assigns in insertion order, so it is stable across frames.
function splitmix32(i) {
  let x = (i + 0x9E3779B9) >>> 0;
  x = Math.imul(x ^ (x >>> 16), 0x21F0AAAD) >>> 0;
  x = Math.imul(x ^ (x >>> 15), 0x735A2D97) >>> 0;
  return (x ^ (x >>> 15)) >>> 0;
}

function randomColour(out, k) {
  const x = splitmix32(k);
  const h = (x & 0xFFFF) / 65536.0;
  const s = 0.55 + ((x >>> 16) & 0xFF) / 255.0 * 0.35;
  const l = 0.42 + ((x >>> 24) & 0xFF) / 255.0 * 0.22;
  return out.setHSL(h, s, l, THREE.SRGBColorSpace);
}

function colourFor(out, k, a, b, c) {
  if (colorMode === 'uniform') return out.copy(plain);
  if (colorMode === 'random') return randomColour(out, k);
  let v;
  if (colorMode === 'a') v = a;
  else if (colorMode === 'c') v = c;
  else if (colorMode === 'aspect') v = a ? c / a : 0;
  else v = 4 / 3 * Math.PI * a * b * c;
  const [lo, hi] = ranges[colorMode];
  return viridis(out, hi > lo ? (v - lo) / (hi - lo) : 0.5);
}

let current = 0;
let pending = 0;

async function showFrame(i) {
  current = i;
  pending++;
  const token = pending;
  const arr = await frameData(i);
  if (token !== pending) return;                 // a newer request won

  const f = FRAMES[i];
  mesh.count = f.n;
  for (let k = 0; k < f.n; k++) {
    const o = k * STRIDE;
    _p.set(arr[o], arr[o + 1], arr[o + 2]);
    _q.set(arr[o + 3], arr[o + 4], arr[o + 5], arr[o + 6]);
    _s.set(arr[o + 7], arr[o + 8], arr[o + 9]);
    _m.compose(_p, _q, _s);
    mesh.setMatrixAt(k, _m);
    mesh.setColorAt(k, colourFor(_c, k, _s.x, _s.y, _s.z));
  }
  mesh.instanceMatrix.needsUpdate = true;
  if (mesh.instanceColor) mesh.instanceColor.needsUpdate = true;

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
playBtn.textContent = playing ? '\u23F8' : '\u25B6';

slider.addEventListener('input', () => {
  playing = false;
  playBtn.textContent = '\u25B6';
  showFrame(+slider.value);
});
playBtn.addEventListener('click', () => {
  playing = !playing;
  playBtn.textContent = playing ? '\u23F8' : '\u25B6';
});
document.getElementById('color').addEventListener('change', (e) => {
  colorMode = e.target.value;
  showFrame(current);
});document.getElementById('box').addEventListener('change', (e) => {
  cell.visible = e.target.checked;
});
addEventListener('keydown', (e) => {
  if (e.code === 'Space') { e.preventDefault(); playBtn.click(); }
  if (e.code === 'ArrowRight') { slider.value = +slider.value + 1; slider.dispatchEvent(new Event('input')); }
  if (e.code === 'ArrowLeft') { slider.value = +slider.value - 1; slider.dispatchEvent(new Event('input')); }
});

function resize() {
  const w = innerWidth, h = innerHeight;
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
addEventListener('resize', resize);
resize();

let last = performance.now();
let acc = 0;

function animate(now) {
  requestAnimationFrame(animate);
  const dt = (now - last) / 1000;
  last = now;

  if (document.getElementById('spin').checked) {
    const s = 0.35 * dt;
    const p = camera.position.clone().sub(controls.target);
    p.applyAxisAngle(new THREE.Vector3(0, 0, 1), s);
    camera.position.copy(controls.target).add(p);
  }

  if (playing && FRAMES.length > 1) {
    acc += dt;
    const dwell = 1 / DATA.fps;
    if (acc >= dwell) {
      acc = 0;
      showFrame((current + 1) % FRAMES.length);
    }
  }

  controls.update();
  renderer.render(scene, camera);
}

showFrame(0).then(() => requestAnimationFrame(animate));
</script>
</body>
</html>
"""


def _pack(snapshot):
    """(base64 float32 blob, per-frame stats) for one snapshot.

    Layout per particle (little-endian float32):
        px py pz  qx qy qz qw  a b c
    The quaternion is reordered from the file's scalar-first (w,x,y,z) to the
    (x,y,z,w) order that three.js uses."""
    n = len(snapshot)
    buf = bytearray(n * 10 * 4)
    stats = []
    for i in range(n):
        x, y, z = snapshot.positions[i]
        w, qx, qy, qz = snapshot.quats[i]
        a, b, c = snapshot.axes[i]
        struct.pack_into("<10f", buf, i * 40,
                         x, y, z, qx, qy, qz, w, a, b, c)
        stats.append({"a": a, "b": b, "c": c,
                      "aspect": (c / a) if a else 0.0,
                      "volume": 4.0 / 3.0 * 3.141592653589793 * a * b * c})
    return base64.b64encode(bytes(buf)).decode("ascii"), stats


def build(paths, out_path, title, fps=4.0, color="uniform"):
    times = read_times(os.path.dirname(os.path.abspath(paths[0])) or ".")

    frames = []
    lo = [float("inf")] * 3
    hi = [float("-inf")] * 3
    for i, path in enumerate(paths):
        snap = read_snapshot(path)
        blobb64, stats = _pack(snap)
        if i < len(times):
            text = "%s   t = %g" % (os.path.basename(path), times[i])
        else:
            text = os.path.basename(path)
        frames.append({"label": text, "n": len(snap),
                       "b64": blobb64, "stats": stats})
        blo, bhi = bounding_box(snap)
        for k in range(3):
            lo[k] = min(lo[k], blo[k])
            hi[k] = max(hi[k], bhi[k])

    import json
    payload = json.dumps({"box": {"min": lo, "max": hi},
                          "fps": fps, "frames": frames})
    # </script> inside the JSON would terminate the host <script> element
    payload = payload.replace("</", "<\\/")

    html = (PAGE
            .replace("__TITLE__", title)
            .replace("__THREE__", THREE_VERSION)
            .replace("__LAST__", str(max(0, len(frames) - 1)))
            .replace("__COLORMODE__", color)
            .replace("__BOXSIZE__", "&times;".join(
                "%g" % (hi[k] - lo[k]) for k in range(3)))
            .replace("__DATA__", payload))

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
