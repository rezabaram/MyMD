"""Shared pieces of the interactive three.js viewers.

Two pages render the same scene: ``web_viewer.py`` builds a self-contained HTML
file with a whole trajectory embedded in it, and ``live_viewer.py`` serves a
dashboard over HTTP that follows a run as it is produced.  The scene, the
colour maps and the per-frame update live here so the two cannot drift apart.

The JavaScript is emitted as one string meant to be concatenated with a
page-specific script *inside a single* ``<script type="module">``.  ES modules
cannot be split across script elements, so the shared part is a prefix, not a
separate file.
"""

from __future__ import annotations

import json
import struct

# Little-endian float32 per particle: px py pz  qx qy qz qw  a b c
# The quaternion is reordered from the file's scalar-first (w,x,y,z) to the
# (x,y,z,w) order three.js uses.
STRIDE = 10


def pack_snapshot(snapshot):
    """The float32 blob for one snapshot, as raw bytes."""
    buf = bytearray(len(snapshot) * STRIDE * 4)
    for i in range(len(snapshot)):
        x, y, z = snapshot.positions[i]
        w, qx, qy, qz = snapshot.quats[i]
        a, b, c = snapshot.axes[i]
        struct.pack_into("<10f", buf, i * 40, x, y, z, qx, qy, qz, w, a, b, c)
    return bytes(buf)


def colour_ranges():
    return {k: [float("inf"), float("-inf")]
            for k in ("a", "c", "aspect", "volume")}


def note_ranges(rng, a, b, c):
    """Widen the colour-mode ranges with one particle's semi-axes."""
    import math
    vals = {"a": a, "c": c,
            "aspect": (c / a) if a else 0.0,
            "volume": 4.0 / 3.0 * math.pi * a * b * c}
    for k, v in vals.items():
        if v < rng[k][0]:
            rng[k][0] = v
        if v > rng[k][1]:
            rng[k][1] = v

THREE_VERSION = "0.160.0"

CSS = r"""
  :root { color-scheme: dark; }
  * { box-sizing: border-box; }
  html, body { margin: 0; height: 100%; overflow: hidden;
               background: #14161a; color: #e8e8ea;
               font: 13px/1.4 ui-sans-serif, -apple-system, "Segoe UI", sans-serif; }
  #view { position: fixed; inset: 0; }
  .panel { position: fixed; padding: 9px 12px; border-radius: 10px;
           background: rgba(24,26,31,.9); backdrop-filter: blur(8px);
           border: 1px solid rgba(255,255,255,.10); }
  #ui { left: 12px; bottom: 12px; right: 12px;
        display: flex; flex-wrap: wrap; align-items: center; gap: 10px; }
  button { height: 30px; padding: 0 10px; border-radius: 7px; cursor: pointer;
           border: 1px solid rgba(255,255,255,.16);
           background: rgba(255,255,255,.07); color: inherit; font: inherit; }
  button:hover { background: rgba(255,255,255,.14); }
  button:disabled { opacity: .4; cursor: default; }
  #play { width: 38px; padding: 0; }
  #slider { flex: 1 1 200px; min-width: 120px; accent-color: #6aa9ff; }
  select, input[type=text], input[type=number], label {
        background: rgba(255,255,255,.07); color: inherit;
        border: 1px solid rgba(255,255,255,.16); border-radius: 7px;
        padding: 4px 7px; font: inherit; }
  input[type=text], input[type=number] { width: 100%; }
  label { display: inline-flex; align-items: center; gap: 5px; cursor: pointer;
          border: none; background: none; }
  label input { accent-color: #6aa9ff; margin: 0; }
  #label { font-variant-numeric: tabular-nums; min-width: 140px; opacity: .9; }
  #hud { left: 12px; top: 12px; line-height: 1.5; max-width: 300px; }
  #hud b { font-weight: 600; }
  .dim { opacity: .6; }
  #err { position: fixed; inset: 0; display: none; place-items: center;
         padding: 40px; text-align: center; }
  code { background: rgba(255,255,255,.1); padding: 2px 5px; border-radius: 4px; }
"""

ERROR_BLOCK = r"""
<div id="err">
  <div>
    <p><b>Could not load three.js from the CDN.</b></p>
    <p>This page needs network access the first time it is opened.
       If you are offline, download <code>three.module.js</code> and the
       <code>OrbitControls</code> addon, put them next to this file and edit
       the import map below.</p>
  </div>
</div>
"""


def importmap(three_version=THREE_VERSION):
    return (
        '<script type="importmap">\n'
        "{\n"
        '  "imports": {\n'
        '    "three": "https://unpkg.com/three@%s/build/three.module.js",\n'
        '    "three/addons/": "https://unpkg.com/three@%s/examples/jsm/"\n'
        "  }\n"
        "}\n"
        "</script>\n"
        '<script>\n'
        "  addEventListener('error', (e) => {\n"
        "    if (String(e.message || '').includes('three') ||\n"
        "        String((e.filename || '')).includes('unpkg'))\n"
        "      document.getElementById('err').style.display = 'grid';\n"
        "  }, true);\n"
        "</script>\n" % (three_version, three_version)
    )


def scene_js(box, ranges, color_mode="uniform", max_n=1, extra_globals=""):
    """The shared scene: renderer, camera, lights, cell and particle mesh.

    Defines, for the page-specific script that follows:
        applyFrame(blob, n, label)  -- upload one frame
        setColorMode(mode), colorMode
        ensureCapacity(n)
        spinCamera(dt), drawScene()
        cell, camera, controls, renderer, scene
    """
    return r"""
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';

let BOX = __BOX__;
let ranges = __RANGES__;
const STRIDE = 10;                     // px,py,pz, qx,qy,qz,qw, a,b,c
let colorMode = '__COLORMODE__';
__EXTRA__
// ---------------------------------------------------------------- decoding
// Little-endian float32, STRIDE floats per particle, base64 from Python.
function decodeB64(b64) {
  const bin = atob(b64);
  const u8 = new Uint8Array(bin.length);
  for (let i = 0; i < bin.length; i++) u8[i] = bin.charCodeAt(i);
  return new Float32Array(u8.buffer);
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

const boxMin = new THREE.Vector3(...BOX.min);
const boxMax = new THREE.Vector3(...BOX.max);
const boxCentre = new THREE.Vector3();
const boxSize = new THREE.Vector3();
let diagonal = 1;
function recomputeBox() {
  boxMin.set(...BOX.min);
  boxMax.set(...BOX.max);
  boxCentre.copy(boxMin).add(boxMax).multiplyScalar(0.5);
  boxSize.copy(boxMax).sub(boxMin);
  diagonal = Math.max(boxSize.length(), 1e-6);
}

recomputeBox();
const camera = new THREE.PerspectiveCamera(38, 1, 0.01, 1000);
camera.up.set(0, 0, 1);                       // gravity is along -z here
function homeCamera() {
  camera.position.copy(boxCentre).add(
    new THREE.Vector3(1.15, -1.75, 1.0).normalize().multiplyScalar(diagonal * 0.95));
  camera.near = diagonal / 500;
  camera.far = diagonal * 40;
  camera.updateProjectionMatrix();
  controls.target.copy(boxCentre);
}
homeCamera();

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

scene.add(new THREE.HemisphereLight(0xdfe8ff, 0x2a2c33, 1.6));
const key = new THREE.DirectionalLight(0xffffff, 2.6);
key.position.set(1, -1.3, 2).multiplyScalar(diagonal);
scene.add(key);
const fill = new THREE.DirectionalLight(0x9fb8ff, 1.1);
fill.position.set(-1.4, 1.0, -0.6).multiplyScalar(diagonal);
scene.add(fill);

// simulation cell -- rebuilt when the box changes, which it does in the live
// viewer when a run with different boxsize is started
const cellMaterial = new THREE.LineBasicMaterial({ color: 0x6f7683 });
let cell = null;
function buildCell() {
  if (cell) { scene.remove(cell); cell.geometry.dispose(); }
  cell = new THREE.LineSegments(
    new THREE.EdgesGeometry(new THREE.BoxGeometry(
      boxSize.x || 1, boxSize.y || 1, boxSize.z || 1)), cellMaterial);
  cell.position.copy(boxCentre);
  scene.add(cell);
}
buildCell();

// Adopt a new bounding box (the live viewer learns it from the config it is
// about to run).  The camera is only reset the first time, so a box change
// does not throw away wherever the user has orbited to.
function setBox(b) {
  const first = !boxReady;
  BOX = b;
  recomputeBox();
  buildCell();
  if (first) { boxReady = true; homeCamera(); }
}
let boxReady = false;

// -------------------------------------------------------------- particles
const geometry = new THREE.SphereGeometry(1, 24, 16);
const material = new THREE.MeshStandardMaterial({
  roughness: 0.45, metalness: 0.12, flatShading: false,
});
let mesh = null;
let capacity = 0;

// The live viewer does not know how many particles a run will end with, so
// the instanced mesh is grown on demand rather than sized up front.
function ensureCapacity(n) {
  if (mesh && n <= capacity) return;
  const want = Math.max(n, Math.ceil(capacity * 1.5), 256);
  if (mesh) { scene.remove(mesh); mesh.dispose(); }
  mesh = new THREE.InstancedMesh(geometry, material, want);
  mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);
  mesh.frustumCulled = false;
  mesh.count = 0;
  scene.add(mesh);
  capacity = want;
}

const _m = new THREE.Matrix4();
const _p = new THREE.Vector3();
const _q = new THREE.Quaternion();
const _s = new THREE.Vector3();
const _c = new THREE.Color();
const plain = new THREE.Color(0x7fb0e8);

// Deterministic pseudo-random colour per particle, matching particle_color()
// in tools/ellipmd_io.py (splitmix32) so the same particle gets the same
// colour here and in OVITO.  Keyed on the particle index, which the solver
// assigns in insertion order, so it is stable across frames.
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
  const r = ranges[colorMode] || [0, 1];
  return viridis(out, r[1] > r[0] ? (v - r[0]) / (r[1] - r[0]) : 0.5);
}
function setColorMode(m) { colorMode = m; }

// Upload one frame.  `arr` is a Float32Array of n*STRIDE values.
function applyFrame(arr, n) {
  ensureCapacity(n);
  mesh.count = n;
  for (let k = 0; k < n; k++) {
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
}

function spinCamera(dt) {
  const s = 0.35 * dt;
  const p = camera.position.clone().sub(controls.target);
  p.applyAxisAngle(new THREE.Vector3(0, 0, 1), s);
  camera.position.copy(controls.target).add(p);
}

function resizeRenderer() {
  const w = innerWidth, h = innerHeight;
  renderer.setSize(w, h);
  camera.aspect = w / h;
  camera.updateProjectionMatrix();
}
function drawScene() {
  controls.update();
  renderer.render(scene, camera);
}
addEventListener('resize', resizeRenderer);
resizeRenderer();
""".replace("__BOX__", json.dumps(box)) \
   .replace("__RANGES__", json.dumps(ranges)) \
   .replace("__COLORMODE__", color_mode) \
   .replace("__EXTRA__", extra_globals)
