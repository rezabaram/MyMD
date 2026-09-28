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
    <p><b>The viewer did not start.</b></p>
    <p id="errmsg" class="dim"></p>
    <p>Usually this means three.js could not be fetched from the CDN, which the
       page needs the first time it is opened.  If you are offline, download
       <code>three.module.js</code> and the <code>OrbitControls</code> addon,
       put them next to this file and edit the import map.</p>
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
        "  function showErr(msg) {\n"
        "    const el = document.getElementById('err');\n"
        "    if (!el) return;\n"
        "    const m = document.getElementById('errmsg');\n"
        "    if (m && msg) m.textContent = String(msg);\n"
        "    el.style.display = 'grid';\n"
        "  }\n"
        "  addEventListener('error', (e) => showErr(e.message || e.error || ''), true);\n"
        "  addEventListener('unhandledrejection', (e) => showErr(e.reason || ''));\n"
        "  setTimeout(() => {\n"
        "    if (!window.__mymd_started) showErr('the module did not finish loading');\n"
        "  }, 4000);\n"
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
// preserveDrawingBuffer keeps the rendered frame readable after the
// compositor has taken it, which canvas.toBlob() during a video export needs.
const renderer = new THREE.WebGLRenderer({ antialias: true,
                                           preserveDrawingBuffer: true });
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

const controls = new OrbitControls(camera, renderer.domElement);
controls.enableDamping = true;
controls.dampingFactor = 0.08;

// Framing the camera touches controls.target, so it has to be defined after
// controls -- calling it earlier is a temporal-dead-zone ReferenceError that
// takes the whole module down, which shows up as a black canvas and a Start
// button that does nothing.  It is called at the end of this file, once the
// renderer has a real aspect ratio to fit against.
//
// The box is fitted by its bounding sphere.  The old code put the camera at
// 0.95 * the box diagonal, which framed roughly the middle 60% of the box:
// fitting the sphere needs radius / sin(fov/2), about 1.6x further out for
// these proportions.  camera.fov is the *vertical* field of view, so a window
// wider than it is tall is limited by height and a tall one by width.
function fitDistance() {
  const radius = Math.max(boxSize.length() * 0.5, 1e-6);
  const vFov = camera.fov * Math.PI / 180;
  const hFov = 2 * Math.atan(Math.tan(vFov / 2) * Math.max(camera.aspect, 1e-6));
  return Math.max(radius / Math.sin(vFov / 2),
                  radius / Math.sin(hFov / 2)) * 1.06;
}

function homeCamera() {
  const dist = fitDistance();
  camera.position.copy(boxCentre).add(
    new THREE.Vector3(1.15, -1.75, 1.0).normalize().multiplyScalar(dist));
  camera.near = Math.max(dist / 200, 1e-6);
  camera.far = dist * 200;
  camera.updateProjectionMatrix();
  controls.target.copy(boxCentre);
  controls.update();
}

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

// Set while a video is being recorded, so the wall-clock spin in animate()
// does not fight the fixed per-frame spin the export applies.
let exporting = false;
/// True when the page is served by live_viewer and ffmpeg is available there.
let serverEncode = false;

function spinCamera(dt) {
  if (exporting) return;
  const s = 0.35 * dt;
  const p = camera.position.clone().sub(controls.target);
  p.applyAxisAngle(new THREE.Vector3(0, 0, 1), s);
  camera.position.copy(controls.target).add(p);
}



// Switching the backing store to the export resolution while leaving the CSS
// size alone stretches the on-screen image -- the canvas keeps its displayed
// size and its contents are rescaled, which changes the apparent aspect ratio
// for the whole of the recording.  Fit the canvas into its container at the
// export aspect instead, letterboxed, and put it back afterwards.
function fitCanvasForExport(width, height) {
  const canvas = renderer.domElement;
  const prev = { w: canvas.style.width, h: canvas.style.height,
                 pos: canvas.style.position, left: canvas.style.left,
                 top: canvas.style.top };
  const host = view.getBoundingClientRect();
  if (host.width > 0 && host.height > 0) {
    const scale = Math.min(host.width / width, host.height / height);
    canvas.style.position = 'absolute';
    canvas.style.width = Math.round(width * scale) + 'px';
    canvas.style.height = Math.round(height * scale) + 'px';
    canvas.style.left = Math.round((host.width - width * scale) / 2) + 'px';
    canvas.style.top = Math.round((host.height - height * scale) / 2) + 'px';
  }
  return prev;
}

function unfitCanvasAfterExport(prev) {
  const canvas = renderer.domElement;
  canvas.style.width = prev.w;
  canvas.style.height = prev.h;
  canvas.style.position = prev.pos;
  canvas.style.left = prev.left;
  canvas.style.top = prev.top;
  resizeRenderer();            // backing store back to the window, and CSS with it
}

// ------------------------------------------------------------ video export
//
// Records the 3D canvas to a video file.  Only the canvas is captured, so none
// of the panels, buttons or sliders appear -- they are separate DOM elements,
// not part of the drawing surface.  Whatever the page is drawing *is* in the
// video: the simulation box is drawn when the "box" toggle is on, and the
// camera spin is applied one fixed step per exported frame when "spin" is on,
// so the rotation is smooth and reproducible rather than tied to however long
// each frame took to render.
//
// Recording is done by the browser's own MediaRecorder against a stream taken
// from the canvas, so nothing is uploaded anywhere.  MP4 is preferred; where
// the browser cannot produce it (Chrome has historically offered only WebM)
// the file is written as WebM and named .webm rather than being called an mp4
// it is not.
//
// The page supplies where the frames come from, because the two pages keep
// them differently -- the static viewer has a fixed FRAMES array, the live one
// a list that grows as the run goes.


// Encodes through the server's ffmpeg instead of the browser's MediaRecorder.
//
// GitHub plays H.264 in an MP4 container and nothing else inline, and a
// browser will not necessarily produce that -- Chrome offers WebM.  When the
// page is served by tools/live_viewer.py and ffmpeg is on PATH, the frames are
// rendered to PNG here and piped to ffmpeg there, which gives a file GitHub
// will play, with yuv420p for compatibility and +faststart so it starts
// playing before it has finished downloading.
async function exportMp4ViaServer(opts, source) {
  const list = source.frames();
  if (!list || !list.length) {
    alert('There are no frames to record yet.');
    return false;
  }
  const canvas = renderer.domElement;
  const size = new THREE.Vector2();
  renderer.getSize(size);
  const oldAspect = camera.aspect;
  const oldCellVisible = cell.visible;
  cell.visible = true;
  renderer.setSize(opts.width, opts.height, false);
  const prevLayout = fitCanvasForExport(opts.width, opts.height);
  camera.aspect = opts.width / opts.height;
  camera.updateProjectionMatrix();
  exporting = true;

  let token = null;
  try {
    const r0 = await fetch('/api/encode/start', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ fps: opts.fps, crf: opts.crf }),
    });
    if (!r0.ok) {
      throw new Error('the server has no encoder (HTTP ' + r0.status + ' from '
        + '/api/encode/start).  If "make live" was started before the video '
        + 'export was added, restart it.');
    }
    let st;
    try {
      st = await r0.json();
    } catch (e) {
      throw new Error('the server did not return JSON from /api/encode/start '
        + '(HTTP ' + r0.status + ')');
    }
    token = st.token;
    if (!token) throw new Error(st.error || 'could not start the encoder');

    for (let i = 0; i < list.length; i++) {
      await source.show(i);
      spinCamera(1.0 / opts.fps);      // one fixed step, as in the browser path
      drawScene();
      const blob = await new Promise((res) => canvas.toBlob(res, 'image/png'));
      if (!blob) throw new Error('could not read the canvas');
      const r = await fetch('/api/encode/frame?token=' + encodeURIComponent(token),
                            { method: 'POST', body: blob });
      if (!r.ok) {
        throw new Error('HTTP ' + r.status + ' from /api/encode/frame: '
                        + (await r.text()));
      }
      if (opts.onProgress) opts.onProgress(i + 1, list.length);
    }

    const resp = await fetch('/api/encode/finish?token=' + encodeURIComponent(token),
                             { method: 'POST' });
    if (!resp.ok) {
      throw new Error('HTTP ' + resp.status + ' from /api/encode/finish: '
                      + (await resp.text()));
    }
    token = null;
    const video = await resp.blob();
    const url = URL.createObjectURL(video);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'ellipmd-' + list.length + 'frames.mp4';
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    setTimeout(() => URL.revokeObjectURL(url), 10000);
    return video.size;
  } catch (e) {
    if (token) {
      try { await fetch('/api/encode/abort?token=' + encodeURIComponent(token),
                        { method: 'POST' }); } catch (e2) { }
    }
    alert('mp4 export failed: ' + (e && e.message ? e.message : e));
    return false;
  } finally {
    exporting = false;
    unfitCanvasAfterExport(prevLayout);
    camera.aspect = oldAspect;
    camera.updateProjectionMatrix();
    cell.visible = oldCellVisible;
    if (typeof userMovedCamera !== 'undefined' && !userMovedCamera) homeCamera();
  }
}

function pickVideoMime() {
  if (typeof MediaRecorder === 'undefined') return null;
  const candidates = [
    'video/mp4;codecs=avc1.4d002a',     // mp4, main profile
    'video/mp4;codecs=avc1.42E01E',     // mp4, baseline
    'video/mp4',
    'video/webm;codecs=vp9',
    'video/webm;codecs=vp8',
    'video/webm',
  ];
  for (const m of candidates) {
    try { if (MediaRecorder.isTypeSupported(m)) return m; } catch (e) { }
  }
  return null;
}

function videoFormatLabel() {
  const m = pickVideoMime();
  if (!m) return 'unavailable';
  return m.indexOf('mp4') >= 0 ? 'mp4' : 'webm';
}

const waitMs = (ms) => new Promise((r) => setTimeout(r, ms));

// `source` = { frames: () => array, show: (i) => Promise }
async function exportVideo(opts, source) {
  opts = Object.assign({ fps: 30, width: 1280, height: 720,
                         bitrate: 12000000, onProgress: null }, opts || {});
  const mime = pickVideoMime();
  if (!mime) {
    alert('This browser cannot record video: MediaRecorder is unavailable.');
    return false;
  }
  const list = source.frames();
  if (!list || !list.length) {
    alert('There are no frames to record yet.');
    return false;
  }

  const canvas = renderer.domElement;

  // Record at a fixed size rather than whatever the window happens to be.
  const size = new THREE.Vector2();
  renderer.getSize(size);
  const oldAspect = camera.aspect;
  // The simulation box and the camera spin are part of the render, not view
  // furniture, so the video always has both regardless of the toggles.
  const oldCellVisible = cell.visible;
  cell.visible = true;
  renderer.setSize(opts.width, opts.height, false);
  const prevLayout = fitCanvasForExport(opts.width, opts.height);
  camera.aspect = opts.width / opts.height;
  camera.updateProjectionMatrix();

  const stream = canvas.captureStream(0);          // 0 = only on requestFrame
  const track = stream.getVideoTracks()[0];
  const chunks = [];
  const rec = new MediaRecorder(stream,
                                { mimeType: mime, videoBitsPerSecond: opts.bitrate });
  rec.ondataavailable = (e) => { if (e.data && e.data.size) chunks.push(e.data); };
  const stopped = new Promise((res) => { rec.onstop = res; });

  const period = 1000 / opts.fps;
  exporting = true;                                // animate() stops spinning
  rec.start();
  try {
    for (let i = 0; i < list.length; i++) {
      const t0 = performance.now();
      await source.show(i);
      // One fixed step of rotation per frame: the wall-clock spin in animate()
      // would give a different angle per frame depending on render time.
      spinCamera(1.0 / opts.fps);
      drawScene();                       // applyFrame only uploads; this draws
      if (track.requestFrame) track.requestFrame();
      if (opts.onProgress) opts.onProgress(i + 1, list.length);
      // MediaRecorder timestamps by wall clock, so pace the frames or a slow
      // render makes the video play at the wrong speed.
      const spent = performance.now() - t0;
      if (spent < period) await waitMs(period - spent);
    }
  } finally {
    exporting = false;
    rec.stop();
    await stopped;
    unfitCanvasAfterExport(prevLayout);
    camera.aspect = oldAspect;
    camera.updateProjectionMatrix();
    cell.visible = oldCellVisible;
    if (typeof userMovedCamera !== 'undefined' && !userMovedCamera) homeCamera();
  }

  const ext = mime.indexOf('mp4') >= 0 ? 'mp4' : 'webm';
  const blob = new Blob(chunks, { type: mime.split(';')[0] });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url;
  a.download = 'ellipmd-' + list.length + 'frames.' + ext;
  document.body.appendChild(a);
  a.click();
  document.body.removeChild(a);
  setTimeout(() => URL.revokeObjectURL(url), 10000);
  return true;
}

/// Wire up an export button and a resolution select, if the page has them.
function initVideoExport(btnId, resId, source) {
  const btn = document.getElementById(btnId);
  if (!btn) return;
  const res = resId ? document.getElementById(resId) : null;
  let idle = 'record';

  function refresh() {
    if (serverEncode) {
      btn.textContent = idle;
      btn.title = 'record the 3D view (without the interface) to an H.264 MP4 '
                + 'via ffmpeg -- the format GitHub plays inline';
      btn.disabled = false;
    } else {
      const fmt = videoFormatLabel();
      btn.textContent = fmt === 'unavailable' ? 'no video' : 'record ' + fmt;
      btn.title = fmt === 'unavailable'
        ? 'this browser cannot record video and the page has no server-side encoder'
        : 'record the 3D view (without the interface) to a ' + fmt
          + ' file.  GitHub plays MP4 inline; WebM it will not.';
      btn.disabled = (fmt === 'unavailable');
    }
  }
  refresh();
  initVideoExport.refresh = refresh;

  btn.addEventListener('click', async () => {
    let w = 1280, h = 720;
    if (res && res.value) { const p = res.value.split('x'); w = +p[0]; h = +p[1]; }
    btn.disabled = true;
    const progress = (i, n) => { btn.textContent = Math.round(100 * i / n) + '%'; };
    let size = 0;
    if (serverEncode) {
      size = await exportMp4ViaServer({ width: w, height: h, fps: 30, crf: 23,
                                        onProgress: progress }, source);
    } else {
      const ok = await exportVideo({ width: w, height: h, fps: 30,
                                     onProgress: progress }, source);
      size = ok ? 1 : 0;
    }
    refresh();
    if (size) {
      btn.textContent = size > 1 ? (size / 1e6).toFixed(1) + ' MB' : idle;
      setTimeout(refresh, 4000);
    }
  });
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
resizeRenderer();
homeCamera();

// Re-fit when the window changes shape, but only until the user has orbited:
// after that the camera is theirs and snapping it back would be rude.
let userMovedCamera = false;
controls.addEventListener('start', () => { userMovedCamera = true; });
addEventListener('resize', () => {
  resizeRenderer();
  if (!userMovedCamera) homeCamera();
});
window.__mymd_started = true;
""".replace("__BOX__", json.dumps(box)) \
   .replace("__RANGES__", json.dumps(ranges)) \
   .replace("__COLORMODE__", color_mode) \
   .replace("__EXTRA__", extra_globals)
