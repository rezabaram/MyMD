#!/usr/bin/env python3
"""A local web dashboard for running ellipmd and watching it as it goes.

    python3 tools/live_viewer.py                # then open the printed URL
    python3 tools/live_viewer.py --port 8765 --config config_rain

The page has a parameter form, Start and Stop, a progress bar, live statistics
and a 3D view that picks up each new snapshot as the solver writes it.  Nothing
is installed and nothing leaves the machine: the server is the Python standard
library's http.server and the page is served from it.

Endpoints
    GET  /                 the dashboard
    GET  /api/status       everything the page polls for, as JSON
    GET  /api/frame/<name> one snapshot as a float32 blob (see viewer_common)
    POST /api/run          {"config": "...", "seed": 1} -> start a run
    POST /api/stop         terminate the run
    GET  /api/config       the config text currently loaded in the form
"""

from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ellipmd_io import read_snapshot  # noqa: E402
from viewer_common import (  # noqa: E402
    CSS, ERROR_BLOCK, THREE_VERSION, colour_ranges, importmap, note_ranges,
    pack_snapshot, scene_js,
)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ELF = os.path.join(ROOT, "ellipmd")

# How many log_energy rows to keep for the page.  At 20 frames/s and half a
# minute of run that is plenty for a sparkline and costs nothing.
ENERGY_KEEP = 600


def _config_params(text):
    """The handful of parameters the page needs to show progress against.

    Deliberately forgiving: this parses the text the user typed, before the
    solver has validated it, so anything unparseable is simply absent.
    """
    want = {
        "maxTime": float, "nParticle": int, "outDt": float, "timeStep": float,
        "particleSize": float, "strainRate": float, "scaling": float,
        "boxcorner": (lambda s: [float(v) for v in s.split()]),
        "boxsize": (lambda s: [float(v) for v in s.split()]),
        "method": str, "boundary": str, "rainRate": float,
    }
    out = {}
    for line in text.splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        key, value = parts[0], parts[1].strip()
        fn = want.get(key)
        if fn is None:
            continue
        try:
            out[key] = fn(value)
        except ValueError:
            pass
    return out


class Run:
    """One solver process and the directory it is producing."""

    def __init__(self, rundir):
        self.rundir = rundir
        self.proc = None
        self.seed = 0
        self.config_text = ""
        self.log = ""
        self.started = None
        self.exit_code = None
        self._lock = threading.Lock()

    # -- process -----------------------------------------------------------
    def start(self, config_text, seed):
        self.stop()
        if os.path.isdir(self.rundir):
            shutil.rmtree(self.rundir)
        os.makedirs(self.rundir)
        with open(os.path.join(self.rundir, "config"), "w") as fh:
            fh.write(config_text)

        self.config_text = config_text
        self.seed = seed
        self.exit_code = None
        self.log = ""
        self.started = time.time()
        with open(os.path.join(self.rundir, "solver.log"), "wb") as logfh:
            self.proc = subprocess.Popen(
                [ELF, str(seed), "config"], cwd=self.rundir,
                stdout=logfh, stderr=subprocess.STDOUT)

    def stop(self):
        p = self.proc
        if p is not None and p.poll() is None:
            p.terminate()
            try:
                p.wait(timeout=5)
            except subprocess.TimeoutExpired:
                p.kill()
        self.proc = None

    def running(self):
        return self.proc is not None and self.proc.poll() is None

    def note_exit(self):
        if self.proc is not None and self.exit_code is None:
            code = self.proc.poll()
            if code is not None:
                self.exit_code = code


class FrameCache:
    """Reads each snapshot once and remembers its time, size and float blob."""

    def __init__(self):
        self.frames = []          # [{"name","t","n"}]
        self.blobs = {}           # name -> bytes
        self.ranges = colour_ranges()
        self._seen = set()
        self.lock = threading.Lock()

    def scan(self, rundir):
        with self.lock:
            self._scan(rundir)

    def _scan(self, rundir):
        if not os.path.isdir(rundir):
            return
        names = []
        for entry in os.listdir(rundir):
            if entry == "outend":
                names.append((10 ** 9, entry))
            else:
                m = re.fullmatch(r"out(\d+)", entry)
                if m:
                    names.append((int(m.group(1)), entry))
        names.sort()

        for _, name in names:
            if name in self._seen:
                continue
            path = os.path.join(rundir, name)
            try:
                snap = read_snapshot(path)
            except OSError:
                continue
            t = _snapshot_time(path)
            self.blobs[name] = pack_snapshot(snap)
            for a, b, c in snap.axes:
                note_ranges(self.ranges, a, b, c)
            self.frames.append({"name": name, "t": t, "n": len(snap)})
            self._seen.add(name)

    def reset(self):
        with self.lock:
            self.frames = []
            self.blobs = {}
            self.ranges = colour_ranges()
            self._seen = set()

    def blob(self, name):
        with self.lock:
            return self.blobs.get(name)


def _snapshot_time(path):
    """The t= from a snapshot header, or 0 for a file written before them."""
    try:
        with open(path, "r", errors="replace") as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                m = re.search(r"t=([0-9.eE+-]+)", line)
                if m:
                    try:
                        return float(m.group(1))
                    except ValueError:
                        pass
    except OSError:
        pass
    return 0.0


def _energy_tail(rundir):
    path = os.path.join(rundir, "log_energy")
    try:
        with open(path, "r", errors="replace") as fh:
            lines = fh.readlines()[-ENERGY_KEEP:]
    except OSError:
        return []
    rows = []
    for line in lines:
        parts = line.split()
        if len(parts) < 5:
            continue
        try:
            rows.append([float(v) for v in parts[:5]])
        except ValueError:
            pass
    return rows


def _bounds(params):
    corner = params.get("boxcorner") or [0.0, 0.0, 0.0]
    size = params.get("boxsize") or [1.0, 1.0, 1.0]
    if len(corner) < 3 or len(size) < 3:
        return {"min": [0.0, 0.0, 0.0], "max": [1.0, 1.0, 1.0]}
    return {"min": [float(corner[k]) for k in range(3)],
            "max": [float(corner[k]) + float(size[k]) for k in range(3)]}


class Viewer:
    def __init__(self, config_path, rundir):
        self.rundir = rundir
        self.config_path = config_path
        self.run = Run(rundir)
        self.cache = FrameCache()
        with open(config_path) as fh:
            self.initial_config = fh.read()

    def status(self):
        self.run.note_exit()
        self.cache.scan(self.rundir)
        params = _config_params(self.run.config_text or self.initial_config)
        energy = _energy_tail(self.rundir)
        t = energy[-1][0] if energy else (
            self.cache.frames[-1]["t"] if self.cache.frames else 0.0)
        n = self.cache.frames[-1]["n"] if self.cache.frames else 0
        box = _bounds(params)
        rng = {}
        for k, (lo, hi) in self.cache.ranges.items():
            if lo == float("inf"):
                rng[k] = [0.0, 1.0]
            elif hi <= lo:
                rng[k] = [lo, lo + 1e-12]
            else:
                rng[k] = [lo, hi]
        return {
            "running": self.run.running(),
            "exit_code": self.run.exit_code,
            "seed": self.run.seed,
            "rundir": self.rundir,
            "box": box,
            "ranges": rng,
            "params": params,
            "t": t,
            "n": n,
            "elapsed": (time.time() - self.run.started)
                       if self.run.started and self.run.running() else None,
            "frames": list(self.cache.frames),
            "energy": energy,
            "log": self._log_tail(),
        }

    def _log_tail(self):
        path = os.path.join(self.rundir, "solver.log")
        try:
            with open(path, "rb") as fh:
                fh.seek(0, os.SEEK_END)
                size = fh.tell()
                fh.seek(max(0, size - 4000))
                return fh.read().decode("utf-8", "replace")
        except OSError:
            return ""


LIVE_JS = r"""
// ------------------------------------------------------------------ state
let frames = [];
let current = -1;
let follow = true;
let playing = false;
let seeded = false;
const cache = new Map();

const $ = (id) => document.getElementById(id);
const slider = $('slider'), label = $('label'), count = $('count');
const playBtn = $('play'), statusEl = $('status'), statsEl = $('stats');
const barEl = $('bar'), barText = $('bartext'), logEl = $('log');
const configEl = $('config'), seedEl = $('seed');
const startBtn = $('start'), stopBtn = $('stop');

async function getFrame(name) {
  if (cache.has(name)) return cache.get(name);
  const r = await fetch('/api/frame/' + encodeURIComponent(name));
  if (!r.ok) throw new Error('frame ' + name + ': HTTP ' + r.status);
  const arr = new Float32Array(await r.arrayBuffer());
  cache.set(name, arr);
  return arr;
}

async function showFrame(i) {
  if (i < 0 || i >= frames.length) return;
  current = i;
  const f = frames[i];
  const arr = await getFrame(f.name);
  if (current !== i) return;                 // a newer request won
  applyFrame(arr, f.n);
  slider.value = i;
  label.textContent = f.name + (f.t ? '   t = ' + f.t.toFixed(4) : '');
  count.textContent = f.n + (f.n === 1 ? ' particle' : ' particles');
}

// ------------------------------------------------------------------- stats
function fmt(x) {
  if (x === null || x === undefined) return '\u2014';
  if (typeof x !== 'number') return String(x);
  if (x !== 0 && (Math.abs(x) < 1e-3 || Math.abs(x) >= 1e5))
    return x.toExponential(3);
  return x.toPrecision(5);
}

function updateStats(st) {
  const p = st.params || {};
  const running = st.running;

  statusEl.textContent = running ? 'running' : (st.exit_code === null ? 'idle'
                          : (st.exit_code === 0 ? 'finished' : 'exited ' + st.exit_code));
  statusEl.className = running ? 'run' : (st.exit_code === 0 ? 'ok' : 'idle');
  startBtn.disabled = running;
  stopBtn.disabled = !running;
  configEl.readOnly = running;

  // progress: simulated time against maxTime, particles against nParticle
  let frac = 0, text = '';
  if (p.maxTime) {
    frac = Math.max(0, Math.min(1, (st.t || 0) / p.maxTime));
    text = 't = ' + fmt(st.t) + ' / ' + fmt(p.maxTime);
  }
  barEl.style.width = (frac * 100).toFixed(1) + '%';
  if (p.nParticle) text += '     ' + st.n + ' / ' + p.nParticle + ' particles';
  barText.textContent = text;

  const rows = [
    ['simulated time', fmt(st.t)],
    ['particles', st.n],
    ['frames', st.frames.length],
    ['wall time', st.elapsed ? fmt(st.elapsed) + ' s' : '\u2014'],
    ['seed', st.seed],
  ];
  const last = st.energy && st.energy.length ? st.energy[st.energy.length - 1] : null;
  if (last) {
    rows.push(['total energy', fmt(last[1])]);
    rows.push(['kinetic', fmt(last[2])]);
    rows.push(['potential', fmt(last[3])]);
    rows.push(['repulsive', fmt(last[4])]);
  }
  statsEl.innerHTML = rows.map(
    ([k, v]) => '<tr><td class="k">' + k + '</td><td class="v">' + v + '</td></tr>'
  ).join('');

  if (st.log) {
    logEl.textContent = st.log;
    logEl.scrollTop = logEl.scrollHeight;
  }
}

// ------------------------------------------------------------------ polling
let polling = false;
async function poll() {
  let st;
  try {
    st = await (await fetch('/api/status')).json();
  } catch (e) {
    statusEl.textContent = 'server gone';
    return;
  }
  if (st.box && !boxReady) setBox(st.box);
  for (const k in st.ranges) ranges[k] = st.ranges[k];
  frames = st.frames;

  const grew = frames.length && (current >= frames.length - 1);
  slider.max = Math.max(0, frames.length - 1);
  updateStats(st);

  if (frames.length && current < 0) await showFrame(0);
  if (follow && frames.length) await showFrame(frames.length - 1);
  else if (grew && playing) await showFrame(frames.length - 1);
}

async function loop() {
  if (polling) return;
  polling = true;
  while (true) {
    await poll();
    await new Promise(r => setTimeout(r, 800));
  }
}

// ----------------------------------------------------------------------- UI
startBtn.addEventListener('click', async () => {
  startBtn.disabled = true;
  cache.clear();
  current = -1;
  follow = true;
  $('follow').checked = true;
  const body = JSON.stringify({
    config: configEl.value,
    seed: parseInt(seedEl.value || '0', 10),
  });
  const r = await fetch('/api/run', { method: 'POST', body });
  if (!r.ok) {
    const t = await r.text();
    alert('could not start: ' + t);
  }
});

stopBtn.addEventListener('click', async () => {
  await fetch('/api/stop', { method: 'POST' });
});

$('reload').addEventListener('click', async () => {
  const st = await (await fetch('/api/config')).json();
  configEl.value = st.config;
});

playBtn.addEventListener('click', () => {
  playing = !playing;
  follow = false;
  $('follow').checked = false;
  playBtn.textContent = playing ? '\u23F8' : '\u25B6';
});

slider.addEventListener('input', () => {
  playing = false;
  playBtn.textContent = '\u25B6';
  follow = false;
  $('follow').checked = false;
  showFrame(+slider.value);
});

$('follow').addEventListener('change', (e) => {
  follow = e.target.checked;
  if (follow) { playing = false; playBtn.textContent = '\u25B6'; }
});

$('color').addEventListener('change', (e) => {
  setColorMode(e.target.value);
  if (current >= 0) showFrame(current);
});
$('box').addEventListener('change', (e) => { cell.visible = e.target.checked; });
addEventListener('keydown', (e) => {
  if (e.code === 'Space') { e.preventDefault(); playBtn.click(); }
});

let last = performance.now(), acc = 0;
function animate(now) {
  requestAnimationFrame(animate);
  const dt = (now - last) / 1000;
  last = now;
  if ($('spin').checked) spinCamera(dt);
  if (playing && frames.length > 1) {
    acc += dt;
    if (acc >= 1 / 20) { acc = 0; showFrame((current + 1) % frames.length); }
  }
  drawScene();
}
requestAnimationFrame(animate);
loop();
"""


def page(viewer):
    params = _config_params(viewer.initial_config)
    box = _bounds(params)
    return (PAGE_HEAD
            .replace("__CSS__", CSS)
            .replace("__CONFIG__", _escape(viewer.initial_config))
            + importmap(THREE_VERSION)
            + '<script type="module">\n'
            + scene_js(box, {k: [0.0, 1.0] for k in
                             ("a", "c", "aspect", "volume")}, "random")
            + LIVE_JS
            + PAGE_TAIL)


def _escape(text):
    return (text.replace("&", "&amp;").replace("<", "&lt;")
                .replace(">", "&gt;"))


PAGE_HEAD = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>ellipmd live</title>
<style>__CSS__
  #side { position: fixed; right: 0; top: 0; bottom: 0; width: 360px;
          overflow-y: auto; padding: 12px; display: flex; flex-direction: column;
          gap: 10px; background: rgba(20,22,26,.96);
          border-left: 1px solid rgba(255,255,255,.10); }
  #view { right: 360px; }
  #ui { right: 372px; }
  h2 { font-size: 13px; margin: 2px 0 4px; font-weight: 600; opacity: .85; }
  textarea { width: 100%; height: 260px; resize: vertical;
             background: rgba(255,255,255,.05); color: inherit;
             border: 1px solid rgba(255,255,255,.16); border-radius: 7px;
             padding: 7px; font: 11px/1.45 ui-monospace, Menlo, Consolas, monospace; }
  .row { display: flex; gap: 8px; align-items: center; }
  .row input[type=number] { width: 80px; }
  #start { background: #2c6e3f; border-color: #3c8f52; }
  #stop  { background: #6e2c2c; border-color: #8f3c3c; }
  #barwrap { position: relative; height: 18px; border-radius: 9px;
             background: rgba(255,255,255,.08); overflow: hidden; }
  #bar { height: 100%; width: 0%; background: #3f7fd0; transition: width .3s; }
  #bartext { position: absolute; inset: 0; text-align: center;
             font-size: 11px; line-height: 18px; text-shadow: 0 1px 2px #000; }
  #status { font-weight: 600; }
  #status.run { color: #7fd08f; }
  #status.ok { color: #9fb8ff; }
  #status.idle { color: #c9a06a; }
  table { width: 100%; border-collapse: collapse; font-variant-numeric: tabular-nums; }
  td { padding: 2px 0; }
  td.k { opacity: .65; }
  td.v { text-align: right; }
  #log { height: 110px; overflow: auto; font: 10px/1.4 ui-monospace, Menlo, monospace;
         background: rgba(0,0,0,.35); border-radius: 6px; padding: 6px;
         white-space: pre-wrap; opacity: .85; }
  .hint { opacity: .55; font-size: 11px; }
</style>
</head>
<body>
<div id="view"></div>
<div id="hud" class="panel">
  <b>ellipmd live</b><br>
  <span id="count"></span><br>
  <span class="dim">drag = orbit &middot; shift+drag = pan &middot; wheel = zoom &middot; z is up</span>
</div>
<div id="ui" class="panel">
  <button id="play" title="play/pause (space)">&#9654;</button>
  <input id="slider" type="range" min="0" max="0" value="0" step="1">
  <span id="label"></span>
  <label><input id="follow" type="checkbox" checked> follow</label>
  <label>colour
    <select id="color">
      <option value="uniform">uniform</option>
      <option value="random" selected>random per particle</option>
      <option value="c">long semi-axis c</option>
      <option value="a">short semi-axis a</option>
      <option value="aspect">aspect ratio c/a</option>
      <option value="volume">volume</option>
    </select>
  </label>
  <label><input id="box" type="checkbox" checked> box</label>
  <label><input id="spin" type="checkbox"> spin</label>
</div>

<div id="side">
  <div>
    <h2>run</h2>
    <div class="row">
      <button id="start">Start</button>
      <button id="stop" disabled>Stop</button>
      <span id="status" class="idle">idle</span>
    </div>
    <div class="row" style="margin-top:6px">
      <label>seed <input id="seed" type="number" value="1" min="0" step="1"></label>
      <button id="reload" title="reload the config file from disk">reload file</button>
    </div>
  </div>

  <div>
    <h2>progress</h2>
    <div id="barwrap"><div id="bar"></div><div id="bartext"></div></div>
    <table id="stats" style="margin-top:6px"></table>
  </div>

  <div>
    <h2>config</h2>
    <textarea id="config" spellcheck="false">__CONFIG__</textarea>
    <p class="hint">Edit and press Start.  The solver validates it; anything it
       rejects appears in the log below.</p>
  </div>

  <div>
    <h2>solver output</h2>
    <div id="log"></div>
  </div>
</div>
""" + ERROR_BLOCK + r"""
"""

PAGE_TAIL = r"""
</script>
</body>
</html>
"""


class Handler(BaseHTTPRequestHandler):
    viewer = None

    def log_message(self, *a):
        pass                                  # keep the console for the banner

    def _send(self, code, body, ctype="text/plain; charset=utf-8"):
        if isinstance(body, str):
            body = body.encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
        try:
            self.wfile.write(body)
        except BrokenPipeError:
            pass

    def _json(self, obj, code=200):
        self._send(code, json.dumps(obj), "application/json")

    def do_GET(self):
        path = self.path.split("?", 1)[0]
        v = self.viewer
        if path in ("/", "/index.html"):
            return self._send(200, page(v), "text/html; charset=utf-8")
        if path == "/api/status":
            return self._json(v.status())
        if path == "/api/config":
            return self._json({"config": v.initial_config})
        if path.startswith("/api/frame/"):
            name = path[len("/api/frame/"):]
            blob = v.cache.blob(name)
            if blob is None:
                return self._send(404, "no such frame")
            return self._send(200, blob, "application/octet-stream")
        return self._send(404, "not found")

    def do_POST(self):
        v = self.viewer
        path = self.path.split("?", 1)[0]
        length = int(self.headers.get("Content-Length") or 0)
        body = self.rfile.read(length) if length else b""
        if path == "/api/run":
            try:
                payload = json.loads(body or b"{}")
            except ValueError:
                return self._send(400, "bad JSON")
            config = payload.get("config") or v.initial_config
            seed = int(payload.get("seed", 0))
            if not os.path.exists(ELF):
                return self._send(500, "solver not built; run 'make ellipmd'")
            try:
                v.cache.reset()
                v.run.start(config, seed)
            except OSError as exc:
                return self._send(500, str(exc))
            return self._json({"ok": True, "rundir": v.run.rundir})
        if path == "/api/stop":
            v.run.stop()
            return self._json({"ok": True})
        return self._send(404, "not found")


def _free_port(preferred):
    with socket.socket() as s:
        try:
            s.bind(("127.0.0.1", preferred))
            return preferred
        except OSError:
            s.bind(("127.0.0.1", 0))
            return s.getsockname()[1]


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Run ellipmd and watch it live in the browser.")
    ap.add_argument("--port", type=int, default=8770)
    ap.add_argument("--config", default=os.path.join(ROOT, "config_rain"),
                    help="config file loaded into the form (default config_rain)")
    ap.add_argument("--rundir", default=os.path.join(ROOT, "viz", "live"),
                    help="where the run writes its snapshots "
                         "(default viz/live)")
    args = ap.parse_args(argv)

    if not os.path.exists(args.config):
        ap.error("no such config file: %s" % args.config)
    if not os.path.exists(ELF):
        print("warning: %s does not exist; run 'make ellipmd' first" % ELF,
              file=sys.stderr)

    viewer = Viewer(args.config, args.rundir)
    Handler.viewer = viewer

    port = _free_port(args.port)
    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    url = "http://127.0.0.1:%d/" % port
    print("ellipmd live viewer on %s" % url)
    print("  config:  %s" % args.config)
    print("  rundir:  %s" % args.rundir)
    print("  ctrl-c to stop")
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print()
    finally:
        viewer.run.stop()
        httpd.server_close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
