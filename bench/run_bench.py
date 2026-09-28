#!/usr/bin/env python3
"""Timing benchmark for ellipmd -- the driver behind BENCHMARK.md.

Runs four configurations *strictly one after another* (never concurrently, or
the timings are meaningless) and rewrites ``results.json`` after each run, so
partial results survive an interruption.  ``configs/B1..B4`` are the exact
configurations used for the numbers in BENCHMARK.md:

  B1  250 particles, dt=1e-4    2750 steps   baseline
  B2 2500 particles, dt=1e-4    2750 steps   isolates the 10x particle count
  B3  250 particles, dt=1e-5   27500 steps   isolates the 10x smaller dt
  B4 2500 particles, dt=1e-5   27500 steps   the combined run

Usage::

    python3 bench/run_bench.py                 # writes bench/results.json
    python3 bench/run_bench.py --only B4       # single config

Re-run this after an optimisation and diff ``results.json`` to see what moved.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
ELF = os.path.join(ROOT, "ellipmd")
SEED = "3"

RUNS = [
    ("B1", "250 particles, dt=1e-4", 2750),
    ("B2", "2500 particles, dt=1e-4", 2750),
    ("B3", "250 particles, dt=1e-5", 27500),
    ("B4", "2500 particles, dt=1e-5", 27500),
]


def particle_count(directory):
    files = sorted(glob.glob(os.path.join(directory, "out0*")))
    if not files:
        return 0
    n = 0
    with open(files[-1], errors="replace") as fh:
        for line in fh:
            if line.split()[:1] == ["14"]:
                n += 1
    return n


def internal_cpu(path):
    """The program prints (clock()-start)/CLOCKS_PER_SEC once per output step."""
    last = None
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        float(parts[0])
                        float(parts[1])
                    except ValueError:
                        continue
                    last = float(parts[0])
    except OSError:
        pass
    return last


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", action="append", default=None,
                    help="run only this configuration (repeatable)")
    ap.add_argument("--fast", action="store_true",
                    help="skip the longest configuration (B4, ~9 min of the "
                         "~16 min total).  Use this while iterating; do a full "
                         "sweep before recording a baseline.")
    ap.add_argument("--workdir", default=os.path.join(HERE, "runs"),
                    help="where to put the per-run output directories")
    ap.add_argument("--repeat", type=int, default=1, metavar="N",
                    help="run each configuration N times and report the "
                         "median.  A single run varies by up to ~20%% on a "
                         "laptop, which is far more than most optimisations "
                         "are worth, so use this before believing a delta.")
    args = ap.parse_args(argv)

    if not os.path.exists(ELF):
        print("error: %s not found -- run 'make ellipmd' first" % ELF,
              file=sys.stderr)
        return 2

    os.makedirs(args.workdir, exist_ok=True)
    results = []
    t_start = time.time()

    selected = RUNS
    if args.fast:
        longest = max(RUNS, key=lambda r: r[2])
        selected = [r for r in RUNS if r is not longest]
        print("--fast: skipping %s (%s, %d steps)" % (longest[0], longest[1],
                                                      longest[2]), flush=True)

    for name, label, steps in selected:
        if args.only and name not in args.only:
            continue
        cfg = os.path.join(HERE, "configs", name)
        d = os.path.join(args.workdir, name)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)

        repeat = max(1, args.repeat)
        walls = []
        rc = 1
        for rep in range(repeat):
            t0 = time.time()
            with open(os.path.join(d, "stdout.log"), "w") as log:
                rc = subprocess.call([ELF, SEED, cfg], cwd=d,
                                     stdout=log, stderr=subprocess.STDOUT)
            walls.append(time.time() - t0)
            if rc != 0:
                break
        walls.sort()
        wall = walls[len(walls) // 2]

        rec = {
            "name": name,
            "label": label,
            "steps": steps,
            "wall_s": round(wall, 2),
            "wall_s_all": [round(w, 2) for w in walls],
            "repeats": repeat,
            "cpu_s": internal_cpu(os.path.join(d, "stdout.log")),
            "ms_per_step": round(1000.0 * wall / steps, 3),
            "particles": particle_count(d),
            "rc": rc,
        }
        results.append(rec)
        with open(os.path.join(HERE, "results.json"), "w") as fh:
            json.dump(results, fh, indent=1)

        spread = ""
        if repeat > 1:
            spread = "  (of %s)" % ", ".join("%.2f" % w for w in walls)
        print("%-3s %-26s steps=%-6d wall=%8.2f s  %7.3f ms/step  "
              "cpu=%s  N=%d  rc=%d%s"
              % (name, label, steps, wall, rec["ms_per_step"],
                 rec["cpu_s"], rec["particles"], rc, spread), flush=True)

    print("total %.1f s" % (time.time() - t_start), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
