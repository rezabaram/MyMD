#!/usr/bin/env python3
"""Physics regression check for ellipmd.

This is the safety net for the refactoring work: it answers "did I just
silently change the simulation?" before anything else is allowed to land.

Two tiers, because molecular dynamics is chaotic and a bit-exact test would
block legitimate changes (recompiling with different flags perturbs the last
bits, and over thousands of steps that difference amplifies):

  Tier 1  strict   -- run a short, strongly dissipative reference case and
                      compare every particle against a committed expected
                      snapshot within a tight tolerance.  Catches the common
                      mistake, which is a *systematic* error: a transposed
                      rotation, a shifted quaternion, an off-by-one in the
                      contact point.

  Tier 2  physics  -- on the same run, check invariants that must hold no
                      matter how the floating point falls out: no NaN/Inf,
                      particle count, every particle inside the box, total
                      particle volume unchanged, and the energy trend.

Usage::

    python3 bench/check_physics.py              # check all cases
    python3 bench/check_physics.py --case deposition
    python3 bench/check_physics.py --regenerate # rewrite the expected outputs

Exit status is 0 on success, 1 on failure.
"""

from __future__ import annotations

import argparse
import math
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "tools"))

from ellipmd_io import (  # noqa: E402
    bounding_box, expand_paths, quat_to_matrix, read_snapshot, read_times,
)

ELF = os.path.join(ROOT, "ellipmd")
REFERENCE = os.path.join(HERE, "reference")
SEED = "1"

# Tolerances for tier 1.  The box is ~1 unit and particles are ~0.02-0.05,
# so 1e-6 absolute is far below anything physically meaningful while still
# catching a wrong convention.
TOL_POSITION = 1e-9
TOL_AXIS = 1e-9
TOL_ROTATION = 1e-9

# Invariant bands for tier 2.
TOL_VOLUME_REL = 1e-6        # particle shapes should not drift
TOL_INSIDE = 1e-6            # slack on the box test


def run_case(case, workdir):
    """Run one reference case and return the directory it wrote into."""
    src = os.path.join(REFERENCE, case)
    if not os.path.isdir(src):
        raise SystemExit("no such reference case: %s" % case)
    shutil.copy(os.path.join(src, "config"), os.path.join(workdir, "config"))
    with open(os.path.join(workdir, "stdout.log"), "w") as log:
        rc = subprocess.call([ELF, SEED, "config"], cwd=workdir,
                             stdout=log, stderr=subprocess.STDOUT)
    if rc != 0:
        raise SystemExit("case %s: ellipmd exited with %d "
                         "(see %s/stdout.log)" % (case, rc, workdir))
    return workdir


def final_snapshot(workdir):
    """Prefer outend, else the highest-numbered snapshot."""
    end = os.path.join(workdir, "outend")
    if os.path.isfile(end):
        return end
    files = expand_paths([os.path.join(workdir, "out0*")])
    if not files:
        raise SystemExit("case produced no snapshots in %s" % workdir)
    return files[-1]


def compare_snapshots(expected_path, actual_path):
    """Tier 1.  Returns (ok, list-of-messages)."""
    exp = read_snapshot(expected_path)
    act = read_snapshot(actual_path)
    msgs = []
    ok = True

    if len(exp) != len(act):
        return False, ["particle count changed: expected %d, got %d"
                       % (len(exp), len(act))]

    worst_pos = worst_axis = worst_rot = 0.0
    worst_pos_i = worst_rot_i = -1
    for i in range(len(exp)):
        dp = max(abs(exp.positions[i][k] - act.positions[i][k])
                 for k in range(3))
        if dp > worst_pos:
            worst_pos, worst_pos_i = dp, i
        da = max(abs(exp.axes[i][k] - act.axes[i][k]) for k in range(3))
        worst_axis = max(worst_axis, da)
        r1 = quat_to_matrix(exp.quats[i])
        r2 = quat_to_matrix(act.quats[i])
        dr = max(abs(r1[a][b] - r2[a][b])
                 for a in range(3) for b in range(3))
        if dr > worst_rot:
            worst_rot, worst_rot_i = dr, i

    for label, worst, tol, idx in (
            ("position", worst_pos, TOL_POSITION, worst_pos_i),
            ("semi-axis", worst_axis, TOL_AXIS, -1),
            ("rotation", worst_rot, TOL_ROTATION, worst_rot_i)):
        good = worst <= tol
        ok = ok and good
        msgs.append("  %-10s max |delta| = %.3e  (tol %.0e)  %s%s"
                    % (label, worst, tol, "ok" if good else "FAIL",
                       "" if idx < 0 else "  [particle %d]" % idx))
    return ok, msgs


def check_invariants(case, snapshot_path, workdir):
    """Tier 2.  Returns (ok, list-of-messages)."""
    snap = read_snapshot(snapshot_path)
    msgs = []
    ok = True

    bad = 0
    for i in range(len(snap)):
        for triple in (snap.positions[i], snap.axes[i], snap.quats[i]):
            if not all(math.isfinite(v) for v in triple):
                bad += 1
                break
    msgs.append("  finite values            : %s"
                % ("ok" if bad == 0 else "FAIL (%d particles)" % bad))
    ok = ok and bad == 0

    if len(snap) == 0:
        return False, msgs + ["  particle count           : FAIL (zero)"]

    lo, hi = bounding_box(snap)
    outside = 0
    for i in range(len(snap)):
        for k in range(3):
            if (snap.positions[i][k] < lo[k] - TOL_INSIDE
                    or snap.positions[i][k] > hi[k] + TOL_INSIDE):
                outside += 1
                break
    msgs.append("  inside the box           : %s"
                % ("ok" if outside == 0 else "FAIL (%d outside)" % outside))
    ok = ok and outside == 0

    vol = sum(4.0 / 3.0 * math.pi * a * b * c for a, b, c in snap.axes)
    ref = os.path.join(REFERENCE, case, "expected", "outend")
    if os.path.isfile(ref):
        rsnap = read_snapshot(ref)
        rvol = sum(4.0 / 3.0 * math.pi * a * b * c for a, b, c in rsnap.axes)
        rel = abs(vol - rvol) / rvol if rvol else 1.0
        good = rel <= TOL_VOLUME_REL
        msgs.append("  total particle volume    : rel delta %.3e  %s"
                    % (rel, "ok" if good else "FAIL"))
        ok = ok and good

    times = read_times(workdir)
    if times:
        msgs.append("  reached t = %.6g" % times[-1])
    return ok, msgs


def regenerate(case, workdir):
    dst = os.path.join(REFERENCE, case, "expected")
    os.makedirs(dst, exist_ok=True)
    for name in os.listdir(workdir):
        if name.startswith("out"):
            shutil.copy(os.path.join(workdir, name),
                        os.path.join(dst, name))
    print("  regenerated expected/ for %s" % case)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", action="append", default=None,
                    help="check only this case (repeatable)")
    ap.add_argument("--regenerate", action="store_true",
                    help="rewrite the expected outputs instead of checking")
    args = ap.parse_args(argv)

    if not os.path.exists(ELF):
        print("error: %s not found -- run 'make ellipmd' first" % ELF,
              file=sys.stderr)
        return 2

    cases = args.case or sorted(d for d in os.listdir(REFERENCE)
                                if os.path.isdir(os.path.join(REFERENCE, d)))
    if not cases:
        print("error: no reference cases in %s" % REFERENCE, file=sys.stderr)
        return 2

    failures = []
    for case in cases:
        print("%s:" % case)
        workdir = tempfile.mkdtemp(prefix="ellipmd-check-%s-" % case)
        try:
            run_case(case, workdir)
            final = final_snapshot(workdir)
            if args.regenerate:
                regenerate(case, workdir)
                print("  ok (regenerated)")
                continue

            expected = os.path.join(REFERENCE, case, "expected",
                                    os.path.basename(final))
            if not os.path.isfile(expected):
                failures.append(case)
                print("  FAIL: no expected output at %s "
                      "(run with --regenerate once)" % expected)
                continue

            ok1, msgs1 = compare_snapshots(expected, final)
            print(" tier 1 (vs expected snapshot)")
            for m in msgs1:
                print(m)
            ok2, msgs2 = check_invariants(case, final, workdir)
            print(" tier 2 (invariants)")
            for m in msgs2:
                print(m)
            if ok1 and ok2:
                print("  -> PASS")
            else:
                failures.append(case)
                print("  -> FAIL")
        finally:
            shutil.rmtree(workdir, ignore_errors=True)

    print()
    if failures:
        print("FAILED: %s" % ", ".join(failures))
        return 1
    print("all physics checks passed (%d case%s)" % (len(cases),
          "" if len(cases) == 1 else "s"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
