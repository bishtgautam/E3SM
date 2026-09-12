#!/usr/bin/env python3
"""Fail if a 12-month `1x1_glc` run has regressed against a stored reference.

Why this exists. On 2026-09-11 a change that was correct in isolation -- the
`sabg_lyr` reflection -- improved 2 fields and regressed 33, because it removed
a second bug that had been propping it up. That was caught by a 12-minute run
and a hand-read table of numbers. It should have been one command, and the next
such regression should not depend on someone thinking to look.

The metric is each field's worst-month error scaled to its own annual peak in
ELM, matching tools/compare_h0.py. Never raw relative error: that just finds
whichever month put a small number in the denominator, which is what made
"int_snow is 100x too small" look like the top defect in G24 when int_snow
agreed in eleven months of twelve.

Usage:
    check_h0_regression.py <elmxx_h0_glob> <elm_h0_glob> [--ref FILE]
                           [--update] [--tol 0.10]

    --update  rewrite the reference from this run. Do this ONLY after
              confirming the run is an improvement -- the reference is the
              thing that makes a regression visible, so a careless update
              silently raises the floor.

Exit status is 1 if any field regressed by more than --tol (fractional, so
0.10 = 10% worse than the reference), else 0.
"""
import argparse
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from compare_h0 import dump, CONSERVATION  # noqa: E402

SKIP = {"time", "time_bounds", "lon", "lat", "area", "mcdate", "mcsec",
        "mdcur", "mscur", "nstep", "levgrnd", "date_written", "time_written"}

# Fields whose ELM annual peak is so small that a scaled error carries no
# information -- aerosol deposition at 1e-12, canopy evaporation at 1e-07.
# Reported, never enforced. Keeping them enforceable would mean the check
# fires on noise and gets ignored, which is how test_integration went red.
NEGLIGIBLE_PEAK = 1.0e-6


def load_year(pattern):
    """{field: [12 monthly values]} from one glob."""
    out = {}
    paths = sorted(glob.glob(pattern))[:12]
    if len(paths) < 12:
        sys.exit(f"expected 12 monthly tapes, found {len(paths)}: {pattern}")
    for p in paths:
        for k, v in dump(p).items():
            out.setdefault(k, []).append(v[0])
    return out


def scaled_errors(elmxx, elm):
    """{field: (scaled worst-month error, ELM annual peak)}."""
    fields = sorted((set(elmxx) & set(elm)) - SKIP - CONSERVATION)
    out = {}
    for f in fields:
        if len(elm[f]) != 12 or len(elmxx[f]) != 12:
            continue
        peak = max(abs(v) for v in elm[f])
        scale = peak or 1.0e-30
        out[f] = (max(abs(a - b) for a, b in zip(elmxx[f], elm[f])) / scale,
                  peak)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx")
    ap.add_argument("elm")
    ap.add_argument("--ref", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "h0_reference.json"))
    ap.add_argument("--update", action="store_true")
    ap.add_argument("--tol", type=float, default=0.10)
    a = ap.parse_args()

    cur = scaled_errors(load_year(a.elmxx), load_year(a.elm))

    if a.update or not os.path.exists(a.ref):
        with open(a.ref, "w") as fh:
            json.dump({f: e for f, (e, _) in sorted(cur.items())}, fh,
                      indent=1, sort_keys=True)
        print(f"wrote reference for {len(cur)} fields -> {a.ref}")
        if not a.update:
            print("NOTE: no reference existed, so this run became the "
                  "reference. It has not been checked against anything.")
        return 0

    with open(a.ref) as fh:
        ref = json.load(fh)

    regressed, improved, enforced = [], [], 0
    for f, (err, peak) in sorted(cur.items()):
        if f not in ref:
            print(f"  NEW   {f:<16}{err:>10.3e}")
            continue
        prev = ref[f]
        informative = peak > NEGLIGIBLE_PEAK
        enforced += informative
        if err > prev * (1.0 + a.tol) and err > 1.0e-12:
            (regressed if informative else improved).append((f, prev, err, peak))
        elif err < prev * (1.0 - a.tol):
            improved.append((f, prev, err, peak))

    print(f"{len(cur)} fields compared, {enforced} enforced "
          f"({len(cur) - enforced} have an ELM peak below {NEGLIGIBLE_PEAK:g} "
          f"and are reported only)")

    if improved:
        print("\nimproved (or moved on a negligible field):")
        for f, prev, err, peak in improved:
            print(f"  {f:<16}{prev:>10.3e} -> {err:>10.3e}   (peak {peak:.4g})")

    if regressed:
        print("\nREGRESSED:")
        for f, prev, err, peak in regressed:
            print(f"  {f:<16}{prev:>10.3e} -> {err:>10.3e}   "
                  f"({err / prev:.2f}x, peak {peak:.4g})")
        print(f"\n{len(regressed)} field(s) regressed by more than "
              f"{a.tol:.0%}. If this is intended, re-run with --update and "
              f"say in the commit message why the trade is worth it.")
        return 1

    print("\nno regressions.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
