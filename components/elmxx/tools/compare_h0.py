#!/usr/bin/env python3
"""Compare monthly ELMxx h0 tapes against the ELM twin's, field by field.

This is the whole-model acceptance check behind STATUS §0: not a kernel replay,
but "does a free-running year agree with ELM on every history field". It exists
because nothing else in tools/ reads h0 -- the numbers in §0 were previously
made by hand.

Only `ncdump` is available on this machine (no netCDF4/h5py), so the data
section is parsed as text. The 1x1 tapes carry one value per field per month,
which makes that cheap and exact.

Usage:
    compare_h0.py <elmxx_h0_glob> <elm_h0_glob> [--rtol 1e-3]

Both globs must expand to one file per month, sorted by name.
"""
import argparse, glob, re, subprocess, sys
from collections import defaultdict

# Fields that are supposed to be ~0: grade on MAGNITUDE, never relative error.
CONSERVATION = {"ERRH2O", "ERRSEB", "ERRSOL", "ERRSOI", "ERRH2OSNO"}


def dump(path):
    """{varname: [values]} for every variable in the data section."""
    txt = subprocess.run(["ncdump", path], capture_output=True, text=True,
                         check=True).stdout
    body = txt.split("\ndata:\n", 1)[1]
    out = {}
    for m in re.finditer(r"^\s*(\w+)\s*=\s*(.*?);", body, re.S | re.M):
        name, blob = m.group(1), m.group(2)
        vals = []
        for tok in blob.replace("\n", " ").split(","):
            tok = tok.strip()
            if not tok or tok.startswith('"'):
                continue
            try:
                vals.append(float(tok))
            except ValueError:
                pass
        if vals:
            out[name] = vals
    return out


def relerr(a, b, scale):
    """Error relative to the field's own ANNUAL scale, not to this month's value.

    A plain |a-b|/|b| ranks a field by whichever month happens to put a small
    number in the denominator -- that is what made "int_snow is 100x too small"
    (G24) look like the top defect when int_snow agreed in eleven months of
    twelve. Dividing by the year's peak magnitude asks the question that
    matters: is this month's disagreement large FOR THIS FIELD.
    """
    return abs(a - b) / max(scale, 1e-30)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx"); ap.add_argument("elm")
    ap.add_argument("--rtol", type=float, default=1e-3)
    ap.add_argument("--top", type=int, default=12)
    ap.add_argument("--fields", default="", help="comma-separated: print month-by-month")
    a = ap.parse_args()

    xs, es = sorted(glob.glob(a.elmxx)), sorted(glob.glob(a.elm))
    if len(xs) != len(es):
        # ELMxx often carries a trailing month the twin does not.
        n = min(len(xs), len(es)); xs, es = xs[:n], es[:n]
    print(f"comparing {len(xs)} months\n")

    worst = defaultdict(lambda: (0.0, None))   # field -> (max relerr, month)
    mag   = defaultdict(lambda: (0.0, None))   # conservation -> (max |v|, month)
    series = defaultdict(list)                 # field -> [(month, elmxx, elm)]
    skip = {"time", "time_bounds", "lon", "lat", "area", "mcdate", "mcsec",
            "mdcur", "mscur", "nstep", "date_written", "time_written", "levgrnd"}

    for k, (xf, ef) in enumerate(zip(xs, es), start=1):
        X, E = dump(xf), dump(ef)
        for f in sorted(set(X) & set(E) - skip):
            xv, ev = X[f], E[f]
            if len(xv) != len(ev):
                continue
            if f in CONSERVATION:
                m = max(abs(v) for v in xv)
                if m > mag[f][0]:
                    mag[f] = (m, k)
                continue
            series[f].append((k, xv[0], ev[0]))

    # Scale each field by its own annual peak in ELM before ranking.
    for f, rows in series.items():
        scale = max(abs(e) for _, _, e in rows)
        for k, x, e in rows:
            r = relerr(x, e, scale)
            if r > worst[f][0]:
                worst[f] = (r, k)

    good = [f for f, (r, _) in worst.items() if r < a.rtol]
    print(f"fields under {a.rtol:g} rel in EVERY month: "
          f"{len(good)} of {len(worst)}\n")

    print("conservation fields, worst magnitude over the year (want ~0):")
    for f in sorted(mag):
        m, k = mag[f]
        # k stays None for a field that is exactly zero all year (brazil's
        # snow-balance fields), which crashed the report.
        print(f"  {f:<12} {m:.3e}   " + (f"(month {k:02d})" if k is not None else "(all months)"))

    print(f"\nworst {a.top} fields, error scaled by the field's annual peak:")
    print(f"  {'field':<14}{'err/scale':>10}  {'mon':>3}  {'ELM':>12}{'ELMxx':>12}")
    for f, (r, k) in sorted(worst.items(), key=lambda kv: -kv[1][0])[:a.top]:
        if k is None:
            # Exact in every month: no worst month was ever recorded.
            print(f"  {f:<14}{r:>10.3e}  all  (identical every month)")
            continue
        row = next(rw for rw in series[f] if rw[0] == k)
        print(f"  {f:<14}{r:>10.3e}  {k:>3d}  {row[2]:>12.5g}{row[1]:>12.5g}")

    if a.fields:
        for f in a.fields.split(","):
            if f not in series:
                print(f"\n{f}: not present"); continue
            print(f"\n{f} by month:")
            print(f"  {'mon':>3}{'ELM':>13}{'ELMxx':>13}{'diff':>13}")
            for k, x, e in series[f]:
                print(f"  {k:>3}{e:>13.5g}{x:>13.5g}{x-e:>13.5g}")


if __name__ == "__main__":
    main()
