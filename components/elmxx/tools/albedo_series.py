#!/usr/bin/env python3
"""Daily ELMxx-vs-ELM series for the snow-albedo chain over a timestep window.

`compare_elmxx_trajectory.py` answers "which variable departs first"; this
answers "and what does the departure look like". Written for G26: the April
tapes say ELMxx reflects 27 W/m2 more than ELM on identical incoming, and the
open question is whether that albedo gap leads the snow-mass gap (albedo is
the cause) or trails it (albedo is a feedback on having more snow).

Usage:
    albedo_series.py <elmxx.bin> <elm.bin> [--from TS] [--to TS] [--stride N]
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from compare_elmxx_trajectory import read_bin  # noqa: E402

# (label on the ELMxx side, label on the ELM side, band or None)
PAIRS = [
    ("h2osno",     "elmxx_in:h2osno",     "canhydro_in:h2osno",     None),
    ("frac_sno",   "elmxx_in:frac_sno",   "canhydro_in:frac_sno",   None),
    ("snow_depth", "elmxx_in:snow_depth", "canhydro_in:snow_depth", None),
    # The two inputs to CanopyHydrology's fsca MELT branch. That branch is
    # what collapses frac_sno once a pack starts ablating, and n_melt = 0.2498
    # makes it violently sensitive (G24). If ELM's fires and ELMxx's does not,
    # the albedo gap is a consequence of the melt gap, not its cause.
    ("int_snow",   "elmxx_in:int_snow",   "canhydro_in:int_snow",   None),
    ("q_snow_melt", "elmxx_snow_pre:qflx_snow_melt",
                    "canhydro_in:qflx_snow_melt", None),
    ("albgrd_vis", "elmxx_out:albgrd",    "surfrad_in:albgrd",      0),
    ("albgrd_nir", "elmxx_out:albgrd",    "surfrad_in:albgrd",      1),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx")
    ap.add_argument("elm")
    ap.add_argument("--from", dest="lo", type=int, default=1)
    ap.add_argument("--to", dest="hi", type=int, default=10**9)
    ap.add_argument("--stride", type=int, default=48)  # one sample per day
    args = ap.parse_args()

    xx_keep = {"elmxxmap:col_of_kcol"} | {p[1] for p in PAIRS}
    em_keep = {p[2] for p in PAIRS}
    xx = read_bin(args.elmxx, xx_keep)
    em = read_bin(args.elm, em_keep)

    col = int(xx.get(0, {})["elmxxmap:col_of_kcol"][0]) - 1

    steps = [t for t in sorted(xx) if t > 0 and args.lo <= t <= args.hi and t in em]
    steps = steps[:: args.stride]
    if not steps:
        sys.exit("no overlapping timesteps in that window")

    names = [p[0] for p in PAIRS]
    print(f"column {col}, {len(steps)} samples, stride {args.stride}\n")
    hdr = f"{'ts':>6} " + " ".join(f"{n:>22}" for n in names)
    print(hdr)
    print(f"{'':>6} " + " ".join(f"{'ELMxx / ELM':>22}" for _ in names))
    print("-" * len(hdr))
    for ts in steps:
        cells = []
        for _, xl, el, band in PAIRS:
            a, b = xx[ts].get(xl), em[ts].get(el)
            if a is None or b is None:
                cells.append(f"{'-':>22}")
                continue
            av = a[0, band] if band is not None else a[0]
            bv = b[col, band] if band is not None else b[col]
            cells.append(f"{av:>10.4f} /{bv:>10.4f}")
        print(f"{ts:>6} " + " ".join(cells))


if __name__ == "__main__":
    main()
