#!/usr/bin/env python3
"""Daily ELMxx-vs-ELM series for snow aerosol mass, column-summed per species.

Written for STATUS §4 item 1: April snow is brighter in ELMxx than in ELM while
its grains are coarser, so grain size cannot be the cause and aerosol is the
suspect. Melt-season scavenging (S3) was graded only by construction; this
grades the aerosol burden itself against ELM.

Alignment. ELMxx's elmxx_in:mss_* is the state at the START of step N; ELM's
snowlayer_out:mss_* at step N-1 is the pack at the end of the previous step
(Combine/Divide conserve mass), so that is the pair compared.

Usage:
    aerosol_series.py <elmxx.bin> <elm.bin> [--from TS] [--to TS] [--stride N]
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from compare_elmxx_trajectory import read_bin  # noqa: E402

SPECIES = ["bcphi", "bcpho", "ocphi", "ocpho", "dst1", "dst2", "dst3", "dst4"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx")
    ap.add_argument("elm")
    ap.add_argument("--from", dest="lo", type=int, default=1)
    ap.add_argument("--to", dest="hi", type=int, default=10**9)
    ap.add_argument("--stride", type=int, default=48)  # one sample per day
    args = ap.parse_args()

    xx_keep = {"elmxxmap:col_of_kcol", "elmxx_in:h2osno"} | \
              {f"elmxx_in:mss_{s}" for s in SPECIES}
    em_keep = {"canhydro_in:h2osno"} | {f"snowlayer_out:mss_{s}" for s in SPECIES}
    xx = read_bin(args.elmxx, xx_keep)
    em = read_bin(args.elm, em_keep)

    col = int(xx.get(0, {})["elmxxmap:col_of_kcol"][0]) - 1
    steps = [t for t in sorted(xx)
             if t > 1 and args.lo <= t <= args.hi and (t - 1) in em
             and f"snowlayer_out:mss_{SPECIES[0]}" in em[t - 1]]
    steps = steps[:: args.stride]
    if not steps:
        sys.exit("no overlapping timesteps in that window")

    print(f"column {col}, {len(steps)} samples, stride {args.stride}")
    print("column-summed aerosol mass [kg/m2], ELMxx / ELM (ratio)\n")
    hdr = f"{'ts':>6} {'h2osno xx/elm':>17} " + " ".join(f"{s:>21}" for s in SPECIES)
    print(hdr)
    print("-" * len(hdr))
    for ts in steps:
        hx = xx[ts].get("elmxx_in:h2osno")
        he = em.get(ts, {}).get("canhydro_in:h2osno")
        hcell = (f"{hx[0]:>7.2f}/{he[col]:<7.2f}"
                 if hx is not None and he is not None else "-")
        cells = []
        for s in SPECIES:
            a = xx[ts].get(f"elmxx_in:mss_{s}")
            b = em[ts - 1].get(f"snowlayer_out:mss_{s}")
            if a is None or b is None:
                cells.append(f"{'-':>21}")
                continue
            av, bv = float(np.sum(a[0, :])), float(np.sum(b[col, :]))
            ratio = av / bv if bv != 0 else float("nan")
            cells.append(f"{av:>9.2e}/{bv:>9.2e}"[:21] if bv == 0
                         else f"{av:>8.2e} {ratio:>6.3f}x".rjust(21))
        print(f"{ts:>6} {hcell:>17} " + " ".join(cells))


if __name__ == "__main__":
    main()
