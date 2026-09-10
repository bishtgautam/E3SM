#!/usr/bin/env python3
"""Compare an ELMxx finidat-restart run against ELM step by step from step 1.

The point of a finidat restart is that step 1 starts from ELM's own state, so
accumulated drift is removed and whatever differs was caused in that step.
That makes "how big is the error at step 1" a different and more useful
question than "when does it first exceed a tolerance".

Usage:
    finidat_step_compare.py <elmxx.bin> <elm.bin> --offset N [--steps M]

`offset` maps ELMxx step 1 onto ELM step 1+offset; verify it by checking that
a conserved state variable matches to round-off at step 1.
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from compare_elmxx_trajectory import read_bin  # noqa: E402

# name -> (ELMxx label, ELM label, entity)
PAIRS = [
    ("qflx_ev_snow",  "elmxx_st:qflx_ev_snow",  "soilflx_in:qflx_ev_snow",  "col"),
    ("t_grnd",        "elmxx_in:t_grnd",        "canhydro_in:t_grnd",       "col"),
    ("h2osno",        "elmxx_in:h2osno",        "canhydro_in:h2osno",       "col"),
    ("frac_sno",      "elmxx_in:frac_sno",      "canhydro_in:frac_sno",     "col"),
    ("qflx_evap_soi", "elmxx_in:qflx_evap_soi", "soilflx_in:qflx_evap_soi", "col"),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx")
    ap.add_argument("elm")
    ap.add_argument("--offset", type=int, required=True)
    ap.add_argument("--steps", type=int, default=288)
    args = ap.parse_args()

    xx = read_bin(args.elmxx, {"elmxxmap:col_of_kcol"} | {p[1] for p in PAIRS})
    em = read_bin(args.elm, {p[2] for p in PAIRS})
    col = int(xx[0]["elmxxmap:col_of_kcol"][0]) - 1

    def get(ts, name):
        for n, xl, el, _ in PAIRS:
            if n != name:
                continue
            a = xx.get(ts, {}).get(xl)
            b = em.get(ts + args.offset, {}).get(el)
            if a is None or b is None:
                return None, None
            return float(a[0]), float(b[col])
        return None, None

    print(f"ELMxx step 1 == ELM step {1 + args.offset}\n")
    hdr = f"{'step':>5} " + " ".join(f"{n:>26}" for n, *_ in PAIRS)
    print(hdr)
    print(f"{'':>5} " + " ".join(f"{'ELMxx / ELM':>26}" for _ in PAIRS))
    print("-" * len(hdr))
    for s in list(range(1, 13)) + list(range(24, min(args.steps, 289), 24)):
        cells = []
        for n, *_ in PAIRS:
            a, b = get(s, n)
            cells.append(f"{'-':>26}" if a is None else f"{a:>12.5g} /{b:>12.5g}")
        print(f"{s:>5} " + " ".join(cells))

    # Running means, which is what the monthly QSOIL comparison actually saw.
    print()
    for n, *_ in PAIRS:
        va = [get(s, n) for s in range(1, args.steps + 1)]
        va = [(a, b) for a, b in va if a is not None]
        if not va:
            continue
        A = np.array([a for a, _ in va]); B = np.array([b for _, b in va])
        print(f"{n:>14}  mean over {len(va)} steps: ELMxx {A.mean(): .5e}  "
              f"ELM {B.mean(): .5e}  ratio {A.mean()/B.mean() if B.mean() else float('nan'): .3f}")


if __name__ == "__main__":
    main()
