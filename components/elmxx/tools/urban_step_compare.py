#!/usr/bin/env python3
"""Urban identical-IC probe: ELMxx's urban trace against ELM's, step by step.

Both models start from one ELM restart (finidat; ELMxx seeds urban state in
elmxxUrbanMod :: elmxx_urban_apply_finidat), so the two binaries are aligned
step for step and whatever differs at step 1 was caused in step 1 -- or was
never seeded. Never call a later difference "amplification" while step 1
shows one.

Maps (one-cell domain, natural first, lake last): ELMxx packed urban column
c (0-based) is ELM column c+1 (0-based), urban patch p is ELM patch 17+p,
urban landunit k is ELM landunit k+1. Only ELM's ACTIVE urban columns are
compared (urbanflux_in:filter_urbanc): an inactive landunit is not stepped.

Usage:
    urban_step_compare.py <elmxx.bin> <elm.bin> [--steps N] [--every K]
"""
import argparse
import sys

import numpy as np

sys.path.insert(0, __file__.rsplit("/", 1)[0])
from compare_elmxx_trajectory import read_bin  # noqa: E402

SURF = ["roof", "sunwall", "shdwall", "improad", "perroad"]
NLEVSNO, NLEVURB = 5, 5

# name -> (ELMxx label, ELM label, kind). kind: col / patch / lun / soil.
PAIRS = [
    ("t_grnd",      "elmxx_urbin:t_grnd",      "canhydro_in:t_grnd",        "col"),
    ("h2osno",      "elmxx_urbin:h2osno",      "canhydro_in:h2osno",        "col"),
    ("snl",         "elmxx_urbin:snl",         "canhydro_in:snl",           "col"),
    ("snow_depth",  "elmxx_urbin:snow_depth",  "canhydro_in:snow_depth",    "col"),
    ("frac_sno",    "elmxx_urbin:frac_sno",    "canhydro_in:frac_sno",      "col"),
    ("sabg",        "elmxx_urbflx:sabg",       "urbanrad_out:sabg",         "patch"),
    ("lwrad_net",   "elmxx_urbflx:eflx_lwrad_net", "urbanrad_out:eflx_lwrad_net", "patch"),
    ("sh_grnd",     "elmxx_urbflx:eflx_sh_grnd",   "urbanflux_out:eflx_sh_grnd",  "patch"),
    ("evap_soi",    "elmxx_urbflx:qflx_evap_soi",  "urbanflux_out:qflx_evap_soi", "patch"),
    ("t_ref2m",     "elmxx_urbflx:t_ref2m",    "urbanflux_out:t_ref2m",     "patch"),
    ("taf",         "elmxx_urbflx:taf",        "urbanflux_out:taf",         "lun"),
    ("qaf",         "elmxx_urbflx:qaf",        "urbanflux_out:qaf",         "lun"),
    ("t_building",  "elmxx_urbflx:t_building", "urbanflux_out:t_building", "lun"),
    ("t_soi",       "elmxx_urbout:t_soisno_soi",   "soiltemp_out:t_soisno",     "soil"),
    # Top of step on both sides: ELM's first full-column water record.
    ("liq_soi",     "elmxx_urbin:h2osoi_liq_soi",  "cantemp_in:h2osoi_liq",     "soil"),
    ("ice_soi",     "elmxx_urbin:h2osoi_ice_soi",  "cantemp_in:h2osoi_ice",     "soil"),
]
FILTER = "urbanflux_in:filter_urbanc"


def pair(xx, em, ts, kind, xl, el, active):
    """(ELMxx, ELM) arrays over the active urban columns, or (None, None)."""
    a = xx.get(ts, {}).get(xl)
    b = em.get(ts, {}).get(el)
    if a is None or b is None:
        return None, None
    if kind == "col":
        return a[active], b[active + 1]
    if kind == "patch":
        return a[active], b[active + 17]
    if kind == "lun":
        lk = np.unique(active // 5)
        return a[lk], b[lk + 1]
    # soil: (15, nlevgrnd) vs ELM (17, nlevtot); the surface's own layers.
    if b.shape[1] > a.shape[1]:
        b = b[:, NLEVSNO:]
    out_a, out_b = [], []
    for c in active:
        n = NLEVURB if c % 5 < 3 else a.shape[1]
        out_a.append(a[c, :n])
        out_b.append(b[c + 1, :n])
    return np.concatenate(out_a), np.concatenate(out_b)


def worst_surface(xx, em, ts, kind, xl, el, active):
    """The surface carrying the worst difference at step ts."""
    if kind == "lun":
        return "landunit"
    best, who = -1.0, "?"
    for c in active:
        x, e = pair(xx, em, ts, kind, xl, el, np.array([c]))
        d = worst(x, e)
        if d == d and d > best:
            best, who = d, SURF[c % 5] + str(c // 5)
    return who


def worst(x, e):
    """max |x - e| over entries both models carry (walls hold no snow or
    water: spval on ELMxx's side), nan if there are none."""
    ok = (np.abs(x) < 1e30) & (np.abs(e) < 1e30)
    return float(np.max(np.abs(x[ok] - e[ok]))) if ok.any() else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx")
    ap.add_argument("elm")
    ap.add_argument("--steps", type=int, default=1488)
    ap.add_argument("--every", type=int, default=48)
    ap.add_argument("--first", type=float, default=1e-9,
                    help="report the first step each field exceeds this")
    a = ap.parse_args()

    xx = read_bin(a.elmxx, {p[1] for p in PAIRS})
    em = read_bin(a.elm, {p[2] for p in PAIRS} | {FILTER})

    filt = next(em[t][FILTER] for t in sorted(em) if FILTER in em[t])
    active = np.array(sorted(int(v) - 2 for v in filt))   # 1-based ELM col -> packed
    print(f"active urban columns (packed): {list(active)}  "
          f"= {[SURF[c % 5] + str(c // 5) for c in active]}\n")

    # Step 1, per surface: the seeding verdict.
    print("STEP 1, |ELMxx - ELM| per active surface (soil: worst layer)")
    print(f"{'field':<11}" + "".join(f"{SURF[c % 5] + str(c // 5):>11}" for c in active))
    for name, xl, el, kind in PAIRS:
        if kind == "lun":
            x, e = pair(xx, em, 1, kind, xl, el, active)
            if x is None:
                continue
            print(f"{name:<11}" + "".join(f"{worst(x[i:i+1], e[i:i+1]):>11.2e}"
                                          for i in range(len(x)))
                  + "   (per landunit)")
            continue
        row = []
        for c in active:
            x, e = pair(xx, em, 1, kind, xl, el, np.array([c]))
            row.append(np.nan if x is None else worst(x, e))
        print(f"{name:<11}" + "".join(f"{d:>11.2e}" for d in row))

    # Evolution: worst over active surfaces.
    names = [p[0] for p in PAIRS]
    print(f"\nworst |ELMxx - ELM| over active surfaces, every {a.every} steps")
    print(f"{'step':>5}" + "".join(f"{n:>11}" for n in names))
    peak = {n: (0.0, 0) for n in names}
    first = {}
    for ts in range(1, a.steps + 1):
        row = []
        for name, xl, el, kind in PAIRS:
            x, e = pair(xx, em, ts, kind, xl, el, active)
            d = np.nan if x is None else worst(x, e)
            row.append(d)
            if d == d and d > a.first and name not in first:
                first[name] = (ts, d, worst_surface(xx, em, ts, kind, xl, el, active))
            if d == d and d > peak[name][0]:
                peak[name] = (d, ts)
        if ts == 1 or ts % a.every == 0:
            print(f"{ts:>5}" + "".join(f"{d:>11.2e}" for d in row))
    print("\npeak (value @ step), and the first step above --first (surface):")
    for n in names:
        f = first.get(n)
        fs = f"   first > {a.first:g} @ {f[0]} ({f[1]:.2e}, {f[2]})" if f else ""
        print(f"  {n:<11}{peak[n][0]:>11.3e} @ {peak[n][1]:<6}{fs}")


if __name__ == "__main__":
    main()
