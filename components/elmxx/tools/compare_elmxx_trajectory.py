#!/usr/bin/env python3
"""Diff a free-running ELMxx run against ELM, timestep by timestep.

Every other check in this project is a *replay*: ELMxx is handed ELM's own
inputs for one kernel at one timestep. That says nothing about error growth.
This compares the state ELMxx carries forward on its own against the state ELM
carries forward, and reports the first timestep at which each variable departs.

ELMxx writes 'elmxx_in:<var>' at the top of each step, before any kernel. ELM's
matching instant is 'canhydro_in:' -- its first kernel -- with 'cantemp_in:' /
'canflx_in:' as the next-best anchor for variables it does not record there.
Those anchors are named per variable below rather than guessed.

Usage:
    compare_elmxx_trajectory.py <elmxx_diagnostics.bin> <elm_diagnostics.bin>
                                [--rtol 1e-6] [--max-ts N] [--var NAME]
"""
import argparse
import struct
import sys

import numpy as np

MAGIC = b"ELMDIAG1"

# ELMxx variable -> the ELM label describing the same quantity at the same
# instant. Chosen by inspecting which labels the ELM binary actually carries;
# do not substitute an "_out" label for an "_in" one, they are different times.
ANCHOR = {
    "t_grnd":         "canhydro_in:t_grnd",
    "t_h2osfc":       "cantemp_in:t_h2osfc",
    "h2osfc":         "canhydro_in:h2osfc",
    "h2osno":         "canhydro_in:h2osno",
    "snow_depth":     "canhydro_in:snow_depth",
    "frac_sno":       "canhydro_in:frac_sno",
    "frac_h2osfc":    "cantemp_in:frac_h2osfc",
    "int_snow":       "canhydro_in:int_snow",
    "snl":            "canhydro_in:snl",
    "t_soisno":       "cantemp_in:t_soisno",
    "h2osoi_liq_soi": "cantemp_in:h2osoi_liq",
    "h2osoi_ice_soi": "cantemp_in:h2osoi_ice",
    "t_veg":          "canflx_in:t_veg",
    "btran":          "canflx_in:btran",
    "h2ocan":         "canhydro_in:h2ocan",
}

# Sampled at the END of the step (elmxx_out:), because ELM records these
# during the step. Comparing a top-of-step value against them shifts the whole
# diurnal cycle by one timestep and looks like a large error.
END_OF_STEP = {"fsa", "fsr", "albgrd", "albgri", "albsod"}
ANCHOR["fsa"] = "surfrad_out:fsa"
ANCHOR["fsr"] = "surfrad_out:fsr"

# Ground albedo, the quantity April says is wrong (FSA 87.9 vs ELM's 114.9 on
# identical incoming). The instants match: elmxx_diag_snapshot_fluxes runs
# BEFORE the phenology/SurfaceAlbedo block, so elmxx_out:albgrd at step N is
# the albedo computed at the end of step N-1 -- exactly the one ELM's
# surfrad_in:albgrd at step N describes.
ANCHOR["albgrd"] = "surfrad_in:albgrd"
ANCHOR["albgri"] = "surfrad_in:albgri"
ANCHOR["albsod"] = "surfrad_in:albsod"

COLUMN_VARS = {"t_grnd", "t_h2osfc", "h2osfc", "h2osno", "snow_depth",
               "frac_sno", "frac_h2osfc", "int_snow", "snl",
               "t_soisno", "h2osoi_liq_soi", "h2osoi_ice_soi",
               "albgrd", "albgri", "albsod"}

# Full-depth column profiles: slots 0..NLEVSNO-1 are snow, the rest soil.
NLEVSNO = 5
PROFILE_NLEVTOT = {"t_soisno"}

# ELM's own active-patch filter. Comparing the inactive patches is meaningless
# -- most of the 17 carry zero weight and are never updated -- and it produced
# relative errors of exactly 1.0 that looked like real divergence.
PATCH_FILTER = "canflx_in:filter_nolakeurbanp"


def read_bin(path, keep=None):
    """Read an ELMDIAG1 file into {timestep: {label: array}}.

    `keep`, when given, is the set of labels to retain. Every record is still
    parsed -- the format has no index to seek by -- but the unwanted arrays
    are dropped instead of held. The Jan-Apr reference binary is 4.0 GB and
    does not fit in memory otherwise.
    """
    out = {}
    with open(path, "rb") as f:
        if f.read(8) != MAGIC:
            sys.exit(f"{path}: not an ELMDIAG1 file")
        (label_len,) = struct.unpack("<i", f.read(4))
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            (ts,) = struct.unpack("<i", head)
            label = f.read(label_len).decode("ascii", "replace").strip()
            (ndims,) = struct.unpack("<i", f.read(4))
            if ndims == 0:
                val = np.frombuffer(f.read(8), "<f8")
            elif ndims == 1:
                (n,) = struct.unpack("<i", f.read(4))
                val = np.frombuffer(f.read(8 * n), "<f8")
            elif ndims == 2:
                n1, n2 = struct.unpack("<ii", f.read(8))
                val = np.frombuffer(f.read(8 * n1 * n2), "<f8").reshape((n1, n2), order="F")
            elif ndims == -1:
                (n,) = struct.unpack("<i", f.read(4))
                val = np.frombuffer(f.read(4 * n), "<i4").astype(float)
            else:
                sys.exit(f"{path}: unknown ndims {ndims}")
            if keep is None or label in keep:
                out.setdefault(ts, {})[label] = np.array(val)
    return out


def relerr(a, b):
    """Relative error, falling back to absolute where the reference is ~0."""
    denom = np.maximum(np.abs(b), 1e-30)
    err = np.abs(a - b) / denom
    tiny = np.abs(b) < 1e-12
    err[tiny] = np.abs(a - b)[tiny]
    return err


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("elmxx")
    ap.add_argument("elm")
    ap.add_argument("--rtol", type=float, default=1e-6)
    ap.add_argument("--max-ts", type=int, default=None)
    ap.add_argument("--var", default=None, help="restrict to one variable")
    args = ap.parse_args()

    wanted = [v for v in ANCHOR if args.var in (None, v)]
    xx_keep = {"elmxxmap:col_of_kcol", "elmxxmap:patch_of_kpatch"}
    for v in wanted:
        xx_keep.add(f"elmxx_out:{v}" if v in END_OF_STEP else f"elmxx_in:{v}")
    em_keep = {PATCH_FILTER} | {ANCHOR[v] for v in wanted}

    xx = read_bin(args.elmxx, xx_keep)
    em = read_bin(args.elm, em_keep)

    maps = xx.get(0, {})
    if "elmxxmap:col_of_kcol" not in maps:
        sys.exit("ELMxx binary has no elmxxmap records; cannot align indices")
    col_map = maps["elmxxmap:col_of_kcol"].astype(int) - 1      # -> 0-based ELM col
    pat_map = maps["elmxxmap:patch_of_kpatch"].astype(int) - 1  # -> 0-based ELM patch

    steps = sorted(t for t in xx if t > 0)
    if args.max_ts:
        steps = [t for t in steps if t <= args.max_ts]
    shared = [t for t in steps if t in em]
    print(f"ELMxx steps: {len(steps)}   overlapping with ELM: {len(shared)}")
    print(f"columns: {len(col_map)}   patches: {len(pat_map)}   rtol: {args.rtol:g}\n")
    if not shared:
        sys.exit("no overlapping timesteps -- did both runs start from the same date?")

    # Restrict patches to the ones ELM actually integrates.
    active_p = None
    for ts in shared:
        if PATCH_FILTER in em[ts]:
            active_p = em[ts][PATCH_FILTER].astype(int) - 1  # -> 0-based ELM patch
            break
    if active_p is not None:
        keep = [k for k, gp in enumerate(pat_map) if gp in set(active_p.tolist())]
        pat_map = pat_map[keep]
        patch_sel = np.array(keep)
        print(f"active patches: {len(pat_map)} of {len(maps['elmxxmap:patch_of_kpatch'])}"
              f"  (ELM patches {(pat_map + 1).tolist()})\n")
    else:
        patch_sel = None
        print("WARNING: no patch filter found; comparing all patches\n")

    variables = wanted
    rows = []
    for var in sorted(variables):
        anchor = ANCHOR[var]
        first_bad, worst, worst_ts = None, 0.0, None
        seen = 0
        # A variable that is identically zero in BOTH models over every
        # compared entity is not agreement, it is an empty comparison. On a
        # tropical site every snow field is like this, and reporting them as
        # "never diverges" overstates what the trace covers.
        nonzero = False
        for ts in shared:
            prefix = "elmxx_out" if var in END_OF_STEP else "elmxx_in"
            a = xx[ts].get(f"{prefix}:{var}")
            b = em[ts].get(anchor)
            if a is None or b is None:
                continue
            is_col = var in COLUMN_VARS
            idx = col_map if is_col else pat_map
            if not is_col and patch_sel is not None:
                a = a[patch_sel]
            try:
                bb = b[idx] if b.ndim == 1 else b[idx, :]
            except IndexError:
                continue
            if var in PROFILE_NLEVTOT:
                # Drop the snow slots: with snl == 0 they are inactive in both
                # models and hold stale values, not physics.
                a, bb = a[:, NLEVSNO:], bb[:, NLEVSNO:]
            if a.shape != bb.shape:
                # ELM carries NLEVTOT profiles; ELMxx's soil-only fields are shorter.
                n = min(a.shape[-1], bb.shape[-1])
                a2, bb = a[..., -n:], bb[..., -n:]
            else:
                a2 = a
            seen += 1
            aa = np.asarray(a2, float)
            bb = np.asarray(bb, float)
            if np.nanmax(np.abs(aa)) > 0.0 or np.nanmax(np.abs(bb)) > 0.0:
                nonzero = True
            e = relerr(aa, bb)
            m = float(np.nanmax(e))
            if m > worst:
                worst, worst_ts = m, ts
            if first_bad is None and m > args.rtol:
                first_bad = ts
        if seen == 0:
            rows.append((var, anchor, "-", "no overlap", "-"))
        elif not nonzero:
            rows.append((var, anchor, "VACUOUS", "both all-zero", "-"))
        else:
            rows.append((var, anchor, str(first_bad) if first_bad else "never",
                         f"{worst:.3e}", str(worst_ts)))

    w = max(len(r[0]) for r in rows) + 1
    print(f"{'variable':<{w}} {'ELM anchor':<26} {'diverges@':>10} {'worst':>12} {'@ts':>6}")
    print("-" * (w + 60))
    for var, anchor, fb, worst, wts in rows:
        print(f"{var:<{w}} {anchor:<26} {fb:>10} {worst:>12} {wts:>6}")

    vac = [r for r in rows if r[2] == "VACUOUS"]
    if vac:
        print(f"\n{len(vac)} of {len(rows)} variables are VACUOUS on this case -- "
              f"identically zero in both models, so they grade nothing:")
        print("   " + ", ".join(r[0] for r in vac))
    bad = [r for r in rows if r[2] not in ("never", "-", "VACUOUS")]
    print()
    if bad:
        earliest = min(int(r[2]) for r in bad)
        print(f"{len(bad)} variable(s) exceed rtol; earliest departure at timestep {earliest}:")
        for r in sorted(bad, key=lambda r: int(r[2])):
            print(f"   ts {r[2]:>5}  {r[0]}  (worst {r[3]})")
    else:
        print(f"No variable exceeded rtol={args.rtol:g} over {len(shared)} timesteps.")


if __name__ == "__main__":
    main()
