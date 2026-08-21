#!/usr/bin/env python3
"""Replay ELMxx's root water stress kernel against ELM's own records.

Drives elmxx_root_stress_kernel -- the same source the coupled run calls --
with ELM's recorded calc_root_moist_stress inputs, and compares rootr, btran
and rresis against ELM's own outputs.

rootr matters more than btran: it is the per-layer sink SoilWater uses
(qflx_rootsoi = rootr*qflx_tran_veg), and it is normalised BY btran, so btran
agreeing tells you little about whether the layers partition correctly.

Build the driver first:

  gfortran -O2 -o root_stress_replay \\
      share/util/shr_kind_mod.F90 \\
      share/util/shr_const_mod.F90 \\
      components/elmxx/src/main/elmxxRootKernelMod.F90 \\
      components/elmxx/tools/root_stress_replay.F90

Usage:

  python3 components/elmxx/tools/validate_root_stress.py \\
      --driver /path/to/root_stress_replay \\
      --diag   /path/to/run/elm_diagnostics.bin \\
      [--steps 1,50,120,240] [--rtol 1e-12]
"""
from __future__ import annotations

import argparse
import struct
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

MAGIC = b"ELMDIAG1"


def read_bin(path: Path):
    out = {}
    with open(path, "rb") as f:
        if f.read(8) != MAGIC:
            raise SystemExit(f"{path}: not an ELMDIAG1 file")
        (label_len,) = struct.unpack("<i", f.read(4))
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            (ts,) = struct.unpack("<i", head)
            label = f.read(label_len).decode("ascii", "replace").strip()
            (nd,) = struct.unpack("<i", f.read(4))
            if nd == 0:
                val = np.frombuffer(f.read(8), "<f8")
            elif nd == 1:
                (n,) = struct.unpack("<i", f.read(4))
                val = np.frombuffer(f.read(8 * n), "<f8")
            elif nd == 2:
                n1, n2 = struct.unpack("<ii", f.read(8))
                val = np.frombuffer(f.read(8 * n1 * n2), "<f8").reshape((n1, n2), order="F")
            elif nd == -1:
                (n,) = struct.unpack("<i", f.read(4))
                val = np.frombuffer(f.read(4 * n), "<i4").astype(float)
            else:
                raise SystemExit(f"unknown ndims {nd}")
            out.setdefault(ts, {})[label] = np.array(val)
    return out


def build_deck(rec, deck: Path):
    def g(name):
        key = f"rootstress_in:{name}"
        if key not in rec:
            raise SystemExit(f"missing {key} -- is the ELM build instrumented?")
        return rec[key]

    itype = g("patch_itype").astype(int)
    pcol = g("patch_column").astype(int)
    rootfr = g("rootfr")
    np_, nlevgrnd = rootfr.shape
    nc = g("watsat").shape[0]
    nlevbed = int(g("nlevbed")[0])
    tc_stress = float(g("tc_stress")[0])

    # ELM resolves smpsc/smpso per patch; the kernel takes them per PFT.
    # Scatter back, and fill unused PFT slots with the first seen value so the
    # array is well defined (only itype values that occur are ever indexed).
    npft = int(itype.max()) + 1
    smpsc_p, smpso_p = g("smpsc"), g("smpso")
    smpsc = np.zeros(npft)
    smpso = np.zeros(npft)
    for p in range(np_):
        smpsc[itype[p]] = smpsc_p[p]
        smpso[itype[p]] = smpso_p[p]

    def rows(a, n):
        return [" ".join(repr(float(v)) for v in a[i, :]) for i in range(n)]

    lines = [f"{np_} {nc} {nlevbed} {nlevgrnd} {npft}",
             f"{tc_stress!r} 0.0 917.0 1000.0",
             " ".join(str(int(v)) for v in itype),
             " ".join(str(int(v)) for v in pcol),
             " ".join(repr(float(v)) for v in smpsc),
             " ".join(repr(float(v)) for v in smpso)]
    lines += rows(rootfr, np_)
    for name in ("h2osoi_liq", "h2osoi_ice", "dz", "t_soisno",
                 "watsat", "bsw", "sucsat"):
        lines += rows(g(name), nc)
    deck.write_text("\n".join(lines) + "\n")
    return np_, nlevgrnd


def run_driver(driver: Path, deck: Path, results: Path, np_, nlevgrnd):
    try:
        subprocess.run([str(driver.resolve()), str(deck), str(results)], check=True)
    except FileNotFoundError as exc:
        raise SystemExit(f"cannot run {driver} -- build it first") from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"{driver} failed: {exc}") from exc
    vals = results.read_text().split()
    it = iter(vals[2:])                      # skip the np/nlevgrnd header
    btran = np.array([float(next(it)) for _ in range(np_)])
    rootr = np.array([[float(next(it)) for _ in range(nlevgrnd)] for _ in range(np_)])
    rresis = np.array([[float(next(it)) for _ in range(nlevgrnd)] for _ in range(np_)])
    return btran, rootr, rresis


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--driver", required=True, type=Path)
    ap.add_argument("--diag", required=True, type=Path)
    ap.add_argument("--steps", default="1,50,120,240")
    ap.add_argument("--rtol", type=float, default=1e-12)
    args = ap.parse_args()

    data = read_bin(args.diag)
    steps = [int(s) for s in args.steps.split(",")]

    worst_all = 0.0
    graded = 0
    with tempfile.TemporaryDirectory() as tmp:
        deck = Path(tmp) / "deck.txt"
        results = Path(tmp) / "results.txt"
        print(f"{'ts':>5} {'active':>7} {'btran':>12} {'rootr':>12} {'rresis':>12}")
        print("-" * 54)
        for ts in steps:
            if ts not in data or "rootstress_out:rootr" not in data[ts]:
                print(f"{ts:5d}   no rootstress records")
                continue
            np_, nlevgrnd = build_deck(data[ts], deck)
            btran, rootr, rresis = run_driver(args.driver, deck, results, np_, nlevgrnd)

            eb = data[ts]["rootstress_out:btran"]
            er = data[ts]["rootstress_out:rootr"]
            ee = data[ts]["rootstress_out:rresis"]

            # Grade only the patches ELM actually integrates. calc_root_moist_stress
            # runs over filter(nc)%num_nolakep; every other patch is never
            # touched and keeps a zero, which against a computed value reads as
            # a relative error of exactly 1.0 and looks like total disagreement.
            filt = None
            for cand in ("cantemp_in:filter_nolakep", "canflx_in:filter_nolakeurbanp"):
                if cand in data[ts]:
                    filt = data[ts][cand].astype(int) - 1   # -> 0-based
                    break
            if filt is None:
                raise SystemExit("no patch filter in the binary; cannot select "
                                 "the patches ELM integrates")
            itype = data[ts]["rootstress_in:patch_itype"].astype(int)
            act = np.array([p for p in filt if itype[p] != 0])
            if act.size == 0:
                print(f"{ts:5d}   no vegetated patches")
                continue

            def err(a, b):
                d = np.abs(a - b)
                return float((d / np.maximum(np.abs(b), 1.0)).max())

            wb = err(btran[act], eb[act])
            wr = err(rootr[act, :], er[act, :])
            we = err(rresis[act, :], ee[act, :])
            worst_all = max(worst_all, wb, wr, we)
            graded += 1
            print(f"{ts:5d} {act.size:7d} {wb:12.3e} {wr:12.3e} {we:12.3e}")

    print()
    if graded == 0:
        print("NOTHING GRADED -- no timestep carried rootstress records")
        return 1
    verdict = "PASS" if worst_all < args.rtol else "FAIL"
    print(f"{verdict}: worst relative error {worst_all:.3e} over {graded} timesteps "
          f"(rtol {args.rtol:g})")
    return 0 if worst_all < args.rtol else 1


if __name__ == "__main__":
    sys.exit(main())
