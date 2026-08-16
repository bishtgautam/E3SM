#!/usr/bin/env python3
"""Validate ELMxx's ground surface energy balance against ELM's own snapshots.

The kernel is replayed on ELM's inputs and compared against ELM's outputs, so
the only thing under test is the port. Nothing here depends on ELMxx having
run, and nothing depends on the two models' states agreeing -- which they do
not, and are not required to (see STATUS section J).

This is the same method as `validate_surface_albedo.py`, with one difference
that matters: the ground heat fluxes are **not** in a restart file. They are
intermediate quantities, computed and consumed inside SoilTemperature within a
single step. So the reference is ELM's instrumented diagnostic binary, which
records them under `soiltemp_mid:` -- and, being per-timestep, gives hundreds
of independent samples rather than one per restart.

  python3 components/elmxx/tools/validate_ground_heat_flux.py \
    --driver /path/to/ground_heat_flux_replay \
    --diag   /path/to/twin/run/elm_diagnostics.bin

Build the driver first:

  gfortran -O2 -o ground_heat_flux_replay \\
      share/util/shr_kind_mod.F90 \\
      components/elmxx/src/main/elmxxGroundHeatFluxKernelMod.F90 \\
      components/elmxx/tools/ground_heat_flux_replay.F90

Only the twin whose case enables the full diagnostic set carries the
`soiltemp_*` records -- on this machine that is `2x1_brazil`, not
`1x1_brazil`, whose binary holds forcing and canopy hydrology only.
"""

from __future__ import annotations

import argparse
import math
import re
import struct
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

# ELM landunit types this port covers: istsoil and istcrop. The urban branch
# of ComputeGroundHeatFluxAndDeriv uses eflx_lwrad_net and adds wasteheat, air
# conditioning and traffic; it is deliberately not ported, so comparing an
# urban column would grade a known gap rather than the code.
SOIL_LANDUNITS = (1, 2)

# Stefan-Boltzmann is READ from share/util/shr_const_mod.F90, not written here.
# Hardcoding the textbook 5.67e-8 instead of E3SM's 5.670374419e-8 is a 6.6e-5
# relative error, which lands as ~0.02 W/m2 on hs_soil -- small enough to look
# like a rounding difference and large enough to hide a real one. The port uses
# SHR_CONST_STEBOL, so the checker must use the same number, from the same file.
SHR_CONST = Path(__file__).resolve().parents[3] / "share" / "util" / "shr_const_mod.F90"

# Everything the kernel needs, by diagnostic label.
COLUMN_IN = ("soiltemp_in:emg", "soiltemp_in:htvp", "soiltemp_in:t_grnd",
             "soiltemp_in:t_h2osfc", "forcing:forc_lwrad")
PATCH_IN = ("surfrad_out:sabg_soil", "canflx_out:dlrad", "canflx_out:cgrnd",
            "canflx_out:eflx_sh_snow", "canflx_out:eflx_sh_soil",
            "canflx_out:eflx_sh_h2osfc", "canflx_out:qflx_ev_snow",
            "canflx_out:qflx_ev_soil", "canflx_out:qflx_ev_h2osfc")
TOPOLOGY = ("colmeta:lunit_itype", "soilflx_in:col_pfti", "soilflx_in:col_npfts",
            "soilflx_in:wtcol", "soilflx_in:frac_veg_nosno")
ARRAYS_2D = ("soiltemp_in:t_soisno", "surfrad_out:sabg_lyr")
INT_1D = ("soiltemp_in:snl",)
OUTPUTS = ("soiltemp_mid:hs_soil", "soiltemp_mid:hs_top_snow",
           "soiltemp_mid:hs_h2osfc", "soiltemp_mid:dhsdT")

WANTED = frozenset(COLUMN_IN + PATCH_IN + TOPOLOGY + ARRAYS_2D + INT_1D + OUTPUTS)


def stefan_boltzmann(source: Path) -> float:
    """SHR_CONST_STEBOL, out of E3SM's own constants module."""
    text = source.read_text(errors="replace")
    match = re.search(r"SHR_CONST_STEBOL\s*=\s*([0-9.eE+_\-]+?)_R8", text)
    if not match:
        raise SystemExit(f"{source}: cannot find SHR_CONST_STEBOL")
    return float(match.group(1).replace("_", ""))


def read_diagnostics(path: Path, wanted: frozenset):
    """One pass over the diagnostic binary, keeping only the labels asked for.

    Format is documented at the top of ELM's ElmDiagnostics.F90: an 8-byte
    magic and the label length, then records of
    [timestep][label][ndims][dims...][data]. Everything is little-endian
    stream, no Fortran record markers.

    The FIRST occurrence of a label within a timestep wins, matching the
    convention the Python harness in elmpy uses: the writer emits some fields
    twice within one routine, post-operation value first.
    """
    out: Dict[Tuple[int, str], object] = {}
    with path.open("rb") as handle:
        magic = handle.read(8)
        if magic != b"ELMDIAG1":
            raise SystemExit(f"{path}: not an ELM diagnostic binary")
        (label_len,) = struct.unpack("i", handle.read(4))
        while True:
            head = handle.read(4)
            if len(head) < 4:
                break
            (step,) = struct.unpack("i", head)
            label = handle.read(label_len).decode().strip()
            (ndims,) = struct.unpack("i", handle.read(4))
            keep = label in wanted
            if ndims == 0:
                payload = handle.read(8)
                value = struct.unpack("d", payload)[0] if keep else None
            elif ndims == 1:
                (n1,) = struct.unpack("i", handle.read(4))
                payload = handle.read(8 * n1)
                value = list(struct.unpack(f"{n1}d", payload)) if keep else None
            elif ndims == 2:
                n1, n2 = struct.unpack("ii", handle.read(8))
                payload = handle.read(8 * n1 * n2)
                # Fortran order: first index fastest.
                flat = struct.unpack(f"{n1 * n2}d", payload) if keep else None
                value = (n1, n2, flat) if keep else None
            elif ndims == -1:
                (n1,) = struct.unpack("i", handle.read(4))
                payload = handle.read(4 * n1)
                value = list(struct.unpack(f"{n1}i", payload)) if keep else None
            else:
                raise SystemExit(f"{path}: unknown ndims {ndims} for {label}")
            if keep and (step, label) not in out:
                out[(step, label)] = value
    return out


def element(block, i: int, j: int) -> float:
    """(i, j) of a 2-D record, 0-based, stored first-index-fastest."""
    n1, _n2, flat = block
    return flat[j * n1 + i]


def build_deck(rec, step: int, columns, patches_of, nlevsno, nlevtot, sb, deck: Path):
    """One ComputeGroundHeatFluxAndDeriv input set, for the soil columns."""
    def col(label):
        return rec[(step, label)]

    patches = [p for c in columns for p in patches_of[c]]
    place = {c: n + 1 for n, c in enumerate(columns)}    # 1-based for Fortran

    t_soisno = col("soiltemp_in:t_soisno")
    sabg_lyr = col("surfrad_out:sabg_lyr")
    nsnw = sabg_lyr[1]

    patch_col = [place[c] for c in columns for _ in patches_of[c]]
    wtcol = col("soilflx_in:wtcol")
    fvn = col("soilflx_in:frac_veg_nosno")

    def prow(label):
        values = col(label)
        return " ".join(repr(values[p]) for p in patches)

    lines = [f"{len(columns)} {len(patches)} {nlevsno} {nlevtot} {nsnw}",
             repr(sb),
             " ".join(str(v) for v in patch_col),
             " ".join(repr(wtcol[p]) for p in patches),
             " ".join(str(fvn[p]) for p in patches)]
    for label in ("soiltemp_in:emg", "soiltemp_in:htvp", "soiltemp_in:t_grnd",
                  "soiltemp_in:t_h2osfc"):
        lines.append(" ".join(repr(col(label)[c]) for c in columns))
    lines.append(" ".join(str(col("soiltemp_in:snl")[c]) for c in columns))
    lines.append(" ".join(repr(col("forcing:forc_lwrad")[c]) for c in columns))
    for j in range(nlevtot):
        lines.append(" ".join(repr(element(t_soisno, c, j)) for c in columns))
    for label in ("surfrad_out:sabg_soil", "canflx_out:dlrad", "canflx_out:cgrnd",
                  "canflx_out:eflx_sh_snow", "canflx_out:eflx_sh_soil",
                  "canflx_out:eflx_sh_h2osfc", "canflx_out:qflx_ev_snow",
                  "canflx_out:qflx_ev_soil", "canflx_out:qflx_ev_h2osfc"):
        lines.append(prow(label))
    for j in range(nsnw):
        lines.append(" ".join(repr(element(sabg_lyr, p, j)) for p in patches))
    deck.write_text("\n".join(lines) + "\n")
    return patches


def run_driver(driver: Path, deck: Path, results: Path):
    try:
        subprocess.run([str(driver.resolve()), str(deck), str(results)], check=True)
    except FileNotFoundError as exc:
        raise SystemExit(f"cannot run {driver} -- build it first") from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"{driver} failed: {exc}") from exc
    out: Dict[str, List[float]] = {}
    for line in results.read_text().splitlines():
        fields = line.split()
        if fields:
            out[fields[0]] = [float(v) for v in fields[1:]]
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--driver", required=True, type=Path)
    parser.add_argument("--diag", required=True, type=Path)
    parser.add_argument("--nlevsno", type=int, default=5)
    parser.add_argument("--tolerance", type=float, default=1.0e-10,
                        help="max absolute difference accepted, W/m2 (default 1e-10)")
    parser.add_argument("--steps", type=int, default=0,
                        help="compare only the first N timesteps (0 = all)")
    parser.add_argument("--show", type=int, default=6,
                        help="failing timesteps to print in full (default 6)")
    parser.add_argument("--shr-const", type=Path, default=SHR_CONST,
                        help="shr_const_mod.F90 to read SHR_CONST_STEBOL from")
    args = parser.parse_args()

    sb = stefan_boltzmann(args.shr_const)
    print(f"Stefan-Boltzmann {sb!r} from {args.shr_const.name}")
    print(f"reading {args.diag} ...")
    rec = read_diagnostics(args.diag, WANTED)
    if not rec:
        raise SystemExit(f"{args.diag}: none of the needed records are in this "
                         "binary -- is this the twin with full diagnostics on?")

    steps = sorted({step for (step, label) in rec
                    if label == "soiltemp_mid:hs_soil"})
    usable = [s for s in steps
              if all((s, label) in rec for label in WANTED - {"colmeta:lunit_itype"})]
    if args.steps:
        usable = usable[:args.steps]
    if not usable:
        raise SystemExit("no timestep carries the full input set")

    ltype = rec[(min(steps), "colmeta:lunit_itype")]
    columns = [i for i, t in enumerate(ltype) if t in SOIL_LANDUNITS]
    pfti = rec[(usable[0], "soilflx_in:col_pfti")]
    npfts = rec[(usable[0], "soilflx_in:col_npfts")]
    patches_of = {c: list(range(pfti[c] - 1, pfti[c] - 1 + npfts[c])) for c in columns}
    nlevtot = rec[(usable[0], "soiltemp_in:t_soisno")][1]

    npatch = sum(len(patches_of[c]) for c in columns)
    print(f"{len(steps)} timesteps recorded, {len(usable)} with a full input set")
    print(f"{len(columns)} soil columns of {len(ltype)}, {npatch} patches, "
          f"nlevtot {nlevtot}")

    # What the comparison actually exercised. A pass over a twin that never
    # grows snow says nothing about the snow branch, and a pass over bare
    # patches says nothing about the canopy term -- so measure both and print
    # them next to the result, where they cannot be overlooked.
    snl_seen = sorted({rec[(s, "soiltemp_in:snl")][c] for s in usable for c in columns})
    fvn_seen = sorted({rec[(s, "soilflx_in:frac_veg_nosno")][p]
                       for s in usable for c in columns for p in patches_of[c]
                       if rec[(s, "soilflx_in:wtcol")][p] != 0.0})
    print(f"coverage: snl {snl_seen}, frac_veg_nosno {fvn_seen} "
          f"(weighted patches only)")
    if snl_seen == [0]:
        print("          NO SNOW at any step -- lyr_top is always the first soil "
              "layer, so\n          hs_top_snow is identically hs_soil and the "
              "snow branch is untested.")
    if fvn_seen == [0]:
        print("          NO CANOPY on any weighted patch -- the (1-frac_veg_nosno) "
              "longwave\n          term is never tested at anything but 1.")
    elif fvn_seen == [1]:
        print("          ALL patches canopied -- the (1-frac_veg_nosno) longwave "
              "term is\n          never tested at anything but 0.")
    print()

    fields = {"soiltemp_mid:hs_soil": "hs_soil",
              "soiltemp_mid:hs_top_snow": "hs_top_snow",
              "soiltemp_mid:hs_h2osfc": "hs_h2osfc",
              "soiltemp_mid:dhsdT": "dhsdT"}
    worst = {name: (0.0, -1, -1, 0.0, 0.0) for name in fields.values()}
    failures, nonfinite = [], 0

    with tempfile.TemporaryDirectory() as scratch:
        deck = Path(scratch) / "deck.txt"
        results = Path(scratch) / "results.txt"
        for step in usable:
            build_deck(rec, step, columns, patches_of,
                       args.nlevsno, nlevtot, sb, deck)
            ported = run_driver(args.driver, deck, results)
            nonfinite += int(ported["n_nonfinite"][0])
            bad = 0.0
            for label, name in fields.items():
                reference = rec[(step, label)]
                for n, c in enumerate(columns):
                    diff = abs(reference[c] - ported[name][n])
                    bad = max(bad, diff)
                    if diff > worst[name][0]:
                        worst[name] = (diff, step, c, reference[c], ported[name][n])
            if bad > args.tolerance:
                failures.append((bad, step))

    print()
    print(f"{'field':<14s} {'worst |ELM-ELMxx|':>18s}   where")
    for name, (diff, step, c, ref, got) in sorted(worst.items(),
                                                  key=lambda kv: -kv[1][0]):
        print(f"{name:<14s} {diff:18.4e}   step {step:4d} col {c:3d}: "
              f"ELM {ref: .10f}  ELMxx {got: .10f}")
    if nonfinite:
        print(f"\n{nonfinite} patch-steps produced a non-finite flux")

    print()
    if failures:
        failures.sort(reverse=True)
        print(f"{len(failures)} of {len(usable)} timesteps exceed "
              f"{args.tolerance:g}; worst first:")
        for diff, step in failures[:args.show]:
            print(f"    step {step:4d}   {diff:.4e}")
        return 1
    print(f"{len(usable)} timesteps compared, all within {args.tolerance:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
