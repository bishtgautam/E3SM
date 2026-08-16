#!/usr/bin/env python3
"""Validate ELMxx's surface albedo kernel against ELM's own restarts.

The kernel is replayed on ELM's inputs and compared against ELM's outputs, so
the only thing under test is the port. Nothing here depends on ELMxx having
run, and nothing depends on the two models' states agreeing -- which they do
not, and are not required to (see STATUS section J).

  python3 components/elmxx/tools/validate_surface_albedo.py \
    --driver      /path/to/surface_albedo_replay \
    --elm-log     /path/to/twin/run/lnd.log.NNNNNN \
    --surfdata    /path/to/surfdata.nc \
    --params      /path/to/clm_params.nc \
    --elm-restart /path/to/twin/run/'*.elm.r.*.nc'

Build the driver first:

  gfortran -O2 -o surface_albedo_replay \\
      share/util/shr_kind_mod.F90 \\
      components/elmxx/src/main/elmxxSurfaceAlbedoKernelMod.F90 \\
      components/elmxx/tools/surface_albedo_replay.F90

Night restarts are graded too, and deliberately. At night ELM does no solar
calculation, but what it leaves behind is not arbitrary -- albd = albi = 1 with
zero transmittance is its marker for "nothing was computed here", tlai_z still
carries elai because canopy layering is not gated on the sun, and all of it is
what the coupler and the next step read. A port that invents its own night
convention diverges silently. Restarts are reported one per line with their
sun angle, so a failure names the geometry it happened at.
"""

from __future__ import annotations

import argparse
import glob
import math
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Sequence

NUMBER = re.compile(r"[-+]?(?:\d+\.?(?:\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")

# ELM landunit types this port covers: istsoil and istcrop. Everything else --
# lake, wetland, glacier, urban -- goes down a SoilAlbedo branch that is
# deliberately not ported, so comparing it would grade a known gap.
SOIL_LANDUNITS = (1, 2)

DENH2O = 1000.0
DENICE = 917.0
NLEVSNO = 5

# ELM restart name -> replay output name. Column fields first, then patch.
COLUMN_FIELDS = ("albsod", "albsoi", "albgrd", "albgri")
PATCH_BANDED = ("albd", "albi", "fabd", "fabi", "ftdd", "ftid", "ftii")
PATCH_SCALAR = ("tlai_z", "fsun_z", "fabd_sun_z", "fabd_sha_z",
                "fabi_sun_z", "fabi_sha_z", "vcmaxcintsun", "vcmaxcintsha")


def cdl_values(filename: Path, names: Sequence[str]) -> Dict[str, List[float]]:
    """Read selected variables via the standard ncdump utility."""
    try:
        output = subprocess.check_output(
            ["ncdump", "-v", ",".join(names), str(filename)], text=True
        )
    except FileNotFoundError as exc:
        raise SystemExit("ncdump is required to read a NetCDF file") from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"cannot read NetCDF file {filename}: {exc}") from exc

    data = output.split("data:", 1)
    if len(data) != 2:
        raise SystemExit(f"ncdump produced no data section for {filename}")
    result: Dict[str, List[float]] = {}
    for name in names:
        match = re.search(rf"^\s*{re.escape(name)}\s*=\s*(.*?);", data[1],
                          re.S | re.M)
        if not match:
            raise SystemExit(f"{filename} is missing {name}")
        result[name] = cdl_numbers(match.group(1))
    return result


def cdl_numbers(text: str) -> List[float]:
    """Split one CDL value list, keeping fill values in place.

    ncdump prints an unwritten element as a bare underscore, and an ELM
    restart is full of them: H2OSOI_LIQ is (column, levtot) and the five snow
    slots of a snow-free column are never written. Scanning for numbers and
    ignoring everything else therefore does not lose one value, it SHIFTS
    every value after it -- soil layer 1 reads as layer 6, silently, with a
    perfectly plausible number. Split on the separator instead and let a fill
    become NaN, which is loud.
    """
    values: List[float] = []
    for token in text.split(","):
        token = token.strip()
        if not token:
            continue
        if token == "_":
            values.append(float("nan"))
        else:
            values.append(float(token.replace("D", "E").replace("d", "e")))
    return values


def soil_layer_thickness(elm_log: Path) -> float:
    """ELM prints the vertical grid; layer 1 thickness is all this needs.

    Read rather than hardcoded: the grid is a function of nlevgrnd and of
    initVerticalMod's formula, and a tool that silently assumed one would stop
    being true the first time either changed.
    """
    for line in elm_log.read_text(errors="replace").splitlines():
        if "dzsoi:" in line:
            values = NUMBER.findall(line.split("dzsoi:", 1)[1])
            if values:
                return float(values[0].replace("D", "E"))
    raise SystemExit(f"{elm_log}: no 'dzsoi:' line -- is this an ELM lnd.log?")


def pft_parameters(params: Path):
    """Canopy optics, per PFT and band, from ELM's parameter file."""
    names = ("rholvis", "rholnir", "rhosvis", "rhosnir",
             "taulvis", "taulnir", "tausvis", "tausnir", "xl")
    raw = cdl_values(params, names)
    npft = len(raw["xl"])
    for name, values in raw.items():
        if len(values) != npft:
            raise SystemExit(f"{params}: {name} has {len(values)}, expected {npft}")
    return npft, raw


def build_deck(restart: Path, surfdata: Path, npft: int, params,
               dz1: float, deck: Path):
    """Extract one SurfaceAlbedo input set from an ELM restart."""
    topo = cdl_values(restart, (
        "cols1d_ityplun", "cols1d_gridcell_index",
        "pfts1d_ityplun", "pfts1d_column_index", "pfts1d_itypveg",
        "coszen", "frac_sno", "H2OSOI_LIQ", "H2OSOI_ICE",
        "elai", "esai", "T_VEG", "FWET",
    ))
    colour = cdl_values(surfdata, ("SOIL_COLOR",))["SOIL_COLOR"]

    col_ltype = [int(v) for v in topo["cols1d_ityplun"]]
    col_grid = [int(v) for v in topo["cols1d_gridcell_index"]]
    ncol_all = len(col_ltype)
    nlevtot = len(topo["H2OSOI_LIQ"]) // ncol_all

    # Columns this port covers, in restart order, and their 1-based position.
    columns = [i for i, lt in enumerate(col_ltype) if lt in SOIL_LANDUNITS]
    position = {c: n + 1 for n, c in enumerate(columns)}

    patches = [i for i, lt in enumerate(int(v) for v in topo["pfts1d_ityplun"])
               if lt in SOIL_LANDUNITS]

    nc, np_ = len(columns), len(patches)
    if nc == 0 or np_ == 0:
        raise SystemExit(f"{restart.name}: no soil column or patch to compare")

    coszen = [topo["coszen"][c] for c in columns]
    fracsno = [topo["frac_sno"][c] for c in columns]
    # SOIL_COLOR is per gridcell; cols1d_gridcell_index is 1-based.
    soilcolor = [int(colour[col_grid[c] - 1]) for c in columns]
    # h2osoi_vol is diagnostic and not on the restart, so rebuild it exactly
    # as HydrologyNoDrainage does, for soil layer 1 (slot nlevsno+1 of levtot).
    top = NLEVSNO
    h2osoi = [topo["H2OSOI_LIQ"][c * nlevtot + top] / (dz1 * DENH2O)
              + topo["H2OSOI_ICE"][c * nlevtot + top] / (dz1 * DENICE)
              for c in columns]

    patch_col = [position[int(topo["pfts1d_column_index"][p]) - 1]
                 for p in patches]
    patch_ivt = [int(topo["pfts1d_itypveg"][p]) for p in patches]
    elai = [topo["elai"][p] for p in patches]
    esai = [topo["esai"][p] for p in patches]
    tveg = [topo["T_VEG"][p] for p in patches]
    fwet = [topo["FWET"][p] for p in patches]

    def row(values):
        return " ".join(repr(v) for v in values)

    lines = [f"{nc} {np_} {npft}", row(patch_col), row(coszen), row(soilcolor),
             row(h2osoi), row(fracsno), row(patch_ivt), row(elai), row(esai),
             row(tveg), row(fwet)]
    for name in ("rhol", "rhos", "taul", "taus"):
        lines.append(row(params[name + "vis"]))
        lines.append(row(params[name + "nir"]))
    lines.append(row(params["xl"]))
    deck.write_text("\n".join(lines) + "\n")

    return columns, patches, coszen


def elm_reference(restart: Path, columns, patches):
    """ELM's own SurfaceAlbedo outputs, on the same selection."""
    names = list(COLUMN_FIELDS) + list(PATCH_BANDED) + list(PATCH_SCALAR)
    raw = cdl_values(restart, names)
    out: Dict[str, List[float]] = {}
    for name in COLUMN_FIELDS:
        for band, offset in (("vis", 0), ("nir", 1)):
            out[f"{name}_{band}"] = [raw[name][c * 2 + offset] for c in columns]
    for name in PATCH_BANDED:
        for band, offset in (("vis", 0), ("nir", 1)):
            out[f"{name}_{band}"] = [raw[name][p * 2 + offset] for p in patches]
    for name in PATCH_SCALAR:          # levcan = 1, so one value per patch
        out[name] = [raw[name][p] for p in patches]
    return out


def run_driver(driver: Path, deck: Path, results: Path):
    try:
        subprocess.run([str(driver), str(deck), str(results)], check=True)
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


def compare(elm, elmxx, coszen, patch_col):
    """Largest absolute difference per field, and where it happened."""
    report = []
    for name, reference in sorted(elm.items()):
        if name not in elmxx:
            continue
        ported = elmxx[name]
        if len(ported) != len(reference):
            raise SystemExit(f"{name}: {len(ported)} values, expected {len(reference)}")
        worst, where = 0.0, -1
        for i, (a, b) in enumerate(zip(reference, ported)):
            if math.isnan(b) or math.isinf(b):
                worst, where = float("inf"), i
                break
            if abs(a - b) > worst:
                worst, where = abs(a - b), i
        report.append((worst, name, where,
                       reference[where] if where >= 0 else 0.0,
                       ported[where] if where >= 0 else 0.0))
    report.sort(reverse=True, key=lambda row: row[0])
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--driver", required=True, type=Path)
    parser.add_argument("--elm-log", required=True, type=Path)
    parser.add_argument("--surfdata", required=True, type=Path)
    parser.add_argument("--params", required=True, type=Path)
    parser.add_argument("--elm-restart", required=True, nargs="+")
    parser.add_argument("--tolerance", type=float, default=1.0e-12,
                        help="max absolute difference accepted (default 1e-12)")
    parser.add_argument("--top", type=int, default=8,
                        help="fields to print per restart (default 8)")
    args = parser.parse_args()

    restarts = sorted({path for pattern in args.elm_restart
                       for path in glob.glob(pattern)})
    if not restarts:
        raise SystemExit("no ELM restart matched")

    dz1 = soil_layer_thickness(args.elm_log)
    npft, params = pft_parameters(args.params)
    print(f"soil layer 1 thickness {dz1:.10e} m, {npft} PFTs on the parameter file")

    failures, daytime = 0, 0
    with tempfile.TemporaryDirectory() as scratch:
        deck = Path(scratch) / "deck.txt"
        results = Path(scratch) / "results.txt"
        for path in restarts:
            restart = Path(path)
            columns, patches, coszen = build_deck(restart, args.surfdata, npft,
                                                  params, dz1, deck)
            if max(coszen) > 0.0:
                daytime += 1
            elmxx = run_driver(args.driver, deck, results)
            elm = elm_reference(restart, columns, patches)
            report = compare(elm, elmxx, coszen, None)
            stamp = restart.name.split(".elm.r.")[-1].removesuffix(".nc")
            worst = report[0][0] if report else 0.0
            verdict = "PASS" if worst <= args.tolerance else "FAIL"
            print(f"\n{stamp}  coszen {max(coszen):.6f}  "
                  f"{len(columns)} columns, {len(patches)} patches  "
                  f"worst {worst:.3e}  {verdict}")
            for value, name, index, reference, ported in report[:args.top]:
                print(f"    {name:<14s} {value:12.4e}  "
                      f"at {index:3d}: ELM {reference: .10f}  ELMxx {ported: .10f}")
            if verdict == "FAIL":
                failures += 1

    if daytime == 0:
        print("\nNo daytime restart in the set -- the two-stream was never "
              "exercised, so this run graded almost nothing.", file=sys.stderr)
        return 2
    print(f"\n{len(restarts)} restarts compared ({daytime} daytime), "
          f"{failures} failing at tolerance {args.tolerance:g}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
