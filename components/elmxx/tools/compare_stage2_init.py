#!/usr/bin/env python3
"""Compare ELMxx's Stage 2 snapshot against an ELM twin restart.

Run the ELM and ELMxx twins for the same short interval, then invoke this
script from the ELMxx run directory (or pass its rank-local snapshots):

  python3 components/elmxx/tools/compare_stage2_init.py \
    --elm-restart /path/to/twin.elm.r.YYYY-MM-DD-SSSSS.nc \
    --surfdata /path/to/shared-surfdata.nc \
    --elmxx-snapshot /path/to/elmxx_stage2_rank*.dat

The ELM restart is the ELM-side extraction: it contains the resolved subgrid
and the monthly-interpolated TLai/TSai/height fields. The ELMxx snapshots also
carry the full per-column soil state, making the state used by this comparison
auditable in the same artifact.
"""

from __future__ import annotations

import argparse
import glob
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


NUMBER = re.compile(r"[-+]?(?:\d+\.?(?:\d*)?|\.\d+)(?:[EeDd][-+]?\d+)?")


def cdl_values(filename: Path, names: Sequence[str]) -> Dict[str, List[float]]:
    """Read selected variables via the standard ncdump utility."""
    try:
        output = subprocess.check_output(
            ["ncdump", "-v", ",".join(names), str(filename)], text=True
        )
    except FileNotFoundError as exc:
        raise SystemExit("ncdump is required to read an ELM restart") from exc
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"cannot read NetCDF file {filename}: {exc}") from exc

    data = output.split("data:", 1)
    if len(data) != 2:
        raise SystemExit(f"ncdump produced no data section for {filename}")
    result: Dict[str, List[float]] = {}
    for name in names:
        match = re.search(rf"\b{re.escape(name)}\s*=\s*(.*?);", data[1], re.S)
        if not match:
            raise SystemExit(f"{filename} is missing {name}")
        result[name] = [float(value.replace("D", "E").replace("d", "e"))
                        for value in NUMBER.findall(match.group(1))]
    return result


def records(paths: Iterable[str]):
    """Load all rank-local ELMxx snapshot records."""
    land: Dict[Tuple[int, int], Tuple[float, ...]] = {}
    column: Dict[Tuple[int, int, int], Tuple[float, ...]] = {}
    patch: Dict[Tuple[int, int, int, int], Tuple[float, ...]] = {}
    dates = set()
    for pattern in paths:
        matched = glob.glob(pattern)
        if not matched:
            raise SystemExit(f"no ELMxx snapshot matches {pattern}")
        for filename in matched:
            for line in Path(filename).read_text().splitlines():
                fields = line.split()
                if not fields:
                    continue
                if fields[0] == "ELMXX_STAGE2":
                    dates.add(tuple(int(value) for value in fields[2:4]))
                elif fields[0] == "L":
                    key = (int(fields[1]), int(fields[2]))
                    land[key] = (float(fields[3]),)
                elif fields[0] == "C":
                    key = tuple(int(value) for value in fields[1:4])
                    column[key] = tuple(float(value) for value in fields[4:])
                elif fields[0] == "P":
                    key = tuple(int(value) for value in fields[1:5])
                    patch[key] = tuple(float(value) for value in fields[5:])
    if len(dates) != 1:
        raise SystemExit(f"ELMxx snapshots do not have one common month/day: {sorted(dates)}")
    return land, column, patch, dates.pop()


def assert_same_keys(label: str, actual: Dict, expected: Dict) -> None:
    if actual.keys() != expected.keys():
        only_actual = sorted(actual.keys() - expected.keys())[:5]
        only_expected = sorted(expected.keys() - actual.keys())[:5]
        raise AssertionError(
            f"{label} topology differs: only ELMxx={only_actual}; only ELM={only_expected}"
        )


def assert_values(label: str, actual: Dict, expected: Dict, atol: float) -> int:
    checked = 0
    for key in sorted(expected):
        avalue, evalue = actual[key], expected[key]
        if len(avalue) != len(evalue):
            raise AssertionError(f"{label} {key}: field count differs")
        for position, (a, e) in enumerate(zip(avalue, evalue), start=1):
            if abs(a - e) > atol:
                raise AssertionError(
                    f"{label} {key}, value {position}: ELMxx={a:.17g}, ELM={e:.17g}, "
                    f"abs diff={abs(a - e):.3e}"
                )
            checked += 1
    return checked


def surfdata_dimensions(surfdata: Path) -> Tuple[int, int]:
    """Return (nlevsoi, flattened global gridcell count) from a surfdata file."""
    try:
        header = subprocess.check_output(["ncdump", "-h", str(surfdata)], text=True)
    except subprocess.CalledProcessError as exc:
        raise SystemExit(f"cannot read surfdata header {surfdata}: {exc}") from exc
    dims = {name: int(value) for name, value in re.findall(r"^\s*(\w+)\s*=\s*(\d+)\s*;", header, re.M)}
    if "nlevsoi" not in dims:
        raise SystemExit(f"{surfdata} has no nlevsoi dimension")
    if "gridcell" in dims:
        ngrid = dims["gridcell"]
    elif "lsmlon" in dims and "lsmlat" in dims:
        ngrid = dims["lsmlon"] * dims["lsmlat"]
    else:
        raise SystemExit(f"{surfdata} has neither gridcell nor lsmlon/lsmlat dimensions")
    return dims["nlevsoi"], ngrid


def check_source_soil(columns: Dict, surfdata: Path, atol: float) -> int:
    """Check that each ELMxx column retains the shared surfdata soil row."""
    nlevsoi, ngrid = surfdata_dimensions(surfdata)
    source = cdl_values(surfdata, ("PCT_SAND", "PCT_CLAY", "ORGANIC", "SOIL_COLOR"))
    expected_size = nlevsoi * ngrid
    for name in ("PCT_SAND", "PCT_CLAY", "ORGANIC"):
        if len(source[name]) != expected_size:
            raise SystemExit(f"{surfdata}: unexpected {name} size {len(source[name])}")
    if len(source["SOIL_COLOR"]) != ngrid:
        raise SystemExit(f"{surfdata}: unexpected SOIL_COLOR size {len(source['SOIL_COLOR'])}")

    checked = 0
    for (gridcell, _ltype, _ctype), values in columns.items():
        if gridcell < 1 or gridcell > ngrid:
            raise AssertionError(f"column gridcell {gridcell} is outside surfdata")
        if len(values) != 3 * nlevsoi + 2:
            raise AssertionError(f"column {gridcell}: unexpected soil payload length")
        offset = gridcell - 1
        actual = values[1:]
        expected = (
            tuple(source["PCT_SAND"][lev * ngrid + offset] for lev in range(nlevsoi))
            + tuple(source["PCT_CLAY"][lev * ngrid + offset] for lev in range(nlevsoi))
            + tuple(source["ORGANIC"][lev * ngrid + offset] for lev in range(nlevsoi))
            + (source["SOIL_COLOR"][offset],)
        )
        for position, (a, e) in enumerate(zip(actual, expected), start=1):
            if abs(a - e) > atol:
                raise AssertionError(
                    f"source soil {(gridcell, _ltype, _ctype)}, value {position}: "
                    f"ELMxx={a:.17g}, surfdata={e:.17g}, abs diff={abs(a - e):.3e}"
                )
            checked += 1
    return checked


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--elm-restart", required=True, type=Path)
    parser.add_argument("--surfdata", required=True, type=Path)
    parser.add_argument("--elmxx-snapshot", required=True, action="append",
                        help="snapshot path or glob; repeatable")
    parser.add_argument("--atol", type=float, default=1.0e-10,
                        help="absolute tolerance for ELM floating values")
    args = parser.parse_args()

    elmxx_land, elmxx_col, elmxx_patch, date = records(args.elmxx_snapshot)
    names = (
        "land1d_gridcell_index", "land1d_ityplun", "land1d_wtxy",
        "cols1d_gridcell_index", "cols1d_ityplun", "cols1d_ityp", "cols1d_wtlnd",
        "pfts1d_gridcell_index", "pfts1d_ityplun", "pfts1d_itypcol",
        "pfts1d_itypveg", "pfts1d_wtcol", "tlai", "tsai", "htop", "hbot",
    )
    elm = cdl_values(args.elm_restart, names)

    elm_land = {
        (int(g), int(t)): (weight,)
        for g, t, weight in zip(elm["land1d_gridcell_index"], elm["land1d_ityplun"],
                                elm["land1d_wtxy"])
    }
    elm_col = {
        (int(g), int(lt), int(ct)): (weight,)
        for g, lt, ct, weight in zip(elm["cols1d_gridcell_index"], elm["cols1d_ityplun"],
                                     elm["cols1d_ityp"], elm["cols1d_wtlnd"])
    }
    elm_patch = {
        (int(g), int(lt), int(ct), int(pt)): (weight, lai, sai, top, bot)
        for g, lt, ct, pt, weight, lai, sai, top, bot in zip(
            elm["pfts1d_gridcell_index"], elm["pfts1d_ityplun"],
            elm["pfts1d_itypcol"], elm["pfts1d_itypveg"], elm["pfts1d_wtcol"],
            elm["tlai"], elm["tsai"], elm["htop"], elm["hbot"]
        )
    }

    try:
        assert_same_keys("landunit", elmxx_land, elm_land)
        assert_same_keys("column", elmxx_col, elm_col)
        assert_same_keys("patch", elmxx_patch, elm_patch)
        checked = assert_values("landunit weight", elmxx_land, elm_land, args.atol)
        checked += assert_values(
            "column weight", {key: (values[0],) for key, values in elmxx_col.items()},
            elm_col, args.atol
        )
        checked += assert_values("patch state", elmxx_patch, elm_patch, args.atol)
        checked += check_source_soil(elmxx_col, args.surfdata, args.atol)
    except AssertionError as exc:
        print(f"Stage 2 initialization comparison FAILED: {exc}", file=sys.stderr)
        return 1

    print(
        "Stage 2 initialization comparison PASSED: "
        f"{len(elm_land)} landunits, {len(elm_col)} columns, {len(elm_patch)} patches; "
        f"{checked} exact-to-tolerance values; ELMxx month/day={date[0]}/{date[1]}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
