#!/usr/bin/env python3
"""Compare kELM h0 against the ELM twin's, field by field AND cell by cell.

compare_h0.py reads one value per field (the 1x1 twins). Multi-cell twins on
a global grid with a land mask (top4, nhcell4) need this one: ELM writes its
h0 either on the full (lat, lon) grid, fill outside the active cells, or --
with hist_dov2xy = .false., hist_type1d_pertape = 'GRID', which the top4 twin
uses to keep files small -- as a 1-D `gridcell` vector. kELM writes a compact
`lndgrid`. Cells are matched by their lat/lon.

Error is scaled by the field's annual peak IN THAT CELL (as compare_h0.py
does per field), so a field is only as good as its worst cell.

Two modes:

  compare_h0_cells.py KELM_GLOB ELM_GLOB
      kELM vs ELM.

  compare_h0_cells.py KELM_GLOB ELM_GLOB --base KELM_BASE_GLOB ELM_BASE_GLOB
      The RESPONSE to a switch: (kELM - kELM_base) vs (ELM - ELM_base),
      scaled by the ELM response's annual peak in that cell. This is how a
      parameterization (e.g. TOP) is graded independently of the twins'
      plane-parallel residual.

Needs netCDF4 (not in the system python here; use a conda env that has it).
"""
import argparse
import glob
from collections import defaultdict

import netCDF4
import numpy as np

CONSERVATION = {"ERRH2O", "ERRSEB", "ERRSOL", "ERRSOI", "ERRH2OSNO"}
SKIP = {"time", "time_bounds", "lon", "lat", "area", "mcdate", "mcsec",
        "mdcur", "mscur", "nstep", "date_written", "time_written", "levgrnd",
        "landfrac", "landmask", "pftmask", "topo", "edgen", "edgee", "edges",
        "edgew", "lonatm", "latatm", "lonrof", "latrof", "levlak", "levdcmp"}


def read_elm(path, cells):
    """{field: array(ncell)} for ELM's single-level fields, gridded or 1-D.

    `cells` holds (j, i) on a gridded tape, gridcell indices on a 1-D one."""
    out = {}
    with netCDF4.Dataset(path) as f:
        for name, v in f.variables.items():
            if name in SKIP:
                continue
            if v.dimensions == ("time", "lat", "lon"):
                a = np.ma.filled(v[0], np.nan)
                out[name] = np.array([a[j, i] for j, i in cells])
            elif v.dimensions == ("time", "gridcell"):
                a = np.ma.filled(v[0], np.nan)
                out[name] = a[np.array(cells)]
    return out


def read_kelm(path, order):
    out = {}
    with netCDF4.Dataset(path) as f:
        for name, v in f.variables.items():
            if name in SKIP or v.dimensions != ("time", "lndgrid"):
                continue
            a = np.ma.filled(v[0], np.nan)
            out[name] = a[order]
    return out


def cell_map(kelm_file, elm_file):
    """ELM (j, i) and kELM lndgrid index of each kELM cell, by lat/lon."""
    with netCDF4.Dataset(elm_file) as e:
        one_d = "grid1d_lat" in e.variables
        if one_d:
            glat, glon = e["grid1d_lat"][:], e["grid1d_lon"][:] % 360.0
        else:
            elat, elon = e["lat"][:], e["lon"][:]
    with netCDF4.Dataset(kelm_file) as k:
        klat, klon = k["lat"][:], k["lon"][:]
    cells, order, labels = [], [], []
    for n in range(len(klat)):
        lon = klon[n] % 360.0
        if one_d:
            g = int(np.argmin(np.abs(glat - klat[n]) + np.abs(glon - lon)))
            cells.append(g)
        else:
            j = int(np.argmin(np.abs(elat - klat[n])))
            i = int(np.argmin(np.abs((elon % 360.0) - lon)))
            cells.append((j, i))
        order.append(n)
        labels.append(f"{klat[n]:.2f}N {klon[n]:.2f}E")
    return cells, np.array(order), labels


def series(kglob, eglob, cells, order):
    ks, es = sorted(glob.glob(kglob)), sorted(glob.glob(eglob))
    n = min(len(ks), len(es))
    K = [read_kelm(f, order) for f in ks[:n]]
    E = [read_elm(f, cells) for f in es[:n]]
    return K, E


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("kelm"); ap.add_argument("elm")
    ap.add_argument("--base", nargs=2, metavar=("KELM_BASE", "ELM_BASE"))
    ap.add_argument("--rtol", type=float, default=1e-3)
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--fields", default="", help="comma-separated: print month by month")
    a = ap.parse_args()

    k0 = sorted(glob.glob(a.kelm))[0]
    e0 = sorted(glob.glob(a.elm))[0]
    cells, order, labels = cell_map(k0, e0)
    K, E = series(a.kelm, a.elm, cells, order)
    if a.base:
        KB, EB = series(a.base[0], a.base[1], cells, order)
        n = min(len(K), len(KB), len(E), len(EB))
        K = [{f: K[m][f] - KB[m][f] for f in K[m] if f in KB[m]} for m in range(n)]
        E = [{f: E[m][f] - EB[m][f] for f in E[m] if f in EB[m]} for m in range(n)]
    n = min(len(K), len(E))
    print(f"comparing {n} months, {len(cells)} cells: " + "; ".join(labels))
    if a.base:
        print("mode: RESPONSE (run - base), kELM vs ELM")
    print()

    fields = sorted(set.intersection(*[set(K[m]) & set(E[m]) for m in range(n)]))
    worst = {}
    mag = defaultdict(float)
    for f in fields:
        x = np.array([K[m][f] for m in range(n)])   # (month, cell)
        e = np.array([E[m][f] for m in range(n)])
        if f in CONSERVATION:
            mag[f] = np.nanmax(np.abs(x)) if np.isfinite(x).any() else 0.0
            continue
        both = np.isfinite(x) & np.isfinite(e)
        if not both.any():
            continue
        scale = np.nanmax(np.where(both, np.abs(e), np.nan), axis=0)   # per cell
        scale = np.where(scale > 0, scale, 1e-30)
        r = np.where(both, np.abs(x - e) / scale, 0.0)
        m, c = np.unravel_index(int(np.argmax(r)), r.shape)
        worst[f] = (float(r[m, c]), m + 1, c, e[m, c], x[m, c])

    good = [f for f, w in worst.items() if w[0] < a.rtol]
    print(f"fields under {a.rtol:g} (scaled, every month, every cell): {len(good)} of {len(worst)}\n")
    if mag:
        print("conservation fields, worst kELM magnitude (want ~0):")
        for f in sorted(mag):
            print(f"  {f:<12} {mag[f]:.3e}")
        print()
    print(f"worst {a.top} fields, error / that cell's annual peak:")
    print(f"  {'field':<16}{'err/scale':>10}  {'mon':>3}  {'cell':<18}"
          f"{'ELM':>13}{'kELM':>13}{'kELM-ELM':>13}{'rel diff':>11}")
    for f, (r, m, c, ev, xv) in sorted(worst.items(), key=lambda kv: -kv[1][0])[:a.top]:
        print(f"  {f:<16}{r:>10.3e}  {m:>3d}  {labels[c]:<18}"
              f"{ev:>13.5g}{xv:>13.5g}{xv - ev:>13.4e}{reldiff(xv, ev):>11}")

    for f in [s for s in a.fields.split(",") if s]:
        if f not in fields:
            print(f"\n{f}: not present"); continue
        print(f"\n{f} by month" + (" (response: run - base)" if a.base else "")
              + "; rel diff = (kELM - ELM) / |ELM|")
        for c in range(len(cells)):
            print(f"  cell {labels[c]}")
            print(f"    {'mon':>3}{'ELM':>14}{'kELM':>14}{'kELM-ELM':>13}{'rel diff':>11}")
            for m in range(n):
                ev, xv = E[m][f][c], K[m][f][c]
                print(f"    {m+1:>3}{ev:>14.6g}{xv:>14.6g}{xv - ev:>13.4e}{reldiff(xv, ev):>11}")


def reldiff(x, e):
    """(x - e) / |e| as text; '-' where e is zero and x is too, 'inf' otherwise."""
    if e == 0.0:
        return "-" if x == 0.0 else "inf"
    return f"{(x - e) / abs(e):.3e}"

if __name__ == "__main__":
    main()
