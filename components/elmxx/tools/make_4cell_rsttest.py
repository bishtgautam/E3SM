#!/usr/bin/env python3
"""Build the 4-cell domain and surface dataset for the multi-rank tests.

A grid of identical cells cannot catch a cell permutation -- every reorder
compares equal -- so the four cells sit at four coordinates and get different
DATM forcing. Built from the 2x1_brazil files by CDL round trip (ncdump/ncgen;
no NCO or netCDF4 needed):

  * cells 3, 4 repeat cells 1, 2 in every variable whose last dim is the cell
    dim; the domain then gets distinct xc/yc (and LONGXY/LATIXY to match);
  * PCT_URBAN = 0 and PCT_NATVEG = 100 everywhere -- h0 requires a natural
    landunit weight of 1 (elmxx_hist_check_natveg_weight);
  * URBAN_REGION_ID = 0,3,3,3: urban landunits are allocated wherever the
    region is valid, whatever their weight, so one rank of a 2-rank run holds
    none and the packed (density class, cell) urban order differs from the
    restart file's cell-major order.

Usage:
    make_4cell_rsttest.py <2x1_brazil dir> <output dir>
writes domain_4x1_rsttest.nc and surfdata_4x1_rsttest.nc. See CLAUDE.md,
"Restarts", for the test itself.
"""
import os
import re
import subprocess
import sys

LONS = ['-55', '-50', '-60', '-45']
LATS = ['-7', '-10', '-3', '-15']


def vals(body):
    return [v.strip() for v in body.replace('\n', ' ').split(',') if v.strip()]


def fmt(v):
    return '\n  ' + ', '.join(v)


def setv(s, name, newvals):
    s, n = re.subn(r'^ %s =.*?;' % re.escape(name), ' %s =%s ;' % (name, fmt(newvals)),
                   s, flags=re.S | re.M)
    if n != 1:
        sys.exit(f'{name}: expected one data block, found {n}')
    return s


def grow(cdl, celldim):
    """2 cells -> 4, repeating (1, 2) for every variable ending in celldim."""
    hdr, data = cdl.split('\ndata:\n', 1)
    hdr, n = re.subn(r'(\n\s*%s = )2 ;' % celldim, r'\g<1>4 ;', hdr)
    if n != 1:
        sys.exit(f'{celldim} = 2 not found in the header')
    dims = {m.group(2): m.group(3)
            for m in re.finditer(r'\n\s*(\w+) (\w+)\(([^)]*)\) ;', hdr)}

    def repl(m):
        name, body = m.group(1), m.group(2)
        d = [x.strip() for x in dims.get(name.replace('\\', ''), '').split(',')]
        if d[-1] != celldim or body.strip().startswith('"'):
            return m.group(0)
        v = vals(body)
        out = []
        for i in range(0, len(v), 2):
            out += [v[i], v[i + 1], v[i], v[i + 1]]
        return ' %s =%s ;' % (name, fmt(out))

    data = re.sub(r'^ (\S+) =(.*?);\s*$', lambda m: repl(m) + '\n', data,
                  flags=re.S | re.M)
    return hdr + '\ndata:\n' + data


def ncdump(path):
    return subprocess.run(['ncdump', path], check=True, capture_output=True,
                          text=True).stdout


def ncgen(cdl, path):
    tmp = path + '.cdl'
    open(tmp, 'w').write(cdl)
    subprocess.run(['ncgen', '-k', '2', '-o', path, tmp], check=True)
    os.remove(tmp)


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = sys.argv[1:]

    dom = grow(ncdump(os.path.join(src, 'domain_2x1_brazil_c260222.nc')), 'ni')
    dom = setv(dom, 'xc', LONS)
    dom = setv(dom, 'yc', LATS)
    xv, yv = [], []
    for lo, la in zip(map(float, LONS), map(float, LATS)):
        xv += [str(lo - 0.5), str(lo + 0.5), str(lo + 0.5), str(lo - 0.5)]
        yv += [str(la - 0.5), str(la - 0.5), str(la + 0.5), str(la + 0.5)]
    dom = setv(dom, 'xv', xv)
    dom = setv(dom, 'yv', yv)
    ncgen(dom, os.path.join(out, 'domain_4x1_rsttest.nc'))

    sd = grow(ncdump(os.path.join(src, 'surfdata_2x1_brazil_c260222.nc')), 'gridcell')
    sd = setv(sd, 'PCT_NATVEG', ['100'] * 4)
    sd = setv(sd, 'PCT_URBAN', ['0'] * 12)
    sd = setv(sd, 'URBAN_REGION_ID', ['0', '3', '3', '3'])
    sd = setv(sd, 'LONGXY', [str(float(x) % 360) for x in LONS])
    sd = setv(sd, 'LATIXY', LATS)
    ncgen(sd, os.path.join(out, 'surfdata_4x1_rsttest.nc'))
    print(f'wrote domain_4x1_rsttest.nc and surfdata_4x1_rsttest.nc in {out}')


if __name__ == '__main__':
    main()
