#!/usr/bin/env python3
"""Build the 4-cell glacier domain for the glacier twins and layout test
(glacier plan, G-D3 items 2 and 3).

The nhcell_urban domain (full f19) with four f19 cells active, read from the
standard global surfdata, chosen so the glacier path meets real cells:

  (83,126)  67.3N 315E    Greenland: 100% glacier (urban region 12, 0 weight)
  (79,87)   59.7N 217.5E  Alaska: 67% glacier, 32% natural, 1% lake
  (68,29)   38.8N 72.5E   Pamir: 17% glacier, 82% natural, 1% lake
  (72,114)  46.4N 285E    nhcell: natural, 7% lake, urban; no glacier

None has crop or wetland (kELM aborts on either). Two ranks split them two
and two, so the packed glacier order on a rank differs from the restart's.
mask and frac are set to 1 on the four cells; everything else is the
nhcell file's. CDL round trip, no netCDF4 needed.

Usage:
    make_glacier4_test.py <nhcell domain.nc> <output.nc>
"""
import re
import subprocess
import sys

NI = 144
CELLS = [(83, 126), (79, 87), (68, 29), (72, 114)]   # (j, i), 0-based


def setcells(cdl, name, value):
    hdr, data = cdl.split('\ndata:\n', 1)
    m = re.search(r'^ %s =(.*?);' % name, data, re.S | re.M)
    v = [x.strip() for x in m.group(1).replace('\n', ' ').split(',')]
    for j, i in CELLS:
        v[j * NI + i] = value
    body = ',\n  '.join(', '.join(v[k:k + 12]) for k in range(0, len(v), 12))
    data = data[:m.start(1)] + '\n  ' + body + ' ' + data[m.end(1):]
    return hdr + '\ndata:\n' + data


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = sys.argv[1:]
    cdl = subprocess.run(['ncdump', src], check=True, capture_output=True,
                         text=True).stdout
    cdl = setcells(cdl, 'mask', '1')
    cdl = setcells(cdl, 'frac', '1')
    cdl = re.sub(r'^netcdf \S+ \{', 'netcdf glacier4_test {', cdl, count=1,
                 flags=re.M)
    tmp = out + '.cdl'
    open(tmp, 'w').write(cdl)
    subprocess.run(['ncgen', '-k', '2', '-o', out, tmp], check=True)
    subprocess.run(['rm', '-f', tmp], check=True)
    print(f'wrote {out}: {len(CELLS)} active cells {CELLS}')


if __name__ == '__main__':
    main()
