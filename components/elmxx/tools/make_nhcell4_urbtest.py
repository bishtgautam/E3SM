#!/usr/bin/env python3
"""Build the 4-cell urban domain for the multi-rank urban layout test (U8).

The 4x1 restart grid (make_4cell_rsttest.py) carries zero urban weight, so it
never steps an urban surface. This one does: it is the nhcell_urban domain
(full f19, URBANxx's cell) with three more f19 cells switched on, read from
the same global surfdata, so the cells differ in everything that matters to
the urban path:

  (72,114)  46.4N 285E  nhcell itself: region 6, TBD inactive, 7% lake
  (71,112)  44.5N 280E  all three density classes active, 27% lake
  (69,113)  40.7N 282.5E region 18 (other urban parameters), no lake
  (72,115)  46.4N 287.5E TBD active, snowy

Every one has zero crop, glacier and wetland (ELMxx has none, and h0 needs
natural + lake + urban weight 1). Two ranks split them two and two, so the
packed (density class, cell) urban order differs from the restart file's.

mask and frac are set to 1 on the four cells (frac is 1 there in the
standard f19 domain); everything else is the nhcell file's. CDL round trip,
no netCDF4 needed.

Usage:
    make_nhcell4_urbtest.py <nhcell domain.nc> <output.nc>
"""
import re
import subprocess
import sys

NI = 144
CELLS = [(72, 114), (71, 112), (69, 113), (72, 115)]   # (j, i), 0-based


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
    cdl = re.sub(r'^netcdf \S+ \{', 'netcdf nhcell4_urbtest {', cdl, count=1,
                 flags=re.M)
    tmp = out + '.cdl'
    open(tmp, 'w').write(cdl)
    subprocess.run(['ncgen', '-k', '2', '-o', out, tmp], check=True)
    subprocess.run(['rm', '-f', tmp], check=True)
    print(f'wrote {out}: {len(CELLS)} active cells {CELLS}')


if __name__ == '__main__':
    main()
