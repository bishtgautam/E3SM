#!/usr/bin/env python3
"""Build the 4-cell r05 domain for the TOP solar-radiation twins.

The domain is the full r05 grid (domain.lnd.r05_oEC60to30v3.190418.nc) with
mask = 0 everywhere except four high-relief cells. The surface dataset is
the global surfdata_0.5x0.5_simyr1850_c200609_with_TOP.nc, which carries
the five TOP inputs (STDEV_ELEV, SKY_VIEW, TERRAIN_CONFIG, SINSL_COSAS,
SINSL_SINAS). DATM's global Qian forcing covers any land cell.

  (253,123)  36.75N 118.25W  Sierra Nevada: stdev_elev 978 m, natural only
  (252,122)  36.25N 118.75W  Sierra Nevada: 923 m, 2% lake
  (234,537)  27.25N  88.75E  Himalaya: 1149 m, 2% urban (TOP skips urban)
  (242,536)  31.25N  88.25E  Tibetan plateau: 320 m, 14% lake, snowy

Every cell has zero crop, glacier and wetland (kELM has none of those
landunits), |lat| < 60 deg (TOP acts only there), stdev_elev > 0 and frac 1
in the standard r05 domain. Only mask and frac are changed. CDL round trip,
no netCDF4 needed.

Usage:
    make_top4_domain.py <r05 domain.nc> <output.nc>
"""
import re
import subprocess
import sys

NI = 720
CELLS = [(253, 123), (252, 122), (234, 537), (242, 536)]   # (j, i), 0-based


def setall(cdl, name, value, cells_value):
    """Set every value of `name` to `value`, then the CELLS to `cells_value`."""
    hdr, data = cdl.split('\ndata:\n', 1)
    m = re.search(r'^ %s =(.*?);' % name, data, re.S | re.M)
    v = [x.strip() for x in m.group(1).replace('\n', ' ').split(',')]
    v = [value] * len(v)
    for j, i in CELLS:
        v[j * NI + i] = cells_value
    body = ',\n  '.join(', '.join(v[k:k + 12]) for k in range(0, len(v), 12))
    data = data[:m.start(1)] + '\n  ' + body + ' ' + data[m.end(1):]
    return hdr + '\ndata:\n' + data


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    src, out = sys.argv[1:]
    cdl = subprocess.run(['ncdump', src], check=True, capture_output=True,
                         text=True).stdout
    cdl = setall(cdl, 'mask', '0', '1')
    cdl = setall(cdl, 'frac', '0', '1')
    cdl = re.sub(r'^netcdf \S+ \{', 'netcdf top4 {', cdl, count=1, flags=re.M)
    tmp = out + '.cdl'
    open(tmp, 'w').write(cdl)
    subprocess.run(['ncgen', '-k', '2', '-o', out, tmp], check=True)
    subprocess.run(['rm', '-f', tmp], check=True)
    print(f'wrote {out}: {len(CELLS)} active cells {CELLS}')


if __name__ == '__main__':
    main()
