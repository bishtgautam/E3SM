#!/usr/bin/env python3
"""Make a glacier variant of a single-cell surface dataset (glacier plan, G0).

Rewrites the landunit percentages of a 1-cell surfdata file and leaves every
other variable alone. Built by CDL round trip (ncdump/ncgen; no NCO or netCDF4
needed), the same way as make_lake_surfdata.py.

PCT_URBAN and PCT_LAKE are zeroed. PCT_CROP and PCT_WETLAND must already be
zero; the script aborts otherwise rather than renormalize them away (kELM
supports neither, and the twins must see the same land surface).

Usage:
    make_glacier_surfdata.py <input surfdata> <output surfdata> <pct_glacier>
PCT_NATVEG becomes 100 - pct_glacier.
"""
import re
import subprocess
import sys


def set_scalar(cdl, name, value):
    # One gridcell: the data section holds a single value, " NAME = v ;".
    pat = re.compile(r'(\n %s =\s*)[^;]*;' % name)
    if len(pat.findall(cdl)) != 1:
        sys.exit('ERROR: expected exactly one data entry for %s' % name)
    return pat.sub(lambda m: m.group(1) + value + ' ;', cdl)


def get_values(cdl, name):
    m = re.search(r'\n %s =\s*([^;]*);' % name, cdl)
    if m is None:
        sys.exit('ERROR: %s not found in the data section' % name)
    return [float(v) for v in m.group(1).replace(',', ' ').split()]


def main():
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    src, dst, pct = sys.argv[1], sys.argv[2], float(sys.argv[3])
    if not 0.0 < pct <= 100.0:
        sys.exit('ERROR: pct_glacier must be in (0, 100]')

    cdl = subprocess.run(['ncdump', src], check=True, capture_output=True,
                         text=True).stdout
    data = cdl.index('\ndata:')
    head, body = cdl[:data], cdl[data:]

    for name in ('PCT_CROP', 'PCT_WETLAND'):
        if any(v != 0.0 for v in get_values(body, name)):
            sys.exit('ERROR: %s is nonzero on %s' % (name, src))
    nurb = len(get_values(body, 'PCT_URBAN'))

    body = set_scalar(body, 'PCT_GLACIER', '%g' % pct)
    body = set_scalar(body, 'PCT_LAKE', '0')
    body = set_scalar(body, 'PCT_NATVEG', '%g' % (100.0 - pct))
    body = set_scalar(body, 'PCT_URBAN', ', '.join(['0'] * nurb))

    note = ('\t\t:kelm_glacier_variant = "from %s: PCT_GLACIER=%g, '
            'PCT_NATVEG=%g, PCT_LAKE=0, PCT_URBAN=0 '
            '(components/kelm/tools/make_glacier_surfdata.py)" ;\n'
            % (src.split('/')[-1], pct, 100.0 - pct))
    head = head.rstrip('\n') + '\n' + note

    subprocess.run(['ncgen', '-o', dst], input=head + body, check=True,
                   text=True)


if __name__ == '__main__':
    main()
