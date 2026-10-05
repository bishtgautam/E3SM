#!/usr/bin/env python3
"""Make a lake variant of a single-cell surface dataset (plan, Stage 6.5 L0).

Rewrites the landunit percentages of a 1-cell surfdata file and leaves every
other variable alone. Built by CDL round trip (ncdump/ncgen; no NCO or netCDF4
needed), the same way as make_4cell_rsttest.py.

PCT_URBAN is always zeroed. The older with_lakes.nc carries 2/3/5 urban next to
its lake, and ELMxx runs no urban physics, so that area sat idle in every ELMxx
run of it. PCT_CROP, PCT_WETLAND and PCT_GLACIER must already be zero; the
script aborts otherwise rather than renormalize them away.

Usage:
    make_lake_surfdata.py <input surfdata> <output surfdata> <pct_lake>
PCT_NATVEG becomes 100 - pct_lake. LAKEDEPTH must already be on the input.
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
    src, dst, pct_lake = sys.argv[1], sys.argv[2], float(sys.argv[3])
    if not 0.0 < pct_lake <= 100.0:
        sys.exit('ERROR: pct_lake must be in (0, 100]')

    cdl = subprocess.run(['ncdump', src], check=True, capture_output=True,
                         text=True).stdout
    data = cdl.index('\ndata:')
    head, body = cdl[:data], cdl[data:]

    for name in ('PCT_CROP', 'PCT_WETLAND', 'PCT_GLACIER'):
        if any(v != 0.0 for v in get_values(body, name)):
            sys.exit('ERROR: %s is nonzero on %s' % (name, src))
    depth = get_values(body, 'LAKEDEPTH')
    if len(depth) != 1 or not depth[0] > 0.0:
        sys.exit('ERROR: need one positive LAKEDEPTH, got %s' % depth)
    nurb = len(get_values(body, 'PCT_URBAN'))

    body = set_scalar(body, 'PCT_LAKE', '%g' % pct_lake)
    body = set_scalar(body, 'PCT_NATVEG', '%g' % (100.0 - pct_lake))
    body = set_scalar(body, 'PCT_URBAN', ', '.join(['0'] * nurb))

    # Record provenance as a global attribute (ends the header's attribute list).
    note = ('\t\t:elmxx_lake_variant = "from %s: PCT_LAKE=%g, PCT_NATVEG=%g, '
            'PCT_URBAN=0 (components/elmxx/tools/make_lake_surfdata.py)" ;\n'
            % (src.split('/')[-1], pct_lake, 100.0 - pct_lake))
    head = head.rstrip('\n') + '\n' + note

    subprocess.run(['ncgen', '-o', dst], input=head + body, check=True,
                   text=True)


if __name__ == '__main__':
    main()
