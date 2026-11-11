#!/usr/bin/env python3
"""Replicate a one-cell domain and surface dataset onto an N-cell grid.

The scaling cases need a grid of any size whose cells all behave alike. Every
cell here is a copy of the source cell, at the same location, on a grid of
ni = N, nj = 1. DATM maps its forcing onto the domain as it does for any grid.

Usage:
    make_replicated_inputs.py <ncells> <output dir> [--domain FILE] [--surfdata FILE]

Writes domain.lnd.brazil_<N>.nc and surfdata_brazil_<N>.nc, and leaves files
that already exist alone. Needs netCDF4: module load python.
"""
import argparse
import os
import sys

import netCDF4
import numpy as np

DOMAIN = ("/global/cfs/cdirs/e3sm/inputdata/share/domains/domain.clm/"
          "domain.lnd.1x1pt-brazil_navy.090715.nc")
SURFDATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..",
                        "cime", "scripts", "surfdata_1x1_brazil_simyr1850_c140610.nc")


def define_var(src, dst, name):
    """Define a variable and its attributes on dst, without writing its data."""
    var = src.variables[name]
    fill = var.getncattr("_FillValue") if "_FillValue" in var.ncattrs() else None
    out = dst.createVariable(name, var.dtype, var.dimensions, fill_value=fill)
    out.setncatts({a: var.getncattr(a) for a in var.ncattrs() if a != "_FillValue"})


def write_var(src, dst, name, grid_dim):
    """Write a variable, tiling it along grid_dim when it has that dimension."""
    var = src.variables[name]
    out = dst.variables[name]

    if grid_dim not in var.dimensions:
        out[...] = var[...]
        return

    # The grid dimension has size 1 on the source, so it broadcasts wherever it
    # sits. A record at a time, so memory stays near one record's size.
    if var.dimensions[0] == "time":
        for t in range(var.shape[0]):
            out[t, ...] = np.broadcast_to(var[t, ...], out.shape[1:])
    else:
        out[...] = np.broadcast_to(var[...], out.shape)


def replicate(src_path, dst_path, grid_dim, ncells):
    """Write dst_path, a copy of src_path with grid_dim set to ncells."""
    if os.path.exists(dst_path):
        print(f"exists: {dst_path}")
        return
    tmp_path = dst_path + ".tmp"
    with netCDF4.Dataset(src_path) as src, \
            netCDF4.Dataset(tmp_path, "w", format="NETCDF3_64BIT_OFFSET") as dst:
        src.set_auto_mask(False)
        dst.set_auto_mask(False)
        # A classic file moves all its data whenever its header grows, so define
        # every variable before writing any, with nothing prefilled to move.
        dst.set_fill_off()
        if src.dimensions[grid_dim].size != 1:
            sys.exit(f"ERROR: {src_path} has more than one cell")
        for name, dim in src.dimensions.items():
            size = None if dim.isunlimited() else len(dim)
            dst.createDimension(name, ncells if name == grid_dim else size)
        dst.setncatts({a: src.getncattr(a) for a in src.ncattrs()})
        dst.setncattr("replicated_from", src_path)
        for name in src.variables:
            define_var(src, dst, name)
        for name in src.variables:
            write_var(src, dst, name, grid_dim)
    os.rename(tmp_path, dst_path)
    print(f"wrote: {dst_path} ({os.path.getsize(dst_path) / 1e9:.3f} GB)")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("ncells", type=int)
    parser.add_argument("outdir")
    parser.add_argument("--domain", default=DOMAIN)
    parser.add_argument("--surfdata", default=os.path.normpath(SURFDATA))
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    replicate(args.domain, os.path.join(args.outdir, f"domain.lnd.brazil_{args.ncells}.nc"),
              "ni", args.ncells)
    replicate(args.surfdata, os.path.join(args.outdir, f"surfdata_brazil_{args.ncells}.nc"),
              "lsmlon", args.ncells)


if __name__ == "__main__":
    main()
