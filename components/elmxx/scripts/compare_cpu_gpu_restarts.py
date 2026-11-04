#!/usr/bin/env python3
"""Compare the answers of two elmxx_perf runs, field by field.

CPU and GPU builds cannot match bit for bit, so elmxx_perf only checks each
machine against its own reference. This measures how far the two machines are
apart afterwards, from the restarts that submit_perf_scaling.sh saves in ref/.

Two measures are reported for each field that differs:
  rel    worst |a - b| / max(|a|, |b|), the measure elmxx_perf --tolerance uses.
         A residual near zero, such as a balance error, can look large here.
  scale  worst |a - b| / max(|field|), the difference against the field's own
         size, which is not inflated by values near zero.

Usage:
    compare_cpu_gpu_restarts.py <cpu run> <gpu run> [--top 20] [--csv FILE]

Each run is a restart file, or a results directory from
submit_perf_scaling.sh. The two runs must have the same number of steps.
"""
import argparse
import csv
import glob
import os
import sys

import netCDF4
import numpy as np


def find_restart(path):
    """Return the restart file for a restart path or a results directory."""
    if os.path.isfile(path):
        return path
    found = sorted(glob.glob(os.path.join(path, "ref", "*.elmxx.r.*.nc")))
    if len(found) != 1:
        sys.exit(f"ERROR: expected one restart in {path}/ref, found {len(found)}")
    return found[0]


def compare_field(a, b):
    """Return the differences between two arrays, or None if they match."""
    a = np.asarray(a, dtype=np.float64).ravel()
    b = np.asarray(b, dtype=np.float64).ravel()
    nan_a, nan_b = np.isnan(a), np.isnan(b)
    nan_mismatch = int(np.count_nonzero(nan_a != nan_b))
    both = ~(nan_a | nan_b)
    a, b = a[both], b[both]

    differs = a != b
    ndiff = int(np.count_nonzero(differs))
    if ndiff == 0 and nan_mismatch == 0:
        return None

    diff = np.abs(a - b)
    largest = np.maximum(np.abs(a), np.abs(b))
    rel = np.divide(diff, largest, out=np.zeros_like(diff), where=largest > 0)
    field_size = largest.max() if largest.size else 0.0
    return {
        "ndiff": ndiff,
        "nvalues": a.size,
        "nan_mismatch": nan_mismatch,
        "max_abs": float(diff.max()) if diff.size else 0.0,
        "max_rel": float(rel.max()) if rel.size else 0.0,
        "max_scale": float(diff.max() / field_size) if field_size > 0 else 0.0,
        "field_size": float(field_size),
    }


def compare_files(cpu_file, gpu_file):
    """Compare every numeric field, and return (results, skipped, identical)."""
    results = {}
    skipped = []
    identical = 0
    with netCDF4.Dataset(cpu_file) as cpu, netCDF4.Dataset(gpu_file) as gpu:
        cpu.set_auto_mask(False)
        gpu.set_auto_mask(False)
        for name in sorted(set(cpu.variables) | set(gpu.variables)):
            if name not in cpu.variables or name not in gpu.variables:
                skipped.append(f"{name}: in only one file")
                continue
            a, b = cpu.variables[name], gpu.variables[name]
            if a.dtype.kind not in "fiu" or b.dtype.kind not in "fiu":
                continue
            if a.shape != b.shape:
                skipped.append(f"{name}: shapes {a.shape} and {b.shape}")
                continue
            result = compare_field(a[...], b[...])
            if result is None:
                identical += 1
            else:
                results[name] = result
    return results, skipped, identical


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("cpu", help="CPU restart or results directory")
    parser.add_argument("gpu", help="GPU restart or results directory")
    parser.add_argument("--top", type=int, default=20, help="fields to print")
    parser.add_argument("--csv", help="write every differing field to this file")
    args = parser.parse_args()

    cpu_file, gpu_file = find_restart(args.cpu), find_restart(args.gpu)
    results, skipped, identical = compare_files(cpu_file, gpu_file)
    ranked = sorted(results.items(), key=lambda r: r[1]["max_scale"], reverse=True)

    print(f"cpu: {cpu_file}\ngpu: {gpu_file}")
    print(f"fields: {identical} identical, {len(results)} differ, "
          f"{len(skipped)} not compared")
    for line in skipped:
        print(f"  not compared: {line}")
    if ranked:
        worst_rel = max(r["max_rel"] for r in results.values())
        worst_scale = ranked[0][1]["max_scale"]
        nan_fields = sum(1 for r in results.values() if r["nan_mismatch"])
        print(f"worst rel {worst_rel:.3e}, worst scale {worst_scale:.3e}, "
              f"fields with a NaN mismatch: {nan_fields}")
        print(f"\n{'field':<36} {'differ':>11} {'max abs':>10} "
              f"{'rel':>10} {'scale':>10} {'size':>10}")
    for name, r in ranked[:args.top]:
        print(f"{name:<36} {r['ndiff']:>5}/{r['nvalues']:<5} {r['max_abs']:>10.3e} "
              f"{r['max_rel']:>10.3e} {r['max_scale']:>10.3e} {r['field_size']:>10.3e}")

    if args.csv:
        with open(args.csv, "w", newline="") as f:
            writer = csv.writer(f)
            keys = ["ndiff", "nvalues", "nan_mismatch", "max_abs", "max_rel",
                    "max_scale", "field_size"]
            writer.writerow(["field"] + keys)
            for name, r in ranked:
                writer.writerow([name] + [r[k] for k in keys])


if __name__ == "__main__":
    main()
