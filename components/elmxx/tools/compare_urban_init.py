#!/usr/bin/env python3
"""Compare ELMxx's urban parameters with ELM's own, landunit by landunit.

Plan Stage 6.7, U1 exit: every time-constant urban parameter ELMxx seeds
must equal ELM's exactly. ELM's values are its step-1 records in the
instrumented twin's diagnostics binary (urbanalb_in, urbanrad_in,
urbanflux_in); ELMxx's are the "elmxxurb:" records elmxxUrbanMod writes at
init into the ELMXX_DIAG trace, in packed order, with the packed -> subgrid
landunit map. ELMxx's subgrid is built to match ELM's exactly (Stage 2), so
packed landunit k is ELM landunit lun_of_kurb(k). The script also checks the
landunit types agree, so a permuted map cannot pass.

  python3 components/elmxx/tools/compare_urban_init.py \\
      --elm   <ELM twin>/elm_diagnostics.bin \\
      --elmxx <ELMXX twin>/run/elmxx_diag.bin

Exit status 0 only if every field matches bit for bit.
"""

import argparse
import sys

import numpy as np

from compare_elmxx_trajectory import read_bin

# ELMxx label (after "elmxxurb:") -> ELM label holding the same landunit array.
ELM_OF = {
    "canyon_hwr":          "urbanalb_in:canyon_hwr",
    "wtroad_perv":         "urbanalb_in:wtroad_perv",
    "ht_roof":             "urbanflux_in:ht_roof",
    "wtlunit_roof":        "urbanflux_in:wtlunit_roof",
    "wind_hgt_canyon":     "urbanflux_in:wind_hgt_canyon",
    "z_d_town":            "urbanflux_in:z_d_town",
    "z_0_town":            "urbanflux_in:z_0_town",
    "eflx_traffic_factor": "urbanflux_in:eflx_traffic_factor",
    "vf_sr": "urbanalb_in:vf_sr",
    "vf_wr": "urbanalb_in:vf_wr",
    "vf_sw": "urbanalb_in:vf_sw",
    "vf_rw": "urbanalb_in:vf_rw",
    "vf_ww": "urbanalb_in:vf_ww",
    "em_roof":    "urbanrad_in:em_roof",
    "em_wall":    "urbanrad_in:em_wall",
    "em_improad": "urbanrad_in:em_improad",
    "em_perroad": "urbanrad_in:em_perroad",
    "alb_roof_dir":    "urbanalb_in:alb_roof_dir",
    "alb_roof_dif":    "urbanalb_in:alb_roof_dif",
    "alb_wall_dir":    "urbanalb_in:alb_wall_dir",
    "alb_wall_dif":    "urbanalb_in:alb_wall_dif",
    "alb_improad_dir": "urbanalb_in:alb_improad_dir",
    "alb_improad_dif": "urbanalb_in:alb_improad_dif",
    "alb_perroad_dir": "urbanalb_in:alb_perroad_dir",
    "alb_perroad_dif": "urbanalb_in:alb_perroad_dif",
}


def first_step_with(data, label):
    for ts in sorted(data):
        if label in data[ts]:
            return data[ts][label]
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--elm", required=True, help="ELM twin diagnostics binary")
    ap.add_argument("--elmxx", required=True, help="ELMxx ELMXX_DIAG trace")
    args = ap.parse_args()

    elm = read_bin(args.elm, keep=set(ELM_OF.values()) | {"urbanalb_in:lun_itype"})
    xx = read_bin(args.elmxx, keep={"elmxxurb:" + k for k in ELM_OF} |
                  {"elmxxurb:lun_of_kurb"})

    kmap = first_step_with(xx, "elmxxurb:lun_of_kurb")
    if kmap is None:
        sys.exit("ELMxx trace has no elmxxurb: records (no urban landunits, or "
                 "ELMXX_DIAG unset at init)")
    lidx = kmap.astype(int) - 1          # 0-based ELM landunit per packed k
    print(f"{len(lidx)} urban landunits; ELM landunits {list(lidx + 1)}")

    itype = first_step_with(elm, "urbanalb_in:lun_itype")
    if itype is not None:
        types = itype[lidx].astype(int)
        print(f"  ELM landunit types at those indices: {list(types)}")
        if not np.all((types >= 7) & (types <= 9)):
            sys.exit("FAIL: the packed map does not point at ELM's urban landunits")

    nbad = 0
    for k, elm_label in ELM_OF.items():
        a = first_step_with(xx, "elmxxurb:" + k)
        b = first_step_with(elm, elm_label)
        if a is None or b is None:
            print(f"  {k:20s} MISSING ({'ELMxx' if a is None else 'ELM'})")
            nbad += 1
            continue
        b = b[lidx]
        diff = np.abs(a - b)
        worst = float(np.max(diff))
        ok = worst == 0.0
        nbad += 0 if ok else 1
        print(f"  {k:20s} {'exact' if ok else 'DIFF'}  max |d| {worst:.3e}"
              + ("" if ok else f"   ELMxx {a.ravel()[:6]}  ELM {b.ravel()[:6]}"))

    print("PASS" if nbad == 0 else f"FAIL: {nbad} field(s)")
    sys.exit(0 if nbad == 0 else 1)


if __name__ == "__main__":
    main()
