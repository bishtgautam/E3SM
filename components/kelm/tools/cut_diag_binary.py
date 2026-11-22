#!/usr/bin/env python3
"""Cut an ELMDIAG1 diagnostic binary down to a range of timesteps.

DiagnosticReader (the C++ test reader) loads a whole binary into memory, so a
month-long recording is cut before a test is pointed at it. Records are
copied byte for byte; the header is kept. Timesteps 0 and 1 are always kept:
the tests' topology discovery reads colmeta and several fixtures read their
layout from timestep 1.

Usage:
    cut_diag_binary.py <in.bin> <out.bin> <first_ts> <last_ts>
"""
import struct
import sys

MAGIC = b"ELMDIAG1"


def main():
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    src, dst = sys.argv[1], sys.argv[2]
    first, last = int(sys.argv[3]), int(sys.argv[4])
    kept = dropped = 0
    with open(src, "rb") as f, open(dst, "wb") as g:
        if f.read(8) != MAGIC:
            sys.exit(f"{src}: not an ELMDIAG1 file")
        lab = f.read(4)
        (label_len,) = struct.unpack("<i", lab)
        g.write(MAGIC + lab)
        while True:
            head = f.read(4)
            if len(head) < 4:
                break
            (ts,) = struct.unpack("<i", head)
            rec = head + f.read(label_len)
            nd_b = f.read(4)
            (nd,) = struct.unpack("<i", nd_b)
            rec += nd_b
            if nd == 0:
                rec += f.read(8)
            elif nd == 1:
                nb = f.read(4)
                rec += nb + f.read(8 * struct.unpack("<i", nb)[0])
            elif nd == 2:
                nb = f.read(8)
                n1, n2 = struct.unpack("<ii", nb)
                rec += nb + f.read(8 * n1 * n2)
            elif nd == -1:
                nb = f.read(4)
                rec += nb + f.read(4 * struct.unpack("<i", nb)[0])
            else:
                sys.exit(f"{src}: unknown ndims {nd} at timestep {ts}")
            if ts in (0, 1) or first <= ts <= last:
                g.write(rec)
                kept += 1
            else:
                dropped += 1
            if ts > last and ts > 1:
                break
    print(f"kept {kept} records, dropped {dropped} (stopped after ts {last})")


if __name__ == "__main__":
    main()
