#!/usr/bin/env python3
"""Joins the per worker files written by star-fletcher into a single .rsf@ volume.

Usage: congeal.py <folder>/out-<name>.meta [output]

The index has a header (volume_width, esize, align, output) followed by one
"t i j k file offset bytes" line per block. Each block is a cube of side c stored
with x as the fastest axis, the same layout as the volume (see src/includes/macros.h).
The output defaults to the one in the index, which is the `in=` of out-<name>.rsf.
"""

import os
import sys


def read_index(path):
    header = {}
    entries = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) == 2:
                header[fields[0]] = fields[1]
            elif len(fields) == 7:
                t, i, j, k = (int(x) for x in fields[:4])
                entries.append((t, i, j, k, fields[4], int(fields[5]), int(fields[6])))
            else:
                raise ValueError(f"malformed line in {path}: {line!r}")
    return header, entries


def cube_side(n):
    c = round(n ** (1 / 3))
    if c * c * c != n:
        raise ValueError(f"block with {n} elements is not a cube")
    return c


def main():
    if len(sys.argv) not in (2, 3):
        print(__doc__, file=sys.stderr)
        return 1

    meta_path = sys.argv[1]
    folder = os.path.dirname(meta_path) or "."
    header, entries = read_index(meta_path)

    width = int(header["volume_width"])
    esize = int(header["esize"])
    out_path = sys.argv[2] if len(sys.argv) == 3 else os.path.join(folder, header["output"])

    if not entries:
        print("index has no blocks, nothing to do", file=sys.stderr)
        return 1

    moments = max(e[0] for e in entries) + 1
    moment_bytes = width ** 3 * esize

    seen = set()
    files = {}
    with open(out_path, "wb") as out:
        out.truncate(moment_bytes * moments)

        for t, i, j, k, name, offset, nbytes in entries:
            if (t, i, j, k) in seen:
                print(f"warning: block {(i, j, k)} of moment {t} written twice", file=sys.stderr)
            seen.add((t, i, j, k))

            if name not in files:
                files[name] = open(os.path.join(folder, name), "rb")
            src = files[name]
            src.seek(offset)
            block = src.read(nbytes)
            if len(block) != nbytes:
                raise ValueError(f"{name} is too short for the block at offset {offset}")

            c = cube_side(nbytes // esize)
            row = c * esize
            # each row (y, z) of the block is contiguous in the volume
            for z in range(c):
                for y in range(c):
                    vx, vy, vz = i * c, j * c + y, k * c + z
                    out.seek(t * moment_bytes + (vx + width * (vy + vz * width)) * esize)
                    src_off = (y + z * c) * row
                    out.write(block[src_off:src_off + row])

    for f in files.values():
        f.close()

    blocks_per_dim = width // cube_side(entries[0][6] // esize)
    expected = moments * blocks_per_dim ** 3
    if len(seen) != expected:
        print(f"warning: {expected - len(seen)} blocks missing from the index", file=sys.stderr)

    print(f"wrote {out_path}: {moments} moments of {width}^3, {len(seen)} blocks")
    return 0


if __name__ == "__main__":
    sys.exit(main())
