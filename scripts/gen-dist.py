#!/usr/bin/env python3
# Writes a block distribution file for DIST_FILE.
# Format: "N N N", then N*N lines (k outer, j inner) of N ranks (along i), from 1 to the number of nodes.
# The idx-th number is the owner of the block idx = i + N * (j + N * k).
#
# usage: gen-dist.py <mode> <blocks per dim> <nodes> [--grid PXxPYxPZ] [--tile T|full] [--speeds s1,s2,...] [--stats]
#   slab-z, cyclic, single: fixed distributions
#   bc:   3D block-cyclic over a PX x PY x PZ grid of nodes, `tile` consecutive blocks per node on each axis
#         (tile 1 is pure cyclic, tile full gives one contiguous cuboid per node)
#   1d1d: 1D x 1D (x 1D) of Beaumont et al., as in Nesi et al. (ICPADS 2020): the shuffle of
#         Algorithm 2 on each axis, over the virtual slabs of the grid weighted by the node speeds
# In bc and 1d1d the node at grid coordinate (gx, gy, gz) is the rank gx + PX * (gy + PY * gz).
import argparse
import math
import sys

parser = argparse.ArgumentParser()
parser.add_argument("mode", choices=["slab-z", "cyclic", "single", "bc", "1d1d"])
parser.add_argument("n", type=int, help="blocks per dimension")
parser.add_argument("nodes", type=int)
parser.add_argument("--grid", help="grid of nodes PXxPYxPZ, by default the most cubic one (z split first)")
parser.add_argument("--tile", default="1", help="bc: consecutive blocks per node on each axis, or full")
parser.add_argument("--speeds", help="1d1d: relative speed of each rank, by default all equal")
parser.add_argument("--stats", action="store_true", help="print the load and communication of the distribution")
args = parser.parse_args()
n, nodes = args.n, args.nodes


# most cubic factorization of the nodes that fits in the blocks, split z first, then y
def default_grid():
    grids = [(px, py, nodes // (px * py))
             for px in range(1, nodes + 1) for py in range(1, nodes + 1)
             if nodes % (px * py) == 0 and max(px, py, nodes // (px * py)) <= n]
    if not grids:
        sys.exit(f"{nodes} nodes do not fit in {n} blocks per dimension")
    return min(grids, key=lambda g: (sum(g), -g[2], -g[1]))


# Algorithm 2: gives the indices to the slab that would finish first with one more
def shuffle(slab_speeds):
    c = [0] * len(slab_speeds)
    a = [0] * n
    for k in reversed(range(n)):
        p = min(range(len(c)), key=lambda s: (c[s] + 1) / slab_speeds[s])
        a[k] = p
        c[p] += 1
    return a


if args.mode in ("bc", "1d1d"):
    grid = tuple(int(p) for p in args.grid.split("x")) if args.grid else default_grid()
    if len(grid) != 3 or math.prod(grid) != nodes or max(grid) > n:
        sys.exit(f"invalid grid {grid} for {nodes} nodes and {n} blocks per dimension")
    px, py, pz = grid
    coords = [(r % px, (r // px) % py, r // (px * py)) for r in range(nodes)]

    if args.mode == "bc":
        tiles = [math.ceil(n / p) if args.tile == "full" else int(args.tile) for p in grid]
        axes = [[(i // t) % p for i in range(n)] for t, p in zip(tiles, grid)]
    else:
        speeds = [float(s) for s in args.speeds.split(",")] if args.speeds else [1.0] * nodes
        if len(speeds) != nodes:
            sys.exit(f"{len(speeds)} speeds for {nodes} nodes")
        # the speed of a virtual slab is the sum of the speeds of its nodes
        axes = [shuffle([sum(s for s, c in zip(speeds, coords) if c[axis] == slab) for slab in range(p)])
                for axis, p in enumerate(grid)]
    ax, ay, az = axes


def owner(i, j, k):
    if args.mode == "slab-z":
        return k * nodes // n
    if args.mode == "cyclic":
        return (i + j + k) % nodes
    if args.mode == "single":
        return 0
    return ax[i] + px * (ay[j] + py * az[k])


owners = {(i, j, k): owner(i, j, k) for k in range(n) for j in range(n) for i in range(n)}

print(n, n, n)
for k in range(n):
    for j in range(n):
        print(" ".join(str(owners[i, j, k] + 1) for i in range(n)))

if args.stats:
    blocks = [0] * nodes
    for o in owners.values():
        blocks[o] += 1
    # the 18 neighbours of t - 1 read by each wave_propagation task (src/main.c)
    plane = [(di, dj) for dj in (-1, 0, 1) for di in (-1, 0, 1) if (di, dj) != (0, 0)]
    cross = [(0, 0), (0, -1), (-1, 0), (1, 0), (0, 1)]
    offsets = [(di, dj, -1) for di, dj in cross] + [(di, dj, 0) for di, dj in plane] + [(di, dj, 1) for di, dj in cross]
    # starpu-mpi sends each block once per step to every other rank that reads it
    transfers = set()
    for (i, j, k), o in owners.items():
        for di, dj, dk in offsets:
            nb = (i + di, j + dj, k + dk)
            if nb in owners and owners[nb] != o:
                transfers.add((nb, o))
    # the share of each rank is proportional to its speed (all equal without --speeds)
    weights = [float(s) for s in args.speeds.split(",")] if args.speeds else [1.0] * nodes
    share = [n ** 3 * w / sum(weights) for w in weights]
    if args.mode in ("bc", "1d1d"):
        print(f"grid {px}x{py}x{pz}, x: {ax}, y: {ay}, z: {az}", file=sys.stderr)
    print(f"blocks per rank: {blocks}, imbalance (max blocks/share): {max(b / s for b, s in zip(blocks, share)):.3f}", file=sys.stderr)
    print(f"block transfers per step and wave: {len(transfers)} "
          f"({len(transfers) / n ** 3:.2f} per block)", file=sys.stderr)
