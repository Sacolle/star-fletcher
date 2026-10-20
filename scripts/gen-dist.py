#!/usr/bin/env python3
# Writes a block distribution file for DIST_FILE.
# Format: "N N N", then N*N lines (k outer, j inner) of N ranks (along i), from 1 to the number of nodes.
# usage: gen-dist.py <slab-z|cyclic|single> <blocks per dim> <nodes>
import sys

kind, n, nodes = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])

def owner(i, j, k):
    if kind == "slab-z":
        return k * nodes // n
    if kind == "cyclic":
        return (i + j + k) % nodes
    if kind == "single":
        return 0
    sys.exit(f"unknown distribution {kind}")

print(n, n, n)
for k in range(n):
    for j in range(n):
        print(" ".join(str(owner(i, j, k) + 1) for i in range(n)))
