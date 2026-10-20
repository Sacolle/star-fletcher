# Local CPU benchmark

Reference numbers for comparing later executions on the same machine.

## Setup
- Machine: Intel Core i7-9700 @ 3.00GHz, 8 cores (no SMT), 31 GB RAM, Linux 6.18
- StarPU 1.4.12 from `nix develop .#no-cuda` (gcc 13.4), default StarPU env (`STARPU_NCPU` unset, default scheduler)
- Build: `make RELEASE_MODE=1` (CPU only)
- Run: `ENABLE_IO=0 ./main TTI 200 200 200 8 12.5 12.5 12.5 0.0001 0.01 <blocks> 0.0005`
  - volume 224 x 224 x 224 (200 + 2*8 absorb + 2*4 border), 100 time steps
  - `main` (279a40f) takes the number of segments per axis (`4` = 56³ cubes),
    `block-shapes` (7b52d01) takes the block size per axis (`bx by bz`)
- Each configuration: 1 warm-up run (calibrates the StarPU history perfmodels, discarded) + 3 timed runs.
  Time is "Computation Elapsed time", throughput is "Msamples/s"; the table shows the median of the 3 runs.
- Script: for each run, `grep "Computation Elapsed"` and `grep "Msamples/s"` of the program output.

## Results (2026-10-20)

| Version | Blocks (bx×by×bz) | # blocks | Block surface (points²) | 56³ surface / block surface | Median time (s) | Median Msamples/s | Runs (s) |
|---|---|---|---|---|---|---|---|
| main | 56×56×56 | 64 | 18816 | 1.000 | 48.96 | 20.58 | 47.22, 48.96, 50.16 |
| main (rerun at the end) | 56×56×56 | 64 | 18816 | 1.000 | 49.60 | 20.32 | 49.60, 49.37, 49.77 |
| block-shapes | 56×56×56 | 64 | 18816 | 1.000 | 50.65 | 19.90 | 49.65, 50.65, 51.22 |
| block-shapes | 112×56×28 | 64 | 21952 | 0.857 | 50.45 | 19.97 | 50.76, 50.45, 50.37 |
| block-shapes | 28×56×112 | 64 | 21952 | 0.857 | 51.78 | 19.46 | 51.20, 51.78, 51.89 |
| block-shapes | 224×56×14 | 64 | 32928 | 0.571 | 65.16 | 15.47 | 59.68, 65.16, 69.71 (*) |
| block-shapes | 224×28×28 | 64 | 26656 | 0.706 | 50.49 | 19.96 | 50.47, 50.49, 50.40 |
| block-shapes | 32×32×32 | 343 | 6144 | 3.063 | 54.98 | 18.33 | 54.92, 55.47, 54.98 |
| block-shapes | 112×112×112 | 8 | 75264 | 0.250 | 50.22 | 20.07 | 50.22, 50.22, 50.21 |

(*) Runs got slower one after the other (warm-up 54.2 s). This is the only configuration that drifted
like this, so it should be rerun before drawing conclusions about it.

## Notes
- Noise: the first main series went 42.7 s (warm-up) -> 50.2 s, the rerun stayed at 49.4-49.8 s.
  Differences under ~5% between configurations are within run-to-run variation on this machine.
- block-shapes with 56³ cubes produces bit-identical output to main (same MD5 of the `.rsf@`),
  so any time difference between them is noise or code generation, not a different computation.
- Block surface is the surface area of one block, 2·(bx·by + by·bz + bx·bz). The ratio column is
  18816 (the surface of a 56³ block) divided by that value.
