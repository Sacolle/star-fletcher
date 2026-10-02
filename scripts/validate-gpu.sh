#!/usr/bin/env bash
# Validates the CUDA kernel against the CPU one on the same machine. Run it from the repository
# root inside the CUDA + StarPU master shell, on a machine with a GPU:
#
#   nix develop .#cuda-lattest -c ./scripts/validate-gpu.sh
#
# 1. CPU reference, built with EXACT_FP=1 (no fused multiply-add). Also compared with the MD5s
#    recorded on x86-64, for information only: other architectures (e.g. aarch64 libm) may differ.
# 2. GPU bit-exact build (EXACT_FP=1, RTM tasks only on the GPU): its MD5 must equal step 1's.
# 3. GPU release build (with FMA): each output frame compared with step 1 (max abs/rel diff).
#
# ARCH (nvcc -arch=sm_XX) is detected by the Makefile from nvidia-smi, or can be given here.
# GPU_RUN wraps the GPU runs; by default nixglhost (host GPU driver on non-NixOS machines) when
# available. On NixOS use e.g. GPU_RUN="env LD_LIBRARY_PATH=/run/opengl-driver/lib".
set -u

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

if ! nvidia-smi -L > /dev/null 2>&1; then
    echo "No GPU visible (nvidia-smi -L failed): run this on a machine with an NVIDIA GPU."
    exit 1
fi
nvidia-smi -L

if [ -z "${GPU_RUN+set}" ]; then
    GPU_RUN=$(command -v nixglhost || true)
fi
echo "GPU runs wrapped with: ${GPU_RUN:-<nothing>}"

MAKE_ARCH=()
[ -n "${ARCH:-}" ] && MAKE_ARCH=("ARCH=$ARCH")

# name | MD5 recorded on x86-64 | arguments
CASES=(
    "32x3|6e61419977ee51f5a65ea0f3ae268688|TTI 32 32 32 8 12.5 12.5 12.5 0.0001 0.0003 4 0.0001"
    "64|09a09b089a06291093473ae21d0bd464|TTI 64 64 64 8 12.5 12.5 12.5 0.0005 0.0025 4 0.0025"
    "16|010ca5f86c493cb2a8c478bfd0d49aca|TTI 16 16 16 4 12.5 12.5 12.5 0.001 0.01 2 0.01"
    "32x50|fba1436f256c4d1825d742d4a2549876|TTI 32 32 32 8 12.5 12.5 12.5 0.0001 0.005 4 0.001"
)

build(){
    make clean > /dev/null
    if ! make "$@" > "$WORK/build.log" 2>&1; then
        echo "build failed: make $*"; tail -20 "$WORK/build.log"; exit 1
    fi
}

# run <label> <case name> <args...>: prints the output MD5 (or "failed");
# output goes to $WORK/out-<label>-<case name>.rsf@
run(){
    local label=$1 name=$2; shift 2
    OUTPUT_FOLDER="$WORK" OUTPUT_FILE="$label-$name" $RUNNER ./main "$@" > "$WORK/$label-$name.log" 2>&1
    local rc=$?
    if [ $rc -ne 0 ]; then
        echo "  $label $name: exit code $rc, last lines of its output:" >&2
        tail -3 "$WORK/$label-$name.log" | sed 's/^/      /' >&2; echo >&2
        if grep -q "No such device" "$WORK/$label-$name.log"; then
            echo "      -> StarPU found no worker for the RTM task: it does not see the GPU (driver not found? check GPU_RUN)" >&2
        fi
        echo "failed"
        return 1
    fi
    md5sum < "$WORK/out-$label-$name.rsf@" | cut -c1-32
}

declare -A CPU_MD5

echo "== 1. CPU reference (release, EXACT_FP=1)"
build RELEASE_MODE=1 EXACT_FP=1
RUNNER=""
for c in "${CASES[@]}"; do
    IFS='|' read -r name recorded args <<< "$c"
    got=$(STARPU_NCUDA=0 run cpu "$name" $args)
    CPU_MD5[$name]=$got
    if [ "$got" = "$recorded" ]; then note="matches the x86-64 baseline"; else note="differs from the x86-64 baseline (expected off x86-64)"; fi
    echo "  $name: $got ($note)"
done

echo "== 2. GPU bit-exact (EXACT_FP=1, RTM tasks only on the GPU) vs step 1"
build CUDA_BACKEND=1 NO_CPU_KERNEL=1 RELEASE_MODE=1 EXACT_FP=1 "${MAKE_ARCH[@]}"
grep -o -- "-arch=[a-z_0-9]*" "$WORK/build.log" | head -1 | sed 's/^/  nvcc /'
RUNNER=$GPU_RUN
failed=0
for c in "${CASES[@]}"; do
    IFS='|' read -r name recorded args <<< "$c"
    # the perturbation and output tasks only have CPU implementations: keep one CPU worker
    got=$(STARPU_NCPU=1 STARPU_NCUDA=1 run gpuexact "$name" $args)
    if [ "$got" = "${CPU_MD5[$name]}" ]; then echo "  OK   $name"; else echo "  FAIL $name: GPU $got, CPU ${CPU_MD5[$name]}"; failed=1; fi
done

echo "== 3. GPU release (with FMA) vs step 1, frame by frame"
build CUDA_BACKEND=1 NO_CPU_KERNEL=1 RELEASE_MODE=1 "${MAKE_ARCH[@]}"
for c in "${CASES[@]}"; do
    IFS='|' read -r name recorded args <<< "$c"
    STARPU_NCPU=1 STARPU_NCUDA=1 run gpu "$name" $args > /dev/null || continue
    width=$(grep n1= "$WORK/out-gpu-$name.rsf" | cut -d= -f2)
    python3 - "$WORK/out-cpu-$name.rsf@" "$WORK/out-gpu-$name.rsf@" "$width" "$name" <<'EOF'
import sys
import numpy as np
cpu_file, gpu_file, width, name = sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4]
cpu = np.fromfile(cpu_file, dtype=np.float32).reshape(-1, width, width, width)
gpu = np.fromfile(gpu_file, dtype=np.float32).reshape(-1, width, width, width)
diff = np.abs(cpu.astype(np.float64) - gpu.astype(np.float64))
scale = np.abs(cpu).max()
print(f"  {name}: frames={cpu.shape[0]} max|cpu|={scale:.3e} max|diff|={diff.max():.3e} "
      f"max|diff|/max|cpu|={diff.max() / scale if scale else 0.0:.3e} nan={np.isnan(gpu).sum()}")
EOF
done

make clean > /dev/null
exit $failed
