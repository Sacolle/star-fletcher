
#include <starpu.h>
#include "kernel.h"

extern size_t g_cuda_thread_x, g_cuda_thread_y, g_cuda_thread_z;

#define DESCR_COUNT 52

#define CODE_IMPL
#define CUDA_CODE
#include "../derivatives/derivatives-impl.h"
#undef CUDA_CODE
#undef CODE_IMPL

#include "wave-neighborhood.h"


struct rtm_kernel_params {
    // Spatial bounds & dimensions (32 bits: cheaper index math on the GPU)
    int32_t x_start, y_start, z_start;
    int32_t x_end, y_end, z_end;
    int32_t cube_width_x, cube_width_y, cube_width_z;
    int32_t stride_y, stride_z;

    // Finite difference coefficients
    FP dt, dxxinv, dyyinv, dzzinv, dxyinv, dxzinv, dyzinv;

    // Buffer pointers, in task order (the medium, write and t - 2 buffers are read from here)
    FP* ptrs[DESCR_COUNT];

    // t - 1 central cube and neighbor face/edge views of each wave (see wave-neighborhood.h).
    // Each view keeps its own strides: on the GPU the faces are compact separate allocations.
    block_view_t p_wave_neighborhood[NEIGHBORHOOD_SIZE];
    block_view_t q_wave_neighborhood[NEIGHBORHOOD_SIZE];
};

// __grid_constant__: the neighborhoods are passed by address to the derivatives straight from
// the kernel parameters, without the compiler copying the struct to per-thread local memory.
__global__ void rtm_cuda_kernel_impl(const __grid_constant__ struct rtm_kernel_params p) {
    #define ch1dxx 	(p.ptrs[0])
    #define ch1dyy	(p.ptrs[1])
    #define ch1dzz	(p.ptrs[2])
    #define ch1dxy	(p.ptrs[3])
    #define ch1dyz	(p.ptrs[4])
    #define ch1dxz	(p.ptrs[5])
    #define v2px	(p.ptrs[6])
    #define v2pz	(p.ptrs[7])
    #define v2sz  	(p.ptrs[8])
    #define v2pn	(p.ptrs[9])
    #define pwwrite 	(p.ptrs[10])
    #define pwcentralt1 (p.p_wave_neighborhood[NEIGHBOR_IDX(0, 0, 0)].ptr)
    #define pwcentralt2 (p.ptrs[30])
    #define qwwrite  	(p.ptrs[31])
    #define qwcentralt1 (p.q_wave_neighborhood[NEIGHBOR_IDX(0, 0, 0)].ptr)
    #define qwcentralt2 (p.ptrs[51])


    // 1. Calculate thread's 3D coordinate
    const int32_t x = blockIdx.x * blockDim.x + threadIdx.x;
    const int32_t y = blockIdx.y * blockDim.y + threadIdx.y;
    const int32_t z = blockIdx.z * blockDim.z + threadIdx.z;

    // 2. Bounds check against the total cube dimension
    if (x >= p.cube_width_x || y >= p.cube_width_y || z >= p.cube_width_z) {
        return;
    }

    const int32_t idx = x + y * p.stride_y + z * p.stride_z;

    // 3. Apply the internal boundary logic from your original code
    if (
        (z < p.z_start || z >= p.z_end) ||
        (y < p.y_start || y >= p.y_end) ||
        (x < p.x_start || x >= p.x_end)
    ) {
        pwwrite[idx] = FP_LIT(0.0);
        qwwrite[idx] = FP_LIT(0.0);
        return;
    }

    // p derivatives, H1(p) and H2(p)
    const FP pxx = snd_deriv_x(p.p_wave_neighborhood, x, y, z, p.cube_width_x, p.dxxinv);
    const FP pyy = snd_deriv_y(p.p_wave_neighborhood, x, y, z, p.cube_width_y, p.dyyinv);
    const FP pzz = snd_deriv_z(p.p_wave_neighborhood, x, y, z, p.cube_width_z, p.dzzinv);
    const FP pxy = cross_deriv_xy(p.p_wave_neighborhood, x, y, z, p.cube_width_x, p.dxyinv);
    const FP pyz = cross_deriv_yz(p.p_wave_neighborhood, x, y, z, p.cube_width_y, p.dyzinv);
    const FP pxz = cross_deriv_xz(p.p_wave_neighborhood, x, y, z, p.cube_width_x, p.dxzinv);

    const FP cpxx = ch1dxx[idx] * pxx;
    const FP cpyy = ch1dyy[idx] * pyy;
    const FP cpzz = ch1dzz[idx] * pzz;
    const FP cpxy = ch1dxy[idx] * pxy;
    const FP cpxz = ch1dxz[idx] * pxz;
    const FP cpyz = ch1dyz[idx] * pyz;
    const FP h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
    const FP h2p = pxx + pyy + pzz - h1p;

    // q derivatives, H1(q) and H2(q)
    const FP qxx = snd_deriv_x(p.q_wave_neighborhood, x, y, z, p.cube_width_x, p.dxxinv);
    const FP qyy = snd_deriv_y(p.q_wave_neighborhood, x, y, z, p.cube_width_y, p.dyyinv);
    const FP qzz = snd_deriv_z(p.q_wave_neighborhood, x, y, z, p.cube_width_z, p.dzzinv);
    const FP qxy = cross_deriv_xy(p.q_wave_neighborhood, x, y, z, p.cube_width_x, p.dxyinv);
    const FP qyz = cross_deriv_yz(p.q_wave_neighborhood, x, y, z, p.cube_width_y, p.dyzinv);
    const FP qxz = cross_deriv_xz(p.q_wave_neighborhood, x, y, z, p.cube_width_x, p.dxzinv);

    const FP cqxx = ch1dxx[idx] * qxx;
    const FP cqyy = ch1dyy[idx] * qyy;
    const FP cqzz = ch1dzz[idx] * qzz;
    const FP cqxy = ch1dxy[idx] * qxy;
    const FP cqxz = ch1dxz[idx] * qxz;
    const FP cqyz = ch1dyz[idx] * qyz;
    const FP h1q = cqxx + cqyy + cqzz + cqxy + cqxz + cqyz;
    const FP h2q = qxx + qyy + qzz - h1q;

    // p-q derivatives, H1(p-q) and H2(p-q)
    const FP h1pmq = h1p - h1q;
    const FP h2pmq = h2p - h2q;

    // rhs of p and q equations
    const FP rhsp = v2px[idx] * h2p + v2pz[idx] * h1q + v2sz[idx] * h1pmq;
    const FP rhsq = v2pn[idx] * h2p + v2pz[idx] * h1q - v2sz[idx] * h2pmq;

    // new p and q
    pwwrite[idx] = FP_LIT(2.0) * pwcentralt1[idx] - pwcentralt2[idx] + rhsp * p.dt * p.dt;
    qwwrite[idx] = FP_LIT(2.0) * qwcentralt1[idx] - qwcentralt2[idx] + rhsq * p.dt * p.dt;
}



extern "C" void rtm_kernel_cuda(void *descr[], void *cl_args) {
    const rtm_args_t* args = (rtm_args_t*) cl_args;

    const FP dx = args->dx;
    const FP dy = args->dy;
    const FP dz = args->dz;
    const FP dt = args->dt;

    struct rtm_kernel_params params;

    params.x_start = (int32_t) args->x_start;
    params.y_start = (int32_t) args->y_start;
    params.z_start = (int32_t) args->z_start;
    params.x_end = (int32_t) args->x_end;
    params.y_end = (int32_t) args->y_end;
    params.z_end = (int32_t) args->z_end;

    params.cube_width_x = (int32_t) STARPU_BLOCK_GET_NX(descr[0]);
    params.cube_width_y = (int32_t) STARPU_BLOCK_GET_NY(descr[0]);
    params.cube_width_z = (int32_t) STARPU_BLOCK_GET_NZ(descr[0]);

    // the medium, write and t - 2 buffers are whole cubes with the same layout
    params.stride_y = (int32_t) STARPU_BLOCK_GET_LDY(descr[0]);
    params.stride_z = (int32_t) STARPU_BLOCK_GET_LDZ(descr[0]);

    params.dt = dt;
    params.dxxinv = FP_LIT(1.0) / (dx * dx);
    params.dyyinv = FP_LIT(1.0) / (dy * dy);
    params.dzzinv = FP_LIT(1.0) / (dz * dz);
    params.dxyinv = FP_LIT(1.0) / (dx * dy);
    params.dxzinv = FP_LIT(1.0) / (dx * dz);
    params.dyzinv = FP_LIT(1.0) / (dy * dz);

    for(size_t i = 0; i < DESCR_COUNT; i++){
        params.ptrs[i] = (FP*) STARPU_BLOCK_GET_PTR(descr[i]);
    }

    fill_wave_neighborhood(params.p_wave_neighborhood, descr, 11);
    fill_wave_neighborhood(params.q_wave_neighborhood, descr, 32);

    dim3 threads_per_block(g_cuda_thread_x, g_cuda_thread_y, g_cuda_thread_z);
    dim3 num_blocks(
        (params.cube_width_x + threads_per_block.x - 1) / threads_per_block.x,
        (params.cube_width_y + threads_per_block.y - 1) / threads_per_block.y,
        (params.cube_width_z + threads_per_block.z - 1) / threads_per_block.z
    );

    // Launch the kernel asynchronously on StarPU's managed stream
    rtm_cuda_kernel_impl<<<num_blocks, threads_per_block, 0, starpu_cuda_get_local_stream()>>>(params);

    // Standard error checking
    cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        STARPU_CUDA_REPORT_ERROR(status);
    }
}
