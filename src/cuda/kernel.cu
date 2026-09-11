
#include <starpu.h>
#include "kernel.h"

__device__ static inline size_t cuda_idx(
    const size_t x, const size_t y, const size_t z, 
    const size_t ldy, const size_t ldz
){
    return x + ldy * y + z * ldz;
}

#ifndef THREAD_X
#define THREAD_X 8
#endif

#ifndef THREAD_Y
#define THREAD_Y 8
#endif

#ifndef THREAD_Z
#define THREAD_Z 8
#endif

#define DESCR_COUNT 52

#define CODE_IMPL
#define CUDA_CODE
#include "../derivatives/derivatives-impl.h"
#undef CUDA_CODE
#undef CODE_IMPL


struct rtm_kernel_params {
    // Spatial bounds & dimensions
    size_t x_start, y_start, z_start;
    size_t x_end, y_end, z_end;
    size_t cube_width_x, cube_width_y, cube_width_z;
    size_t stride_x, stride_y, stride_z;

    // Finite difference coefficients
    FP dt, dxxinv, dyyinv, dzzinv, dxyinv, dxzinv, dyzinv;

    // Buffer pointers
    FP* ptrs[DESCR_COUNT];
};

__global__ void rtm_cuda_kernel_impl(struct rtm_kernel_params p) {
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
    #define pwcentralt1 (p.ptrs[11])
    #define pwip0jp0km1 (p.ptrs[12])
    #define pwip0jm1km1 (p.ptrs[13])
    #define pwim1jp0km1 (p.ptrs[14])
    #define pwip1jp0km1 (p.ptrs[15])
    #define pwip0jp1km1 (p.ptrs[16])
    #define pwim1jm1kp0 (p.ptrs[17])
    #define pwip0jm1kp0 (p.ptrs[18])
    #define pwip1jm1kp0 (p.ptrs[19])
    #define pwim1jp0kp0 (p.ptrs[20])
    #define pwip1jp0kp0 (p.ptrs[21])
    #define pwim1jp1kp0 (p.ptrs[22])
    #define pwip0jp1kp0 (p.ptrs[23])
    #define pwip1jp1kp0 (p.ptrs[24])
    #define pwip0jp0kp1 (p.ptrs[25])
    #define pwip0jm1kp1 (p.ptrs[26])
    #define pwim1jp0kp1 (p.ptrs[27])
    #define pwip1jp0kp1 (p.ptrs[28])
    #define pwip0jp1kp1 (p.ptrs[29])
    #define pwcentralt2 (p.ptrs[30])
    #define qwwrite  	(p.ptrs[31])
    #define qwcentralt1 (p.ptrs[32])
    #define qwip0jp0km1 (p.ptrs[33])
    #define qwip0jm1km1 (p.ptrs[34])
    #define qwim1jp0km1 (p.ptrs[35])
    #define qwip1jp0km1 (p.ptrs[36])
    #define qwip0jp1km1 (p.ptrs[37])
    #define qwim1jm1kp0 (p.ptrs[38])
    #define qwip0jm1kp0 (p.ptrs[39])
    #define qwip1jm1kp0 (p.ptrs[40])
    #define qwim1jp0kp0 (p.ptrs[41])
    #define qwip1jp0kp0 (p.ptrs[42])
    #define qwim1jp1kp0 (p.ptrs[43])
    #define qwip0jp1kp0 (p.ptrs[44])
    #define qwip1jp1kp0 (p.ptrs[45])
    #define qwip0jp0kp1 (p.ptrs[46])
    #define qwip0jm1kp1 (p.ptrs[47])
    #define qwim1jp0kp1 (p.ptrs[48])
    #define qwip1jp0kp1 (p.ptrs[49])
    #define qwip0jp1kp1 (p.ptrs[50])
    #define qwcentralt2 (p.ptrs[51])


    // 1. Calculate thread's 3D coordinate
    const size_t x = blockIdx.x * blockDim.x + threadIdx.x;
    const size_t y = blockIdx.y * blockDim.y + threadIdx.y;
    const size_t z = blockIdx.z * blockDim.z + threadIdx.z;

    // 2. Bounds check against the total cube dimension
    if (x >= p.cube_width_x || y >= p.cube_width_y || z >= p.cube_width_z) {
        return;
    }

    const size_t idx = cuda_idx(x, y, z, p.stride_y, p.stride_z);

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

    const FP pxx = snd_deriv_dir(pwcentralt1, pwim1jp0kp0, pwip1jp0kp0, x, idx, p.stride_x, p.dxxinv, p.cube_width_x);
    const FP pyy = snd_deriv_dir(pwcentralt1, pwip0jm1kp0, pwip0jp1kp0, y, idx, p.stride_y, p.dyyinv, p.cube_width_y);
    const FP pzz = snd_deriv_dir(pwcentralt1, pwip0jp0km1, pwip0jp0kp1, z, idx, p.stride_z, p.dzzinv, p.cube_width_z);
    const FP pxy = cross_deriv_ddir(
        pwcentralt1, idx, 
        x, pwim1jp0kp0, pwip1jp0kp0, p.stride_x, 
        y, pwip0jm1kp0, pwip0jp1kp0, p.stride_y, 
        pwip1jp1kp0, pwip1jm1kp0, pwim1jp1kp0, pwim1jm1kp0,
        p.cube_width_x, p.dxyinv
    ); 
    const FP pyz = cross_deriv_ddir(
        pwcentralt1, idx, 
        y, pwip0jm1kp0, pwip0jp1kp0, p.stride_y, 
        z, pwip0jp0km1, pwip0jp0kp1, p.stride_z, 
        pwip0jp1kp1, pwip0jp1km1, pwip0jm1kp1, pwip0jm1km1,
        p.cube_width_y, p.dyzinv
    ); 
    const FP pxz = cross_deriv_ddir(
        pwcentralt1, idx, 
        x, pwim1jp0kp0, pwip1jp0kp0, p.stride_x, 
        z, pwip0jp0km1, pwip0jp0kp1, p.stride_z, 
        pwip1jp0kp1, pwip1jp0km1, pwim1jp0kp1, pwim1jp0km1,
        p.cube_width_x, p.dxzinv
    ); 

    const FP cpxx = ch1dxx[idx] * pxx;
    const FP cpyy = ch1dyy[idx] * pyy;
    const FP cpzz = ch1dzz[idx] * pzz;
    const FP cpxy = ch1dxy[idx] * pxy;
    const FP cpxz = ch1dxz[idx] * pxz;
    const FP cpyz = ch1dyz[idx] * pyz;
    const FP h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
    const FP h2p = pxx + pyy + pzz - h1p;

    // q derivatives, H1(q) and H2(q)
    const FP qxx = snd_deriv_dir(qwcentralt1, qwim1jp0kp0, qwip1jp0kp0, x, idx, p.stride_x, p.dxxinv, p.cube_width_x);
    const FP qyy = snd_deriv_dir(qwcentralt1, qwip0jm1kp0, qwip0jp1kp0, y, idx, p.stride_y, p.dyyinv, p.cube_width_y);
    const FP qzz = snd_deriv_dir(qwcentralt1, qwip0jp0km1, qwip0jp0kp1, z, idx, p.stride_z, p.dzzinv, p.cube_width_z);
    const FP qxy = cross_deriv_ddir(
        qwcentralt1, idx, 
        x, qwim1jp0kp0, qwip1jp0kp0, p.stride_x, 
        y, qwip0jm1kp0, qwip0jp1kp0, p.stride_y, 
        qwip1jp1kp0, qwip1jm1kp0, qwim1jp1kp0, qwim1jm1kp0,
        p.cube_width_x, p.dxyinv
    ); 
    const FP qyz = cross_deriv_ddir(
        qwcentralt1, idx, 
        y, qwip0jm1kp0, qwip0jp1kp0, p.stride_y, 
        z, qwip0jp0km1, qwip0jp0kp1, p.stride_z, 
        qwip0jp1kp1, qwip0jp1km1, qwip0jm1kp1, qwip0jm1km1,
        p.cube_width_y, p.dyzinv
    ); 
    const FP qxz = cross_deriv_ddir(
        qwcentralt1, idx, 
        x, qwim1jp0kp0, qwip1jp0kp0, p.stride_x, 
        z, qwip0jp0km1, qwip0jp0kp1, p.stride_z, 
        qwip1jp0kp1, qwip1jp0km1, qwim1jp0kp1, qwim1jp0km1,
        p.cube_width_x, p.dxzinv
    ); 

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

    const FP dxxinv = FP_LIT(1.0) / (dx * dx);
    const FP dyyinv = FP_LIT(1.0) / (dy * dy);
    const FP dzzinv = FP_LIT(1.0) / (dz * dz);
    const FP dxyinv = FP_LIT(1.0) / (dx * dy);
    const FP dxzinv = FP_LIT(1.0) / (dx * dz);
    const FP dyzinv = FP_LIT(1.0) / (dy * dz);

    const size_t cube_width_x = STARPU_BLOCK_GET_NX(descr[0]);
    const size_t cube_width_y = STARPU_BLOCK_GET_NY(descr[0]);
    const size_t cube_width_z = STARPU_BLOCK_GET_NZ(descr[0]);

    const size_t stride_x = 1;
    const size_t stride_y = STARPU_BLOCK_GET_LDY(descr[0]);
    const size_t stride_z = STARPU_BLOCK_GET_LDZ(descr[0]);

    struct rtm_kernel_params params = {
	.x_start = args->x_start, .y_start = args->y_start, .z_start = args->z_start,
	.x_end = args->x_end, .y_end = args->y_end, .z_end = args->z_end,
	.cube_width_x = cube_width_x, .cube_width_y = cube_width_y, .cube_width_z = cube_width_z,
	.stride_x = stride_x, .stride_y = stride_y, .stride_z = stride_z, .dt = dt,
	.dxxinv = dxxinv, .dyyinv = dyyinv, .dzzinv = dzzinv,
	.dxyinv = dxyinv, .dxzinv = dxzinv, .dyzinv = dyzinv 
    };

    // pass only the pointer to the struct
    for(size_t i = 0; i < DESCR_COUNT; i++){
	params.ptrs[i] = (FP*) STARPU_BLOCK_GET_PTR(descr[i]);
    } 

    dim3 threads_per_block(THREAD_X, THREAD_Y, THREAD_Z); 
    dim3 num_blocks(
        (cube_width_x + threads_per_block.x - 1) / threads_per_block.x,
        (cube_width_y + threads_per_block.y - 1) / threads_per_block.y,
        (cube_width_z + threads_per_block.z - 1) / threads_per_block.z
    );

    // Launch the kernel asynchronously on StarPU's managed stream
    rtm_cuda_kernel_impl<<<num_blocks, threads_per_block, 0, starpu_cuda_get_local_stream()>>>(params);

    // Standard error checking
    cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        STARPU_CUDA_REPORT_ERROR(status);
    }
}


    
