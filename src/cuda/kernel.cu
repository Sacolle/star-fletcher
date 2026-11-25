#include <starpu.h>
#include <math.h>
#include "kernel.h"
#include "global.h"
#include "indexing.h"
#include "halo.h"
#include "stencil.h"

// value written in the scratch where there is no slab (domain edge or 3D corner), the stencil never reads it
#ifdef RELEASE
#define MISSING_VALUE FP_LIT(0.0)
#else
#define MISSING_VALUE NAN
#endif


struct gather_params {
    struct halo_view view;
    FP* scratch;
    uint32_t scratch_ldy, scratch_ldz, scratch_ldt;
    ivec3_t scratch_sizes;
    // halo width in each axis, scratch coordinate = block coordinate + halo width
    ivec3_t halo_widths;
};

// per point in the scratch, compute the correct value comming from the descr
__global__ void gather_kernel(struct gather_params p){
    const int scratch_x = blockIdx.x * blockDim.x + threadIdx.x;
    const int scratch_y = blockIdx.y * blockDim.y + threadIdx.y;
    const int scratch_z = blockIdx.z * blockDim.z + threadIdx.z;
    if(scratch_x >= p.scratch_sizes.x || scratch_y >= p.scratch_sizes.y || scratch_z >= p.scratch_sizes.z){
        return;
    }
    // block coordinates, up to the halo width outside of the block
    const int x = scratch_x - p.halo_widths.x;
    const int y = scratch_y - p.halo_widths.y;
    const int z = scratch_z - p.halo_widths.z;
    const uint32_t scratch_idx = idx3_u32(scratch_x, scratch_y, scratch_z, p.scratch_ldy, p.scratch_ldz);
    p.scratch[scratch_idx] = halo_fetch(&p.view, x, y, z, 0, MISSING_VALUE);
    p.scratch[scratch_idx + p.scratch_ldt] = halo_fetch(&p.view, x, y, z, 1, MISSING_VALUE);
}


struct stencil_params {
    // the scratch, p at field 0 and q at field 1
    const FP* scratch;
    uint32_t scratch_ldy, scratch_ldz, scratch_ldt;
    // halo width in each axis, scratch coordinate = block coordinate + halo width
    uvec3_t halo_widths;

    // the medium tensor, every block has the same layout (ldy, ldz)
    const FP* medium;
    uint32_t medium_ldt;
    const FP *p_t1, *q_t1, *p_t2, *q_t2;
    FP *p_t, *q_t;
    uint32_t ldy, ldz;
    uvec3_t block_sizes;

    // computed points, [start, end) in each axis
    uvec3_t start, end;

    // finite difference coefficients
    FP dt, dxxinv, dyyinv, dzzinv, dxyinv, dxzinv, dyzinv;

    int has_source;
    uvec3_t source_point;
    FP source_value;
};

// computes the stencil in the scratch memory
__global__ void stencil_kernel(struct stencil_params p){
    #define ch1dxx (p.medium + MEDIUM_CH1DXX * p.medium_ldt)
    #define ch1dyy (p.medium + MEDIUM_CH1DYY * p.medium_ldt)
    #define ch1dzz (p.medium + MEDIUM_CH1DZZ * p.medium_ldt)
    #define ch1dxy (p.medium + MEDIUM_CH1DXY * p.medium_ldt)
    #define ch1dyz (p.medium + MEDIUM_CH1DYZ * p.medium_ldt)
    #define ch1dxz (p.medium + MEDIUM_CH1DXZ * p.medium_ldt)
    #define v2px   (p.medium + MEDIUM_V2PX * p.medium_ldt)
    #define v2pz   (p.medium + MEDIUM_V2PZ * p.medium_ldt)
    #define v2sz   (p.medium + MEDIUM_V2SZ * p.medium_ldt)
    #define v2pn   (p.medium + MEDIUM_V2PN * p.medium_ldt)

    const uint32_t x = blockIdx.x * blockDim.x + threadIdx.x;
    const uint32_t y = blockIdx.y * blockDim.y + threadIdx.y;
    const uint32_t z = blockIdx.z * blockDim.z + threadIdx.z;
    if (x >= p.block_sizes.x || y >= p.block_sizes.y || z >= p.block_sizes.z) {
        return;
    }

    const uint32_t idx = idx3_u32(x, y, z, p.ldy, p.ldz);

    if (
        (z < p.start.z || z >= p.end.z) ||
        (y < p.start.y || y >= p.end.y) ||
        (x < p.start.x || x >= p.end.x)
    ) {
        p.p_t[idx] = FP_LIT(0.0);
        p.q_t[idx] = FP_LIT(0.0);
        return;
    }

    const uint32_t scratch_idx = idx3_u32(x + p.halo_widths.x, y + p.halo_widths.y, z + p.halo_widths.z, p.scratch_ldy, p.scratch_ldz);
    const FP* p_center = p.scratch + scratch_idx;
    const FP* q_center = p.scratch + p.scratch_ldt + scratch_idx;
    // x is the contiguous axis
    const uint32_t scratch_ldx = 1;

    // p derivatives, H1(p) and H2(p)
    const FP pxx = stencil_second_derivative(p_center, scratch_ldx, p.dxxinv);
    const FP pyy = stencil_second_derivative(p_center, p.scratch_ldy, p.dyyinv);
    const FP pzz = stencil_second_derivative(p_center, p.scratch_ldz, p.dzzinv);
    const FP pxy = stencil_cross_derivative(p_center, scratch_ldx, p.scratch_ldy, p.dxyinv);
    const FP pyz = stencil_cross_derivative(p_center, p.scratch_ldy, p.scratch_ldz, p.dyzinv);
    const FP pxz = stencil_cross_derivative(p_center, scratch_ldx, p.scratch_ldz, p.dxzinv);

    const FP cpxx = ch1dxx[idx] * pxx;
    const FP cpyy = ch1dyy[idx] * pyy;
    const FP cpzz = ch1dzz[idx] * pzz;
    const FP cpxy = ch1dxy[idx] * pxy;
    const FP cpxz = ch1dxz[idx] * pxz;
    const FP cpyz = ch1dyz[idx] * pyz;
    const FP h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
    const FP h2p = pxx + pyy + pzz - h1p;

    // q derivatives, H1(q) and H2(q)
    const FP qxx = stencil_second_derivative(q_center, scratch_ldx, p.dxxinv);
    const FP qyy = stencil_second_derivative(q_center, p.scratch_ldy, p.dyyinv);
    const FP qzz = stencil_second_derivative(q_center, p.scratch_ldz, p.dzzinv);
    const FP qxy = stencil_cross_derivative(q_center, scratch_ldx, p.scratch_ldy, p.dxyinv);
    const FP qyz = stencil_cross_derivative(q_center, p.scratch_ldy, p.scratch_ldz, p.dyzinv);
    const FP qxz = stencil_cross_derivative(q_center, scratch_ldx, p.scratch_ldz, p.dxzinv);

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

    // new p and q, the thread that owns the source point also adds it
    FP new_p = FP_LIT(2.0) * p.p_t1[idx] - p.p_t2[idx] + rhsp * p.dt * p.dt;
    FP new_q = FP_LIT(2.0) * p.q_t1[idx] - p.q_t2[idx] + rhsq * p.dt * p.dt;
    if(p.has_source && x == p.source_point.x && y == p.source_point.y && z == p.source_point.z){
        new_p += p.source_value;
        new_q += p.source_value;
    }
    p.p_t[idx] = new_p;
    p.q_t[idx] = new_q;

    #undef ch1dxx
    #undef ch1dyy
    #undef ch1dzz
    #undef ch1dxy
    #undef ch1dyz
    #undef ch1dxz
    #undef v2px
    #undef v2pz
    #undef v2sz
    #undef v2pn
}

struct pack_params {
    const FP *p_wave, *q_wave;
    uint32_t ldy, ldz;
    // one entry per sent slab, blockIdx.y picks the slab
    FP* slabs[NEIGHBOR_COUNT];
    uint32_t slab_ldy[NEIGHBOR_COUNT], slab_ldz[NEIGHBOR_COUNT], slab_ldt[NEIGHBOR_COUNT];
    copy_region_t regions[NEIGHBOR_COUNT];
};

// copy the scratch data into the respective slabs
// considers the slab as a linear index, that's why the use of unindex
__global__ void pack_slabs_kernel(struct pack_params p){
    const int slab = blockIdx.y;
    const copy_region_t region = p.regions[slab];
    const uint32_t point_count = region.extent.x * region.extent.y * region.extent.z;
    // point is the starting point of this thread,
    // point += gridDim.x * blockDim.x jumps to this threads next respective point
    for(uint32_t point = blockIdx.x * blockDim.x + threadIdx.x; point < point_count; point += gridDim.x * blockDim.x){
        const uvec3_t coords = unindex_u32(point, region.extent.x, region.extent.y);
        const uint32_t source = idx3_u32(region.origin.x + coords.x, region.origin.y + coords.y, region.origin.z + coords.z, p.ldy, p.ldz);
        const uint32_t destination = idx3_u32(coords.x, coords.y, coords.z, p.slab_ldy[slab], p.slab_ldz[slab]);
        p.slabs[slab][destination] = p.p_wave[source];
        p.slabs[slab][destination + p.slab_ldt[slab]] = p.q_wave[source];
    }
}

#define PACK_THREADS 256
#define PACK_MAX_BLOCKS_PER_SLAB 64

// launches the pack of the slabs in the buffers first_buffer, first_buffer + 1, ...
static void launch_pack(
    void *descr[], int first_buffer, uint32_t neighbor_mask,
    const FP* p_wave, const FP* q_wave, size_t ldy, size_t ldz, cudaStream_t stream
){
    struct pack_params params;
    params.p_wave = p_wave;
    params.q_wave = q_wave;
    params.ldy = ldy;
    params.ldz = ldz;

    int slab_count = 0;
    size_t largest_slab = 0;
    // for every slab neighbor, set its data
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        if(!(neighbor_mask & NBIT(table_index))) continue;
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        void* slab = descr[first_buffer + slab_count];
        params.slabs[slab_count] = (FP*) STARPU_TENSOR_GET_PTR(slab);
        params.slab_ldy[slab_count] = STARPU_TENSOR_GET_LDY(slab);
        params.slab_ldz[slab_count] = STARPU_TENSOR_GET_LDZ(slab);
        params.slab_ldt[slab_count] = STARPU_TENSOR_GET_LDT(slab);
        params.regions[slab_count] = slab_source_region(offset);
        const size_t slab_points = uvec3_prod(params.regions[slab_count].extent);
        if(slab_points > largest_slab) largest_slab = slab_points;
        slab_count++;
    }
    if(slab_count == 0){
        return;
    }

    // code to determine the dimensions of the spawn kernel
    // amount of blocks of PACK_THREADS size, rounded up to a multiple of PACK_THREADS
    size_t blocks_per_slab = (largest_slab + PACK_THREADS - 1) / PACK_THREADS;
    // cap it at PACK_MAX_BLOCKS_PER_SLAB
    if(blocks_per_slab > PACK_MAX_BLOCKS_PER_SLAB) blocks_per_slab = PACK_MAX_BLOCKS_PER_SLAB;
    // set the grid, where the X param is the blocks per slab and the Y is the amount of slabs,
    //so the different slabs compute in parallel.
    const dim3 grid((unsigned) blocks_per_slab, (unsigned) slab_count);
    pack_slabs_kernel<<<grid, PACK_THREADS, 0, stream>>>(params);
}

static void check_launch(void){
    cudaError_t status = cudaGetLastError();
    if (status != cudaSuccess) {
        STARPU_CUDA_REPORT_ERROR(status);
    }
}

extern "C" int stencil_max_threads_per_block(void){
    struct cudaFuncAttributes attributes;
    if(cudaFuncGetAttributes(&attributes, stencil_kernel) != cudaSuccess){
        cudaGetLastError(); // clear the error, so it is not reported by the next launch check
        return -1;
    }
    return attributes.maxThreadsPerBlock;
}

// amount of blocks of size threads rounded up to a multiple of threads
static dim3 grid_for(size_t size_x, size_t size_y, size_t size_z, dim3 threads){
    return dim3(
        (size_x + threads.x - 1) / threads.x,
        (size_y + threads.y - 1) / threads.y,
        (size_z + threads.z - 1) / threads.z
    );
}

extern "C" void rtm_kernel_cuda(void *descr[], void *cl_args) {
    const rtm_args_t* args = (rtm_args_t*) cl_args;
    const uint32_t neighbor_mask = args->neighbor_mask;
    const int neighbor_count = __builtin_popcount(neighbor_mask);
    cudaStream_t stream = starpu_cuda_get_local_stream();
    const dim3 threads_per_block(g.cuda_threads.x, g.cuda_threads.y, g.cuda_threads.z);

    void* p_t1_interface = descr[RTM_HANDLE_POS_P_T1];
    const FP* p_t1 = (FP*) STARPU_BLOCK_GET_PTR(p_t1_interface);
    const FP* q_t1 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_Q_T1]);
    const size_t block_size_x = STARPU_BLOCK_GET_NX(p_t1_interface);
    const size_t block_size_y = STARPU_BLOCK_GET_NY(p_t1_interface);
    const size_t block_size_z = STARPU_BLOCK_GET_NZ(p_t1_interface);
    const size_t ldy = STARPU_BLOCK_GET_LDY(p_t1_interface);
    const size_t ldz = STARPU_BLOCK_GET_LDZ(p_t1_interface);

    void* scratch_interface = descr[RTM_HANDLE_POS_SCRATCH];
    FP* scratch = (FP*) STARPU_TENSOR_GET_PTR(scratch_interface);
    const size_t scratch_size_x = STARPU_TENSOR_GET_NX(scratch_interface);
    const size_t scratch_size_y = STARPU_TENSOR_GET_NY(scratch_interface);
    const size_t scratch_size_z = STARPU_TENSOR_GET_NZ(scratch_interface);

    // gather the slab data in the neighbors

    struct gather_params gather;
    const uvec3_t block_sizes = {(uint32_t) block_size_x, (uint32_t) block_size_y, (uint32_t) block_size_z};
    halo_view_init(&gather.view, block_sizes);
    halo_view_set_self(&gather.view, p_t1, q_t1, ldy, ldz);
    int slab_position = 0;
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        if(!(neighbor_mask & NBIT(table_index))) continue;
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        void* slab = descr[rtm_handle_pos_received_slab(slab_position++)];
        halo_view_set_slab(&gather.view, offset, (FP*) STARPU_TENSOR_GET_PTR(slab),
            STARPU_TENSOR_GET_LDY(slab), STARPU_TENSOR_GET_LDZ(slab), STARPU_TENSOR_GET_LDT(slab));
    }
    gather.scratch = scratch;
    gather.scratch_ldy = STARPU_TENSOR_GET_LDY(scratch_interface);
    gather.scratch_ldz = STARPU_TENSOR_GET_LDZ(scratch_interface);
    gather.scratch_ldt = STARPU_TENSOR_GET_LDT(scratch_interface);
    gather.scratch_sizes = {(int32_t) scratch_size_x, (int32_t) scratch_size_y, (int32_t) scratch_size_z};
    // the scratch is the block plus one halo on each side
    const uvec3_t halo_widths = {
        (uint32_t) (scratch_size_x - block_size_x) / 2,
        (uint32_t) (scratch_size_y - block_size_y) / 2,
        (uint32_t) (scratch_size_z - block_size_z) / 2,
    };
    gather.halo_widths = {(int32_t) halo_widths.x, (int32_t) halo_widths.y, (int32_t) halo_widths.z};

    gather_kernel<<<grid_for(scratch_size_x, scratch_size_y, scratch_size_z, threads_per_block), threads_per_block, 0, stream>>>(gather);
    check_launch();

    // set the parameters to call the stencil kernel
    
    void* medium_interface = descr[RTM_HANDLE_POS_MEDIUM];
    const FP dx = args->dx, dy = args->dy, dz = args->dz;
    struct stencil_params stencil;
    stencil.scratch = scratch;
    stencil.scratch_ldy = gather.scratch_ldy;
    stencil.scratch_ldz = gather.scratch_ldz;
    stencil.scratch_ldt = gather.scratch_ldt;
    stencil.halo_widths = halo_widths;
    stencil.medium = (FP*) STARPU_TENSOR_GET_PTR(medium_interface);
    stencil.medium_ldt = STARPU_TENSOR_GET_LDT(medium_interface);
    stencil.p_t1 = p_t1;
    stencil.q_t1 = q_t1;
    stencil.p_t2 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_P_T2]);
    stencil.q_t2 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_Q_T2]);
    stencil.p_t = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_P_T]);
    stencil.q_t = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_Q_T]);
    stencil.ldy = ldy;
    stencil.ldz = ldz;
    stencil.block_sizes = block_sizes;
    stencil.start = {(uint32_t) args->x_start, (uint32_t) args->y_start, (uint32_t) args->z_start};
    stencil.end = {(uint32_t) args->x_end, (uint32_t) args->y_end, (uint32_t) args->z_end};
    stencil.dt = args->dt;
    stencil.dxxinv = FP_LIT(1.0) / (dx * dx);
    stencil.dyyinv = FP_LIT(1.0) / (dy * dy);
    stencil.dzzinv = FP_LIT(1.0) / (dz * dz);
    stencil.dxyinv = FP_LIT(1.0) / (dx * dy);
    stencil.dxzinv = FP_LIT(1.0) / (dx * dz);
    stencil.dyzinv = FP_LIT(1.0) / (dy * dz);
    stencil.has_source = args->has_source;
    stencil.source_point = {(uint32_t) args->source_x, (uint32_t) args->source_y, (uint32_t) args->source_z};
    stencil.source_value = args->source_value;

    stencil_kernel<<<grid_for(block_size_x, block_size_y, block_size_z, threads_per_block), threads_per_block, 0, stream>>>(stencil);
    check_launch();

    // 3. pack
    launch_pack(descr, rtm_handle_pos_sent_slab(neighbor_count, 0), neighbor_mask,
        stencil.p_t, stencil.q_t, ldy, ldz, stream);
    check_launch();
}

extern "C" void pack_kernel_cuda(void *descr[], void *cl_args) {
    const pack_args_t* args = (pack_args_t*) cl_args;
    void* p_interface = descr[PACK_HANDLE_POS_P];
    launch_pack(descr, pack_handle_pos_sent_slab(0), args->neighbor_mask,
        (FP*) STARPU_BLOCK_GET_PTR(p_interface), (FP*) STARPU_BLOCK_GET_PTR(descr[PACK_HANDLE_POS_Q]),
        STARPU_BLOCK_GET_LDY(p_interface), STARPU_BLOCK_GET_LDZ(p_interface), starpu_cuda_get_local_stream());
    check_launch();
}
