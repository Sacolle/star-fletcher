#include <starpu.h>

#include "kernel.h"
#include "floatingpoint.h"
#include "indexing.h"
#include "halo.h"
#include "stencil.h"
#include <errno.h>
#include <math.h>


err_t make_rtm_args(rtm_args_t** rtm_args, const rtm_args_t* values){
    *rtm_args = (rtm_args_t*) malloc(sizeof(rtm_args_t));
    if(*rtm_args == NULL){
        return errno;
    }
    **rtm_args = *values;
    return 0;
}

err_t make_pack_args(pack_args_t** pack_args, uint32_t neighbor_mask){
    *pack_args = (pack_args_t*) malloc(sizeof(pack_args_t));
    if(*pack_args == NULL){
        return errno;
    }
    (*pack_args)->neighbor_mask = neighbor_mask;
    return 0;
}

// copies the borders of the p and q blocks into the slabs sent to each existing neighbor,
// the slabs are the buffers first_buffer, first_buffer + 1, ... in the order of NEIGHBOR_OFFSETS
static void pack_sent_slabs(
    void *descr[], int first_buffer, uint32_t neighbor_mask,
    const FP* p_wave, const FP* q_wave, size_t ldy, size_t ldz
){
    int slab_position = 0;
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        // skip it doesnt have the neighbor
        if(!(neighbor_mask & NBIT(table_index))) continue;
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        void* slab = descr[first_buffer + slab_position++];
	// put the p and q wave parts into the slab
        halo_pack_slab(p_wave, q_wave, ldy, ldz, offset,
            (FP*) STARPU_TENSOR_GET_PTR(slab),
            STARPU_TENSOR_GET_LDY(slab), STARPU_TENSOR_GET_LDZ(slab), STARPU_TENSOR_GET_LDT(slab));
    }
}

void rtm_kernel(void *descr[], void *cl_args){
    const rtm_args_t* args = (rtm_args_t*) cl_args;
    const uint32_t neighbor_mask = args->neighbor_mask;
    const int neighbor_count = __builtin_popcount(neighbor_mask);

    const size_t x_start = args->x_start;
    const size_t y_start = args->y_start;
    const size_t z_start = args->z_start;
    const size_t x_end = args->x_end;
    const size_t y_end = args->y_end;
    const size_t z_end = args->z_end;

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

    // the waves at t - 1 and t - 2, and the new level t. every block of the task has the same layout
    const FP* pwcentralt1 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_P_T1]);
    const FP* qwcentralt1 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_Q_T1]);
    const FP* pwcentralt2 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_P_T2]);
    const FP* qwcentralt2 = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_Q_T2]);
    FP *const pwwrite = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_P_T]);
    FP *const qwwrite = (FP*) STARPU_BLOCK_GET_PTR(descr[RTM_HANDLE_POS_Q_T]);

    void* center_block = descr[RTM_HANDLE_POS_P_T1];
    const size_t block_size_x = STARPU_BLOCK_GET_NX(center_block);
    const size_t block_size_y = STARPU_BLOCK_GET_NY(center_block);
    const size_t block_size_z = STARPU_BLOCK_GET_NZ(center_block);
    const size_t ldy = STARPU_BLOCK_GET_LDY(center_block);
    const size_t ldz = STARPU_BLOCK_GET_LDZ(center_block);

    // unpack the precomputed vlaues which existed as 10 separate handles and
    // now are used as a tensor
    void* medium_interface = descr[RTM_HANDLE_POS_MEDIUM];
    const FP* medium = (FP*) STARPU_TENSOR_GET_PTR(medium_interface);
    const size_t medium_ldt = STARPU_TENSOR_GET_LDT(medium_interface);
    const FP* ch1dxx = medium + MEDIUM_CH1DXX * medium_ldt;
    const FP* ch1dyy = medium + MEDIUM_CH1DYY * medium_ldt;
    const FP* ch1dzz = medium + MEDIUM_CH1DZZ * medium_ldt;
    const FP* ch1dxy = medium + MEDIUM_CH1DXY * medium_ldt;
    const FP* ch1dyz = medium + MEDIUM_CH1DYZ * medium_ldt;
    const FP* ch1dxz = medium + MEDIUM_CH1DXZ * medium_ldt;
    const FP* v2px = medium + MEDIUM_V2PX * medium_ldt;
    const FP* v2pz = medium + MEDIUM_V2PZ * medium_ldt;
    const FP* v2sz = medium + MEDIUM_V2SZ * medium_ldt;
    const FP* v2pn = medium + MEDIUM_V2PN * medium_ldt;

    // the padded scratch buffer, p at field 0 and q at field 1
    void* scratch_interface = descr[RTM_HANDLE_POS_SCRATCH];
    FP *const scratch = (FP*) STARPU_TENSOR_GET_PTR(scratch_interface);
    const size_t scratch_ldx = 1;
    const size_t scratch_ldy = STARPU_TENSOR_GET_LDY(scratch_interface);
    const size_t scratch_ldz = STARPU_TENSOR_GET_LDZ(scratch_interface);
    const size_t scratch_ldt = STARPU_TENSOR_GET_LDT(scratch_interface);
    // halo width in each axis, the scratch is the block plus one halo on each side
    const uvec3_t halo_widths = {
        (STARPU_TENSOR_GET_NX(scratch_interface) - block_size_x) / 2,
        (STARPU_TENSOR_GET_NY(scratch_interface) - block_size_y) / 2,
        (STARPU_TENSOR_GET_NZ(scratch_interface) - block_size_z) / 2,
    };

    // 1. gather: the level t - 1 of the block and the slabs it received into the scratch
#ifndef RELEASE
    // the halo of the sides without a neighbor is never read, NaN makes any read of it show in the output
    for(size_t i = 0; i < 2 * scratch_ldt; i++){
        scratch[i] = NAN;
    }
#endif

    // this is the pack section of the code, where it loads the slabs into the
    // halo and the scratch
    
    const copy_region_t whole_block = {
        .origin = {0, 0, 0},
        .extent = {block_size_x, block_size_y, block_size_z}};

    // copy the center block into the scratch

    // point where the center volume starts
    const size_t block_in_scratch = idx3(
        halo_widths.x, halo_widths.y, halo_widths.z, scratch_ldy, scratch_ldz);

    extract_region(pwcentralt1, ldy, ldz, whole_block,
                   scratch + block_in_scratch, scratch_ldy, scratch_ldz);

    extract_region(qwcentralt1, ldy, ldz, whole_block,
                   scratch + scratch_ldt + block_in_scratch, scratch_ldy,
                   scratch_ldz);

    int slab_position = 0;
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        // if i don't have that neighbor, skip
        if (!(neighbor_mask & NBIT(table_index)))
          continue;
        
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        void* slab_interface = descr[rtm_handle_pos_received_slab(slab_position++)];
        const FP* slab = (FP*) STARPU_TENSOR_GET_PTR(slab_interface);
        const size_t slab_ldy = STARPU_TENSOR_GET_LDY(slab_interface);
        const size_t slab_ldz = STARPU_TENSOR_GET_LDZ(slab_interface);
        const size_t slab_ldt = STARPU_TENSOR_GET_LDT(slab_interface);
        const copy_region_t whole_slab = {
            .origin = {0, 0, 0},
            .extent = {
                STARPU_TENSOR_GET_NX(slab_interface),
                STARPU_TENSOR_GET_NY(slab_interface),
                STARPU_TENSOR_GET_NZ(slab_interface),
            },
        };
        const copy_region_t landing = slab_halo_region(offset);
        const size_t landing_in_scratch =
            idx3(landing.origin.x, landing.origin.y, landing.origin.z,
                 scratch_ldy, scratch_ldz);
	// copy slab p wave into scratch
        extract_region(slab, slab_ldy, slab_ldz, whole_slab,
                       scratch + landing_in_scratch, scratch_ldy, scratch_ldz);
	// copy slab q wave into scratch
        extract_region(slab + slab_ldt, slab_ldy, slab_ldz, whole_slab,
                       scratch + scratch_ldt + landing_in_scratch, scratch_ldy,
                       scratch_ldz);
    }

    // 2. stencil on the scratch, writes every point of the level t
    const FP* p_scratch = scratch;
    const FP* q_scratch = scratch + scratch_ldt;

    for(size_t z = 0; z < block_size_z; z++)
    for(size_t y = 0; y < block_size_y; y++)
    for(size_t x = 0; x < block_size_x; x++){
        const size_t idx = idx3(x, y, z, ldy, ldz);

        if(
            (z < z_start || z >= z_end) ||
            (y < y_start || y >= y_end) ||
            (x < x_start || x >= x_end)
        ){
            pwwrite[idx] = FP_LIT(0.0);
            qwwrite[idx] = FP_LIT(0.0);
            continue;
        }

        const size_t scratch_idx = idx3(x + halo_widths.x, y + halo_widths.y, z + halo_widths.z, scratch_ldy, scratch_ldz);
        const FP* p_center = p_scratch + scratch_idx;
        const FP* q_center = q_scratch + scratch_idx;

        // p derivatives, H1(p) and H2(p)
        const FP pxx = stencil_second_derivative(p_center, scratch_ldx, dxxinv);
        const FP pyy = stencil_second_derivative(p_center, scratch_ldy, dyyinv);
        const FP pzz = stencil_second_derivative(p_center, scratch_ldz, dzzinv);
        const FP pxy = stencil_cross_derivative(p_center, scratch_ldx, scratch_ldy, dxyinv);
        const FP pyz = stencil_cross_derivative(p_center, scratch_ldy, scratch_ldz, dyzinv);
        const FP pxz = stencil_cross_derivative(p_center, scratch_ldx, scratch_ldz, dxzinv);

        const FP cpxx = ch1dxx[idx] * pxx;
        const FP cpyy = ch1dyy[idx] * pyy;
        const FP cpzz = ch1dzz[idx] * pzz;
        const FP cpxy = ch1dxy[idx] * pxy;
        const FP cpxz = ch1dxz[idx] * pxz;
        const FP cpyz = ch1dyz[idx] * pyz;
        const FP h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
        const FP h2p = pxx + pyy + pzz - h1p;

        // q derivatives, H1(q) and H2(q)
        const FP qxx = stencil_second_derivative(q_center, scratch_ldx, dxxinv);
        const FP qyy = stencil_second_derivative(q_center, scratch_ldy, dyyinv);
        const FP qzz = stencil_second_derivative(q_center, scratch_ldz, dzzinv);
        const FP qxy = stencil_cross_derivative(q_center, scratch_ldx, scratch_ldy, dxyinv);
        const FP qyz = stencil_cross_derivative(q_center, scratch_ldy, scratch_ldz, dyzinv);
        const FP qxz = stencil_cross_derivative(q_center, scratch_ldx, scratch_ldz, dxzinv);

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
        pwwrite[idx] = FP_LIT(2.0) * pwcentralt1[idx] - pwcentralt2[idx] + rhsp * dt * dt;
        qwwrite[idx] = FP_LIT(2.0) * qwcentralt1[idx] - qwcentralt2[idx] + rhsq * dt * dt;
    }

    // the source is added to the level t before packing, so the neighbors see it
    if(args->has_source){
        const size_t source_idx = idx3(args->source_x, args->source_y, args->source_z, ldy, ldz);
        //print usado no teste sourcediff.sh
        //printf("%ld: %.9f + %.9f\n", args->t, pwwrite[source_idx], args->source_value);
        pwwrite[source_idx] += args->source_value;
        qwwrite[source_idx] += args->source_value;
    }

    // 3. pack: the borders of the level t into the slabs sent to the neighbors
    pack_sent_slabs(descr, rtm_handle_pos_sent_slab(neighbor_count, 0), neighbor_mask, pwwrite, qwwrite, ldy, ldz);
}

void pack_kernel(void *descr[], void *cl_args){
    const pack_args_t* args = (pack_args_t*) cl_args;
    void* p_interface = descr[PACK_HANDLE_POS_P];
    pack_sent_slabs(descr, pack_handle_pos_sent_slab(0), args->neighbor_mask,
        (FP*) STARPU_BLOCK_GET_PTR(p_interface), (FP*) STARPU_BLOCK_GET_PTR(descr[PACK_HANDLE_POS_Q]),
        STARPU_BLOCK_GET_LDY(p_interface), STARPU_BLOCK_GET_LDZ(p_interface));
}
