#include <criterion/criterion.h>
#include <criterion/theories.h>
#include <criterion/new/assert.h>
#include <criterion/logging.h>

#include "global.h"
#include "indexing.h"
#include "halo.h"
#include "stencil.h"
#include "medium.h"
#include "mem.h"
#include "kernel.h"

#include <starpu.h>

#include <stdio.h>

#define Der2(p, i, s, d2inv) ((K0*p[i]+ K1*(p[i+s]+p[i-s])+ K2*(p[i+2*s]+p[i-2*s]) + K3*(p[i+3*s]+p[i-3*s]) + K4*(p[i+4*s]+p[i-4*s]))*(d2inv))
#define DerCross(p, i, s11, s21, dinv) ((L11*(p[i+s21+s11]-p[i+s21-s11]-p[i-s21+s11]+p[i-s21-s11])+                                                                                    \
       L12*(p[i+s21+(2*s11)]-p[i+s21-(2*s11)]-p[i-s21+(2*s11)]+p[i-s21-(2*s11)]+p[i+(2*s21)+s11]-p[i+(2*s21)-s11]-p[i-(2*s21)+s11]+p[i-(2*s21)-s11])+                                  \
       L13*(p[i+s21+(3*s11)]-p[i+s21-(3*s11)]-p[i-s21+(3*s11)]+p[i-s21-(3*s11)]+p[i+(3*s21)+s11]-p[i+(3*s21)-s11]-p[i-(3*s21)+s11]+p[i-(3*s21)-s11])+                                  \
       L14*(p[i+s21+(4*s11)]-p[i+s21-(4*s11)]-p[i-s21+(4*s11)]+p[i-s21-(4*s11)]+p[i+(4*s21)+s11]-p[i+(4*s21)-s11]-p[i-(4*s21)+s11]+p[i-(4*s21)-s11])+                                  \
       L22*(p[i+(2*s21)+(2*s11)]-p[i+(2*s21)-(2*s11)]-p[i-(2*s21)+(2*s11)]+p[i-(2*s21)-(2*s11)])+                                                                                      \
       L23*(p[i+(2*s21)+(3*s11)]-p[i+(2*s21)-(3*s11)]-p[i-(2*s21)+(3*s11)]+p[i-(2*s21)-(3*s11)]+p[i+(3*s21)+(2*s11)]-p[i+(3*s21)-(2*s11)]-p[i-(3*s21)+(2*s11)]+p[i-(3*s21)-(2*s11)])+  \
       L24*(p[i+(2*s21)+(4*s11)]-p[i+(2*s21)-(4*s11)]-p[i-(2*s21)+(4*s11)]+p[i-(2*s21)-(4*s11)]+p[i+(4*s21)+(2*s11)]-p[i+(4*s21)-(2*s11)]-p[i-(4*s21)+(2*s11)]+p[i-(4*s21)-(2*s11)])+  \
       L33*(p[i+(3*s21)+(3*s11)]-p[i+(3*s21)-(3*s11)]-p[i-(3*s21)+(3*s11)]+p[i-(3*s21)-(3*s11)])+                                                                                      \
       L34*(p[i+(3*s21)+(4*s11)]-p[i+(3*s21)-(4*s11)]-p[i-(3*s21)+(4*s11)]+p[i-(3*s21)-(4*s11)]+p[i+(4*s21)+(3*s11)]-p[i+(4*s21)-(3*s11)]-p[i-(4*s21)+(3*s11)]+p[i-(4*s21)-(3*s11)])+  \
       L44*(p[i+(4*s21)+(4*s11)]-p[i+(4*s21)-(4*s11)]-p[i-(4*s21)+(4*s11)]+p[i-(4*s21)-(4*s11)]))*(dinv))

#define DerCrossPrint(block, i, s11, s21, dinv) \
	printf("Computed operation Macro:\n %.9f *  (\n        %.9f - %.9f - %.9f + %.9f\n    ) +\n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f + \n        %.9f - %.9f - %.9f + %.9f\n    ) +\n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f +\n        %.9f - %.9f - %.9f + %.9f\n    ) +\n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f +\n        %.9f - %.9f - %.9f + %.9f\n    ) +\n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f\n    ) +\n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f + \n        %.9f - %.9f - %.9f + %.9f\n	) +        \n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f + \n        %.9f - %.9f - %.9f + %.9f\n	) +      \n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f\n\n	) + \n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f + \n        %.9f - %.9f - %.9f + %.9f\n	) + \n    %.9f *  (\n        %.9f - %.9f - %.9f + %.9f\n    )) * %.9f\n", \
       L11, block[i+s21+s11], block[i+s21-s11], block[i-s21+s11], block[i-s21-s11],                                                                                    \
       L12, block[i+s21+(2*s11)], block[i+s21-(2*s11)], block[i-s21+(2*s11)], block[i-s21-(2*s11)], block[i+(2*s21)+s11], block[i+(2*s21)-s11], block[i-(2*s21)+s11], block[i-(2*s21)-s11],                                  \
       L13, block[i+s21+(3*s11)], block[i+s21-(3*s11)], block[i-s21+(3*s11)], block[i-s21-(3*s11)], block[i+(3*s21)+s11], block[i+(3*s21)-s11], block[i-(3*s21)+s11], block[i-(3*s21)-s11],                                  \
       L14, block[i+s21+(4*s11)], block[i+s21-(4*s11)], block[i-s21+(4*s11)], block[i-s21-(4*s11)], block[i+(4*s21)+s11], block[i+(4*s21)-s11], block[i-(4*s21)+s11], block[i-(4*s21)-s11],                                  \
       L22, block[i+(2*s21)+(2*s11)], block[i+(2*s21)-(2*s11)], block[i-(2*s21)+(2*s11)], block[i-(2*s21)-(2*s11)],                                                                                      \
       L23, block[i+(2*s21)+(3*s11)], block[i+(2*s21)-(3*s11)], block[i-(2*s21)+(3*s11)], block[i-(2*s21)-(3*s11)], block[i+(3*s21)+(2*s11)], block[i+(3*s21)-(2*s11)], block[i-(3*s21)+(2*s11)], block[i-(3*s21)-(2*s11)],  \
       L24, block[i+(2*s21)+(4*s11)], block[i+(2*s21)-(4*s11)], block[i-(2*s21)+(4*s11)], block[i-(2*s21)-(4*s11)], block[i+(4*s21)+(2*s11)], block[i+(4*s21)-(2*s11)], block[i-(4*s21)+(2*s11)], block[i-(4*s21)-(2*s11)],  \
       L33, block[i+(3*s21)+(3*s11)], block[i+(3*s21)-(3*s11)], block[i-(3*s21)+(3*s11)], block[i-(3*s21)-(3*s11)],                                                                                      \
       L34, block[i+(3*s21)+(4*s11)], block[i+(3*s21)-(4*s11)], block[i-(3*s21)+(4*s11)], block[i-(3*s21)-(4*s11)], block[i+(4*s21)+(3*s11)], block[i+(4*s21)-(3*s11)], block[i-(4*s21)+(3*s11)], block[i-(4*s21)-(3*s11)],  \
       L44, block[i+(4*s21)+(4*s11)], block[i+(4*s21)-(4*s11)], block[i-(4*s21)+(4*s11)], block[i-(4*s21)-(4*s11)], dinv)



#define SEED 42

#define EPSILON FLT_EPSILON
// the block kernel and the fletcher-base implementation run the same operations in the same order
#define KERNEL_EPSILON 1e-6
#define FP_CRIT flt

#define BORDER_WIDTH 4

#define TRY(x) TRYTO(x, {}, failed_setup)

void setup_seed(void){
    srand(SEED);
}

struct global_data g;

const size_t absorb_width = 4;

FP* g_volume_matrix_pp;
FP* g_volume_matrix_pc;
FP* g_volume_matrix_qp;
FP* g_volume_matrix_qc;
FP** g_segment_matrix_p[3];
FP** g_segment_matrix_q[3];

// one buffer per block with the MEDIUM_FIELD_COUNT fields one after the other
FP **g_medium;
#define MEDIUM_AT(field, block, point) (g_medium[block][(field) * g.block_size + (point)])

mem_vec_t allocs = NULL;

void base_computation_implementation(
    float dt, float dx, float dy, float dz, 
    size_t i, size_t med_i, size_t med_j, 
    float* pc, float* pp, float* qc, float* qp
){
    const int strideX = 1;
    const int strideY = g.volume_dims.x;
    const int strideZ = g.volume_dims.x * g.volume_dims.x;

    const float dxxinv=1.0f/(dx*dx);
    const float dyyinv=1.0f/(dy*dy);
    const float dzzinv=1.0f/(dz*dz);
    const float dxyinv=1.0f/(dx*dy);
    const float dxzinv=1.0f/(dx*dz);
    const float dyzinv=1.0f/(dy*dz);

    // p derivatives, H1(p) and H2(p)

    const float pxx= Der2(pc, i, strideX, dxxinv);
    const float pyy= Der2(pc, i, strideY, dyyinv);
    const float pzz= Der2(pc, i, strideZ, dzzinv);
    const float pxy= DerCross(pc, i, strideX, strideY, dxyinv);
    const float pyz= DerCross(pc, i, strideY, strideZ, dyzinv);
    const float pxz= DerCross(pc, i, strideX, strideZ, dxzinv);

    const float cpxx=MEDIUM_AT(MEDIUM_CH1DXX, med_i, med_j)*pxx;
    const float cpyy=MEDIUM_AT(MEDIUM_CH1DYY, med_i, med_j)*pyy;
    const float cpzz=MEDIUM_AT(MEDIUM_CH1DZZ, med_i, med_j)*pzz;
    const float cpxy=MEDIUM_AT(MEDIUM_CH1DXY, med_i, med_j)*pxy;
    const float cpxz=MEDIUM_AT(MEDIUM_CH1DXZ, med_i, med_j)*pxz;
    const float cpyz=MEDIUM_AT(MEDIUM_CH1DYZ, med_i, med_j)*pyz;
    const float h1p=cpxx+cpyy+cpzz+cpxy+cpxz+cpyz;
    const float h2p=pxx+pyy+pzz-h1p;

    // q derivatives, H1(q) and H2(q)

    const float qxx= Der2(qc, i, strideX, dxxinv);
    const float qyy= Der2(qc, i, strideY, dyyinv);
    const float qzz= Der2(qc, i, strideZ, dzzinv);
    const float qxy= DerCross(qc, i, strideX,  strideY, dxyinv);
    const float qyz= DerCross(qc, i, strideY,  strideZ, dyzinv);
    const float qxz= DerCross(qc, i, strideX,  strideZ, dxzinv);

    const float cqxx=MEDIUM_AT(MEDIUM_CH1DXX, med_i, med_j)*qxx;
    const float cqyy=MEDIUM_AT(MEDIUM_CH1DYY, med_i, med_j)*qyy;
    const float cqzz=MEDIUM_AT(MEDIUM_CH1DZZ, med_i, med_j)*qzz;
    const float cqxy=MEDIUM_AT(MEDIUM_CH1DXY, med_i, med_j)*qxy;
    const float cqxz=MEDIUM_AT(MEDIUM_CH1DXZ, med_i, med_j)*qxz;
    const float cqyz=MEDIUM_AT(MEDIUM_CH1DYZ, med_i, med_j)*qyz;
    const float h1q=cqxx+cqyy+cqzz+cqxy+cqxz+cqyz;
    const float h2q=qxx+qyy+qzz-h1q;

    // p-q derivatives, H1(p-q) and H2(p-q)

    const float h1pmq=h1p-h1q;
    const float h2pmq=h2p-h2q;

    // rhs of p and q equations

    const float rhsp=MEDIUM_AT(MEDIUM_V2PX, med_i, med_j)*h2p + MEDIUM_AT(MEDIUM_V2PZ, med_i, med_j)*h1q + MEDIUM_AT(MEDIUM_V2SZ, med_i, med_j)*h1pmq;
    const float rhsq=MEDIUM_AT(MEDIUM_V2PN, med_i, med_j)*h2p + MEDIUM_AT(MEDIUM_V2PZ, med_i, med_j)*h1q - MEDIUM_AT(MEDIUM_V2SZ, med_i, med_j)*h2pmq;

    // new p and q

    pp[i]=2.0f*pc[i] - pp[i] + rhsp*dt*dt;
    qp[i]=2.0f*qc[i] - qp[i] + rhsq*dt*dt;
}

void build_matricies(){
    setup_seed();

    // use 3 because it simplifies
    g.num_blocks = (uvec3_t){3, 3, 3};

    //can change to see the effect of diferent values
    g.block_sizes = (uvec3_t){16, 16, 16};
    g.volume_dims = (uvec3_t){16 * 3, 16 * 3, 16 * 3};
    g.block_size = uvec3_prod(g.block_sizes);
    g.total_blocks = uvec3_prod(g.num_blocks);
    TRY(halo_init(g.num_blocks, g.block_sizes, &g.halo_widths));

    const enum Form form = TTI;
    
    mem_vec_t medium_allocs = NULL;
    FP *vpz, *vsv, *epsilon, *delta, *phi, *theta;
    const size_t medium_size = sizeof(FP) * uvec3_prod(g.volume_dims);

    TRY(mem_allocate_local(&medium_allocs, (void**) &vpz, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &vsv, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &epsilon, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &delta, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &phi, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &theta, medium_size));

    // inicialize the buffers above based on the type of medium
    medium_initialize(form, uvec3_prod(g.volume_dims), vpz, vsv, epsilon, delta, phi, theta);
    
    // set the absorption zone for vpz and vsv
    medium_random_velocity_boundary(BORDER_WIDTH, absorb_width, vpz, vsv);

    // buffer_count buffers of buffer_size points each
    #define ALLOCATE_NESTED_BUFFER(v, buffer_count, buffer_size) \
        TRY(mem_allocate_local(&allocs, (void**) &v, (buffer_count) * sizeof(FP*))); \
        for(size_t i = 0; i < (buffer_count); i++) \
            TRY(mem_allocate_local(&allocs, (void**)(v + i), (buffer_size) * sizeof(FP)));

    TRY(mem_allocate_local(&allocs, (void**) &g_medium, g.total_blocks * sizeof(FP*)));
    for(size_t i = 0; i < g.total_blocks; i++)
        TRY(mem_allocate_local(&allocs, (void**)(g_medium + i), MEDIUM_FIELD_COUNT * g.block_size * sizeof(FP)));

    medium_calc_intermediary_values(vpz, vsv, epsilon, delta, phi, theta, g_medium);

    //at this point the values for the medium will not be used again
    mem_free_local(&medium_allocs);

    TRY(mem_allocate_local(&allocs, (void**) &g_volume_matrix_pp, uvec3_prod(g.volume_dims) * sizeof(FP))); 
    TRY(mem_allocate_local(&allocs, (void**) &g_volume_matrix_pc, uvec3_prod(g.volume_dims) * sizeof(FP))); 
    TRY(mem_allocate_local(&allocs, (void**) &g_volume_matrix_qp, uvec3_prod(g.volume_dims) * sizeof(FP))); 
    TRY(mem_allocate_local(&allocs, (void**) &g_volume_matrix_qc, uvec3_prod(g.volume_dims) * sizeof(FP))); 

    //TODO: size is diferent due to border
    // I do not acess the inner values of the border ones:
    // so they do not need to be initialized, but i need to be able to acess them
    const size_t padded_block_count = (g.num_blocks.x + 2) * (g.num_blocks.y + 2) * (g.num_blocks.z + 2);
    ALLOCATE_NESTED_BUFFER(g_segment_matrix_p[0], padded_block_count, g.block_size);
    ALLOCATE_NESTED_BUFFER(g_segment_matrix_p[1], padded_block_count, g.block_size);
    ALLOCATE_NESTED_BUFFER(g_segment_matrix_p[2], padded_block_count, g.block_size);

    ALLOCATE_NESTED_BUFFER(g_segment_matrix_q[0], padded_block_count, g.block_size);
    ALLOCATE_NESTED_BUFFER(g_segment_matrix_q[1], padded_block_count, g.block_size);
    ALLOCATE_NESTED_BUFFER(g_segment_matrix_q[2], padded_block_count, g.block_size);

    for(size_t k = 0; k < g.num_blocks.x; k++)
    for(size_t j = 0; j < g.num_blocks.x; j++)
    for(size_t i = 0; i < g.num_blocks.x; i++){
        for(size_t z = 0; z < g.block_sizes.x; z++)
        for(size_t y = 0; y < g.block_sizes.x; y++)
        for(size_t x = 0; x < g.block_sizes.x; x++){
            const FP rand_val = FP_RAND();
            for(size_t d = 0; d < 3; d++){
                g_segment_matrix_p[d][block_idx(i, j, k)][local_idx(x, y, z)] = 0.0f;
                g_segment_matrix_q[d][block_idx(i, j, k)][local_idx(x, y, z)] = 0.0f;
            }

            g_volume_matrix_pp[local_to_volume_idx(x, y, z, i, j, k)] = 0.0f;
            g_volume_matrix_pc[local_to_volume_idx(x, y, z, i, j, k)] = 0.0f;

            g_volume_matrix_qp[local_to_volume_idx(x, y, z, i, j, k)] = 0.0f;
            g_volume_matrix_qc[local_to_volume_idx(x, y, z, i, j, k)] = 0.0f;
        }
    }

  return;

    failed_setup:
    mem_free_local(&medium_allocs);
    mem_free_local(&allocs);
    cr_assert(false);
}

void teardown_values(){
    mem_free_local(&allocs);
}

TestSuite(fletcher_kernel, .init = build_matricies, .fini = teardown_values);

#define ASBLK(_ptr) ((struct starpu_block_interface) { \
    .id = STARPU_BLOCK_INTERFACE_ID, .ptr = (uintptr_t) (_ptr), \
    .nx = g.block_sizes.x, .ny = g.block_sizes.y, .nz = g.block_sizes.z, \
    .ldy = g.block_sizes.x, .ldz = g.block_sizes.x * g.block_sizes.y, .elemsize = sizeof(FP)})

// a contiguous nx * ny * nz * nt tensor
#define ASTENSOR(_ptr, _nx, _ny, _nz, _nt) ((struct starpu_tensor_interface) { \
    .id = STARPU_TENSOR_INTERFACE_ID, .ptr = (uintptr_t) (_ptr), \
    .nx = (_nx), .ny = (_ny), .nz = (_nz), .nt = (_nt), \
    .ldy = (_nx), .ldz = (_nx) * (_ny), .ldt = (_nx) * (_ny) * (_nz), .elemsize = sizeof(FP)})

// copies the block (i, j, k) of a volume array into a contiguous block buffer
static void extract_block(const FP* volume, size_t i, size_t j, size_t k, FP* block){
    const copy_region_t region = {
        .origin = {i * g.block_sizes.x, j * g.block_sizes.y, k * g.block_sizes.z},
        .extent = g.block_sizes,
    };
    extract_region(volume, g.volume_dims.x, g.volume_dims.x * g.volume_dims.y, region,
        block, g.block_sizes.x, g.block_sizes.x * g.block_sizes.y);
}

// runs rtm_kernel on the block (i, j, k), reading the levels t - 1 (pc, qc) and t - 2 (pp, qp) from the volume
// arrays: the slabs it receives are packed from the neighbor blocks, like the update tasks of main do.
// the level t of the block is written into p_t and q_t
static void run_kernel_on_block(size_t i, size_t j, size_t k, const rtm_args_t* args,
    const FP* pc, const FP* pp, const FP* qc, const FP* qp, FP* p_t, FP* q_t
){
    const size_t block_bytes = g.block_size * sizeof(FP);
    const uint32_t neighbor_mask = args->neighbor_mask;
    const int neighbor_count = __builtin_popcount(neighbor_mask);

    FP* p_t1 = malloc(block_bytes); FP* q_t1 = malloc(block_bytes);
    FP* p_t2 = malloc(block_bytes); FP* q_t2 = malloc(block_bytes);
    FP* neighbor_p = malloc(block_bytes); FP* neighbor_q = malloc(block_bytes);
    cr_assert(p_t1 && q_t1 && p_t2 && q_t2 && neighbor_p && neighbor_q);
    extract_block(pc, i, j, k, p_t1);
    extract_block(qc, i, j, k, q_t1);
    extract_block(pp, i, j, k, p_t2);
    extract_block(qp, i, j, k, q_t2);

    const uvec3_t scratch_dims = {
        g.block_sizes.x + 2 * g.halo_widths.x, g.block_sizes.y + 2 * g.halo_widths.y, g.block_sizes.z + 2 * g.halo_widths.z,
    };
    FP* scratch = malloc(2 * uvec3_prod(scratch_dims) * sizeof(FP));
    cr_assert(ne(ptr, scratch, NULL));

    struct starpu_block_interface blocks[6] = {
        ASBLK(p_t1), ASBLK(q_t1), ASBLK(p_t2), ASBLK(q_t2), ASBLK(p_t), ASBLK(q_t),
    };
    struct starpu_tensor_interface medium = ASTENSOR(g_medium[block_idx(i, j, k)],
        g.block_sizes.x, g.block_sizes.y, g.block_sizes.z, MEDIUM_FIELD_COUNT);
    struct starpu_tensor_interface scratch_interface = ASTENSOR(scratch, scratch_dims.x, scratch_dims.y, scratch_dims.z, 2);
    struct starpu_tensor_interface received[NEIGHBOR_COUNT], sent[NEIGHBOR_COUNT];
    FP* slab_buffers[2 * NEIGHBOR_COUNT];

    void* descr[64];
    descr[RTM_HANDLE_POS_MEDIUM] = &medium;
    descr[RTM_HANDLE_POS_P_T1] = &blocks[0];
    descr[RTM_HANDLE_POS_Q_T1] = &blocks[1];
    descr[RTM_HANDLE_POS_P_T2] = &blocks[2];
    descr[RTM_HANDLE_POS_Q_T2] = &blocks[3];
    descr[RTM_HANDLE_POS_SCRATCH] = &scratch_interface;
    descr[RTM_HANDLE_POS_P_T] = &blocks[4];
    descr[RTM_HANDLE_POS_Q_T] = &blocks[5];

    int slab_position = 0;
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        if(!((neighbor_mask >> table_index) & 1)) continue;
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        const uvec3_t dims = slab_dims(offset);
        FP* received_slab = malloc(2 * uvec3_prod(dims) * sizeof(FP));
        FP* sent_slab = malloc(2 * uvec3_prod(dims) * sizeof(FP));
        cr_assert(received_slab && sent_slab);
        slab_buffers[2 * slab_position] = received_slab;
        slab_buffers[2 * slab_position + 1] = sent_slab;

        // the neighbor at `offset` sends the slab towards the opposite offset
        extract_block(pc, i + offset.x, j + offset.y, k + offset.z, neighbor_p);
        extract_block(qc, i + offset.x, j + offset.y, k + offset.z, neighbor_q);
        halo_pack_slab(neighbor_p, neighbor_q, g.block_sizes.x, g.block_sizes.x * g.block_sizes.y,
            NEIGHBOR_OFFSETS[halo_opposite_neighbor(table_index)], received_slab, dims.x, dims.x * dims.y, uvec3_prod(dims));

        received[slab_position] = ASTENSOR(received_slab, dims.x, dims.y, dims.z, 2);
        sent[slab_position] = ASTENSOR(sent_slab, dims.x, dims.y, dims.z, 2);
        descr[rtm_handle_pos_received_slab(slab_position)] = &received[slab_position];
        descr[rtm_handle_pos_sent_slab(neighbor_count, slab_position)] = &sent[slab_position];
        slab_position++;
    }

    rtm_kernel(descr, (void*) args);

    for(int slab = 0; slab < 2 * neighbor_count; slab++) free(slab_buffers[slab]);
    free(scratch);
    free(p_t1); free(q_t1); free(p_t2); free(q_t2); free(neighbor_p); free(neighbor_q);
}

// arguments of the update of block (i, j, k), without the source
static rtm_args_t block_args(size_t i, size_t j, size_t k, FP dt, FP dx, FP dy, FP dz){
    #define AS_BORDER(v) ((v) ? BORDER_WIDTH : 0)
    const rtm_args_t args = {
        .x_start = AS_BORDER(i == 0), .x_end = g.block_sizes.x - AS_BORDER(i == g.num_blocks.x - 1),
        .y_start = AS_BORDER(j == 0), .y_end = g.block_sizes.y - AS_BORDER(j == g.num_blocks.y - 1),
        .z_start = AS_BORDER(k == 0), .z_end = g.block_sizes.z - AS_BORDER(k == g.num_blocks.z - 1),
        .dx = dx, .dy = dy, .dz = dz, .dt = dt,
        .neighbor_mask = halo_neighbor_mask(i, j, k),
    };
    #undef AS_BORDER
    return args;
}

// computes the block (i, j, k) with the fletcher-base implementation, in place in pp and qp
static void base_on_block(size_t i, size_t j, size_t k, const rtm_args_t* args, FP dt, FP dx, FP dy, FP dz){
    for(size_t z = args->z_start; z < args->z_end; z++)
    for(size_t y = args->y_start; y < args->y_end; y++)
    for(size_t x = args->x_start; x < args->x_end; x++){
        base_computation_implementation(dt, dx, dy, dz,
            local_to_volume_idx(x, y, z, i, j, k), block_idx(i, j, k), local_idx(x, y, z),
            g_volume_matrix_pc, g_volume_matrix_pp, g_volume_matrix_qc, g_volume_matrix_qp
        );
    }
}

// compares the computed points of the block (i, j, k) with the volume arrays pp and qp
static void expect_block_equals_base(size_t i, size_t j, size_t k, const rtm_args_t* args, const FP* p_t, const FP* q_t){
    for(size_t z = args->z_start; z < args->z_end; z++)
    for(size_t y = args->y_start; y < args->y_end; y++)
    for(size_t x = args->x_start; x < args->x_end; x++){
        const size_t vol_i = local_to_volume_idx(x, y, z, i, j, k);
        cr_expect(epsilon_eq(flt, g_volume_matrix_pp[vol_i], p_t[local_idx(x, y, z)], KERNEL_EPSILON),
            "p diff at block (%ld, %ld, %ld) point (%ld, %ld, %ld)", i, j, k, x, y, z);
        cr_expect(epsilon_eq(flt, g_volume_matrix_qp[vol_i], q_t[local_idx(x, y, z)], KERNEL_EPSILON),
            "q diff at block (%ld, %ld, %ld) point (%ld, %ld, %ld)", i, j, k, x, y, z);
    }
}

// one step from random levels t - 1 and t - 2, so every term of the stencil and every halo cell is used
Test(fletcher_kernel, compute_one_step) {
    const FP dt = 0.001f;
    const FP dx = 12.5f, dy = 12.5f, dz = 12.5f;

    for(size_t point = 0; point < uvec3_prod(g.volume_dims); point++){
        g_volume_matrix_pc[point] = FP_RAND() - FP_LIT(0.5);
        g_volume_matrix_pp[point] = FP_RAND() - FP_LIT(0.5);
        g_volume_matrix_qc[point] = FP_RAND() - FP_LIT(0.5);
        g_volume_matrix_qp[point] = FP_RAND() - FP_LIT(0.5);
    }

    FP* p_t = malloc(g.block_size * sizeof(FP));
    FP* q_t = malloc(g.block_size * sizeof(FP));

    for(size_t k = 0; k < g.num_blocks.z; k++)
    for(size_t j = 0; j < g.num_blocks.y; j++)
    for(size_t i = 0; i < g.num_blocks.x; i++){
        const rtm_args_t args = block_args(i, j, k, dt, dx, dy, dz);
        run_kernel_on_block(i, j, k, &args, g_volume_matrix_pc, g_volume_matrix_pp, g_volume_matrix_qc, g_volume_matrix_qp, p_t, q_t);
        // only the block (i, j, k) of pp and qp is overwritten, the next blocks read pc, qc and their own part of pp, qp
        base_on_block(i, j, k, &args, dt, dx, dy, dz);
        expect_block_equals_base(i, j, k, &args, p_t, q_t);
    }

    free(p_t);
    free(q_t);
}

// several steps from the source alone, the update of the block of the source adds it to the level t
Test(fletcher_kernel, compute_many_step) {
    const int64_t amount_of_steps = 10;
    const FP dt = 0.001f;
    const FP dx = 12.5f, dy = 12.5f, dz = 12.5f;

    const size_t volume_center_idx = volume_idx(g.volume_dims.x / 2, g.volume_dims.y / 2, g.volume_dims.z / 2);
    const size_t source_block = volume_to_block_idx(volume_center_idx);
    const uvec3_t source_point = unindex(volume_to_local_idx(volume_center_idx), g.block_sizes.x, g.block_sizes.y);

    // the level 0 holds the source value of t = 0
    g_volume_matrix_pc[volume_center_idx] += medium_source_value(dt, 0);
    g_volume_matrix_qc[volume_center_idx] += medium_source_value(dt, 0);

    FP* p_t = malloc(g.block_size * sizeof(FP));
    FP* q_t = malloc(g.block_size * sizeof(FP));

    for(int64_t t = 1; t <= amount_of_steps; t++){
        const FP source = medium_source_value(dt, t);

        for(size_t k = 0; k < g.num_blocks.z; k++)
        for(size_t j = 0; j < g.num_blocks.y; j++)
        for(size_t i = 0; i < g.num_blocks.x; i++){
            rtm_args_t args = block_args(i, j, k, dt, dx, dy, dz);
            const int has_source = block_idx(i, j, k) == source_block;
            args.has_source = has_source;
            args.source_x = source_point.x;
            args.source_y = source_point.y;
            args.source_z = source_point.z;
            args.source_value = source;
            args.t = t;

            run_kernel_on_block(i, j, k, &args, g_volume_matrix_pc, g_volume_matrix_pp, g_volume_matrix_qc, g_volume_matrix_qp, p_t, q_t);
            base_on_block(i, j, k, &args, dt, dx, dy, dz);
            if(has_source){
                g_volume_matrix_pp[volume_center_idx] += source;
                g_volume_matrix_qp[volume_center_idx] += source;
            }
            expect_block_equals_base(i, j, k, &args, p_t, q_t);
        }

        // the level t becomes the level t - 1
        FP* tmp = g_volume_matrix_pc;
        g_volume_matrix_pc = g_volume_matrix_pp;
        g_volume_matrix_pp = tmp;

        tmp = g_volume_matrix_qc;
        g_volume_matrix_qc = g_volume_matrix_qp;
        g_volume_matrix_qp = tmp;
    }

    free(p_t);
    free(q_t);
}
