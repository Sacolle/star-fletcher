#include <criterion/criterion.h>
#include <criterion/theories.h>
#include <criterion/new/assert.h>
#include <criterion/logging.h>

#include "macros.h"

#include "derivatives.h"

#include <stdio.h>

//derivadas definidas dentro do fletcher
#define Der1(p, i, s, dinv) (L1*(p[i+s]-p[i-s])+ L2*(p[i+2*s]-p[i-2*s]) + L3*(p[i+3*s]-p[i-3*s]) + L4*(p[i+4*s]-p[i-4*s]))*(dinv)
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
#define FP_CRIT flt

void setup_seed(void){
    srand(SEED);
}


Test(random, random_retries){
    const size_t size = 1000;
    FP first_list[size];
    FP second_list[size];

    srand(SEED);
    for(int i = 0; i < size; i++){
        first_list[i] = FP_RAND();
    }

    srand(SEED);
    for(int i = 0; i < size; i++){
        second_list[i] = FP_RAND();
    }

    for(int i = 0; i < size; i++){
        cr_assert(epsilon_eq(FP_CRIT, first_list[i], second_list[i], 0.0001));
    }
}

size_t g_volume_width;
size_t g_cube_width;
size_t g_width_in_cubes;

FP* g_volume_matrix;
FP** g_segment_matrix;

void build_matricies(){
    setup_seed();

    // use 3 because it simplifies
    g_width_in_cubes = 3;
    
    //can change to see the effect of diferent values
    g_cube_width = 16;
    g_volume_width = g_cube_width * g_width_in_cubes;

    if((g_volume_matrix = malloc(sizeof(FP) * CUBE(g_volume_width))) == NULL){
        cr_assert(false);
    }

    if((g_segment_matrix = malloc(sizeof(FP*) * CUBE(g_width_in_cubes))) == NULL){
        cr_assert(false);
    }
    for(size_t l = 0; l < CUBE(g_width_in_cubes); l++){
        if((g_segment_matrix[l] = malloc(sizeof(FP) * CUBE(g_cube_width))) == NULL){
            cr_assert(false);
        }
    }

    for(size_t k = 0; k < g_width_in_cubes; k++)
    for(size_t j = 0; j < g_width_in_cubes; j++)
    for(size_t i = 0; i < g_width_in_cubes; i++){
        for(size_t z = 0; z < g_cube_width; z++)
        for(size_t y = 0; y < g_cube_width; y++)
        for(size_t x = 0; x < g_cube_width; x++){
            const FP rand_val = FP_RAND();
            g_segment_matrix[block_idx(i, j, k)][cube_idx(x, y, z)] = rand_val;
            g_volume_matrix[block_cube_to_volume_idx(x, y, z, i, j, k)] = rand_val;
        }
    }
}

void teardown_values(){
    free(g_volume_matrix);
    for(size_t l = 0; l < CUBE(g_width_in_cubes); l++){
        free(g_segment_matrix[l]);
    }
    free(g_segment_matrix);
}

// slice of neighbor (dx, dy, dz) facing the central segment: the low neighbor gives its high
// face, the high neighbor its low face; the central segment (0, 0, 0) is taken whole
static void neighbor_slice(int dx, int dy, int dz, int32_t start[3], int32_t size[3]){
    const int32_t w = (int32_t) g_cube_width;
    const int d[3] = {dx, dy, dz};
    for(int a = 0; a < 3; a++){
        start[a] = d[a] < 0 ? w - STENCIL_RADIUS : 0;
        size[a] = d[a] != 0 ? STENCIL_RADIUS : w;
    }
}

// Neighborhood of the central segment (1, 1, 1) as the RTM task sees it on the CPU: each
// neighbor slice in place, keeping the segment's strides (layout of g_cube_face_filter).
static void build_center_neighborhood(block_view_t neighborhood[NEIGHBORHOOD_SIZE]){
    const int32_t w = (int32_t) g_cube_width;

    for(int dz = -1; dz <= 1; dz++)
    for(int dy = -1; dy <= 1; dy++)
    for(int dx = -1; dx <= 1; dx++){
        int32_t start[3], size[3];
        neighbor_slice(dx, dy, dz, start, size);
        const FP* segment = g_segment_matrix[block_idx(1 + dx, 1 + dy, 1 + dz)];

        neighborhood[NEIGHBOR_IDX(dx, dy, dz)] = (block_view_t){
            .ptr = segment + start[0] + start[1] * w + start[2] * w * w,
            .ldy = w,
            .ldz = w * w,
        };
    }
}

// Same neighborhood as the RTM task sees it on a GPU: each neighbor slice copied into its own
// compact buffer (ldy = slice width, ldz = slice width * height), like StarPU allocates the
// face/edge sub-handles on a device. The buffers go to `compact_buffers`, freed by the caller.
static void build_center_neighborhood_compact(block_view_t neighborhood[NEIGHBORHOOD_SIZE], FP* compact_buffers[NEIGHBORHOOD_SIZE]){
    const int32_t w = (int32_t) g_cube_width;

    for(int dz = -1; dz <= 1; dz++)
    for(int dy = -1; dy <= 1; dy++)
    for(int dx = -1; dx <= 1; dx++){
        int32_t start[3], size[3];
        neighbor_slice(dx, dy, dz, start, size);
        const FP* segment = g_segment_matrix[block_idx(1 + dx, 1 + dy, 1 + dz)];

        FP* buffer = malloc(sizeof(FP) * size[0] * size[1] * size[2]);
        cr_assert(buffer != NULL);
        for(int32_t z = 0; z < size[2]; z++)
        for(int32_t y = 0; y < size[1]; y++)
        for(int32_t x = 0; x < size[0]; x++){
            buffer[x + y * size[0] + z * size[0] * size[1]] =
                segment[(start[0] + x) + (start[1] + y) * w + (start[2] + z) * w * w];
        }

        compact_buffers[NEIGHBOR_IDX(dx, dy, dz)] = buffer;
        neighborhood[NEIGHBOR_IDX(dx, dy, dz)] = (block_view_t){
            .ptr = buffer,
            .ldy = size[0],
            .ldz = size[0] * size[1],
        };
    }
}

static void free_compact_buffers(FP* compact_buffers[NEIGHBORHOOD_SIZE]){
    for(int dz = -1; dz <= 1; dz++)
    for(int dy = -1; dy <= 1; dy++)
    for(int dx = -1; dx <= 1; dx++){
        free(compact_buffers[NEIGHBOR_IDX(dx, dy, dz)]);
    }
}

// Every point of the central segment, including the ones whose stencil reaches the
// neighbor faces and edges, against the fletcher-base macros on the contiguous volume.
static void check_second_derivatives(const block_view_t neighborhood[NEIGHBORHOOD_SIZE]){
    const int32_t w = (int32_t) g_cube_width;

    for(int32_t z = 0; z < w; z++)
    for(int32_t y = 0; y < w; y++)
    for(int32_t x = 0; x < w; x++){
        const size_t volume_idx = block_cube_to_volume_idx(x, y, z, 1, 1, 1);
        const FP baseline_xx = Der2(g_volume_matrix, volume_idx, 1, 1.0);
        const FP baseline_yy = Der2(g_volume_matrix, volume_idx, g_volume_width, 1.0);
        const FP baseline_zz = Der2(g_volume_matrix, volume_idx, SQUARE(g_volume_width), 1.0);

        cr_assert(epsilon_eq(FP_CRIT, baseline_xx, snd_deriv_x(neighborhood, x, y, z, w, 1.0), EPSILON), "xx at (%d, %d, %d)", x, y, z);
        cr_assert(epsilon_eq(FP_CRIT, baseline_yy, snd_deriv_y(neighborhood, x, y, z, w, 1.0), EPSILON), "yy at (%d, %d, %d)", x, y, z);
        cr_assert(epsilon_eq(FP_CRIT, baseline_zz, snd_deriv_z(neighborhood, x, y, z, w, 1.0), EPSILON), "zz at (%d, %d, %d)", x, y, z);
    }
}

static void check_cross_derivatives(const block_view_t neighborhood[NEIGHBORHOOD_SIZE]){
    const int32_t w = (int32_t) g_cube_width;

    for(int32_t z = 0; z < w; z++)
    for(int32_t y = 0; y < w; y++)
    for(int32_t x = 0; x < w; x++){
        const size_t volume_idx = block_cube_to_volume_idx(x, y, z, 1, 1, 1);
        const FP baseline_xy = DerCross(g_volume_matrix, volume_idx, 1, g_volume_width, 1.0);
        const FP baseline_yz = DerCross(g_volume_matrix, volume_idx, g_volume_width, SQUARE(g_volume_width), 1.0);
        const FP baseline_xz = DerCross(g_volume_matrix, volume_idx, 1, SQUARE(g_volume_width), 1.0);

        cr_assert(epsilon_eq(FP_CRIT, baseline_xy, cross_deriv_xy(neighborhood, x, y, z, w, 1.0), EPSILON), "xy at (%d, %d, %d)", x, y, z);
        cr_assert(epsilon_eq(FP_CRIT, baseline_yz, cross_deriv_yz(neighborhood, x, y, z, w, 1.0), EPSILON), "yz at (%d, %d, %d)", x, y, z);
        cr_assert(epsilon_eq(FP_CRIT, baseline_xz, cross_deriv_xz(neighborhood, x, y, z, w, 1.0), EPSILON), "xz at (%d, %d, %d)", x, y, z);
    }
}

TestSuite(derivative, .init = build_matricies, .fini = teardown_values);

Test(derivative, same_random_values) {
    for(size_t k = 0; k < g_width_in_cubes; k++)
    for(size_t j = 0; j < g_width_in_cubes; j++)
    for(size_t i = 0; i < g_width_in_cubes; i++){
        FP* block = g_segment_matrix[block_idx(i, j, k)];
        for(size_t z = 0; z < g_cube_width; z++)
        for(size_t y = 0; y < g_cube_width; y++)
        for(size_t x = 0; x < g_cube_width; x++){
            cr_assert(epsilon_eq(FP_CRIT, 
                block[cube_idx(x,y,z)], 
                g_volume_matrix[block_cube_to_volume_idx(x, y, z, i, j, k)], 
                EPSILON)
            );
        }
    }
}

Test(derivative, second_derivative_all_points) {
    block_view_t neighborhood[NEIGHBORHOOD_SIZE];
    build_center_neighborhood(neighborhood);
    check_second_derivatives(neighborhood);
}

Test(derivative, cross_derivative_all_points) {
    block_view_t neighborhood[NEIGHBORHOOD_SIZE];
    build_center_neighborhood(neighborhood);
    check_cross_derivatives(neighborhood);
}

// device layout: compact face/edge buffers with their own strides
Test(derivative, second_derivative_all_points_compact_faces) {
    block_view_t neighborhood[NEIGHBORHOOD_SIZE];
    FP* compact_buffers[NEIGHBORHOOD_SIZE];
    build_center_neighborhood_compact(neighborhood, compact_buffers);
    check_second_derivatives(neighborhood);
    free_compact_buffers(compact_buffers);
}

Test(derivative, cross_derivative_all_points_compact_faces) {
    block_view_t neighborhood[NEIGHBORHOOD_SIZE];
    FP* compact_buffers[NEIGHBORHOOD_SIZE];
    build_center_neighborhood_compact(neighborhood, compact_buffers);
    check_cross_derivatives(neighborhood);
    free_compact_buffers(compact_buffers);
}
