#include <criterion/criterion.h>
#include <criterion/new/assert.h>
#include <criterion/logging.h>

#include <math.h>
#include <stdlib.h>

#include "global.h"
#include "indexing.h"
#include "halo.h"
#include "stencil.h"

struct global_data g;

// sets up g for a volume split in num_blocks blocks per axis
static err_t configure(uvec3_t volume_dims, uvec3_t num_blocks){
    g.volume_dims = volume_dims;
    g.num_blocks = num_blocks;
    g.block_sizes = (uvec3_t){
        volume_dims.x / num_blocks.x,
        volume_dims.y / num_blocks.y,
        volume_dims.z / num_blocks.z,
    };
    g.block_size = uvec3_prod(g.block_sizes);
    g.total_blocks = uvec3_prod(g.num_blocks);
    return halo_init(g.num_blocks, g.block_sizes, &g.halo_widths);
}

// ---------------------------------------------------------------------------
// neighbor table
// ---------------------------------------------------------------------------

Test(halo, neighbor_table_has_faces_and_edges_only){
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        const int nonzero_count = (offset.x != 0) + (offset.y != 0) + (offset.z != 0);
        cr_expect(eq(i32, nonzero_count >= 1 && nonzero_count <= 2, 1), "entry %d is the block itself or a 3D corner", table_index);
        // every offset is listed once
        for(int other_index = table_index + 1; other_index < NEIGHBOR_COUNT; other_index++){
            const ivec3_t other = NEIGHBOR_OFFSETS[other_index];
            cr_expect(eq(i32, offset.x == other.x && offset.y == other.y && offset.z == other.z, 0),
                "entries %d and %d are the same offset", table_index, other_index);
        }
    }
}

Test(halo, opposite_neighbor_is_negated_offset){
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        const ivec3_t opposite = NEIGHBOR_OFFSETS[halo_opposite_neighbor(table_index)];
        cr_expect(eq(i32, offset.x + opposite.x, 0), "entry %d", table_index);
        cr_expect(eq(i32, offset.y + opposite.y, 0), "entry %d", table_index);
        cr_expect(eq(i32, offset.z + opposite.z, 0), "entry %d", table_index);
    }
}

Test(halo, face_masks_match_table){
    uint32_t low[3] = {0}, high[3] = {0};
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        const int components[3] = {offset.x, offset.y, offset.z};
        for(int axis = 0; axis < 3; axis++){
            if(components[axis] < 0) low[axis] |= NBIT(table_index);
            if(components[axis] > 0) high[axis] |= NBIT(table_index);
        }
    }
    cr_expect(eq(u32, NEIGHBORS_LOW_X, low[0]));
    cr_expect(eq(u32, NEIGHBORS_HIGH_X, high[0]));
    cr_expect(eq(u32, NEIGHBORS_LOW_Y, low[1]));
    cr_expect(eq(u32, NEIGHBORS_HIGH_Y, high[1]));
    cr_expect(eq(u32, NEIGHBORS_LOW_Z, low[2]));
    cr_expect(eq(u32, NEIGHBORS_HIGH_Z, high[2]));
    cr_expect(eq(u32, NEIGHBORS_ALL, (UINT32_C(1) << NEIGHBOR_COUNT) - 1));
}

// if block_coordinate + offset is inside [0, block_count), the direct definition of a neighbor existing
static int in_range(size_t block_coordinate, int offset, size_t block_count){
    const long neighbor_coordinate = (long) block_coordinate + offset;
    return neighbor_coordinate >= 0 && neighbor_coordinate < (long) block_count;
}

Test(halo, mask_matches_direct_check_on_every_block){
    const uvec3_t decompositions[] = {{1, 1, 1}, {1, 1, 4}, {1, 4, 4}, {4, 4, 4}, {2, 3, 4}, {4, 1, 2}};
    for(size_t case_index = 0; case_index < sizeof(decompositions) / sizeof(decompositions[0]); case_index++){
        const uvec3_t num_blocks = decompositions[case_index];
        cr_assert(eq(i64, configure((uvec3_t){48, 48, 48}, num_blocks), 0));
        for(size_t k = 0; k < num_blocks.z; k++)
        for(size_t j = 0; j < num_blocks.y; j++)
        for(size_t i = 0; i < num_blocks.x; i++){
            const uint32_t mask = halo_neighbor_mask(i, j, k);
            for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
                const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
                const int exists = in_range(i, offset.x, num_blocks.x) &&
                                   in_range(j, offset.y, num_blocks.y) &&
                                   in_range(k, offset.z, num_blocks.z);
                cr_expect(eq(u32, (mask >> table_index) & 1, (uint32_t) exists),
                    "case %zu block (%zu, %zu, %zu) entry %d", case_index, i, j, k, table_index);
            }
        }
    }
}

Test(halo, center_block_neighbor_count){
    const struct { uvec3_t num_blocks; int neighbor_count; } cases[] = {
        {{1, 1, 1}, 0},
        {{1, 1, 4}, 2},   // z slabs
        {{4, 1, 1}, 2},
        {{1, 4, 4}, 8},   // pencils
        {{4, 4, 1}, 8},
        {{4, 4, 4}, 18},  // cubes
    };
    for(size_t case_index = 0; case_index < sizeof(cases) / sizeof(cases[0]); case_index++){
        const uvec3_t num_blocks = cases[case_index].num_blocks;
        cr_assert(eq(i64, configure((uvec3_t){32, 32, 32}, num_blocks), 0));
        const uint32_t mask = halo_neighbor_mask(num_blocks.x / 2, num_blocks.y / 2, num_blocks.z / 2);
        cr_expect(eq(i32, __builtin_popcount(mask), cases[case_index].neighbor_count), "case %zu", case_index);
    }
}

Test(halo, no_neighbors_along_axis_not_split){
    cr_assert(eq(i64, configure((uvec3_t){32, 32, 32}, (uvec3_t){1, 2, 4}), 0));
    for(size_t k = 0; k < g.num_blocks.z; k++)
    for(size_t j = 0; j < g.num_blocks.y; j++){
        const uint32_t mask = halo_neighbor_mask(0, j, k);
        for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
            if(NEIGHBOR_OFFSETS[table_index].x != 0){
                cr_expect(eq(u32, (mask >> table_index) & 1, 0), "block (0, %zu, %zu) entry %d", j, k, table_index);
            }
        }
    }
}

Test(halo, halo_width){
    cr_assert(eq(i64, configure((uvec3_t){32, 32, 32}, (uvec3_t){1, 2, 4}), 0));
    cr_expect(eq(sz, g.halo_widths.x, 0));
    cr_expect(eq(sz, g.halo_widths.y, STENCIL_RADIUS));
    cr_expect(eq(sz, g.halo_widths.z, STENCIL_RADIUS));
}

Test(halo, block_thinner_than_radius_fails){
    // 12 / 4 = 3 < R on a split axis
    cr_expect(ne(i64, configure((uvec3_t){12, 12, 12}, (uvec3_t){1, 1, 4}), 0));
    // 16 / 4 = R is the minimum
    cr_expect(eq(i64, configure((uvec3_t){16, 16, 16}, (uvec3_t){1, 1, 4}), 0));
}

Test(halo, neighbor_mask){
    cr_assert(eq(i64, configure((uvec3_t){24, 24, 24}, (uvec3_t){3, 3, 3}), 0));
    // corner block: 3 faces and 3 edges
    cr_expect(eq(i32, __builtin_popcount(halo_neighbor_mask(0, 0, 0)), 6));
    // center block: all of them
    cr_expect(eq(i32, __builtin_popcount(halo_neighbor_mask(1, 1, 1)), 18));
    // block at the center of a face of the volume: 5 faces, 8 edges
    cr_expect(eq(i32, __builtin_popcount(halo_neighbor_mask(1, 1, 0)), 13));

    const uint32_t mask = halo_neighbor_mask(0, 0, 0);
    for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
        const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
        const int expected = offset.x >= 0 && offset.y >= 0 && offset.z >= 0;
        cr_expect(eq(i32, (mask >> table_index) & 1, expected),
            "entry %d (%d, %d, %d)", table_index, offset.x, offset.y, offset.z);
    }
}

Test(halo, slab_regions){
    cr_assert(eq(i64, configure((uvec3_t){24, 32, 40}, (uvec3_t){3, 4, 5}), 0));
    const copy_region_t face = slab_source_region((ivec3_t){1, 0, 0});
    cr_expect(eq(sz, face.origin.x, 8 - STENCIL_RADIUS));
    cr_expect(eq(sz, face.origin.y, 0));
    cr_expect(eq(sz, face.extent.x, STENCIL_RADIUS));
    cr_expect(eq(sz, face.extent.y, 8));
    cr_expect(eq(sz, face.extent.z, 8));

    const copy_region_t edge = slab_source_region((ivec3_t){0, -1, 1});
    cr_expect(eq(sz, edge.origin.y, 0));
    cr_expect(eq(sz, edge.origin.z, 8 - STENCIL_RADIUS));
    cr_expect(eq(sz, edge.extent.x, 8));
    cr_expect(eq(sz, edge.extent.y, STENCIL_RADIUS));
    cr_expect(eq(sz, edge.extent.z, STENCIL_RADIUS));

    const copy_region_t halo = slab_halo_region((ivec3_t){1, -1, 0});
    cr_expect(eq(sz, halo.origin.x, 8 + STENCIL_RADIUS));
    cr_expect(eq(sz, halo.origin.y, 0));
    cr_expect(eq(sz, halo.origin.z, STENCIL_RADIUS));
}

// ---------------------------------------------------------------------------
// pack then fetch round trip
// ---------------------------------------------------------------------------

// a value unique to each point of the volume, exact in float for volumes below 64 per axis
static FP global_value(long x, long y, long z, int field){
    return (FP)(x + 64 * y + 4096 * z) + (field ? FP_LIT(0.5) : FP_LIT(0.0));
}

#define MISSING_VALUE FP_LIT(-1.0)

// packs every block of the volume into the slabs it sends, then checks that every point of the
// padded region of every block, read through halo_fetch, has the value of the volume at that point
static void round_trip(uvec3_t volume_dims, uvec3_t num_blocks){
    cr_assert(eq(i64, configure(volume_dims, num_blocks), 0));
    const uvec3_t block_sizes = g.block_sizes;
    const uvec3_t halo_widths = g.halo_widths;

    FP** p_waves = malloc(g.total_blocks * sizeof(FP*));
    FP** q_waves = malloc(g.total_blocks * sizeof(FP*));
    // sent_slabs[block * NEIGHBOR_COUNT + table_index]: what the block sends towards that neighbor
    FP** sent_slabs = calloc(g.total_blocks * NEIGHBOR_COUNT, sizeof(FP*));
    cr_assert(p_waves && q_waves && sent_slabs);

    for(size_t k = 0; k < num_blocks.z; k++)
    for(size_t j = 0; j < num_blocks.y; j++)
    for(size_t i = 0; i < num_blocks.x; i++){
        const size_t block = block_idx(i, j, k);
        p_waves[block] = malloc(g.block_size * sizeof(FP));
        q_waves[block] = malloc(g.block_size * sizeof(FP));
        for(size_t z = 0; z < block_sizes.z; z++)
        for(size_t y = 0; y < block_sizes.y; y++)
        for(size_t x = 0; x < block_sizes.x; x++){
            const long volume_x = x + i * block_sizes.x;
            const long volume_y = y + j * block_sizes.y;
            const long volume_z = z + k * block_sizes.z;
            p_waves[block][local_idx(x, y, z)] = global_value(volume_x, volume_y, volume_z, 0);
            q_waves[block][local_idx(x, y, z)] = global_value(volume_x, volume_y, volume_z, 1);
        }
    }

    for(size_t k = 0; k < num_blocks.z; k++)
    for(size_t j = 0; j < num_blocks.y; j++)
    for(size_t i = 0; i < num_blocks.x; i++){
        const size_t block = block_idx(i, j, k);
        const uint32_t mask = halo_neighbor_mask(i, j, k);
        for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
            if(!((mask >> table_index) & 1)) continue;
            const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
            const uvec3_t slab_extent = slab_dims(offset);
            FP* slab = malloc(2 * uvec3_prod(slab_extent) * sizeof(FP));
            halo_pack_slab(p_waves[block], q_waves[block], block_sizes.x, block_sizes.x * block_sizes.y,
                offset, slab, slab_extent.x, slab_extent.x * slab_extent.y, uvec3_prod(slab_extent));
            sent_slabs[block * NEIGHBOR_COUNT + table_index] = slab;
        }
    }

    size_t checked_count = 0;
    for(size_t k = 0; k < num_blocks.z; k++)
    for(size_t j = 0; j < num_blocks.y; j++)
    for(size_t i = 0; i < num_blocks.x; i++){
        const size_t block = block_idx(i, j, k);
        struct halo_view view;
        halo_view_init(&view, block_sizes);
        halo_view_set_self(&view, p_waves[block], q_waves[block], block_sizes.x, block_sizes.x * block_sizes.y);
        const uint32_t mask = halo_neighbor_mask(i, j, k);
        // the slab coming from the neighbor at `offset` is what that neighbor sends towards -offset
        for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
            if(!((mask >> table_index) & 1)) continue;
            const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
            const size_t neighbor = block_idx(i + offset.x, j + offset.y, k + offset.z);
            const int opposite_index = halo_opposite_neighbor(table_index);
            const uvec3_t slab_extent = slab_dims(offset);
            halo_view_set_slab(&view, offset,
                sent_slabs[neighbor * NEIGHBOR_COUNT + opposite_index], slab_extent.x, slab_extent.x * slab_extent.y, uvec3_prod(slab_extent));
        }

        for(long z = -(long) halo_widths.z; z < (long)(block_sizes.z + halo_widths.z); z++)
        for(long y = -(long) halo_widths.y; y < (long)(block_sizes.y + halo_widths.y); y++)
        for(long x = -(long) halo_widths.x; x < (long)(block_sizes.x + halo_widths.x); x++){
            const int offset_x = x < 0 ? -1 : x >= (long) block_sizes.x ? 1 : 0;
            const int offset_y = y < 0 ? -1 : y >= (long) block_sizes.y ? 1 : 0;
            const int offset_z = z < 0 ? -1 : z >= (long) block_sizes.z ? 1 : 0;
            const long volume_x = x + (long)(i * block_sizes.x);
            const long volume_y = y + (long)(j * block_sizes.y);
            const long volume_z = z + (long)(k * block_sizes.z);
            const int is_corner = offset_x != 0 && offset_y != 0 && offset_z != 0;
            const int is_outside = volume_x < 0 || volume_y < 0 || volume_z < 0 ||
                volume_x >= (long) volume_dims.x || volume_y >= (long) volume_dims.y || volume_z >= (long) volume_dims.z;
            for(int field = 0; field < 2; field++){
                const FP fetched = halo_fetch(&view, (int) x, (int) y, (int) z, field, MISSING_VALUE);
                const FP expected = (is_corner || is_outside) ? MISSING_VALUE : global_value(volume_x, volume_y, volume_z, field);
                cr_assert(eq(flt, fetched, expected),
                    "block (%zu, %zu, %zu) point (%ld, %ld, %ld) field %d", i, j, k, x, y, z, field);
                checked_count++;
            }
        }
    }
    cr_log_info("checked %zu points", checked_count);

    for(size_t block = 0; block < g.total_blocks; block++){
        free(p_waves[block]);
        free(q_waves[block]);
        for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
            free(sent_slabs[block * NEIGHBOR_COUNT + table_index]);
        }
    }
    free(p_waves);
    free(q_waves);
    free(sent_slabs);
}

Test(halo, round_trip_slabs){
    round_trip((uvec3_t){12, 16, 40}, (uvec3_t){1, 1, 5});
}

Test(halo, round_trip_pencils){
    round_trip((uvec3_t){12, 24, 32}, (uvec3_t){1, 3, 4});
}

Test(halo, round_trip_cubes){
    round_trip((uvec3_t){24, 32, 40}, (uvec3_t){3, 4, 5});
}

Test(halo, round_trip_thinnest_blocks){
    // block size = R, the slab is the whole block
    round_trip((uvec3_t){12, 12, 12}, (uvec3_t){3, 3, 3});
}

// ---------------------------------------------------------------------------
// stencil, exact (up to rounding) on polynomials the 8th order scheme can represent
// ---------------------------------------------------------------------------

#define WIDTH 9
#define PLANE_INDEX(x, y) ((x) + WIDTH * (y))

Test(stencil, second_derivative_of_quadratic){
    FP line[WIDTH];
    for(int x = 0; x < WIDTH; x++){
        const FP position = (FP)(x - 4);
        line[x] = FP_LIT(3.0) * position * position + FP_LIT(2.0) * position + FP_LIT(1.0);
    }
    const FP derivative = stencil_second_derivative(&line[4], 1, FP_LIT(1.0));
    cr_expect(epsilon_eq(flt, derivative, FP_LIT(6.0), 1e-4), "got %f", derivative);
}

Test(stencil, second_derivative_spacing_and_stride){
    // f = x^2 stored every 2 elements with spacing 0.5: f'' = 2
    FP line[2 * WIDTH];
    for(int x = 0; x < WIDTH; x++){
        const FP position = (FP)(x - 4) * FP_LIT(0.5);
        line[2 * x] = position * position;
        line[2 * x + 1] = NAN;
    }
    const FP derivative = stencil_second_derivative(&line[8], 2, FP_LIT(1.0) / (FP_LIT(0.5) * FP_LIT(0.5)));
    cr_expect(epsilon_eq(flt, derivative, FP_LIT(2.0), 1e-4), "got %f", derivative);
}

Test(stencil, cross_derivative_of_bilinear){
    FP plane[WIDTH * WIDTH];
    for(int y = 0; y < WIDTH; y++)
    for(int x = 0; x < WIDTH; x++){
        const FP position_x = (FP)(x - 4), position_y = (FP)(y - 4);
        plane[PLANE_INDEX(x, y)] = FP_LIT(5.0) * position_x * position_y + position_x * position_x + FP_LIT(7.0) * position_y;
    }
    const FP* center = &plane[PLANE_INDEX(4, 4)];
    cr_expect(epsilon_eq(flt, stencil_cross_derivative(center, 1, WIDTH, FP_LIT(1.0)), FP_LIT(5.0), 1e-4));
    // symmetric in the order of the axes
    cr_expect(epsilon_eq(flt, stencil_cross_derivative(center, WIDTH, 1, FP_LIT(1.0)), FP_LIT(5.0), 1e-4));
}

Test(stencil, constant_is_zero){
    FP plane[WIDTH * WIDTH];
    for(int i = 0; i < WIDTH * WIDTH; i++) plane[i] = FP_LIT(42.0);
    const FP* center = &plane[PLANE_INDEX(4, 4)];
    cr_expect(epsilon_eq(flt, stencil_second_derivative(center, 1, FP_LIT(1.0)), FP_LIT(0.0), 1e-3));
    cr_expect(epsilon_eq(flt, stencil_cross_derivative(center, 1, WIDTH, FP_LIT(1.0)), FP_LIT(0.0), 1e-3));
}
