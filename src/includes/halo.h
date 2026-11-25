#ifndef GUARD_HALO
#define GUARD_HALO

#include <stddef.h>
#include <stdint.h>

#include "err.h"
#include "floatingpoint.h"
#include "global.h"
#include "indexing.h"

#ifdef __cplusplus
extern "C" {
#endif

#define STENCIL_RADIUS 4

// the slot of the neighbor at `offset`, each component in {-1, 0, 1}, in a 3x3x3 table
// the center slot, offset (0, 0, 0), is the block itself
#define NEIGHBOR_SLOT_COUNT 27
#define NEIGHBOR_SLOT_SELF 13
HOST_DEVICE static inline int neighbor_slot(ivec3_t offset){
    return (offset.x + 1) + 3 * (offset.y + 1) + 9 * (offset.z + 1);
}

// a rectangular part of a 3D array that is copied between a block, a slab and the scratch buffer
typedef struct copy_region { uvec3_t origin; uvec3_t extent; } copy_region_t;

// neighbors a block can have: 6 faces + 12 edges, the 3D corners are never used by the stencil
#define NEIGHBOR_COUNT 18

// offsets possible from the central value
extern const ivec3_t NEIGHBOR_OFFSETS[NEIGHBOR_COUNT];

// inititalizes the halo_widths with the dimensions of each halo
// 0 if there is no split in that axis and STENCIL_RADIUS if it has
// so a segmentation of (1, 1, 2) (slab) results in a halo_widths of (0, 0, STENCIL_RADIUS)
err_t halo_init(const uvec3_t num_blocks, const uvec3_t block_sizes, uvec3_t* halo_widths);

// table entry of the opposite offset. the table is symmetric: entry n and entry NEIGHBOR_COUNT - 1 - n
// are opposite offsets. the slab a block receives from the neighbor at entry n is the one that
// neighbor sends towards entry halo_opposite_neighbor(n)
HOST_DEVICE static inline int halo_opposite_neighbor(int table_index){
    return NEIGHBOR_COUNT - 1 - table_index;
}

// the bit of entry n of NEIGHBOR_OFFSETS in a neighbor mask
#define NBIT(n) (UINT32_C(1) << (n))

// every entry of NEIGHBOR_OFFSETS
#define NEIGHBORS_ALL    (NBIT(NEIGHBOR_COUNT) - 1)

// masks for where the offsets are acording to their axis
// so NEIGHBORS_LOW_X are the neighbors for which the offsets' x coordinate is
// -1, so all five on the underside
#define NEIGHBORS_LOW_X  (NBIT(1) | NBIT(5) | NBIT(8) | NBIT(10) | NBIT(14))
#define NEIGHBORS_HIGH_X (NBIT(3) | NBIT(7) | NBIT(9) | NBIT(12) | NBIT(16))
#define NEIGHBORS_LOW_Y  (NBIT(0) | NBIT(5) | NBIT(6) | NBIT(7) | NBIT(13))
#define NEIGHBORS_HIGH_Y (NBIT(4) | NBIT(10) | NBIT(11) | NBIT(12) | NBIT(17))
#define NEIGHBORS_LOW_Z  (NBIT(0) | NBIT(1) | NBIT(2) | NBIT(3) | NBIT(4))
#define NEIGHBORS_HIGH_Z (NBIT(13) | NBIT(14) | NBIT(15) | NBIT(16) | NBIT(17))

// bit n is set if the block (i, j, k) has the neighbor at entry n of
// NEIGHBOR_OFFSETS. it is used to see which neighbors a block has, so a block
// at the corner (+x, +y, +z) will not have the neighbors at HIGH_X, HIGH_Y,
// HIGH_Z 
uint32_t halo_neighbor_mask(size_t i, size_t j, size_t k);

// given an offset from the NEIGHBOR_OFFSETS, what are the dimensions of that
// halo slab. If the segmentation is (1, 1, 4) and the halo is at (0, 0, -1),
// then the slab will be dimensioned as (block_sizes.x, block_sizes.y,
// STENCIL_RADIUS)
uvec3_t slab_dims(ivec3_t offset);

// the inner part of a block that has to be copied to the slab. If is at -1,
// then is (0..4), if 0, then is (0..size) if 1 then (size - 4..size)
// Copy region from inside the block
copy_region_t slab_source_region(ivec3_t offset);

// the region where the slab coming from the neighbor at the offset lands, in the coordinates of the
// padded scratch buffer (block coordinate + halo width)
// per axis: offset -1 -> [0, R), offset +1 -> [size + halo, size + halo + R), offset 0 -> [halo, halo + size)
// Copy region from the extended blockm with the halo padding
copy_region_t slab_halo_region(ivec3_t offset);

// copies the points of `region` out of `source` (strides ldy, ldz) into the start of `destination`
void extract_region(
    const FP* source, size_t ldy, size_t ldz, copy_region_t region,
    FP* destination, size_t destination_ldy, size_t destination_ldz
);

// copies the block data into a slab that holds p at field 0 and q at field 1
void halo_pack_slab(
    const FP* p_wave, const FP* q_wave, size_t ldy, size_t ldz,
    ivec3_t offset,
    FP* slab, size_t slab_ldy, size_t slab_ldz, size_t slab_ldt
);


// a block or a slab, with one pointer per field (p wave, q wave), NULL if it does not exist.
// the strides are 32 bits, a block or a slab has less than 2^32 points (main checks it)
struct halo_slab {
    const FP* fields[2];
    uint32_t ldy, ldz;
};

// the block itself in the slot NEIGHBOR_SLOT_SELF, and in the slot neighbor_slot(offset)
// the slab coming from the neighbor at that offset
struct halo_view {
    // signed, halo_fetch compares it with coordinates that can be negative
    ivec3_t block_sizes;
    struct halo_slab slabs[NEIGHBOR_SLOT_COUNT];
};

// clears the view for a block of the given size
void halo_view_init(struct halo_view* view, uvec3_t block_sizes);

// sets the block itself, p and q come from different handles but share the strides
void halo_view_set_self(struct halo_view* view, const FP* p_wave, const FP* q_wave, size_t ldy, size_t ldz);

// sets the slab coming from the neighbor at the offset, with p at field 0 and q at field 1
void halo_view_set_slab(
    struct halo_view* view, ivec3_t offset,
    const FP* slab, size_t ldy, size_t ldz, size_t ldt
);

// the value of `field` (0 = p wave, 1 = q wave) at the point (x, y, z) in block
// coordinates
HOST_DEVICE static inline FP halo_fetch(const struct halo_view* view, int x, int y, int z, int field, FP missing_value){
    // calculates to which block the data has to come from in form of a offset
    const ivec3_t offset = {
        (x < 0)                      ? -1
        : (x >= view->block_sizes.x) ? 1
                                     : 0,
        (y < 0)                      ? -1
        : (y >= view->block_sizes.y) ? 1
                                     : 0,
        (z < 0)                      ? -1
        : (z >= view->block_sizes.z) ? 1
                                     : 0,
    };
    // obtain that slab, which might be the center block
    const struct halo_slab* slab = &view->slabs[neighbor_slot(offset)];
    if(slab->fields[field] == NULL){
        return missing_value;
    }
    // shift into the slab coordinates: -R..-1 -> 0..R-1, size..size+R-1 -> 0..R-1
    const int slab_x = (offset.x < 0) ? x + STENCIL_RADIUS : (offset.x > 0) ? x - view->block_sizes.x : x;
    const int slab_y = (offset.y < 0) ? y + STENCIL_RADIUS : (offset.y > 0) ? y - view->block_sizes.y : y;
    const int slab_z = (offset.z < 0) ? z + STENCIL_RADIUS : (offset.z > 0) ? z - view->block_sizes.z : z;
    //get the value at the computed point
    return slab
        ->fields[field][idx3_u32(slab_x, slab_y, slab_z, slab->ldy, slab->ldz)];
}

#ifdef __cplusplus
}
#endif

#endif
