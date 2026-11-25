#include "halo.h"

#include <string.h>

// the 6 faces and 12 edges, z is the slowest axis and x the fastest.
// the 3D corners are left out, the stencil never needs (+-x, +-y, +-z) at once
const ivec3_t NEIGHBOR_OFFSETS[NEIGHBOR_COUNT] = {
    // offset.z = -1
    { 0, -1, -1}, {-1,  0, -1}, { 0,  0, -1}, { 1,  0, -1}, { 0,  1, -1},
    // offset.z = 0
    {-1, -1,  0}, { 0, -1,  0}, { 1, -1,  0},
    {-1,  0,  0},               { 1,  0,  0},
    {-1,  1,  0}, { 0,  1,  0}, { 1,  1,  0},
    // offset.z = +1
    { 0, -1,  1}, {-1,  0,  1}, { 0,  0,  1}, { 1,  0,  1}, { 0,  1,  1},
};

err_t halo_init(const uvec3_t num_blocks, const uvec3_t block_sizes, uvec3_t* halo_widths){
  //const uvec3_t num_blocks = g.num_blocks;
  //const uvec3_t block_sizes = g.block_sizes;
    const int is_split[3] = {num_blocks.x > 1, num_blocks.y > 1, num_blocks.z > 1};

    // if any of the segmentations fall bellow a width of 4, error
    if((is_split[0] && block_sizes.x < STENCIL_RADIUS) ||
       (is_split[1] && block_sizes.y < STENCIL_RADIUS) ||
       (is_split[2] && block_sizes.z < STENCIL_RADIUS)){
        return ME_COUNT_DONT_MATCH;
    }

    *halo_widths = (uvec3_t){
        is_split[0] ? STENCIL_RADIUS : 0,
        is_split[1] ? STENCIL_RADIUS : 0,
        is_split[2] ? STENCIL_RADIUS : 0,
    };
    return 0;
}

uint32_t halo_neighbor_mask(size_t i, size_t j, size_t k){
    uint32_t mask = NEIGHBORS_ALL;
    if(i == 0)                  mask &= ~NEIGHBORS_LOW_X;
    if(i == g.num_blocks.x - 1) mask &= ~NEIGHBORS_HIGH_X;
    if(j == 0)                  mask &= ~NEIGHBORS_LOW_Y;
    if(j == g.num_blocks.y - 1) mask &= ~NEIGHBORS_HIGH_Y;
    if(k == 0)                  mask &= ~NEIGHBORS_LOW_Z;
    if(k == g.num_blocks.z - 1) mask &= ~NEIGHBORS_HIGH_Z;
    return mask;
}

uvec3_t slab_dims(ivec3_t offset){
    const uvec3_t block_sizes = g.block_sizes;
    return (uvec3_t){
        offset.x ? STENCIL_RADIUS : block_sizes.x,
        offset.y ? STENCIL_RADIUS : block_sizes.y,
        offset.z ? STENCIL_RADIUS : block_sizes.z,
    };
}

// first point, along one axis, of the slab sent towards the offset
// if offset is 1, set to block_size - STENCIL_RADIUS, else set to 0
static inline size_t slab_source_start(int offset, size_t block_size){
    return offset > 0 ? block_size - STENCIL_RADIUS : 0;
}

copy_region_t slab_source_region(ivec3_t offset){
    const uvec3_t block_sizes = g.block_sizes;
    return (copy_region_t){
        .origin = {
            slab_source_start(offset.x, block_sizes.x),
            slab_source_start(offset.y, block_sizes.y),
            slab_source_start(offset.z, block_sizes.z),
        },
        .extent = slab_dims(offset),
    };
}

// first point, along one axis and in scratch coordinates, of the slab coming from the offset
static inline size_t slab_halo_start(int offset, size_t block_size, size_t halo_width){
    return offset < 0 ? 0 : offset > 0 ? block_size + halo_width : halo_width;
}

copy_region_t slab_halo_region(ivec3_t offset){
    const uvec3_t block_sizes = g.block_sizes;
    const uvec3_t halo_widths = g.halo_widths;
    return (copy_region_t){
        .origin = {
            slab_halo_start(offset.x, block_sizes.x, halo_widths.x),
            slab_halo_start(offset.y, block_sizes.y, halo_widths.y),
            slab_halo_start(offset.z, block_sizes.z, halo_widths.z),
        },
        .extent = slab_dims(offset),
    };
}

void extract_region(
    const FP* source, size_t ldy, size_t ldz, copy_region_t region,
    FP* destination, size_t destination_ldy, size_t destination_ldz
){
    for(size_t z = 0; z < region.extent.z; z++)
    for(size_t y = 0; y < region.extent.y; y++){
        const FP* row = source + idx3(region.origin.x, region.origin.y + y, region.origin.z + z, ldy, ldz);
        memcpy(destination + idx3(0, y, z, destination_ldy, destination_ldz), row, region.extent.x * sizeof(FP));
    }
}

// copies the waves to a slab
void halo_pack_slab(
    const FP* p_wave, const FP* q_wave, size_t ldy, size_t ldz,
    ivec3_t offset,
    FP* slab, size_t slab_ldy, size_t slab_ldz, size_t slab_ldt
){
    const copy_region_t region = slab_source_region(offset);
    extract_region(p_wave, ldy, ldz, region, slab, slab_ldy, slab_ldz);
    extract_region(q_wave, ldy, ldz, region, slab + slab_ldt, slab_ldy, slab_ldz);
}

void halo_view_init(struct halo_view* view, uvec3_t block_sizes){
    memset(view, 0, sizeof(*view));
    view->block_sizes = (ivec3_t){ (int32_t) block_sizes.x, (int32_t) block_sizes.y, (int32_t) block_sizes.z };
}

void halo_view_set_self(struct halo_view* view, const FP* p_wave, const FP* q_wave, size_t ldy, size_t ldz){
    view->slabs[NEIGHBOR_SLOT_SELF] = (struct halo_slab){ .fields = {p_wave, q_wave}, .ldy = (uint32_t) ldy, .ldz = (uint32_t) ldz };
}

void halo_view_set_slab(
    struct halo_view* view, ivec3_t offset,
    const FP* slab, size_t ldy, size_t ldz, size_t ldt
){
    view->slabs[neighbor_slot(offset)] =
        (struct halo_slab){ .fields = {slab, slab + ldt}, .ldy = (uint32_t) ldy, .ldz = (uint32_t) ldz };
}
