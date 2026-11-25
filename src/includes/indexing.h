#ifndef GUARD_INDEXING
#define GUARD_INDEXING

#include <stddef.h>

#include "global.h"

// functions marked with it can be called from both the host and CUDA kernels
#ifdef __CUDACC__
#define HOST_DEVICE __host__ __device__
#else
#define HOST_DEVICE
#endif

// linear index of the point (x, y, z) in a 3D array with strides ldy and ldz
HOST_DEVICE static inline size_t idx3(size_t x, size_t y, size_t z, size_t ldy, size_t ldz){
    return x + y * ldy + z * ldz;
}

// inverse of idx3 for a contiguous array of nx * ny * _ points
// x is the fastest axis
HOST_DEVICE static inline uvec3_t unindex(size_t idx, size_t nx, size_t ny){
    uvec3_t coords;
    coords.x = idx % nx;
    coords.y = (idx / nx) % ny;
    coords.z = idx / (nx * ny);
    return coords;
}

// 32 bit versions of idx3 and unindex, for the index math inside the GPU kernels: 64 bit integers take two
// registers each and 64 bit division is emulated on GPUs. only valid while the array has less than 2^32 points,
// which main checks for the blocks, the slabs and the scratch
HOST_DEVICE static inline uint32_t idx3_u32(uint32_t x, uint32_t y, uint32_t z, uint32_t ldy, uint32_t ldz){
    return x + y * ldy + z * ldz;
}

HOST_DEVICE static inline uvec3_t unindex_u32(uint32_t idx, uint32_t nx, uint32_t ny){
    uvec3_t coords;
    coords.x = idx % nx;
    coords.y = (idx / nx) % ny;
    coords.z = idx / (nx * ny);
    return coords;
}

// amount of points of a box of size v, computed in size_t so it does not overflow
HOST_DEVICE static inline size_t uvec3_prod(uvec3_t v){
    return (size_t) v.x * v.y * v.z;
}


// the linear index of the point (x, y, z) inside a block
static inline size_t local_idx(size_t x, size_t y, size_t z){
    return idx3(x, y, z, g.block_sizes.x, (size_t) g.block_sizes.x * g.block_sizes.y);
}

// the linear index of the block (i, j, k) in the volume
static inline size_t block_idx(size_t i, size_t j, size_t k){
    return idx3(i, j, k, g.num_blocks.x, (size_t) g.num_blocks.x * g.num_blocks.y);
}

// the linear index of the point (vx, vy, vz) in the whole volume
static inline size_t volume_idx(size_t vx, size_t vy, size_t vz){
    return idx3(vx, vy, vz, g.volume_dims.x, (size_t) g.volume_dims.x * g.volume_dims.y);
}

// turn the point (x, y, z) of block (i, j, k) into the linear index in the whole volume
static inline size_t local_to_volume_idx(size_t x, size_t y, size_t z, size_t i, size_t j, size_t k){
    return volume_idx(x + i * (size_t) g.block_sizes.x, y + j * (size_t) g.block_sizes.y, z + k * (size_t) g.block_sizes.z);
}

// turn a volume linear index into the index of the block that holds the point
static inline size_t volume_to_block_idx(size_t idx){
    const uvec3_t coords = unindex(idx, g.volume_dims.x, g.volume_dims.y);
    return block_idx(coords.x / g.block_sizes.x, coords.y / g.block_sizes.y, coords.z / g.block_sizes.z);
}

// turn a volume linear index into the index inside the block that holds the point
static inline size_t volume_to_local_idx(size_t idx){
    const uvec3_t coords = unindex(idx, g.volume_dims.x, g.volume_dims.y);
    return local_idx(coords.x % g.block_sizes.x, coords.y % g.block_sizes.y, coords.z % g.block_sizes.z);
}

#endif
