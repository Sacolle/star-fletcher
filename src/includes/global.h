#ifndef GUARD_GLOBAL
#define GUARD_GLOBAL

#include <stddef.h>
#include <stdint.h>

#include "floatingpoint.h"

// a triple of sizes, counts or coordinates. 32 bits per axis is plenty for one axis, but products of
// the components (amount of points of a volume, linear indices) can pass 2^32: compute them in size_t
typedef struct uvec3 { uint32_t x, y, z; } uvec3_t;

// a triple of signed coordinates, for points that can be outside of a block (negative or past its end)
typedef struct ivec3 { int32_t x, y, z; } ivec3_t;

struct global_data {
  // dimensions of the whole volume
  uvec3_t volume_dims;
  // amount of blocks in each axis
  uvec3_t num_blocks;
  // size of a block in each axis
  uvec3_t block_sizes;
  // halo width in each axis: the stencil radius on split axes, 0 otherwise
  uvec3_t halo_widths;
  // amount of points in a block, block_sizes.x * block_sizes.y * block_sizes.z
  size_t block_size;
  // amount of blocks in the volume, num_blocks.x * num_blocks.y * num_blocks.z
  size_t total_blocks;
  // simulation time between two outputs
  FP dt_output;
  // threads per CUDA block
  uvec3_t cuda_threads;
};

extern struct global_data g;

#endif
