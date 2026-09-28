#ifndef GUARD_WAVE_NEIGHBORHOOD
#define GUARD_WAVE_NEIGHBORHOOD

// Translates the buffers of one wave in the RTM task (see the handle order in main.c) into the
// neighborhood table used by the derivatives (block_view_t, NEIGHBOR_IDX in derivatives-impl.h).
// Shared by the CPU kernel (kernel.c) and the CUDA host wrapper (cuda/kernel.cu), so it is
// written in the common subset of C and C++.

#include <stdint.h>
#include <starpu.h>

#include "derivatives.h"

// Offset (dx, dy, dz) of the cube each t-1 buffer of a wave comes from, in submission order.
// The first one is the central cube itself, the others are face/edge sub-handles of the neighbors.
#define WAVE_NEIGHBOR_BUFFERS 19
static const int8_t wave_neighbor_offsets[WAVE_NEIGHBOR_BUFFERS][3] = {
    { 0,  0,  0},
    // layer k - 1
    { 0,  0, -1}, { 0, -1, -1}, {-1,  0, -1}, {+1,  0, -1}, { 0, +1, -1},
    // layer k
    {-1, -1,  0}, { 0, -1,  0}, {+1, -1,  0}, {-1,  0,  0},
    {+1,  0,  0}, {-1, +1,  0}, { 0, +1,  0}, {+1, +1,  0},
    // layer k + 1
    { 0,  0, +1}, { 0, -1, +1}, {-1,  0, +1}, {+1,  0, +1}, { 0, +1, +1},
};

// Builds the neighborhood table of one wave from its WAVE_NEIGHBOR_BUFFERS buffers
// starting at descr[first_buffer]. Unused entries (the corners) are left zeroed.
static inline void fill_wave_neighborhood(block_view_t neighborhood[NEIGHBORHOOD_SIZE], void *descr[], int first_buffer){
    for(int i = 0; i < NEIGHBORHOOD_SIZE; i++){
        neighborhood[i].ptr = NULL;
        neighborhood[i].ldy = 0;
        neighborhood[i].ldz = 0;
    }
    for(int i = 0; i < WAVE_NEIGHBOR_BUFFERS; i++){
        void *buffer = descr[first_buffer + i];
        const int8_t *offset = wave_neighbor_offsets[i];
        block_view_t *view = &neighborhood[NEIGHBOR_IDX(offset[0], offset[1], offset[2])];
        view->ptr = (const FP*) STARPU_BLOCK_GET_PTR(buffer);
        view->ldy = (int32_t) STARPU_BLOCK_GET_LDY(buffer);
        view->ldz = (int32_t) STARPU_BLOCK_GET_LDZ(buffer);
    }
}

#endif
