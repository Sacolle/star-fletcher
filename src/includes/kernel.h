#ifndef GUARD_KERNEL
#define GUARD_KERNEL

#include <stddef.h>
#include <stdint.h>
#include "floatingpoint.h"
#include "err.h"

#ifdef __cplusplus
extern "C" {
#endif

// fields of the medium tensor of a block, field f holds its values from f * ldt to (f + 1) * ldt
enum medium_field {
    MEDIUM_CH1DXX,
    MEDIUM_CH1DYY,
    MEDIUM_CH1DZZ,
    MEDIUM_CH1DXY,
    MEDIUM_CH1DYZ,
    MEDIUM_CH1DXZ,
    MEDIUM_V2PX,
    MEDIUM_V2PZ,
    MEDIUM_V2SZ,
    MEDIUM_V2PN,
    MEDIUM_FIELD_COUNT
};

// arguments of the task that computes the level t of a block
typedef struct rtm_args {
    // the points of the block that are computed, [start, end) in each axis, the others are set to 0
    size_t x_start;
    size_t y_start;
    size_t z_start;
    size_t x_end;
    size_t y_end;
    size_t z_end;
    FP dx;
    FP dy;
    FP dz;
    FP dt;
    // bit n is set if the block has the neighbor at entry n of NEIGHBOR_OFFSETS,
    // the slabs of the task come in the order of that table
    uint32_t neighbor_mask;
    // if the block holds the source: its point inside the block and the value added to the level t
    int has_source;
    size_t source_x;
    size_t source_y;
    size_t source_z;
    FP source_value;
    // the level being computed, only used by debug prints
    int64_t t;
} rtm_args_t;

// mallocs a copy of `values`, to be used as the cl_arg of a task
err_t make_rtm_args(rtm_args_t** rtm_args, const rtm_args_t* values);

// arguments of the task that copies a level of a block into the slabs it sends
typedef struct pack_args {
    uint32_t neighbor_mask;
} pack_args_t;

err_t make_pack_args(pack_args_t** pack_args, uint32_t neighbor_mask);

// the position of all data handles in the handles of a task have to be
// consistent. So the functions are helpers to, based on the amount of neighbors
// a block has, where to put each data handle. The first eight are always there:
// the medium, p and q at t (the level computed), at t - 1 and at t - 2, and the scratch
enum rtm_handle_pos {
    RTM_HANDLE_POS_MEDIUM,
    RTM_HANDLE_POS_P_T,
    RTM_HANDLE_POS_Q_T,
    RTM_HANDLE_POS_P_T1,
    RTM_HANDLE_POS_Q_T1,
    RTM_HANDLE_POS_P_T2,
    RTM_HANDLE_POS_Q_T2,
    RTM_HANDLE_POS_SCRATCH,
    // amount of handles every task has, the slabs come after them
    RTM_HANDLE_POS_FIXED_COUNT
};
// then the slabs, in the order of NEIGHBOR_OFFSETS: the received ones and after them the sent ones
static inline int rtm_handle_pos_received_slab(int slab_position){ return RTM_HANDLE_POS_FIXED_COUNT + slab_position; }
static inline int rtm_handle_pos_sent_slab(int neighbor_count, int slab_position){ return RTM_HANDLE_POS_FIXED_COUNT + neighbor_count + slab_position; }
static inline int rtm_handle_pos_count(int neighbor_count){ return RTM_HANDLE_POS_FIXED_COUNT + 2 * neighbor_count; }

// position of each data handle in the pack task: p and q of the level, then the sent slabs
enum pack_handle_pos {
    PACK_HANDLE_POS_P,
    PACK_HANDLE_POS_Q,
    // amount of handles every pack task has, the slabs come after them
    PACK_HANDLE_POS_FIXED_COUNT
};
static inline int pack_handle_pos_sent_slab(int slab_position){ return PACK_HANDLE_POS_FIXED_COUNT + slab_position; }
static inline int pack_handle_pos_count(int neighbor_count){ return PACK_HANDLE_POS_FIXED_COUNT + neighbor_count; }

// computes the level t of a block: gathers the level t - 1 and the received slabs into the scratch,
// runs the stencil, adds the source and copies the borders of the level t into the sent slabs
void rtm_kernel(void *descr[], void *cl_args);

// copies a level of a block into the slabs it sends
void pack_kernel(void *descr[], void *cl_args);

// the largest thread block the CUDA stencil kernel can be launched with (it is limited by the registers
// the kernel uses), or -1 if the CUDA device can not be queried. only defined in CUDA builds
int stencil_max_threads_per_block(void);

#ifdef __cplusplus
}
#endif

#endif
