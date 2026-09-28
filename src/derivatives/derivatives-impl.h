#ifndef _DERIVATIVES_GUARD_
#define _DERIVATIVES_GUARD_

#include <stddef.h>
#include <stdint.h>

#include "floatingpoint.h"

#ifdef CUDA_CODE
#define ATTRIBUTE __device__
#else
#define ATTRIBUTE
#endif

// Half width of the eighth order stencil. It is also the thickness of the face/edge
// sub-handles produced by g_cube_face_filter (BORDER_WIDTH, checked in main.c).
#define STENCIL_RADIUS 4

// Strided view over a cube, or over a face/edge slice of a neighbor cube.
// Each view keeps its own strides: in place slices share the parent's, but a
// separately allocated slice (e.g. on another device) is compact.
typedef struct {
    const FP* ptr;
    int32_t ldy;
    int32_t ldz;
} block_view_t;

// We are putting the neighboring blocks into a array of size NEIGHBORHOOD_SIZE
// So that we ease indexing. We still want to use notation similar to the main function.
// Therefore the center element of this cube is the (0,0,0) and we index as increments from it
// The block above is (0, -1, 0) (Y axis grows down) and the block to the right is (1, 0, 0)
#define NEIGHBOR_IDX(dx, dy, dz) (((dz) + 1) * 9 + ((dy) + 1) * 3 + ((dx) + 1))
#define NEIGHBORHOOD_SIZE 27

#define L1 FP_LIT(0.8)                    // 4/5
#define L2 FP_LIT(-0.2)                   // -1/5
#define L3 FP_LIT(0.0380952380952381)     // 4/105
#define L4 FP_LIT(-0.0035714285714285713) // -1/280

// eight order finite differences coefficients of the cross second derivative

#define L11 FP_LIT(0.64)                    // L1*L1
#define L12 FP_LIT(-0.16)                   // L1*L2
#define L13 FP_LIT(0.03047619047619047618)  // L1*L3
#define L14 FP_LIT(-0.00285714285714285713) // L1*L4
#define L22 FP_LIT(0.04)                    // L2*L2
#define L23 FP_LIT(-0.00761904761904761904) // L2*L3
#define L24 FP_LIT(0.00071428571428571428)  // L2*L4
#define L33 FP_LIT(0.00145124716553287981)  // L3*L3
#define L34 FP_LIT(-0.00013605442176870748) // L3*L4
#define L44 FP_LIT(0.00001275510204081632)  // L4*L4

// eight order finite differences coefficients of the second derivative
#define K0 FP_LIT(-2.84722222222222222222) // -205/72
#define K1 FP_LIT(1.6)                     // 8/5
#define K2 FP_LIT(-0.2)                    // -1/5
#define K3 FP_LIT(0.02539682539682539682)  // 8/315
#define K4 FP_LIT(-0.00178571428571428571) // -1/560


// Second derivative along one axis at (x, y, z) of the central cube of `neighborhood`.
ATTRIBUTE FP snd_deriv_x(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP d2inv);
ATTRIBUTE FP snd_deriv_y(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP d2inv);
ATTRIBUTE FP snd_deriv_z(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP d2inv);

// Cross derivative on a pair of axes at (x, y, z) of the central cube of `neighborhood`.
ATTRIBUTE FP cross_deriv_xy(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP dinv);
ATTRIBUTE FP cross_deriv_yz(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP dinv);
ATTRIBUTE FP cross_deriv_xz(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP dinv);

// add the impl
#ifdef CODE_IMPL

// Which neighbor a point reaches along one axis: the low side (neg), the high side (pos),
// or none (center). Needs cube_width >= 2 * STENCIL_RADIUS (checked in main.c).
enum { SIDE_CENTER = 0, SIDE_POS = 1, SIDE_NEG = 2 };

ATTRIBUTE static inline int stencil_side(int32_t coord, int32_t cube_width){
    if(coord < STENCIL_RADIUS) return SIDE_NEG;
    if(coord > cube_width - 1 - STENCIL_RADIUS) return SIDE_POS;
    return SIDE_CENTER;
}

#include "./cross-deriv.gen.c"

// `depth` is how many points separate base_idx from the high border of the central cube.
// `border_idx` is the point of `pos` right after that border; both views are walked with
// their own stride along the derivative axis.
ATTRIBUTE static FP snd_deriv_pos_impl(
    const FP* central, const FP* pos,
    int32_t depth, int32_t base_idx, int32_t stride,
    int32_t border_idx, int32_t pos_stride,
    FP d2inv
){
    switch (depth)
    {
    case 0:
        /* right at the border */
        return (
            K0 * central[base_idx] +
            K1 * (pos[border_idx + 0 * pos_stride] + central[base_idx - 1 * stride]) +
            K2 * (pos[border_idx + 1 * pos_stride] + central[base_idx - 2 * stride]) +
            K3 * (pos[border_idx + 2 * pos_stride] + central[base_idx - 3 * stride]) +
            K4 * (pos[border_idx + 3 * pos_stride] + central[base_idx - 4 * stride])
        ) * (d2inv);
    case 1:
        /* right before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (pos[border_idx + 0 * pos_stride] + central[base_idx - 2 * stride]) +
            K3 * (pos[border_idx + 1 * pos_stride] + central[base_idx - 3 * stride]) +
            K4 * (pos[border_idx + 2 * pos_stride] + central[base_idx - 4 * stride])
        ) * (d2inv);
    case 2:
        /* 2 before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + central[base_idx - 2 * stride]) +
            K3 * (pos[border_idx + 0 * pos_stride] + central[base_idx - 3 * stride]) +
            K4 * (pos[border_idx + 1 * pos_stride] + central[base_idx - 4 * stride])
        ) * (d2inv);
    case 3:
        /* 3 before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + central[base_idx - 2 * stride]) +
            K3 * (central[base_idx + 3 * stride] + central[base_idx - 3 * stride]) +
            K4 * (pos[border_idx + 0 * pos_stride] + central[base_idx - 4 * stride])
        ) * (d2inv);
    default:
        /* 4 and less before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + central[base_idx - 2 * stride]) +
            K3 * (central[base_idx + 3 * stride] + central[base_idx - 3 * stride]) +
            K4 * (central[base_idx + 4 * stride] + central[base_idx - 4 * stride])
        ) * (d2inv);
    }
}

// Mirror of snd_deriv_pos_impl for the low border: `depth` points separate base_idx from
// it and `border_idx` is the point of `neg` right before it.
ATTRIBUTE static FP snd_deriv_neg_impl(
    const FP* central, const FP* neg,
    int32_t depth, int32_t base_idx, int32_t stride,
    int32_t border_idx, int32_t neg_stride,
    FP d2inv
){
    switch (depth)
    {
    case 0:
        /* right at the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + neg[border_idx - 0 * neg_stride]) +
            K2 * (central[base_idx + 2 * stride] + neg[border_idx - 1 * neg_stride]) +
            K3 * (central[base_idx + 3 * stride] + neg[border_idx - 2 * neg_stride]) +
            K4 * (central[base_idx + 4 * stride] + neg[border_idx - 3 * neg_stride])
        ) * (d2inv);
    case 1:
        /* right before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + neg[border_idx - 0 * neg_stride]) +
            K3 * (central[base_idx + 3 * stride] + neg[border_idx - 1 * neg_stride]) +
            K4 * (central[base_idx + 4 * stride] + neg[border_idx - 2 * neg_stride])
        ) * (d2inv);
    case 2:
        /* 2 before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + central[base_idx - 2 * stride]) +
            K3 * (central[base_idx + 3 * stride] + neg[border_idx - 0 * neg_stride]) +
            K4 * (central[base_idx + 4 * stride] + neg[border_idx - 1 * neg_stride])
        ) * (d2inv);
    case 3:
        /* 3 before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + central[base_idx - 2 * stride]) +
            K3 * (central[base_idx + 3 * stride] + central[base_idx - 3 * stride]) +
            K4 * (central[base_idx + 4 * stride] + neg[border_idx - 0 * neg_stride])
        ) * (d2inv);
    default:
        /* 4 and less before the border */
        return (
            K0 * central[base_idx] +
            K1 * (central[base_idx + 1 * stride] + central[base_idx - 1 * stride]) +
            K2 * (central[base_idx + 2 * stride] + central[base_idx - 2 * stride]) +
            K3 * (central[base_idx + 3 * stride] + central[base_idx - 3 * stride]) +
            K4 * (central[base_idx + 4 * stride] + central[base_idx - 4 * stride])
        ) * (d2inv);
    }
}

ATTRIBUTE FP snd_deriv_x(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP d2inv){
    const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
    const int32_t base_idx = x + y * center.ldy + z * center.ldz;

    if(x < STENCIL_RADIUS){
        const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
        const int32_t border_idx = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
        return snd_deriv_neg_impl(center.ptr, left.ptr, x, base_idx, 1, border_idx, 1, d2inv);
    }else{
        const block_view_t right = neighborhood[NEIGHBOR_IDX(+1, 0, 0)];
        const int32_t border_idx = y * right.ldy + z * right.ldz;
        return snd_deriv_pos_impl(center.ptr, right.ptr, cube_width - 1 - x, base_idx, 1, border_idx, 1, d2inv);
    }
}

ATTRIBUTE FP snd_deriv_y(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP d2inv){
    const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
    const int32_t base_idx = x + y * center.ldy + z * center.ldz;

    if(y < STENCIL_RADIUS){
        const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
        const int32_t border_idx = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
        return snd_deriv_neg_impl(center.ptr, top.ptr, y, base_idx, center.ldy, border_idx, top.ldy, d2inv);
    }else{
        const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, +1, 0)];
        const int32_t border_idx = x + z * bottom.ldz;
        return snd_deriv_pos_impl(center.ptr, bottom.ptr, cube_width - 1 - y, base_idx, center.ldy, border_idx, bottom.ldy, d2inv);
    }
}

ATTRIBUTE FP snd_deriv_z(const block_view_t *neighborhood, int32_t x, int32_t y, int32_t z, int32_t cube_width, FP d2inv){
    const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
    const int32_t base_idx = x + y * center.ldy + z * center.ldz;

    if(z < STENCIL_RADIUS){
        const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
        const int32_t border_idx = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
        return snd_deriv_neg_impl(center.ptr, front.ptr, z, base_idx, center.ldz, border_idx, front.ldz, d2inv);
    }else{
        const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, +1)];
        const int32_t border_idx = x + y * back.ldy;
        return snd_deriv_pos_impl(center.ptr, back.ptr, cube_width - 1 - z, base_idx, center.ldz, border_idx, back.ldz, d2inv);
    }
}

#endif /* CODE_IMPL */

#endif /* _DERIVATIVES_GUARD_ */
