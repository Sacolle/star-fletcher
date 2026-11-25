#ifndef GUARD_STENCIL
#define GUARD_STENCIL

#include <stddef.h>
#include <stdint.h>

#include "floatingpoint.h"
#include "indexing.h"


// coefficients of the second derivative
#define K0 FP_LIT(-2.84722222222222222222)  // -205/72
#define K1 FP_LIT(1.6)                      // 8/5
#define K2 FP_LIT(-0.2)                     // -1/5
#define K3 FP_LIT(0.02539682539682539682)   // 8/315
#define K4 FP_LIT(-0.00178571428571428571)  // -1/560

// coefficients of the cross second derivative, Lab = La * Lb of the first derivative
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

// second derivative along the axis with the given stride, times 1 / spacing^2
HOST_DEVICE static inline FP stencil_second_derivative(const FP* center, int32_t stride, FP inverse_spacing_squared){
    return (K0 * center[0] +
            K1 * (center[stride] + center[-stride]) +
            K2 * (center[2 * stride] + center[-2 * stride]) +
            K3 * (center[3 * stride] + center[-3 * stride]) +
            K4 * (center[4 * stride] + center[-4 * stride])) * inverse_spacing_squared;
}

#define STENCIL_CROSS(distance_2, distance_1)                     \
      center[(distance_2) * stride_2 + (distance_1) * stride_1]           \
    - center[(distance_2) * stride_2 - (distance_1) * stride_1]           \
    - center[-(distance_2) * stride_2 + (distance_1) * stride_1]          \
    + center[-(distance_2) * stride_2 - (distance_1) * stride_1]

// cross derivative along the axes with strides stride_1 and stride_2, times 1 / (spacing_1 * spacing_2)
HOST_DEVICE static inline FP stencil_cross_derivative(
    const FP* center, int32_t stride_1, int32_t stride_2, FP inverse_spacing_product
){
    return (L11 * (STENCIL_CROSS(1, 1)) +
            L12 * (STENCIL_CROSS(1, 2) + STENCIL_CROSS(2, 1)) +
            L13 * (STENCIL_CROSS(1, 3) + STENCIL_CROSS(3, 1)) +
            L14 * (STENCIL_CROSS(1, 4) + STENCIL_CROSS(4, 1)) +
            L22 * (STENCIL_CROSS(2, 2)) +
            L23 * (STENCIL_CROSS(2, 3) + STENCIL_CROSS(3, 2)) +
            L24 * (STENCIL_CROSS(2, 4) + STENCIL_CROSS(4, 2)) +
            L33 * (STENCIL_CROSS(3, 3)) +
            L34 * (STENCIL_CROSS(3, 4) + STENCIL_CROSS(4, 3)) +
            L44 * (STENCIL_CROSS(4, 4))) * inverse_spacing_product;
}

#undef STENCIL_CROSS

#endif
