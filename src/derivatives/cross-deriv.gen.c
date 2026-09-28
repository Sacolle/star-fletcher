/*
    Este é um arquivo gerado automaticamente pelo script `cross-deriv-gen.py`
    A meta é lidar discernir os casos antes da computação, de forma que na hora do cáculo,
    este possa tomar uma forma idêntica a do fletcher base.

    Na pipeline de compilação, ele é incluido no arquivo derivatives-impl.h. Para cada par de
    eixos (xy, yz, xz) há uma função por quadrante (lado de cada eixo que o ponto toca) e uma
    função cross_deriv_<par> que escolhe o quadrante.
*/
#include <assert.h>
#include <stdio.h>

#ifdef RELEASE
#define UNREACHABLE __builtin_unreachable()
#else
#define UNREACHABLE assert(0 && "Unreachable!"); return 0.0
#endif

#include "floatingpoint.h"
#include "derivatives.h"


ATTRIBUTE static FP cross_deriv_xy_neg_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
	const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
	const block_view_t left_top = neighborhood[NEIGHBOR_IDX(-1, -1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = x;
	const int32_t border_1 = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
	const int32_t depth2 = y;
	const int32_t border_2 = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
	const int32_t border_3 = STENCIL_RADIUS - 1 + (STENCIL_RADIUS - 1) * left_top.ldy + z * left_top.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - top.ptr[border_2 + 1] + left_top.ptr[border_3]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - top.ptr[border_2 + 2] + left_top.ptr[border_3 - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 1] + left_top.ptr[border_3 - 1 * left_top.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 1] + left_top.ptr[border_3 - 2 * left_top.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 3 * top.ldy + 1] + left_top.ptr[border_3 - 3 * left_top.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 2] + left_top.ptr[border_3 - 1 * left_top.ldy - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 2] + left_top.ptr[border_3 - 2 * left_top.ldy - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 2] + left_top.ptr[border_3 - 3 * left_top.ldy - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 3] + left_top.ptr[border_3 - 2 * left_top.ldy - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 3] + left_top.ptr[border_3 - 3 * left_top.ldy - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 3] - top.ptr[border_2 - 3 * top.ldy + 4] + left_top.ptr[border_3 - 3 * left_top.ldy - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 + 1] + left_top.ptr[border_3]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 1] + left_top.ptr[border_3 - 1 * left_top.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 1] + left_top.ptr[border_3 - 2 * left_top.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - top.ptr[border_2 + 2] + left_top.ptr[border_3 - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 2] + left_top.ptr[border_3 - 1 * left_top.ldy - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 2] + left_top.ptr[border_3 - 2 * left_top.ldy - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 3] + left_top.ptr[border_3 - 2 * left_top.ldy - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 + 1] + left_top.ptr[border_3]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 1] + left_top.ptr[border_3 - 1 * left_top.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 + 2] + left_top.ptr[border_3 - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 2] + left_top.ptr[border_3 - 1 * left_top.ldy - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 2] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 3] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 1] + left.ptr[border_1 - 3 * left.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 + 1] + left_top.ptr[border_3]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 + 2] + left_top.ptr[border_3 - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 3] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 3]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - top.ptr[border_2 + 2] + left_top.ptr[border_3] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 2] + left_top.ptr[border_3 - 1 * left_top.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 2] + left_top.ptr[border_3 - 2 * left_top.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 3 * top.ldy + 2] + left_top.ptr[border_3 - 3 * left_top.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 3] + left_top.ptr[border_3 - 2 * left_top.ldy - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 3] + left_top.ptr[border_3 - 3 * left_top.ldy - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 4] + left_top.ptr[border_3 - 3 * left_top.ldy - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 + 2] + left_top.ptr[border_3]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 2] + left_top.ptr[border_3 - 1 * left_top.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 2] + left_top.ptr[border_3 - 2 * left_top.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 3] + left_top.ptr[border_3 - 2 * left_top.ldy - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 + 2] + left_top.ptr[border_3]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 2] + left_top.ptr[border_3 - 1 * left_top.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 2] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 + 2] + left_top.ptr[border_3]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 + 3] + left_top.ptr[border_3 - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 2] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 2]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - top.ptr[border_2 + 3] + left_top.ptr[border_3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 3] + left_top.ptr[border_3 - 2 * left_top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 3 * top.ldy + 3] + left_top.ptr[border_3 - 3 * left_top.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 4] + left_top.ptr[border_3 - 3 * left_top.ldy - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 + 3] + left_top.ptr[border_3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 3] + left_top.ptr[border_3 - 2 * left_top.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 + 3] + left_top.ptr[border_3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 1] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 3] + left_top.ptr[border_3 - 1 * left_top.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 + 3] + left_top.ptr[border_3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 1] - top.ptr[border_2 + 4] + left_top.ptr[border_3 - 1]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - top.ptr[border_2 + 4] + left_top.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 3 * top.ldy + 3] + top.ptr[border_2 - 3 * top.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 3 * top.ldy + 4] + left_top.ptr[border_3 - 3 * left_top.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - top.ptr[border_2 + 4] + left_top.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 2 * top.ldy + 4] + left_top.ptr[border_3 - 2 * left_top.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy] - top.ptr[border_2 + 4] + left_top.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 - 1 * top.ldy + 4] + left_top.ptr[border_3 - 1 * left_top.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy] - top.ptr[border_2 + 4] + left_top.ptr[border_3]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xy_neg_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
	const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, 1, 0)];
	const block_view_t left_bottom = neighborhood[NEIGHBOR_IDX(-1, 1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = x;
	const int32_t border_1 = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
	const int32_t depth2 = cube_width - 1 - y;
	const int32_t border_2 = x + z * bottom.ldz;
	const int32_t border_3 = STENCIL_RADIUS - 1 + z * left_bottom.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - left_bottom.ptr[border_3] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 1] + left.ptr[border_1 - 3 * left.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 1] + left.ptr[border_1 - 4 * left.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy - 1]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 2]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 3] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 2]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 3 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1] - left_bottom.ptr[border_3] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 1] + left.ptr[border_1 - 3 * left.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 1] + left.ptr[border_1 - 4 * left.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy - 1]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 2]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 2]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1] - left_bottom.ptr[border_3] - center.ptr[base_idx - 3 * center.ldy + 1] + left.ptr[border_1 - 3 * left.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 1] + left.ptr[border_1 - 4 * left.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 2]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 2]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 1] + left.ptr[border_1 - 3 * left.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        bottom.ptr[border_2 + 1] - left_bottom.ptr[border_3] - center.ptr[base_idx - 4 * center.ldy + 1] + left.ptr[border_1 - 4 * left.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 3] +
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 2]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 3] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 3]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 1]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 2] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 3 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 1]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 1]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        bottom.ptr[border_2 + 2] - left_bottom.ptr[border_3] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 2] +
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 2] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 2]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 1] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 3 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 1] +
		        bottom.ptr[border_2 + 3] - left_bottom.ptr[border_3] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3 - 1] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 1]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 3] - bottom.ptr[border_2 + 3 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 3 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 3 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 2 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - left_bottom.ptr[border_3 + 1 * left_bottom.ldy] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy] +
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 4] - left_bottom.ptr[border_3] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xy_neg_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = x;
	const int32_t border_1 = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
    switch (depth1){
    case 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 1] + left.ptr[border_1 - 1 * left.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 1] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 1] + left.ptr[border_1 - 3 * left.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - left.ptr[border_1 + 4 * left.ldy] - center.ptr[base_idx - 4 * center.ldy + 1] + left.ptr[border_1 - 4 * left.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 3] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 3]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 2] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 2] + left.ptr[border_1 - 2 * left.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 2] + left.ptr[border_1 - 3 * left.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - left.ptr[border_1 + 4 * left.ldy] - center.ptr[base_idx - 4 * center.ldy + 2] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 2] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 2]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 3] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 3] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 3] + left.ptr[border_1 - 3 * left.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy - 1] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - left.ptr[border_1 + 4 * left.ldy] - center.ptr[base_idx - 4 * center.ldy + 3] + left.ptr[border_1 - 4 * left.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy - 1]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - left.ptr[border_1 + 1 * left.ldy] - center.ptr[base_idx - 1 * center.ldy + 4] + left.ptr[border_1 - 1 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - left.ptr[border_1 + 2 * left.ldy] - center.ptr[base_idx - 2 * center.ldy + 4] + left.ptr[border_1 - 2 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - left.ptr[border_1 + 3 * left.ldy] - center.ptr[base_idx - 3 * center.ldy + 4] + left.ptr[border_1 - 3 * left.ldy] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - left.ptr[border_1 + 4 * left.ldy] - center.ptr[base_idx - 4 * center.ldy + 4] + left.ptr[border_1 - 4 * left.ldy]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xy_pos_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t right = neighborhood[NEIGHBOR_IDX(1, 0, 0)];
	const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
	const block_view_t right_top = neighborhood[NEIGHBOR_IDX(1, -1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - x;
	const int32_t border_1 = y * right.ldy + z * right.ldz;
	const int32_t depth2 = y;
	const int32_t border_2 = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
	const int32_t border_3 = (STENCIL_RADIUS - 1) * right_top.ldy + z * right_top.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right_top.ptr[border_3] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 1] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right_top.ptr[border_3 + 3] + top.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 1] - right_top.ptr[border_3 - 3 * right_top.ldy] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 2] - right_top.ptr[border_3 - 2 * right_top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 - 3 * right_top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 - 2 * right_top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 3 * right_top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 3 * right_top.ldy + 3] + top.ptr[border_2 - 3 * top.ldy - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right_top.ptr[border_3] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 1] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 1] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 3] +
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 2] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 + 3] + top.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 - 2 * right_top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 2 * right_top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 1] - right_top.ptr[border_3] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 1] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 2] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 + 3] + top.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 1] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 1] - right_top.ptr[border_3] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 2] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 + 3] + top.ptr[border_2 - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right_top.ptr[border_3] + top.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 2] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 - 3 * right_top.ldy] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 - 2 * right_top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 3 * right_top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 3 * right_top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right_top.ptr[border_3] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 2] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 2 * right_top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 2] - right_top.ptr[border_3] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 2] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 2] - right_top.ptr[border_3] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 + 2] + top.ptr[border_2 - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right_top.ptr[border_3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 3 * right_top.ldy] + top.ptr[border_2 - 3 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 3 * right_top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right_top.ptr[border_3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 3] - right_top.ptr[border_3] + top.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 3] - right_top.ptr[border_3] + top.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 + 1] + top.ptr[border_2 - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right_top.ptr[border_3] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 3 * top.ldy + 3] + top.ptr[border_2 - 3 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 3 * right_top.ldy] + top.ptr[border_2 - 3 * top.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right_top.ptr[border_3] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 2 * right_top.ldy] + top.ptr[border_2 - 2 * top.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 4] - right_top.ptr[border_3] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3 - 1 * right_top.ldy] + top.ptr[border_2 - 1 * top.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 4] - right_top.ptr[border_3] + top.ptr[border_2 - 4]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xy_pos_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t right = neighborhood[NEIGHBOR_IDX(1, 0, 0)];
	const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, 1, 0)];
	const block_view_t right_bottom = neighborhood[NEIGHBOR_IDX(1, 1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - x;
	const int32_t border_1 = y * right.ldy + z * right.ldz;
	const int32_t depth2 = cube_width - 1 - y;
	const int32_t border_2 = x + z * bottom.ldz;
	const int32_t border_3 = z * right_bottom.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right_bottom.ptr[border_3 + 3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy + 3] - bottom.ptr[border_2 + 3 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 1] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 1] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 1] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 2] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 2] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 2] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 2] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy] - bottom.ptr[border_2 + 3 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 3] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 1] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 3] - bottom.ptr[border_2 + 3 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 3 * right_bottom.ldy] - bottom.ptr[border_2 + 3 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 2 * right_bottom.ldy] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3 + 1 * right_bottom.ldy] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right_bottom.ptr[border_3] - bottom.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xy_pos_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t right = neighborhood[NEIGHBOR_IDX(1, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - x;
	const int32_t border_1 = y * right.ldy + z * right.ldz;
    switch (depth1){
    case 0:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 1] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 1] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 1] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 1] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 2] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 2] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 2] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 2] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 2] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 2] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 3] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 3] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 3] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 3] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 3] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 4] - right.ptr[border_1 - 4 * right.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldy] - center.ptr[base_idx + 1 * center.ldy - 4] - right.ptr[border_1 - 1 * right.ldy] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldy] - center.ptr[base_idx + 2 * center.ldy - 4] - right.ptr[border_1 - 2 * right.ldy] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldy] - center.ptr[base_idx + 3 * center.ldy - 4] - right.ptr[border_1 - 3 * right.ldy] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldy] - center.ptr[base_idx + 4 * center.ldy - 4] - right.ptr[border_1 - 4 * right.ldy] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xy_center_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth2 = y;
	const int32_t border_2 = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
    switch (depth2){
    case 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - top.ptr[border_2 + 4] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 3 * top.ldy + 1] + top.ptr[border_2 - 3 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - top.ptr[border_2 - 1 * top.ldy + 4] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 3 * top.ldy + 2] + top.ptr[border_2 - 3 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - center.ptr[base_idx + 3 * center.ldy - 4] - top.ptr[border_2 - 2 * top.ldy + 4] + top.ptr[border_2 - 2 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 3 * top.ldy + 3] + top.ptr[border_2 - 3 * top.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - center.ptr[base_idx + 4 * center.ldy - 4] - top.ptr[border_2 - 3 * top.ldy + 4] + top.ptr[border_2 - 3 * top.ldy - 4]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 2 * top.ldy + 1] + top.ptr[border_2 - 2 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - top.ptr[border_2 + 4] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 2 * top.ldy + 2] + top.ptr[border_2 - 2 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - center.ptr[base_idx + 3 * center.ldy - 4] - top.ptr[border_2 - 1 * top.ldy + 4] + top.ptr[border_2 - 1 * top.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 2 * top.ldy + 3] + top.ptr[border_2 - 2 * top.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - center.ptr[base_idx + 4 * center.ldy - 4] - top.ptr[border_2 - 2 * top.ldy + 4] + top.ptr[border_2 - 2 * top.ldy - 4]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 - 1 * top.ldy + 1] + top.ptr[border_2 - 1 * top.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 - 1 * top.ldy + 2] + top.ptr[border_2 - 1 * top.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - center.ptr[base_idx + 3 * center.ldy - 4] - top.ptr[border_2 + 4] + top.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 - 1 * top.ldy + 3] + top.ptr[border_2 - 1 * top.ldy - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - center.ptr[base_idx + 4 * center.ldy - 4] - top.ptr[border_2 - 1 * top.ldy + 4] + top.ptr[border_2 - 1 * top.ldy - 4]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - top.ptr[border_2 + 1] + top.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - top.ptr[border_2 + 2] + top.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - center.ptr[base_idx + 3 * center.ldy - 4] - center.ptr[base_idx - 3 * center.ldy + 4] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - top.ptr[border_2 + 3] + top.ptr[border_2 - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldy + 4] - center.ptr[base_idx + 4 * center.ldy - 4] - top.ptr[border_2 + 4] + top.ptr[border_2 - 4]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xy_center_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, 1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth2 = cube_width - 1 - y;
	const int32_t border_2 = x + z * bottom.ldz;
    switch (depth2){
    case 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        bottom.ptr[border_2 + 4] - bottom.ptr[border_2 - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 1] - bottom.ptr[border_2 + 3 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 2] - bottom.ptr[border_2 + 3 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - center.ptr[base_idx - 3 * center.ldy + 4] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3 * bottom.ldy + 3] - bottom.ptr[border_2 + 3 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 3 * bottom.ldy + 4] - bottom.ptr[border_2 + 3 * bottom.ldy - 4] - center.ptr[base_idx - 4 * center.ldy + 4] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 1] - bottom.ptr[border_2 + 2 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        bottom.ptr[border_2 + 4] - bottom.ptr[border_2 - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 2] - bottom.ptr[border_2 + 2 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - center.ptr[base_idx - 3 * center.ldy + 4] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2 * bottom.ldy + 3] - bottom.ptr[border_2 + 2 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 2 * bottom.ldy + 4] - bottom.ptr[border_2 + 2 * bottom.ldy - 4] - center.ptr[base_idx - 4 * center.ldy + 4] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 1] - bottom.ptr[border_2 + 1 * bottom.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 2] - bottom.ptr[border_2 + 1 * bottom.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        bottom.ptr[border_2 + 4] - bottom.ptr[border_2 - 4] - center.ptr[base_idx - 3 * center.ldy + 4] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1 * bottom.ldy + 3] - bottom.ptr[border_2 + 1 * bottom.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 1 * bottom.ldy + 4] - bottom.ptr[border_2 + 1 * bottom.ldy - 4] - center.ptr[base_idx - 4 * center.ldy + 4] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
		        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
		        bottom.ptr[border_2 + 1] - bottom.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
		        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
		        bottom.ptr[border_2 + 2] - bottom.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldy + 4] - center.ptr[base_idx + 3 * center.ldy - 4] - center.ptr[base_idx - 3 * center.ldy + 4] + center.ptr[base_idx - 3 * center.ldy - 4] +
		        bottom.ptr[border_2 + 3] - bottom.ptr[border_2 - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
			) +
		    L44 * (
		        bottom.ptr[border_2 + 4] - bottom.ptr[border_2 - 4] - center.ptr[base_idx - 4 * center.ldy + 4] + center.ptr[base_idx - 4 * center.ldy - 4]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xy_center_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
    	
	
	return ((
	    L11 * (
	        center.ptr[base_idx + 1 * center.ldy + 1] - center.ptr[base_idx + 1 * center.ldy - 1] - center.ptr[base_idx - 1 * center.ldy + 1] + center.ptr[base_idx - 1 * center.ldy - 1]
	    ) +
	    L12 * (
	        center.ptr[base_idx + 1 * center.ldy + 2] - center.ptr[base_idx + 1 * center.ldy - 2] - center.ptr[base_idx - 1 * center.ldy + 2] + center.ptr[base_idx - 1 * center.ldy - 2] +
	        center.ptr[base_idx + 2 * center.ldy + 1] - center.ptr[base_idx + 2 * center.ldy - 1] - center.ptr[base_idx - 2 * center.ldy + 1] + center.ptr[base_idx - 2 * center.ldy - 1]
	    ) +
	    L13 * (
	        center.ptr[base_idx + 1 * center.ldy + 3] - center.ptr[base_idx + 1 * center.ldy - 3] - center.ptr[base_idx - 1 * center.ldy + 3] + center.ptr[base_idx - 1 * center.ldy - 3] +
	        center.ptr[base_idx + 3 * center.ldy + 1] - center.ptr[base_idx + 3 * center.ldy - 1] - center.ptr[base_idx - 3 * center.ldy + 1] + center.ptr[base_idx - 3 * center.ldy - 1]
	    ) +
	    L14 * (
	        center.ptr[base_idx + 1 * center.ldy + 4] - center.ptr[base_idx + 1 * center.ldy - 4] - center.ptr[base_idx - 1 * center.ldy + 4] + center.ptr[base_idx - 1 * center.ldy - 4] +
	        center.ptr[base_idx + 4 * center.ldy + 1] - center.ptr[base_idx + 4 * center.ldy - 1] - center.ptr[base_idx - 4 * center.ldy + 1] + center.ptr[base_idx - 4 * center.ldy - 1]
	    ) +
	    L22 * (
	        center.ptr[base_idx + 2 * center.ldy + 2] - center.ptr[base_idx + 2 * center.ldy - 2] - center.ptr[base_idx - 2 * center.ldy + 2] + center.ptr[base_idx - 2 * center.ldy - 2]
	    ) +
	    L23 * (
	        center.ptr[base_idx + 2 * center.ldy + 3] - center.ptr[base_idx + 2 * center.ldy - 3] - center.ptr[base_idx - 2 * center.ldy + 3] + center.ptr[base_idx - 2 * center.ldy - 3] +
	        center.ptr[base_idx + 3 * center.ldy + 2] - center.ptr[base_idx + 3 * center.ldy - 2] - center.ptr[base_idx - 3 * center.ldy + 2] + center.ptr[base_idx - 3 * center.ldy - 2]
		) +
	    L24 * (
	        center.ptr[base_idx + 2 * center.ldy + 4] - center.ptr[base_idx + 2 * center.ldy - 4] - center.ptr[base_idx - 2 * center.ldy + 4] + center.ptr[base_idx - 2 * center.ldy - 4] +
	        center.ptr[base_idx + 4 * center.ldy + 2] - center.ptr[base_idx + 4 * center.ldy - 2] - center.ptr[base_idx - 4 * center.ldy + 2] + center.ptr[base_idx - 4 * center.ldy - 2]
		) +
	    L33 * (
	        center.ptr[base_idx + 3 * center.ldy + 3] - center.ptr[base_idx + 3 * center.ldy - 3] - center.ptr[base_idx - 3 * center.ldy + 3] + center.ptr[base_idx - 3 * center.ldy - 3]
	
		) +
	    L34 * (
	        center.ptr[base_idx + 3 * center.ldy + 4] - center.ptr[base_idx + 3 * center.ldy - 4] - center.ptr[base_idx - 3 * center.ldy + 4] + center.ptr[base_idx - 3 * center.ldy - 4] +
	        center.ptr[base_idx + 4 * center.ldy + 3] - center.ptr[base_idx + 4 * center.ldy - 3] - center.ptr[base_idx - 4 * center.ldy + 3] + center.ptr[base_idx - 4 * center.ldy - 3]
		) +
	    L44 * (
	        center.ptr[base_idx + 4 * center.ldy + 4] - center.ptr[base_idx + 4 * center.ldy - 4] - center.ptr[base_idx - 4 * center.ldy + 4] + center.ptr[base_idx - 4 * center.ldy - 4]
	    )) * dinv);
}

ATTRIBUTE FP cross_deriv_xy(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    switch (stencil_side(x, cube_width) * 3 + stencil_side(y, cube_width))
    {
    case SIDE_NEG * 3 + SIDE_NEG: return cross_deriv_xy_neg_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_NEG * 3 + SIDE_POS: return cross_deriv_xy_neg_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_NEG * 3 + SIDE_CENTER: return cross_deriv_xy_neg_center(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_NEG: return cross_deriv_xy_pos_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_POS: return cross_deriv_xy_pos_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_CENTER: return cross_deriv_xy_pos_center(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_NEG: return cross_deriv_xy_center_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_POS: return cross_deriv_xy_center_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_CENTER:
    default: return cross_deriv_xy_center_center(neighborhood, x, y, z, cube_width, dinv);
    }
}

ATTRIBUTE static FP cross_deriv_yz_neg_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
	const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
	const block_view_t top_front = neighborhood[NEIGHBOR_IDX(0, -1, -1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = y;
	const int32_t border_1 = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
	const int32_t depth2 = z;
	const int32_t border_2 = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
	const int32_t border_3 = x + (STENCIL_RADIUS - 1) * top_front.ldy + (STENCIL_RADIUS - 1) * top_front.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - front.ptr[border_2 + 1 * front.ldy] + top_front.ptr[border_3]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 2 * top_front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 3 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 2 * top_front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 3 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 3 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 3 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz - 2 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 3 * top.ldy] - front.ptr[border_2 - 3 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz - 3 * top_front.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 + 1 * front.ldy] + top_front.ptr[border_3]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 2 * top_front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 3 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 3 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 2 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 3 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 3 * top_front.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 + 1 * front.ldy] + top_front.ptr[border_3]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 3 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 2 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 3 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 3 * top_front.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 + 1 * front.ldy] + top_front.ptr[border_3]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 3 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 2 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 1 * top_front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 2 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 3 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 3 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz - 2 * top_front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 2 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 2 * top_front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 2 * top_front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 + 2 * front.ldy] + top_front.ptr[border_3]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 1 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 3 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 3 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz - 1 * top_front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz - 1 * top_front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz - 1 * top_front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 + 3 * front.ldy] + top_front.ptr[border_3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 3 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 3 * top_front.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 2 * top_front.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + top_front.ptr[border_3 - 1 * top_front.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - front.ptr[border_2 + 4 * front.ldy] + top_front.ptr[border_3]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_yz_neg_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
	const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, 1)];
	const block_view_t top_back = neighborhood[NEIGHBOR_IDX(0, -1, 1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = y;
	const int32_t border_1 = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
	const int32_t depth2 = cube_width - 1 - z;
	const int32_t border_2 = x + y * back.ldy;
	const int32_t border_3 = x + (STENCIL_RADIUS - 1) * top_back.ldy;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 3 * top_back.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 3 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 3 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz - 3 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 3 * top.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 3 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 3 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 3 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 3 * top.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 3 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 3 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 3 * top.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 3 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 3 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 3 * top.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
			) +
		    L24 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L33 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 2 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy] +
		        back.ptr[border_2 + 3 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3 - 1 * top_back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz] +
		        back.ptr[border_2 + 3 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 3 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 2 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - top_back.ptr[border_3 + 1 * top_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz] +
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 4 * back.ldy] - top_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_yz_neg_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t top = neighborhood[NEIGHBOR_IDX(0, -1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = y;
	const int32_t border_1 = x + (STENCIL_RADIUS - 1) * top.ldy + z * top.ldz;
    switch (depth1){
    case 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 1 * top.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 3 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 3 * top.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 3 * top.ldy]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 2 * top.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 2 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 2 * top.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 2 * top.ldy]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 3 * top.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz - 1 * top.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz - 1 * top.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz - 1 * top.ldy]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 1 * top.ldz] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 1 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 2 * top.ldz] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 2 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 3 * top.ldz] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 3 * top.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - top.ptr[border_1 + 4 * top.ldz] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + top.ptr[border_1 - 4 * top.ldz]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_yz_pos_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, 1, 0)];
	const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
	const block_view_t bottom_front = neighborhood[NEIGHBOR_IDX(0, 1, -1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - y;
	const int32_t border_1 = x + z * bottom.ldz;
	const int32_t depth2 = z;
	const int32_t border_2 = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
	const int32_t border_3 = x + (STENCIL_RADIUS - 1) * bottom_front.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 3 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 3 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 3 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz + 3 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 3 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 3 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 3 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 3 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 3 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 3 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 2 * front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 2 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 2 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 3 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz + 1 * bottom_front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 + 1 * bottom_front.ldy] + front.ptr[border_2 - 4 * front.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 3 * bottom_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 2 * bottom_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3 - 1 * bottom_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom_front.ptr[border_3] + front.ptr[border_2 - 4 * front.ldy]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_yz_pos_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, 1, 0)];
	const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, 1)];
	const block_view_t bottom_back = neighborhood[NEIGHBOR_IDX(0, 1, 1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - y;
	const int32_t border_1 = x + z * bottom.ldz;
	const int32_t depth2 = cube_width - 1 - z;
	const int32_t border_2 = x + y * back.ldy;
	const int32_t border_3 = x;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 1 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 3 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 3 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 3 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz + 3 * bottom_back.ldy] - back.ptr[border_2 + 3 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 1 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 3 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 3 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 1 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 3 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 1 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 3 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 3 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 2 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 2 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 3 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz + 1 * bottom_back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 3 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldy] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 3 * bottom_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 2 * bottom_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3 + 1 * bottom_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom_back.ptr[border_3] - back.ptr[border_2 - 4 * back.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_yz_pos_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t bottom = neighborhood[NEIGHBOR_IDX(0, 1, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - y;
	const int32_t border_1 = x + z * bottom.ldz;
    switch (depth1){
    case 0:
		
		
		return ((
		    L11 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 3 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 3 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 2 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 2 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz + 1 * bottom.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz + 1 * bottom.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        bottom.ptr[border_1 + 1 * bottom.ldz] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 1 * bottom.ldz] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        bottom.ptr[border_1 + 2 * bottom.ldz] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 2 * bottom.ldz] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        bottom.ptr[border_1 + 3 * bottom.ldz] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 3 * bottom.ldz] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        bottom.ptr[border_1 + 4 * bottom.ldz] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - bottom.ptr[border_1 - 4 * bottom.ldz] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_yz_center_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth2 = z;
	const int32_t border_2 = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
    switch (depth2){
    case 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - front.ptr[border_2 + 4 * front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - front.ptr[border_2 - 3 * front.ldz + 4 * front.ldy] + front.ptr[border_2 - 3 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - front.ptr[border_2 + 4 * front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - front.ptr[border_2 - 2 * front.ldz + 4 * front.ldy] + front.ptr[border_2 - 2 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 2 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - front.ptr[border_2 + 4 * front.ldy] + front.ptr[border_2 - 4 * front.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 3 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - front.ptr[border_2 - 1 * front.ldz + 4 * front.ldy] + front.ptr[border_2 - 1 * front.ldz - 4 * front.ldy]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - front.ptr[border_2 + 1 * front.ldy] + front.ptr[border_2 - 1 * front.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - front.ptr[border_2 + 2 * front.ldy] + front.ptr[border_2 - 2 * front.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - front.ptr[border_2 + 3 * front.ldy] + front.ptr[border_2 - 3 * front.ldy]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - front.ptr[border_2 + 4 * front.ldy] + front.ptr[border_2 - 4 * front.ldy]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_yz_center_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, 1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth2 = cube_width - 1 - z;
	const int32_t border_2 = x + y * back.ldy;
    switch (depth2){
    case 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4 * back.ldy] - back.ptr[border_2 - 4 * back.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4 * back.ldy] - back.ptr[border_2 + 3 * back.ldz - 4 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        back.ptr[border_2 + 4 * back.ldy] - back.ptr[border_2 - 4 * back.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4 * back.ldy] - back.ptr[border_2 + 2 * back.ldz - 4 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 1 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 2 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4 * back.ldy] - back.ptr[border_2 - 4 * back.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldz + 3 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4 * back.ldy] - back.ptr[border_2 + 1 * back.ldz - 4 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
		        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 1 * back.ldy] - back.ptr[border_2 - 1 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
		        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 2 * back.ldy] - back.ptr[border_2 - 2 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
		        back.ptr[border_2 + 3 * back.ldy] - back.ptr[border_2 - 3 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
			) +
		    L44 * (
		        back.ptr[border_2 + 4 * back.ldy] - back.ptr[border_2 - 4 * back.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_yz_center_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
    	
	
	return ((
	    L11 * (
	        center.ptr[base_idx + 1 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 1 * center.ldy]
	    ) +
	    L12 * (
	        center.ptr[base_idx + 1 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 2 * center.ldy] +
	        center.ptr[base_idx + 2 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 1 * center.ldy]
	    ) +
	    L13 * (
	        center.ptr[base_idx + 1 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 3 * center.ldy] +
	        center.ptr[base_idx + 3 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 1 * center.ldy]
	    ) +
	    L14 * (
	        center.ptr[base_idx + 1 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 1 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 1 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 1 * center.ldz - 4 * center.ldy] +
	        center.ptr[base_idx + 4 * center.ldz + 1 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 1 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 1 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 1 * center.ldy]
	    ) +
	    L22 * (
	        center.ptr[base_idx + 2 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 2 * center.ldy]
	    ) +
	    L23 * (
	        center.ptr[base_idx + 2 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 3 * center.ldy] +
	        center.ptr[base_idx + 3 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 2 * center.ldy]
		) +
	    L24 * (
	        center.ptr[base_idx + 2 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 2 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 2 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 2 * center.ldz - 4 * center.ldy] +
	        center.ptr[base_idx + 4 * center.ldz + 2 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 2 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 2 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 2 * center.ldy]
		) +
	    L33 * (
	        center.ptr[base_idx + 3 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 3 * center.ldy]
	
		) +
	    L34 * (
	        center.ptr[base_idx + 3 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 3 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 3 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 3 * center.ldz - 4 * center.ldy] +
	        center.ptr[base_idx + 4 * center.ldz + 3 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 3 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 3 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 3 * center.ldy]
		) +
	    L44 * (
	        center.ptr[base_idx + 4 * center.ldz + 4 * center.ldy] - center.ptr[base_idx + 4 * center.ldz - 4 * center.ldy] - center.ptr[base_idx - 4 * center.ldz + 4 * center.ldy] + center.ptr[base_idx - 4 * center.ldz - 4 * center.ldy]
	    )) * dinv);
}

ATTRIBUTE FP cross_deriv_yz(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    switch (stencil_side(y, cube_width) * 3 + stencil_side(z, cube_width))
    {
    case SIDE_NEG * 3 + SIDE_NEG: return cross_deriv_yz_neg_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_NEG * 3 + SIDE_POS: return cross_deriv_yz_neg_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_NEG * 3 + SIDE_CENTER: return cross_deriv_yz_neg_center(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_NEG: return cross_deriv_yz_pos_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_POS: return cross_deriv_yz_pos_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_CENTER: return cross_deriv_yz_pos_center(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_NEG: return cross_deriv_yz_center_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_POS: return cross_deriv_yz_center_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_CENTER:
    default: return cross_deriv_yz_center_center(neighborhood, x, y, z, cube_width, dinv);
    }
}

ATTRIBUTE static FP cross_deriv_xz_neg_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
	const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
	const block_view_t left_front = neighborhood[NEIGHBOR_IDX(-1, 0, -1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = x;
	const int32_t border_1 = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
	const int32_t depth2 = z;
	const int32_t border_2 = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
	const int32_t border_3 = STENCIL_RADIUS - 1 + y * left_front.ldy + (STENCIL_RADIUS - 1) * left_front.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - front.ptr[border_2 + 1] + left_front.ptr[border_3]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - front.ptr[border_2 + 2] + left_front.ptr[border_3 - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 1] + left_front.ptr[border_3 - 1 * left_front.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 1] + left_front.ptr[border_3 - 2 * left_front.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 3 * front.ldz + 1] + left_front.ptr[border_3 - 3 * left_front.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 2] + left_front.ptr[border_3 - 1 * left_front.ldz - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 2] + left_front.ptr[border_3 - 2 * left_front.ldz - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 2] + left_front.ptr[border_3 - 3 * left_front.ldz - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 3] + left_front.ptr[border_3 - 2 * left_front.ldz - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 3] + left_front.ptr[border_3 - 3 * left_front.ldz - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 3] - front.ptr[border_2 - 3 * front.ldz + 4] + left_front.ptr[border_3 - 3 * left_front.ldz - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 + 1] + left_front.ptr[border_3]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 1] + left_front.ptr[border_3 - 1 * left_front.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 1] + left_front.ptr[border_3 - 2 * left_front.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - front.ptr[border_2 + 2] + left_front.ptr[border_3 - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 2] + left_front.ptr[border_3 - 1 * left_front.ldz - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 2] + left_front.ptr[border_3 - 2 * left_front.ldz - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 3] + left_front.ptr[border_3 - 2 * left_front.ldz - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 + 1] + left_front.ptr[border_3]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 1] + left_front.ptr[border_3 - 1 * left_front.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 + 2] + left_front.ptr[border_3 - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 2] + left_front.ptr[border_3 - 1 * left_front.ldz - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 2] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 3] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 1] + left.ptr[border_1 - 3 * left.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 + 1] + left_front.ptr[border_3]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 + 2] + left_front.ptr[border_3 - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 3] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 3]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - front.ptr[border_2 + 2] + left_front.ptr[border_3] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 2] + left_front.ptr[border_3 - 1 * left_front.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 2] + left_front.ptr[border_3 - 2 * left_front.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 3 * front.ldz + 2] + left_front.ptr[border_3 - 3 * left_front.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 3] + left_front.ptr[border_3 - 2 * left_front.ldz - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 3] + left_front.ptr[border_3 - 3 * left_front.ldz - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 4] + left_front.ptr[border_3 - 3 * left_front.ldz - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 + 2] + left_front.ptr[border_3]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 2] + left_front.ptr[border_3 - 1 * left_front.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 2] + left_front.ptr[border_3 - 2 * left_front.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 3] + left_front.ptr[border_3 - 2 * left_front.ldz - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 + 2] + left_front.ptr[border_3]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 2] + left_front.ptr[border_3 - 1 * left_front.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 2] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 + 2] + left_front.ptr[border_3]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 + 3] + left_front.ptr[border_3 - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 2] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 2]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - front.ptr[border_2 + 3] + left_front.ptr[border_3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 3] + left_front.ptr[border_3 - 2 * left_front.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 3 * front.ldz + 3] + left_front.ptr[border_3 - 3 * left_front.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 4] + left_front.ptr[border_3 - 3 * left_front.ldz - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 + 3] + left_front.ptr[border_3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 3] + left_front.ptr[border_3 - 2 * left_front.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 + 3] + left_front.ptr[border_3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 1] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 3] + left_front.ptr[border_3 - 1 * left_front.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 + 3] + left_front.ptr[border_3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 1] - front.ptr[border_2 + 4] + left_front.ptr[border_3 - 1]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - front.ptr[border_2 + 4] + left_front.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 3 * front.ldz + 3] + front.ptr[border_2 - 3 * front.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 3 * front.ldz + 4] + left_front.ptr[border_3 - 3 * left_front.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - front.ptr[border_2 + 4] + left_front.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 2 * front.ldz + 4] + left_front.ptr[border_3 - 2 * left_front.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz] - front.ptr[border_2 + 4] + left_front.ptr[border_3] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 - 1 * front.ldz + 4] + left_front.ptr[border_3 - 1 * left_front.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz] - front.ptr[border_2 + 4] + left_front.ptr[border_3]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xz_neg_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
	const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, 1)];
	const block_view_t left_back = neighborhood[NEIGHBOR_IDX(-1, 0, 1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = x;
	const int32_t border_1 = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
	const int32_t depth2 = cube_width - 1 - z;
	const int32_t border_2 = x + y * back.ldy;
	const int32_t border_3 = STENCIL_RADIUS - 1 + y * left_back.ldy;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - left_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 1] + left.ptr[border_1 - 3 * left.ldz]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - left_back.ptr[border_3 + 3 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 1] + left.ptr[border_1 - 4 * left.ldz]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - left_back.ptr[border_3 + 2 * left_back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz - 1]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - left_back.ptr[border_3 + 3 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3] - left_back.ptr[border_3 + 2 * left_back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 2]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 3] +
		        back.ptr[border_2 + 3 * back.ldz + 3] - left_back.ptr[border_3 + 3 * left_back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 2]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4] - left_back.ptr[border_3 + 3 * left_back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 1] - left_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 1] + left.ptr[border_1 - 3 * left.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 1] + left.ptr[border_1 - 4 * left.ldz]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz - 1]
			) +
		    L24 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - left_back.ptr[border_3 + 2 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 2]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 3] - left_back.ptr[border_3 + 2 * left_back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 2]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 1] - left_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 1] + left.ptr[border_1 - 3 * left.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 1] + left.ptr[border_1 - 4 * left.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 2] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L33 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 2]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 2]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 3]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 1] + left.ptr[border_1 - 3 * left.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        back.ptr[border_2 + 1] - left_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 1] + left.ptr[border_1 - 4 * left.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        back.ptr[border_2 + 2] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 3] +
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 2]
			) +
		    L44 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 3] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 3]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - left_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - left_back.ptr[border_3 + 3 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3] - left_back.ptr[border_3 + 2 * left_back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 1]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 2] +
		        back.ptr[border_2 + 3 * back.ldz + 3] - left_back.ptr[border_3 + 3 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4] - left_back.ptr[border_3 + 3 * left_back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - left_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz]
			) +
		    L24 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 1]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 2] +
		        back.ptr[border_2 + 2 * back.ldz + 3] - left_back.ptr[border_3 + 2 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 2] - left_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L33 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 1]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 2]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        back.ptr[border_2 + 2] - left_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 2] +
		        back.ptr[border_2 + 3] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L44 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 2] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 2]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 1] +
		        back.ptr[border_2 + 3 * back.ldz + 3] - left_back.ptr[border_3 + 3 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4] - left_back.ptr[border_3 + 3 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 1] +
		        back.ptr[border_2 + 2 * back.ldz + 3] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 3] - left_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 1] +
		        back.ptr[border_2 + 1 * back.ldz + 3] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 1]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 1] +
		        back.ptr[border_2 + 3] - left_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L44 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3 - 1] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 1]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz] +
		        back.ptr[border_2 + 3 * back.ldz + 3] - back.ptr[border_2 + 3 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4] - left_back.ptr[border_3 + 3 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz] +
		        back.ptr[border_2 + 2 * back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - left_back.ptr[border_3 + 2 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz] +
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - left_back.ptr[border_3 + 1 * left_back.ldz] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz] +
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 4] - left_back.ptr[border_3] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xz_neg_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t left = neighborhood[NEIGHBOR_IDX(-1, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = x;
	const int32_t border_1 = STENCIL_RADIUS - 1 + y * left.ldy + z * left.ldz;
    switch (depth1){
    case 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 1] + left.ptr[border_1 - 1 * left.ldz]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 1] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 1] + left.ptr[border_1 - 3 * left.ldz]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - left.ptr[border_1 + 4 * left.ldz] - center.ptr[base_idx - 4 * center.ldz + 1] + left.ptr[border_1 - 4 * left.ldz]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz - 1]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz - 1]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 2]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 3] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 2]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 3]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 2] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 2] + left.ptr[border_1 - 2 * left.ldz]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 2] + left.ptr[border_1 - 3 * left.ldz]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - left.ptr[border_1 + 4 * left.ldz] - center.ptr[base_idx - 4 * center.ldz + 2] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz - 1]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 2] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz - 1]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 2]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 3] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 3] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 3] + left.ptr[border_1 - 3 * left.ldz]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz - 1] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - left.ptr[border_1 + 4 * left.ldz] - center.ptr[base_idx - 4 * center.ldz + 3] + left.ptr[border_1 - 4 * left.ldz]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz - 1]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - left.ptr[border_1 + 1 * left.ldz] - center.ptr[base_idx - 1 * center.ldz + 4] + left.ptr[border_1 - 1 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - left.ptr[border_1 + 2 * left.ldz] - center.ptr[base_idx - 2 * center.ldz + 4] + left.ptr[border_1 - 2 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - left.ptr[border_1 + 3 * left.ldz] - center.ptr[base_idx - 3 * center.ldz + 4] + left.ptr[border_1 - 3 * left.ldz] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - left.ptr[border_1 + 4 * left.ldz] - center.ptr[base_idx - 4 * center.ldz + 4] + left.ptr[border_1 - 4 * left.ldz]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xz_pos_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t right = neighborhood[NEIGHBOR_IDX(1, 0, 0)];
	const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
	const block_view_t right_front = neighborhood[NEIGHBOR_IDX(1, 0, -1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - x;
	const int32_t border_1 = y * right.ldy + z * right.ldz;
	const int32_t depth2 = z;
	const int32_t border_2 = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
	const int32_t border_3 = y * right_front.ldy + (STENCIL_RADIUS - 1) * right_front.ldz;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right_front.ptr[border_3] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 1] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right_front.ptr[border_3 + 3] + front.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 1] - right_front.ptr[border_3 - 3 * right_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 2] - right_front.ptr[border_3 - 2 * right_front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 - 3 * right_front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 - 2 * right_front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 3 * right_front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 3 * right_front.ldz + 3] + front.ptr[border_2 - 3 * front.ldz - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right_front.ptr[border_3] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 1] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 1] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 3] +
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 2] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 + 3] + front.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 - 2 * right_front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 2 * right_front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 1] - right_front.ptr[border_3] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 1] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 2] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 + 3] + front.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 1] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 1] - right_front.ptr[border_3] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 2] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 + 3] + front.ptr[border_2 - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right_front.ptr[border_3] + front.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 2] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 - 3 * right_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 - 2 * right_front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 3 * right_front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 3 * right_front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right_front.ptr[border_3] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 2] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 2 * right_front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 2] - right_front.ptr[border_3] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 2] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 2] - right_front.ptr[border_3] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 + 2] + front.ptr[border_2 - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right_front.ptr[border_3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 3 * right_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 3 * right_front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right_front.ptr[border_3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 3] - right_front.ptr[border_3] + front.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 3] - right_front.ptr[border_3] + front.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 + 1] + front.ptr[border_2 - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right_front.ptr[border_3] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 3 * front.ldz + 3] + front.ptr[border_2 - 3 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 3 * right_front.ldz] + front.ptr[border_2 - 3 * front.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right_front.ptr[border_3] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 2 * right_front.ldz] + front.ptr[border_2 - 2 * front.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 4] - right_front.ptr[border_3] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3 - 1 * right_front.ldz] + front.ptr[border_2 - 1 * front.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 4] - right_front.ptr[border_3] + front.ptr[border_2 - 4]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xz_pos_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t right = neighborhood[NEIGHBOR_IDX(1, 0, 0)];
	const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, 1)];
	const block_view_t right_back = neighborhood[NEIGHBOR_IDX(1, 0, 1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - x;
	const int32_t border_1 = y * right.ldy + z * right.ldz;
	const int32_t depth2 = cube_width - 1 - z;
	const int32_t border_2 = x + y * back.ldy;
	const int32_t border_3 = y * right_back.ldy;
    switch (depth1 * STENCIL_RADIUS + depth2){
    case 0 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 1] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right_back.ptr[border_3 + 3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right_back.ptr[border_3 + 3 * right_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 1] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right_back.ptr[border_3 + 2 * right_back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 2] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 3 * right_back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 2] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 3 * right_back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 3 * right_back.ldz + 3] - back.ptr[border_2 + 3 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right_back.ptr[border_3] - back.ptr[border_2 - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 1] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 1] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 2] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2 * right_back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 2] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2 * right_back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right_back.ptr[border_3] - back.ptr[border_2 - 1] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 1] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 2] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 2] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 3] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1 * right_back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 0 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 1] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right_back.ptr[border_3] - back.ptr[border_2 - 1] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 2] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 2] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 3] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 2] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 3 * right_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 2] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 3 * right_back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 3 * right_back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 2] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 2] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2 * right_back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right_back.ptr[border_3] - back.ptr[border_2 - 2] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 2] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 3] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 1 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 2] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right_back.ptr[border_3] - back.ptr[border_2 - 2] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 3] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 2] - back.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 3] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 3 * right_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 3 * right_back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 3] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 3] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 2 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right_back.ptr[border_3] - back.ptr[border_2 - 3] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 1] - back.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 3] - back.ptr[border_2 + 3 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 3 * right_back.ldz] - back.ptr[border_2 + 3 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 2 * right_back.ldz] - back.ptr[border_2 + 2 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3 + 1 * right_back.ldz] - back.ptr[border_2 + 1 * back.ldz - 4] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 3 * STENCIL_RADIUS + 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right_back.ptr[border_3] - back.ptr[border_2 - 4] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
	default: UNREACHABLE; 
	}
}

ATTRIBUTE static FP cross_deriv_xz_pos_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t right = neighborhood[NEIGHBOR_IDX(1, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth1 = cube_width - 1 - x;
	const int32_t border_1 = y * right.ldy + z * right.ldz;
    switch (depth1){
    case 0:
		
		
		return ((
		    L11 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 1] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 1] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 1] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 1] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 2] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 2] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 2] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 2] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 2] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 2] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 3] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 3] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 3] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 3] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 3] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 4] - right.ptr[border_1 - 4 * right.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        right.ptr[border_1 + 1 * right.ldz] - center.ptr[base_idx + 1 * center.ldz - 4] - right.ptr[border_1 - 1 * right.ldz] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        right.ptr[border_1 + 2 * right.ldz] - center.ptr[base_idx + 2 * center.ldz - 4] - right.ptr[border_1 - 2 * right.ldz] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        right.ptr[border_1 + 3 * right.ldz] - center.ptr[base_idx + 3 * center.ldz - 4] - right.ptr[border_1 - 3 * right.ldz] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        right.ptr[border_1 + 4 * right.ldz] - center.ptr[base_idx + 4 * center.ldz - 4] - right.ptr[border_1 - 4 * right.ldz] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xz_center_neg(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t front = neighborhood[NEIGHBOR_IDX(0, 0, -1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth2 = z;
	const int32_t border_2 = x + y * front.ldy + (STENCIL_RADIUS - 1) * front.ldz;
    switch (depth2){
    case 0:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - front.ptr[border_2 + 4] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 3 * front.ldz + 1] + front.ptr[border_2 - 3 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - front.ptr[border_2 - 1 * front.ldz + 4] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 3 * front.ldz + 2] + front.ptr[border_2 - 3 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - center.ptr[base_idx + 3 * center.ldz - 4] - front.ptr[border_2 - 2 * front.ldz + 4] + front.ptr[border_2 - 2 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 3 * front.ldz + 3] + front.ptr[border_2 - 3 * front.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - center.ptr[base_idx + 4 * center.ldz - 4] - front.ptr[border_2 - 3 * front.ldz + 4] + front.ptr[border_2 - 3 * front.ldz - 4]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 2 * front.ldz + 1] + front.ptr[border_2 - 2 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - front.ptr[border_2 + 4] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 2 * front.ldz + 2] + front.ptr[border_2 - 2 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - center.ptr[base_idx + 3 * center.ldz - 4] - front.ptr[border_2 - 1 * front.ldz + 4] + front.ptr[border_2 - 1 * front.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 2 * front.ldz + 3] + front.ptr[border_2 - 2 * front.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - center.ptr[base_idx + 4 * center.ldz - 4] - front.ptr[border_2 - 2 * front.ldz + 4] + front.ptr[border_2 - 2 * front.ldz - 4]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 - 1 * front.ldz + 1] + front.ptr[border_2 - 1 * front.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 - 1 * front.ldz + 2] + front.ptr[border_2 - 1 * front.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - center.ptr[base_idx + 3 * center.ldz - 4] - front.ptr[border_2 + 4] + front.ptr[border_2 - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 - 1 * front.ldz + 3] + front.ptr[border_2 - 1 * front.ldz - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - center.ptr[base_idx + 4 * center.ldz - 4] - front.ptr[border_2 - 1 * front.ldz + 4] + front.ptr[border_2 - 1 * front.ldz - 4]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - front.ptr[border_2 + 1] + front.ptr[border_2 - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - front.ptr[border_2 + 2] + front.ptr[border_2 - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - center.ptr[base_idx + 3 * center.ldz - 4] - center.ptr[base_idx - 3 * center.ldz + 4] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - front.ptr[border_2 + 3] + front.ptr[border_2 - 3]
			) +
		    L44 * (
		        center.ptr[base_idx + 4 * center.ldz + 4] - center.ptr[base_idx + 4 * center.ldz - 4] - front.ptr[border_2 + 4] + front.ptr[border_2 - 4]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xz_center_pos(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const block_view_t back = neighborhood[NEIGHBOR_IDX(0, 0, 1)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
	const int32_t depth2 = cube_width - 1 - z;
	const int32_t border_2 = x + y * back.ldy;
    switch (depth2){
    case 0:
		
		
		return ((
		    L11 * (
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        back.ptr[border_2 + 4] - back.ptr[border_2 - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 1] - back.ptr[border_2 + 3 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - back.ptr[border_2 + 1 * back.ldz - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 2] - back.ptr[border_2 + 3 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 2 * back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - back.ptr[border_2 + 2 * back.ldz - 4] - center.ptr[base_idx - 3 * center.ldz + 4] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 3 * back.ldz + 3] - back.ptr[border_2 + 3 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 3 * back.ldz + 4] - back.ptr[border_2 + 3 * back.ldz - 4] - center.ptr[base_idx - 4 * center.ldz + 4] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 1:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 1] - back.ptr[border_2 + 2 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        back.ptr[border_2 + 4] - back.ptr[border_2 - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 2] - back.ptr[border_2 + 2 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - back.ptr[border_2 + 1 * back.ldz - 4] - center.ptr[base_idx - 3 * center.ldz + 4] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 2 * back.ldz + 3] - back.ptr[border_2 + 2 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 2 * back.ldz + 4] - back.ptr[border_2 + 2 * back.ldz - 4] - center.ptr[base_idx - 4 * center.ldz + 4] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 2:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 1] - back.ptr[border_2 + 1 * back.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 2] - back.ptr[border_2 + 1 * back.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        back.ptr[border_2 + 4] - back.ptr[border_2 - 4] - center.ptr[base_idx - 3 * center.ldz + 4] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 1 * back.ldz + 3] - back.ptr[border_2 + 1 * back.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 1 * back.ldz + 4] - back.ptr[border_2 + 1 * back.ldz - 4] - center.ptr[base_idx - 4 * center.ldz + 4] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
    case 3:
		
		
		return ((
		    L11 * (
		        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
		    ) +
		    L12 * (
		        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
		        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
		    ) +
		    L13 * (
		        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
		    ) +
		    L14 * (
		        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
		        back.ptr[border_2 + 1] - back.ptr[border_2 - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
		    ) +
		    L22 * (
		        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
		    ) +
		    L23 * (
		        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
		        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
			) +
		    L24 * (
		        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
		        back.ptr[border_2 + 2] - back.ptr[border_2 - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
			) +
		    L33 * (
		        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
		
			) +
		    L34 * (
		        center.ptr[base_idx + 3 * center.ldz + 4] - center.ptr[base_idx + 3 * center.ldz - 4] - center.ptr[base_idx - 3 * center.ldz + 4] + center.ptr[base_idx - 3 * center.ldz - 4] +
		        back.ptr[border_2 + 3] - back.ptr[border_2 - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
			) +
		    L44 * (
		        back.ptr[border_2 + 4] - back.ptr[border_2 - 4] - center.ptr[base_idx - 4 * center.ldz + 4] + center.ptr[base_idx - 4 * center.ldz - 4]
		    )) * dinv);
	default: UNREACHABLE;
	}
}

ATTRIBUTE static FP cross_deriv_xz_center_center(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    (void) cube_width;
	const block_view_t center = neighborhood[NEIGHBOR_IDX(0, 0, 0)];
	const int32_t base_idx = x + y * center.ldy + z * center.ldz;
    	
	
	return ((
	    L11 * (
	        center.ptr[base_idx + 1 * center.ldz + 1] - center.ptr[base_idx + 1 * center.ldz - 1] - center.ptr[base_idx - 1 * center.ldz + 1] + center.ptr[base_idx - 1 * center.ldz - 1]
	    ) +
	    L12 * (
	        center.ptr[base_idx + 1 * center.ldz + 2] - center.ptr[base_idx + 1 * center.ldz - 2] - center.ptr[base_idx - 1 * center.ldz + 2] + center.ptr[base_idx - 1 * center.ldz - 2] +
	        center.ptr[base_idx + 2 * center.ldz + 1] - center.ptr[base_idx + 2 * center.ldz - 1] - center.ptr[base_idx - 2 * center.ldz + 1] + center.ptr[base_idx - 2 * center.ldz - 1]
	    ) +
	    L13 * (
	        center.ptr[base_idx + 1 * center.ldz + 3] - center.ptr[base_idx + 1 * center.ldz - 3] - center.ptr[base_idx - 1 * center.ldz + 3] + center.ptr[base_idx - 1 * center.ldz - 3] +
	        center.ptr[base_idx + 3 * center.ldz + 1] - center.ptr[base_idx + 3 * center.ldz - 1] - center.ptr[base_idx - 3 * center.ldz + 1] + center.ptr[base_idx - 3 * center.ldz - 1]
	    ) +
	    L14 * (
	        center.ptr[base_idx + 1 * center.ldz + 4] - center.ptr[base_idx + 1 * center.ldz - 4] - center.ptr[base_idx - 1 * center.ldz + 4] + center.ptr[base_idx - 1 * center.ldz - 4] +
	        center.ptr[base_idx + 4 * center.ldz + 1] - center.ptr[base_idx + 4 * center.ldz - 1] - center.ptr[base_idx - 4 * center.ldz + 1] + center.ptr[base_idx - 4 * center.ldz - 1]
	    ) +
	    L22 * (
	        center.ptr[base_idx + 2 * center.ldz + 2] - center.ptr[base_idx + 2 * center.ldz - 2] - center.ptr[base_idx - 2 * center.ldz + 2] + center.ptr[base_idx - 2 * center.ldz - 2]
	    ) +
	    L23 * (
	        center.ptr[base_idx + 2 * center.ldz + 3] - center.ptr[base_idx + 2 * center.ldz - 3] - center.ptr[base_idx - 2 * center.ldz + 3] + center.ptr[base_idx - 2 * center.ldz - 3] +
	        center.ptr[base_idx + 3 * center.ldz + 2] - center.ptr[base_idx + 3 * center.ldz - 2] - center.ptr[base_idx - 3 * center.ldz + 2] + center.ptr[base_idx - 3 * center.ldz - 2]
		) +
	    L24 * (
	        center.ptr[base_idx + 2 * center.ldz + 4] - center.ptr[base_idx + 2 * center.ldz - 4] - center.ptr[base_idx - 2 * center.ldz + 4] + center.ptr[base_idx - 2 * center.ldz - 4] +
	        center.ptr[base_idx + 4 * center.ldz + 2] - center.ptr[base_idx + 4 * center.ldz - 2] - center.ptr[base_idx - 4 * center.ldz + 2] + center.ptr[base_idx - 4 * center.ldz - 2]
		) +
	    L33 * (
	        center.ptr[base_idx + 3 * center.ldz + 3] - center.ptr[base_idx + 3 * center.ldz - 3] - center.ptr[base_idx - 3 * center.ldz + 3] + center.ptr[base_idx - 3 * center.ldz - 3]
	
		) +
	    L34 * (
	        center.ptr[base_idx + 3 * center.ldz + 4] - center.ptr[base_idx + 3 * center.ldz - 4] - center.ptr[base_idx - 3 * center.ldz + 4] + center.ptr[base_idx - 3 * center.ldz - 4] +
	        center.ptr[base_idx + 4 * center.ldz + 3] - center.ptr[base_idx + 4 * center.ldz - 3] - center.ptr[base_idx - 4 * center.ldz + 3] + center.ptr[base_idx - 4 * center.ldz - 3]
		) +
	    L44 * (
	        center.ptr[base_idx + 4 * center.ldz + 4] - center.ptr[base_idx + 4 * center.ldz - 4] - center.ptr[base_idx - 4 * center.ldz + 4] + center.ptr[base_idx - 4 * center.ldz - 4]
	    )) * dinv);
}

ATTRIBUTE FP cross_deriv_xz(
    const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv
){
    switch (stencil_side(x, cube_width) * 3 + stencil_side(z, cube_width))
    {
    case SIDE_NEG * 3 + SIDE_NEG: return cross_deriv_xz_neg_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_NEG * 3 + SIDE_POS: return cross_deriv_xz_neg_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_NEG * 3 + SIDE_CENTER: return cross_deriv_xz_neg_center(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_NEG: return cross_deriv_xz_pos_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_POS: return cross_deriv_xz_pos_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_POS * 3 + SIDE_CENTER: return cross_deriv_xz_pos_center(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_NEG: return cross_deriv_xz_center_neg(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_POS: return cross_deriv_xz_center_pos(neighborhood, x, y, z, cube_width, dinv);
    case SIDE_CENTER * 3 + SIDE_CENTER:
    default: return cross_deriv_xz_center_center(neighborhood, x, y, z, cube_width, dinv);
    }
}
