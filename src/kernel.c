#include <starpu.h>

#include "kernel.h"
#include "floatingpoint.h"
#include "macros.h"
#include "derivatives.h"
#include "wave-neighborhood.h"
#include <errno.h>


err_t make_rtm_args(rtm_args_t** rtm_args,
    const size_t x_start, const size_t x_end,
    const size_t y_start, const size_t y_end,
    const size_t z_start, const size_t z_end,
    const FP dx, const FP dy, const FP dz, const FP dt
){
    *rtm_args = (rtm_args_t*) malloc(sizeof(rtm_args_t));
    if(*rtm_args == NULL){
        return errno;
    }
    **rtm_args = (rtm_args_t){
        .x_start = x_start, .x_end = x_end,
        .y_start = y_start, .y_end = y_end,
        .z_start = z_start, .z_end = z_end,
        .dx = dx, .dy = dy, .dz = dz, .dt = dt
    };
    return 0;
}


err_t make_perturb_args(perturb_args_t** perturb_args, const size_t idx, const FP value, const size_t t){
    *perturb_args = (perturb_args_t*) malloc(sizeof(perturb_args_t));
    if(*perturb_args == NULL){
        return errno;
    }
    (*perturb_args)->source_idx = idx;
    (*perturb_args)->perturb_value = value;
    (*perturb_args)->t = t;
    return 0;
}

void rtm_kernel(void *descr[], void *cl_args){
    const rtm_args_t* args = (rtm_args_t*) cl_args;

    const int32_t x_start = (int32_t) args->x_start;
    const int32_t y_start = (int32_t) args->y_start;
    const int32_t z_start = (int32_t) args->z_start;
    const int32_t x_end = (int32_t) args->x_end;
    const int32_t y_end = (int32_t) args->y_end;
    const int32_t z_end = (int32_t) args->z_end;

    const FP dx = args->dx;
    const FP dy = args->dy;
    const FP dz = args->dz;
    const FP dt = args->dt;

    const FP dxxinv = FP_LIT(1.0) / (dx * dx);
    const FP dyyinv = FP_LIT(1.0) / (dy * dy);
    const FP dzzinv = FP_LIT(1.0) / (dz * dz);
    const FP dxyinv = FP_LIT(1.0) / (dx * dy);
    const FP dxzinv = FP_LIT(1.0) / (dx * dz);
    const FP dyzinv = FP_LIT(1.0) / (dy * dz);

    const int32_t cube_width_x = (int32_t) STARPU_BLOCK_GET_NX(descr[0]);
    const int32_t cube_width_y = (int32_t) STARPU_BLOCK_GET_NY(descr[0]);
    const int32_t cube_width_z = (int32_t) STARPU_BLOCK_GET_NZ(descr[0]);

    // precomputed values
    const FP* ch1dxx = (FP*) STARPU_BLOCK_GET_PTR(descr[0]);
    const FP* ch1dyy = (FP*) STARPU_BLOCK_GET_PTR(descr[1]);
    const FP* ch1dzz = (FP*) STARPU_BLOCK_GET_PTR(descr[2]);
    const FP* ch1dxy = (FP*) STARPU_BLOCK_GET_PTR(descr[3]);
    const FP* ch1dyz = (FP*) STARPU_BLOCK_GET_PTR(descr[4]);
    const FP* ch1dxz = (FP*) STARPU_BLOCK_GET_PTR(descr[5]);
    const FP* v2px = (FP*) STARPU_BLOCK_GET_PTR(descr[6]);
    const FP* v2pz = (FP*) STARPU_BLOCK_GET_PTR(descr[7]);
    const FP* v2sz = (FP*) STARPU_BLOCK_GET_PTR(descr[8]);
    const FP* v2pn = (FP*) STARPU_BLOCK_GET_PTR(descr[9]);

    // primary wave: t (written), t-1 (central cube + neighbor faces/edges), t-2 (central cube)
    FP *const pwwrite = (FP*) STARPU_BLOCK_GET_PTR(descr[10]);
    block_view_t p_wave_neighborhood[NEIGHBORHOOD_SIZE];
    fill_wave_neighborhood(p_wave_neighborhood, descr, 11);
    const FP* pwcentralt1 = p_wave_neighborhood[NEIGHBOR_IDX(0, 0, 0)].ptr;
    const FP* pwcentralt2 = (FP*) STARPU_BLOCK_GET_PTR(descr[30]);

    // secondary wave, same layout
    FP *const qwwrite = (FP*) STARPU_BLOCK_GET_PTR(descr[31]);
    block_view_t q_wave_neighborhood[NEIGHBORHOOD_SIZE];
    fill_wave_neighborhood(q_wave_neighborhood, descr, 32);
    const FP* qwcentralt1 = q_wave_neighborhood[NEIGHBOR_IDX(0, 0, 0)].ptr;
    const FP* qwcentralt2 = (FP*) STARPU_BLOCK_GET_PTR(descr[51]);

    for(int32_t z = 0; z < cube_width_z; z++)
    for(int32_t y = 0; y < cube_width_y; y++)
    for(int32_t x = 0; x < cube_width_x; x++){
        const size_t idx = cube_idx(x, y, z);

        if(
            (z < z_start || z >= z_end) ||
            (y < y_start || y >= y_end) ||
            (x < x_start || x >= x_end)
        ){
            pwwrite[idx] = FP_LIT(0.0);
            qwwrite[idx] = FP_LIT(0.0);
            continue;
        }

        // p derivatives, H1(p) and H2(p)
        const FP pxx = snd_deriv_x(p_wave_neighborhood, x, y, z, cube_width_x, dxxinv);
        const FP pyy = snd_deriv_y(p_wave_neighborhood, x, y, z, cube_width_y, dyyinv);
        const FP pzz = snd_deriv_z(p_wave_neighborhood, x, y, z, cube_width_z, dzzinv);
        const FP pxy = cross_deriv_xy(p_wave_neighborhood, x, y, z, cube_width_x, dxyinv);
        const FP pyz = cross_deriv_yz(p_wave_neighborhood, x, y, z, cube_width_y, dyzinv);
        const FP pxz = cross_deriv_xz(p_wave_neighborhood, x, y, z, cube_width_x, dxzinv);

        const FP cpxx = ch1dxx[idx] * pxx;
        const FP cpyy = ch1dyy[idx] * pyy;
        const FP cpzz = ch1dzz[idx] * pzz;
        const FP cpxy = ch1dxy[idx] * pxy;
        const FP cpxz = ch1dxz[idx] * pxz;
        const FP cpyz = ch1dyz[idx] * pyz;
        const FP h1p = cpxx + cpyy + cpzz + cpxy + cpxz + cpyz;
        const FP h2p = pxx + pyy + pzz - h1p;

        // q derivatives, H1(q) and H2(q)
        const FP qxx = snd_deriv_x(q_wave_neighborhood, x, y, z, cube_width_x, dxxinv);
        const FP qyy = snd_deriv_y(q_wave_neighborhood, x, y, z, cube_width_y, dyyinv);
        const FP qzz = snd_deriv_z(q_wave_neighborhood, x, y, z, cube_width_z, dzzinv);
        const FP qxy = cross_deriv_xy(q_wave_neighborhood, x, y, z, cube_width_x, dxyinv);
        const FP qyz = cross_deriv_yz(q_wave_neighborhood, x, y, z, cube_width_y, dyzinv);
        const FP qxz = cross_deriv_xz(q_wave_neighborhood, x, y, z, cube_width_x, dxzinv);

        const FP cqxx = ch1dxx[idx] * qxx;
        const FP cqyy = ch1dyy[idx] * qyy;
        const FP cqzz = ch1dzz[idx] * qzz;
        const FP cqxy = ch1dxy[idx] * qxy;
        const FP cqxz = ch1dxz[idx] * qxz;
        const FP cqyz = ch1dyz[idx] * qyz;
        const FP h1q = cqxx + cqyy + cqzz + cqxy + cqxz + cqyz;
        const FP h2q = qxx + qyy + qzz - h1q;

        // p-q derivatives, H1(p-q) and H2(p-q)
        const FP h1pmq = h1p - h1q;
        const FP h2pmq = h2p - h2q;

        // rhs of p and q equations
        const FP rhsp = v2px[idx] * h2p + v2pz[idx] * h1q + v2sz[idx] * h1pmq;
        const FP rhsq = v2pn[idx] * h2p + v2pz[idx] * h1q - v2sz[idx] * h2pmq;

        // new p and q
        pwwrite[idx] = FP_LIT(2.0) * pwcentralt1[idx] - pwcentralt2[idx] + rhsp * dt * dt;
        qwwrite[idx] = FP_LIT(2.0) * qwcentralt1[idx] - qwcentralt2[idx] + rhsq * dt * dt;
    }
}

void perturbation_kernel(void *descr[], void *cl_args){
    const perturb_args_t* args = (perturb_args_t*) cl_args;

    const size_t idx = args->source_idx;
    const FP value = args->perturb_value;

    FP *const p_wave_block = (FP *const) STARPU_BLOCK_GET_PTR(descr[0]);
    FP *const q_wave_block = (FP *const) STARPU_BLOCK_GET_PTR(descr[1]);

    p_wave_block[idx] += value;
    q_wave_block[idx] += value;
}
