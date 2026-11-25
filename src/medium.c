#include "medium.h"

#include <assert.h>
#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>

#include "argparse.h"
#include "kernel.h"
#include "global.h"
#include "indexing.h"

//return != 0 if err
err_t str_to_medium(const char *str, enum Form* form){
    err_t err;
    int ret;
    if((err = str_to_enum(str, &ret, 3, "ISO", ISO, "VTI", VTI, "TTI", TTI)) != 0){
        return err;
    }
    *form = ret;
    return 0;
}

// TODO: test
// stability condition
FP medium_stability_condition(
    const FP dx, const FP dy, const FP dz,
    FP *restrict vpz, FP *restrict epsilon, const size_t size
){
    FP maxvel = vpz[0] * FP_SQRT(FP_LIT(1.0) + FP_LIT(2.0) * epsilon[0]);
    for (size_t i = 1; i < size; i++){
        maxvel = FP_MAX(maxvel, vpz[i] * FP_SQRT(FP_LIT(1.0) + FP_LIT(2.0) * epsilon[i]));
    }
    const FP mindelta = FP_MIN(FP_MIN(dx, dy), dz);

    const FP recdt = (MI * mindelta) / maxvel;
    return recdt;
}


void medium_initialize(
    const enum Form medium, const size_t size,
    FP *restrict vpz, FP *restrict vsv, FP *restrict epsilon,
    FP *restrict delta, FP *restrict phi, FP *restrict theta
){
    switch (medium) 
    {
    case ISO:
    {
        for (size_t i = 0; i < size; i++){
            vpz[i]     = 3000.0;
            epsilon[i] = 0.0;
            delta[i]   = 0.0;
            phi[i]     = 0.0;
            theta[i]   = 0.0;
            vsv[i]     = 0.0;
        }
    }
    break;
    case VTI: 
    {
        if (SIGMA > MAX_SIGMA){
            //printf("Since sigma (%f) is greater that threshold (%f), sigma is considered infinity and vsv is set to zero\n", SIGMA, MAX_SIGMA);
        }
        for (size_t i = 0; i < size; i++){
            vpz[i]     = 3000.0;
            epsilon[i] = 0.24;
            delta[i]   = 0.1;
            phi[i]     = 0.0;
            theta[i]   = 0.0;
            if (SIGMA > MAX_SIGMA){
                vsv[i] = 0.0;
            }else{
                vsv[i] = vpz[i] * FP_SQRT(FP_ABS(epsilon[i] - delta[i]) / SIGMA);
            }
        }
    }
    break;
    case TTI:
    {
        if (SIGMA > MAX_SIGMA)
        {
            // printf("Since sigma (%f) is greater that threshold (%f), sigma is considered infinity and vsv is set to zero\n", SIGMA, MAX_SIGMA);
        }
        for (size_t i = 0; i < size; i++){
            vpz[i]     = 3000.0;
            epsilon[i] = 0.24;
            delta[i]   = 0.1;
            phi[i]     = 1.0; // evitando coeficientes nulos
            theta[i]   = FP_ATAN(1.0);
            if (SIGMA > MAX_SIGMA){
                vsv[i] = 0.0;
            }else{
                vsv[i] = vpz[i] * FP_SQRT(FP_ABS(epsilon[i] - delta[i]) / SIGMA);
            }
	}
    }
    break;
    default: assert(0 && "Unreachable!");
    }
}
//check if value is in the interval (inclusive)[lb, ub]
static inline bool inbounds(size_t val, size_t lb, size_t ub){
    return val >= lb && val <= ub;
}

// max value of v in the box [start, end) of the volume
FP max_speed_in_cube_segment(FP *v, const uvec3_t start, const uvec3_t end){
    FP max_speed = FP_LIT(0.0);

    for(size_t z = start.z; z < end.z; z++){
        for(size_t y = start.y; y < end.y; y++){
            for(size_t x = start.x; x < end.x; x++){
                const size_t i = volume_idx(x, y, z);
                max_speed = FP_MAX(v[i], max_speed);
            }
        }
    }
    return max_speed;
}

void medium_random_velocity_boundary(
    const size_t border_width, const size_t absorb_width,
    FP *restrict vpz, FP *restrict vsv
){
    const uvec3_t volume_dims = g.volume_dims;
    const size_t outer_width = border_width + absorb_width;

    const uvec3_t inside_start = {outer_width, outer_width, outer_width};
    const uvec3_t inside_end = {volume_dims.x - outer_width, volume_dims.y - outer_width, volume_dims.z - outer_width};

    const FP max_speed_p = max_speed_in_cube_segment(vpz, inside_start, inside_end);
    const FP max_speed_s = max_speed_in_cube_segment(vsv, inside_start, inside_end);

    for(size_t z = 0; z < volume_dims.z; z++){
        for(size_t y = 0; y < volume_dims.y; y++){
            for(size_t x = 0; x < volume_dims.x; x++){

                const size_t i = volume_idx(x, y, z);

                // do nothing inside input grid
                if(inbounds(z, inside_start.z, inside_end.z - 1) &&
                   inbounds(y, inside_start.y, inside_end.y - 1) && 
                   inbounds(x, inside_start.x, inside_end.x - 1)){
                    continue;
                }

                // null speed at border
                // which is a point not within the bounds defided by the aborsb - inner - absobr
                if(!(inbounds(z, border_width, volume_dims.z - border_width - 1) &&
                     inbounds(y, border_width, volume_dims.y - border_width - 1) &&
                     inbounds(x, border_width, volume_dims.x - border_width - 1))
                ){
                    vpz[i] = FP_LIT(0.0);
                    vsv[i] = FP_LIT(0.0);
                    continue;
                }


                // find the greates distance from the border for each direction, 
                // if the direction is inside it's absobsion zone, also save the index of the closest
                // inner point to this absorpsion zone
                #define DIR_AMOUNT 3
                size_t ivel[DIR_AMOUNT] = {x, y, z};
                const size_t dir_start[DIR_AMOUNT] = {inside_start.x, inside_start.y, inside_start.z};
                const size_t dir_end[DIR_AMOUNT] = {inside_end.x, inside_end.y, inside_end.z};
                FP dist = FP_LIT(0.0);
                for(size_t dir_i = 0; dir_i < DIR_AMOUNT; dir_i++){
                    const size_t dir = ivel[dir_i];
                    const size_t start = dir_start[dir_i];
                    const size_t end = dir_end[dir_i];
                    if (dir >= end){ 
                        dist = FP_MAX(dist, (FP)(dir - end + 1 )); 
                        ivel[dir_i] = end - 1; 
                    }else if (dir < start){ 
                        dist = FP_MAX(dist, (FP)(start - dir)); 
                        ivel[dir_i] = start; 
                    }
                }

                // random speed inside absortion zone
                const size_t wave_idx = volume_idx(ivel[0], ivel[1], ivel[2]);
                const FP frac = FP_LIT(1.0)/((FP) absorb_width);
                const FP border_distance = dist * frac;
                const FP random_factor = FP_RAND();

                vpz[i] = vpz[wave_idx] * (1.0 - border_distance) +
                            max_speed_p * random_factor * border_distance;
                vsv[i] = vsv[wave_idx] * (1.0 - border_distance) +
                            max_speed_s * random_factor * border_distance;
                /*
                vpz[i] = border_distance;
                vsv[i] = border_distance;
                */
            }
        }
    }
}

void medium_calc_intermediary_values(
    const FP* vpz, const FP* vsv, const FP* epsilon,
    const FP* delta, const FP* phi, const FP* theta,
    FP **restrict medium
){
    // distance between two fields of the medium of a block
    const size_t ldt = g.block_size;

    for(size_t k = 0; k < g.num_blocks.z; k++){
        for(size_t j = 0; j < g.num_blocks.y; j++){
            for(size_t i = 0; i < g.num_blocks.x; i++){
                //index the block
                const size_t b_i = block_idx(i, j , k);
                FP *const block_medium = medium[b_i];

                for(size_t z = 0; z < g.block_sizes.z; z++){
                    for(size_t y = 0; y < g.block_sizes.y; y++){
                        for(size_t x = 0; x < g.block_sizes.x; x++){
                            const size_t c_i = local_idx(x, y, z);
                            const size_t vol_i = local_to_volume_idx(x, y, z, i, j, k);

                            // changing the data type TO be the same as in the original
                            const FP sinTheta = sin(theta[vol_i]);
                            const FP cosTheta = cos(theta[vol_i]);
                            const FP sin2Theta = sin(2.0 * theta[vol_i]);
                            const FP sinPhi = sin(phi[vol_i]);
                            const FP cosPhi = cos(phi[vol_i]);
                            const FP sin2Phi = sin(2.0 * phi[vol_i]);

                            block_medium[MEDIUM_CH1DXX * ldt + c_i] = sinTheta * sinTheta * cosPhi * cosPhi;
                            block_medium[MEDIUM_CH1DYY * ldt + c_i] = sinTheta * sinTheta * sinPhi * sinPhi;
                            block_medium[MEDIUM_CH1DZZ * ldt + c_i] = cosTheta * cosTheta;
                            block_medium[MEDIUM_CH1DXY * ldt + c_i] = sinTheta * sinTheta * sin2Phi;
                            block_medium[MEDIUM_CH1DYZ * ldt + c_i] = sin2Theta * sinPhi;
                            block_medium[MEDIUM_CH1DXZ * ldt + c_i] = sin2Theta * cosPhi;

                            // coeficients of H1 and H2 at PDEs
                            const FP l_v2pz = vpz[vol_i] * vpz[vol_i];
                            block_medium[MEDIUM_V2SZ * ldt + c_i] = vsv[vol_i] * vsv[vol_i];
                            block_medium[MEDIUM_V2PZ * ldt + c_i] = l_v2pz;
                            block_medium[MEDIUM_V2PX * ldt + c_i] = l_v2pz * (1.0 + 2.0 * epsilon[vol_i]);
                            block_medium[MEDIUM_V2PN * ldt + c_i] = l_v2pz * (1.0 + 2.0 * delta[vol_i]);
                        }
                    }
                }
            }
        }
    }
}


// omissão do FP_LIT para obter os mesmos resultados do que o fletcher base
#define FCUT        40.0
#define PICUBE      31.00627668029982017537
#define TWOSQRTPI   3.54490770181103205458
#define THREESQRTPI 5.31736155271654808184

//#include <stdio.h>

FP medium_source_value(const FP dt, const int64_t it){
    const FP tf = TWOSQRTPI / FCUT;
    const FP fc = FCUT / THREESQRTPI;
    const FP fct = fc * (((FP) it) * dt - tf);
    const FP expo = PICUBE * fct * fct;
    const FP result = ((FP_LIT(1.0) - FP_LIT(2.0) * expo) * FP_EXP(-expo));

    //printf("source: %lf\n", result);
    return result;
}
