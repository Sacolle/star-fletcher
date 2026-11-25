#include <criterion/criterion.h>
#include <criterion/theories.h>
#include <criterion/new/assert.h>
#include <criterion/logging.h>

#include <stdint.h>
#include <float.h>

#include "medium.h"
#include "kernel.h"
#include "global.h"
#include "indexing.h"

#define ind(ix,iy,iz) (((iz) * sy + (iy)) * sx + (ix))
// functions extracted from [fletcher-base](https://github.com/gabrielfrtg/fletcher-base/tree/main)
// this one sits in main, so i extratec it to a function
void fletcher_base_medium_initialize(const enum Form prob, size_t sx, size_t sy, size_t sz, 
    FP *restrict vpz, FP *restrict vsv, FP *restrict epsilon,
    FP *restrict delta, FP *restrict phi, FP *restrict theta
){
    size_t i;
    switch (prob){
    case ISO:
        for (i = 0; i < sx * sy * sz; i++)
        {
            vpz[i]     = 3000.0;
            epsilon[i] = 0.0;
            delta[i]   = 0.0;
            phi[i]     = 0.0;
            theta[i]   = 0.0;
            vsv[i]     = 0.0;
        }
        break;
    case VTI:
        for (i = 0; i < sx * sy * sz; i++)
        {
            vpz[i]     = 3000.0;
            epsilon[i] = 0.24;
            delta[i]   = 0.1;
            phi[i]     = 0.0;
            theta[i]   = 0.0;
            if (SIGMA > MAX_SIGMA)
            {
                vsv[i] = 0.0;
            }
            else
            {
                vsv[i] = vpz[i] * sqrtf(fabsf(epsilon[i] - delta[i]) / SIGMA);
            }
        }
        break;

    case TTI:
        for (i = 0; i < sx * sy * sz; i++)
        {
            vpz[i]     = 3000.0;
            epsilon[i] = 0.24;
            delta[i]   = 0.1;
            phi[i]     = 1.0; 
            theta[i] = atanf(1.0);
            if (SIGMA > MAX_SIGMA)
            {
                vsv[i] = 0.0;
            }
            else
            {
                vsv[i] = vpz[i] * sqrtf(fabsf(epsilon[i] - delta[i]) / SIGMA);
            }
        }
    }
}

void RandomVelocityBoundary(int sx, int sy, int sz,
                            int nx, int ny, int nz,
                            int bord, int absorb,
                            float *vpz, float *vsv)
{

    int i, ix, iy, iz;
    int distx, disty, distz, dist;
    int ivelx, ively, ivelz;
    float bordDist;
    float frac, rfac;
    int firstIn, bordLen;
    float maxP, maxS;

    // maximum speed of P and S within bounds
    maxP = 0.0;
    maxS = 0.0;
    for (iz = bord + absorb; iz < nz + bord + absorb; iz++)
    {
        for (iy = bord + absorb; iy < ny + bord + absorb; iy++)
        {
            for (i = ind(bord + absorb, iy, iz); i < ind(nx + bord + absorb, iy, iz); i++)
            {
                maxP = fmaxf(vpz[i], maxP);
                maxS = fmaxf(vsv[i], maxS);
            }
        }
    }

    bordLen = bord + absorb - 1; // last index on low absortion zone
    firstIn = bordLen + 1;       // first index inside input grid
    frac = 1.0/(float)(absorb);

    for (iz = 0; iz < sz; iz++)
    {
        for (iy = 0; iy < sy; iy++)
        {
            for (ix = 0; ix < sx; ix++)
            {
                i = ind(ix, iy, iz);
                // do nothing inside input grid
                if ((iz >= firstIn && iz <= bordLen + nz) &&
                    (iy >= firstIn && iy <= bordLen + ny) &&
                    (ix >= firstIn && ix <= bordLen + nx))
                {
                    continue;
                }
                // random speed inside absortion zone
                else if ((iz >= bord && iz <= (bord + absorb + nz + absorb - 1)) &&
                         (iy >= bord && iy <= (bord + absorb + ny + absorb - 1)) &&
                         (ix >= bord && ix <= (bord + absorb + nx + absorb - 1)))
                {
                    if (iz > bordLen + nz)
                    {
                        distz = iz - bordLen - nz;
                        ivelz = bordLen + nz;
                    }
                    else if (iz < firstIn)
                    {
                        distz = firstIn - iz;
                        ivelz = firstIn;
                    }
                    else
                    {
                        distz = 0;
                        ivelz = iz;
                    }
                    if (iy > bordLen + ny)
                    {
                        disty = iy - bordLen - ny;
                        ively = bordLen + ny;
                    }
                    else if (iy < firstIn)
                    {
                        disty = firstIn - iy;
                        ively = firstIn;
                    }
                    else
                    {
                        disty = 0;
                        ively = iy;
                    }
                    if (ix > bordLen + nx)
                    {
                        distx = ix - bordLen - nx;
                        ivelx = bordLen + nx;
                    }
                    else if (ix < firstIn)
                    {
                        distx = firstIn - ix;
                        ivelx = firstIn;
                    }
                    else
                    {
                        distx = 0;
                        ivelx = ix;
                    }
                    dist = (disty > distz) ? disty : distz;
                    dist = (dist > distx) ? dist : distx;
                    bordDist = (float)(dist)*frac;
                    rfac = (float)rand() / (float)RAND_MAX;

                    vpz[i] = vpz[ind(ivelx, ively, ivelz)] * (1.0 - bordDist) +
                             maxP * rfac * bordDist;
                    vsv[i] = vsv[ind(ivelx, ively, ivelz)] * (1.0 - bordDist) +
                             maxS * rfac * bordDist;

                    /*
                    vpz[i] = bordDist;
                    vsv[i] = bordDist;
                    */
                }
                // null speed at border
                else
                // PPL added {} surrounding vpz and vsv lines below
                {
                    vpz[i] = 0.0;
                    vsv[i] = 0.0;
                }
            }
        }
    }
}


#define FCUT        40.0
#define PICUBE      31.00627668029982017537
#define TWOSQRTPI    3.54490770181103205458
#define THREESQRTPI  5.31736155271654808184

float Source(float dt, int it){
  float tf, fc, fct, expo;
  tf=TWOSQRTPI/FCUT;
  fc=FCUT/THREESQRTPI;
  fct=fc*(((float)it)*dt-tf);
  expo=PICUBE*fct*fct;
  return ((1.0f-2.0f*expo)*expf(-expo));
}

void intermediary_values(
    int sx, int sy , int sz, 
    float *restrict vpz, float *restrict vsv, float *restrict epsilon,
    float *restrict delta, float *restrict phi, float *restrict theta,

    float* ch1dxx, float* ch1dyy, float* ch1dzz, 
    float* ch1dxy, float* ch1dyz, float* ch1dxz, 
    float* v2px, float* v2pz, float* v2sz, float* v2pn
){

    // coeficients of derivatives at H1 operator
    for (int i=0; i<sx*sy*sz; i++) {
        float sinTheta=sin(theta[i]);
        float cosTheta=cos(theta[i]);
        float sin2Theta=sin(2.0*theta[i]);
        float sinPhi=sin(phi[i]);
        float cosPhi=cos(phi[i]);
        float sin2Phi=sin(2.0*phi[i]);
        ch1dxx[i]=sinTheta*sinTheta * cosPhi*cosPhi;
        ch1dyy[i]=sinTheta*sinTheta * sinPhi*sinPhi;
        ch1dzz[i]=cosTheta*cosTheta;
        ch1dxy[i]=sinTheta*sinTheta * sin2Phi;
        ch1dyz[i]=sin2Theta         * sinPhi;
        ch1dxz[i]=sin2Theta         * cosPhi;
    }

    // coeficients of H1 and H2 at PDEs

    for (int i=0; i<sx*sy*sz; i++){
        v2sz[i]=vsv[i]*vsv[i];
        v2pz[i]=vpz[i]*vpz[i];
        v2px[i]=v2pz[i]*(1.0+2.0*epsilon[i]);
        v2pn[i]=v2pz[i]*(1.0+2.0*delta[i]);
    }
}


#define SEED 42

void setup_seed(void)
{
    srand(SEED);
}

TestSuite(medium, .init = setup_seed);

#define EPSILON FLT_EPSILON


struct global_data g;

Test(medium, proper_initialization){

    g.volume_dims = (uvec3_t){10, 10, 10};
    const int sx = 10;
    const int sy = 10;
    const int sz = 10;
    const int size = sx * sy * sz;
    FP vpz_base[size];
    FP vsv_base[size];
    FP epsilon_base[size];
    FP delta_base[size];
    FP phi_base[size];
    FP theta_base[size];

    FP vpz_myimpl[size];
    FP vsv_myimpl[size];
    FP epsilon_myimpl[size];
    FP delta_myimpl[size];
    FP phi_myimpl[size];
    FP theta_myimpl[size];


    for(int i = 0; i < 3; i++){
        const enum Form prob = i;
        fletcher_base_medium_initialize(prob, sx, sy, sz, 
            vpz_base, vsv_base, epsilon_base, delta_base, phi_base, theta_base);
        
        medium_initialize(prob, size, 
            vpz_myimpl, vsv_myimpl, epsilon_myimpl, delta_myimpl, phi_myimpl, theta_myimpl);

        for(int i = 0; i < size; i++){
            cr_assert(epsilon_eq(flt, vpz_base[i], vpz_myimpl[i], EPSILON));
            cr_assert(epsilon_eq(flt, vsv_base[i], vsv_myimpl[i], EPSILON));
            cr_assert(epsilon_eq(flt, epsilon_base[i], epsilon_myimpl[i], EPSILON));
            cr_assert(epsilon_eq(flt, delta_base[i], delta_myimpl[i], EPSILON));
            cr_assert(epsilon_eq(flt, phi_base[i], phi_myimpl[i], EPSILON));
            cr_assert(epsilon_eq(flt, theta_base[i], theta_myimpl[i], EPSILON));
        }
    }

    cr_log_info("Baseline value is %lf and my impl is %lf", vpz_base[157], vpz_myimpl[157]);
    
}

Test(medium, correct_absorb){

    const int bord = 4;
    const int absorb = 4;
    const int nx = 14;
    const int ny = 14;
    const int nz = 14;
    const int sx = nx + 2 * bord + 2 * absorb;
    const int sy = ny + 2 * bord + 2 * absorb;
    const int sz = nz + 2 * bord + 2 * absorb;

    g.volume_dims = (uvec3_t){sx, sy, sz};
    const int size = sx * sy * sz;

    FP vpz_base[size];
    FP vsv_base[size];
    FP epsilon_base[size];
    FP delta_base[size];
    FP phi_base[size];
    FP theta_base[size];

    FP vpz_myimpl[size];
    FP vsv_myimpl[size];
    FP epsilon_myimpl[size];
    FP delta_myimpl[size];
    FP phi_myimpl[size];
    FP theta_myimpl[size];

    fletcher_base_medium_initialize(TTI, sx, sy, sz, 
        vpz_base, vsv_base, epsilon_base, delta_base, phi_base, theta_base);
    
    medium_initialize(TTI, size, 
        vpz_myimpl, vsv_myimpl, epsilon_myimpl, delta_myimpl, phi_myimpl, theta_myimpl);


    setup_seed();
    RandomVelocityBoundary(sx, sy, sz, nx, ny, nz, bord, absorb, vpz_base, vsv_base);

    setup_seed();
    medium_random_velocity_boundary(bord, absorb, vpz_myimpl, vsv_myimpl);


    for(int z = 0; z < sz; z++){
        for(int y = 0; y < sy; y++){
            for(int x = 0; x < sx; x++){
                const size_t i = volume_idx(x, y, z);
                cr_assert(epsilon_eq(flt, vpz_base[i], vpz_myimpl[i], EPSILON), "vpz at (%d, %d, %d)", x, y, z);
                cr_assert(epsilon_eq(flt, vsv_base[i], vsv_myimpl[i], EPSILON), "vsv at (%d, %d, %d)", x, y, z);
            }
        }
    }
}


Test(medium, correct_intermediary_values){
    const size_t seg = 3;
    const size_t c_size = 10;
    g.num_blocks = (uvec3_t){seg, seg, seg};
    g.block_sizes = (uvec3_t){c_size, c_size, c_size};
    g.block_size = c_size * c_size * c_size;

    const int bord = 4;
    const int absorb = 4;
    const int nx = 14;
    const int ny = 14;
    const int nz = 14;
    const int sx = nx + 2 * bord + 2 * absorb;
    const int sy = ny + 2 * bord + 2 * absorb;
    const int sz = nz + 2 * bord + 2 * absorb;

    g.volume_dims = (uvec3_t){sx, sy, sz};
    const int size = sx * sy * sz;

    FP vpz_base[size];
    FP vsv_base[size];
    FP epsilon_base[size];
    FP delta_base[size];
    FP phi_base[size];
    FP theta_base[size];

    FP vpz_myimpl[size];
    FP vsv_myimpl[size];
    FP epsilon_myimpl[size];
    FP delta_myimpl[size];
    FP phi_myimpl[size];
    FP theta_myimpl[size];

    fletcher_base_medium_initialize(TTI, sx, sy, sz, 
        vpz_base, vsv_base, epsilon_base, delta_base, phi_base, theta_base);
    
    medium_initialize(TTI, size, 
        vpz_myimpl, vsv_myimpl, epsilon_myimpl, delta_myimpl, phi_myimpl, theta_myimpl);


    setup_seed();
    RandomVelocityBoundary(sx, sy, sz, nx, ny, nz, bord, absorb, vpz_base, vsv_base);

    setup_seed();
    medium_random_velocity_boundary(bord, absorb, vpz_myimpl, vsv_myimpl);

    const size_t total_seg = seg * seg * seg;
    const size_t total_c_size = c_size * c_size * c_size;
    
    FP* ch1dxx_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,ch1dxx_base, NULL));

    FP* ch1dyy_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,ch1dyy_base, NULL));

    FP* ch1dzz_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,ch1dzz_base, NULL));

    FP* ch1dxy_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,ch1dxy_base, NULL));

    FP* ch1dyz_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,ch1dyz_base, NULL));

    FP* ch1dxz_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,ch1dxz_base, NULL));

    FP* v2px_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,v2px_base, NULL));

    FP* v2pz_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,v2pz_base, NULL));

    FP* v2sz_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,v2sz_base, NULL));

    FP* v2pn_base = (FP*) malloc(sizeof(FP) * size);
    cr_assert(ne(ptr,v2pn_base, NULL));



    // one buffer per block with the MEDIUM_FIELD_COUNT fields one after the other
    FP* medium_myimpl[total_seg];
    for(int seg = 0; seg < total_seg; seg++){
        if((medium_myimpl[seg] = (FP*) malloc(sizeof(FP) * MEDIUM_FIELD_COUNT * total_c_size)) == NULL){
            cr_assert(false);
        }
    }

    intermediary_values(sx, sy, sz, 
        vpz_base, vsv_base, epsilon_base, delta_base, phi_base, theta_base,
        ch1dxx_base, ch1dyy_base, ch1dzz_base, ch1dxy_base, ch1dyz_base, ch1dxz_base, 
        v2px_base, v2pz_base, v2sz_base, v2pn_base
    );

    medium_calc_intermediary_values(
        vpz_myimpl, vsv_myimpl, epsilon_myimpl, delta_myimpl, phi_myimpl, theta_myimpl,
        medium_myimpl
    );

    for(size_t k = 0; k < g.num_blocks.z; k++){
        for(size_t j = 0; j < g.num_blocks.y; j++){
            for(size_t i = 0; i < g.num_blocks.x; i++){
                //index the block
                const size_t b_i = block_idx(i, j , k);

                for(size_t z = 0; z < g.block_sizes.z; z++){
                    for(size_t y = 0; y < g.block_sizes.y; y++){
                        for(size_t x = 0; x < g.block_sizes.x; x++){
                            const size_t c_i = local_idx(x, y, z);
                            const size_t vol_i = local_to_volume_idx(x, y, z, i, j, k);

                            cr_expect(epsilon_eq(flt, ch1dxx_base[vol_i], medium_myimpl[b_i][MEDIUM_CH1DXX * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, ch1dyy_base[vol_i], medium_myimpl[b_i][MEDIUM_CH1DYY * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, ch1dzz_base[vol_i], medium_myimpl[b_i][MEDIUM_CH1DZZ * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, ch1dxy_base[vol_i], medium_myimpl[b_i][MEDIUM_CH1DXY * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, ch1dyz_base[vol_i], medium_myimpl[b_i][MEDIUM_CH1DYZ * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, ch1dxz_base[vol_i], medium_myimpl[b_i][MEDIUM_CH1DXZ * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, v2px_base[vol_i], medium_myimpl[b_i][MEDIUM_V2PX * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, v2pz_base[vol_i], medium_myimpl[b_i][MEDIUM_V2PZ * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, v2sz_base[vol_i], medium_myimpl[b_i][MEDIUM_V2SZ * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);

                            cr_expect(epsilon_eq(flt, v2pn_base[vol_i], medium_myimpl[b_i][MEDIUM_V2PN * total_c_size + c_i], EPSILON), 
                                "vpz at b(%d, %d, %d), c(%d, %d, %d)",i, j, k, x, y, z);
                        }
                    }
                }
            }
        }
    }

    for(int seg = 0; seg < total_seg; seg++){
        free(medium_myimpl[seg]);
    }
    free(ch1dxx_base);
    free(ch1dyy_base);
    free(ch1dzz_base);
    free(ch1dxy_base);
    free(ch1dyz_base);
    free(ch1dxz_base);
    free(v2px_base);
    free(v2pz_base);
    free(v2sz_base);
    free(v2pn_base);
}

Test(medium, kernel_proper_computation){
    const size_t seg = 3;
    const size_t c_size = 10;
    g.num_blocks = (uvec3_t){seg, seg, seg};
    g.block_sizes = (uvec3_t){c_size, c_size, c_size};
    g.block_size = c_size * c_size * c_size;

    const int bord = 4;
    const int absorb = 4;
    const int nx = 14;
    const int ny = 14;
    const int nz = 14;
    const int sx = nx + 2 * bord + 2 * absorb;
    const int sy = ny + 2 * bord + 2 * absorb;
    const int sz = nz + 2 * bord + 2 * absorb;

    g.volume_dims = (uvec3_t){sx, sy, sz};
    const int size = sx * sy * sz;

    FP vpz[size];
    FP vsv[size];
    FP epsilon[size];
    FP delta[size];
    FP phi[size];
    FP theta[size];

    medium_initialize(TTI, size, vpz, vsv, epsilon, delta, phi, theta);

    medium_random_velocity_boundary(bord, absorb, vpz, vsv);

    const size_t total_seg = seg * seg * seg;
    const size_t total_c_size = c_size * c_size * c_size;

    FP* medium[total_seg];
    for(int seg = 0; seg < total_seg; seg++){
        if((medium[seg] = (FP*) malloc(sizeof(FP) * MEDIUM_FIELD_COUNT * total_c_size)) == NULL){
            cr_assert(false);
        }
    }

    medium_calc_intermediary_values(vpz, vsv, epsilon, delta, phi, theta, medium);

    for(size_t k = 0; k < g.num_blocks.z; k++){
        for(size_t j = 0; j < g.num_blocks.y; j++){
            for(size_t i = 0; i < g.num_blocks.x; i++){
                //index the block
                const size_t b_i = block_idx(i, j, k);

            }
        }
    }

    for(int seg = 0; seg < total_seg; seg++){
        free(medium[seg]);
    }
}


Test(medium, source_value){
    const float dt = 1.0;
    for(int i = 0; i < 100; i++){
        const FP proper = Source(dt, i + 1);
        const FP mine = medium_source_value(dt, i + 1);
        cr_expect(epsilon_eq(flt, proper, mine, EPSILON), "i: %d", i);
    }
}


Test(medium, stability_condition){
    g.volume_dims = (uvec3_t){10, 10, 10};
    const float dx = 12.5, dy = 12.5, dz = 12.5;
    const int sx = 10;
    const int sy = 10;
    const int sz = 10;
    const int size = sx * sy * sz;

    FP vpz[size];
    FP vsv[size];
    FP epsilon[size];
    FP delta[size];
    FP phi[size];
    FP theta[size];
 
    for(int prob = 0; prob < 3; prob++){
        medium_initialize(prob, size, vpz, vsv, epsilon, delta, phi, theta);

        medium_random_velocity_boundary(2, 2, vpz, vsv);

        //math done by base fletcher
        float maxvel;
        maxvel=vpz[0]*sqrt(1.0+2*epsilon[0]);
        for(size_t i = 1; i<sx*sy*sz; i++) {
            maxvel = fmaxf(maxvel, vpz[i] * sqrt(1.0 + 2 * epsilon[i]));
        }
        float mindelta = dx;
        if (dy < mindelta)
            mindelta = dy;
        if (dz < mindelta)
            mindelta = dz;
        float recdt;
        recdt = (MI * mindelta) / maxvel;

        const FP my_res = medium_stability_condition(dx, dy, dz, (FP*) vpz, (FP*) epsilon, size);

        cr_expect(epsilon_eq(flt, recdt, my_res, EPSILON));
    }
}