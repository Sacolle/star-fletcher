#include <starpu.h>
#include <stdint.h>
#include <stdbool.h>

#include <stdio.h>
#include <math.h>

#define ERR_NAME_IMPL
#include "err.h"
#include "kernel.h"
#include "global.h"
#include "indexing.h"
#include "floatingpoint.h"
#include "argparse.h"
#include "medium.h"
#include "mem.h"
#include "io.h"
#include "halo.h"

#define DEFAULT_OUTPUT_FOLDER "results"
#define DEFAULT_OUTPUT_NAME "form"

const size_t BORDER_WIDTH = 4;

struct global_data g = {
    .dt_output = FP_LIT(0.01),
    .cuda_threads = {16, 4, 4},
};

// passes the iter[0] to iter[1], iter[1] to iter[2] and iter[2] to iter[0]
static inline void rotate_for_next_iter(starpu_data_handle_t* iter[3]){
    starpu_data_handle_t* tmp = iter[2]; 
    iter[2] = iter[1]; 
    iter[1] = iter[0]; 
    iter[0] = tmp; 
}


typedef struct dump_block_args {
    size_t i;
    size_t j;
    size_t k;
    size_t t;
} dump_block_args_t;

int make_dump_block_args(dump_block_args_t** args, size_t i, size_t j, size_t k, size_t t){
    //NOTE: isso quebra o meu script de reconstrução?
    if((*args = (dump_block_args_t*) malloc(sizeof(dump_block_args_t))) == NULL){
        return errno;
    }
    **args = (dump_block_args_t) {i, j, k, t};
    return 0;
}


void dump_block_kernel(void *descr[], void *cl_args){
    const dump_block_args_t* args = (dump_block_args_t*) cl_args;
    FP* block = (FP*) STARPU_BLOCK_GET_PTR(descr[0]);

    const size_t nx = STARPU_BLOCK_GET_NX(descr[0]);
    const size_t ny = STARPU_BLOCK_GET_NY(descr[0]);
    const size_t nz = STARPU_BLOCK_GET_NZ(descr[0]);

    io_state_write_file(args->i * g.block_sizes.x, args->j * g.block_sizes.y, args->k * g.block_sizes.z, args->t,
        nx, ny, nz, STARPU_BLOCK_GET_LDY(descr[0]), STARPU_BLOCK_GET_LDZ(descr[0]), block);
}

//#include <linux/time.h>
#define NS_PER_SECOND 1000000000ULL
#define SECONDS_PER_NS 1e-9

uint64_t get_timestamp_ns() {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((uint64_t)ts.tv_sec * NS_PER_SECOND) + ts.tv_nsec;
}

double elapsed_seconds(const uint64_t t_end, const uint64_t t_start){
    return ((double)t_end - t_start) * SECONDS_PER_NS;
}


static struct starpu_perfmodel dump_block_perf_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "dump_block_perf_model"
};

struct starpu_codelet dump_block_codelet = {
    .cpu_funcs = {dump_block_kernel},
    .nbuffers = 1,
    .modes = {STARPU_R}, // p wave
    .model = &dump_block_perf_model
};

static struct starpu_perfmodel rtm_perf_model = {
    .type = STARPU_HISTORY_BASED,
    .symbol = "rtm_perf_model"
};

#ifdef CUDA_BACKEND
extern void rtm_kernel_cuda(void *descr[], void *cl_args);
extern void pack_kernel_cuda(void *descr[], void *cl_args);
#endif


#ifdef CUDA_BACKEND
  #ifndef NO_CPU_KERNEL
    #define WHERE STARPU_CPU | STARPU_CUDA
  #else
    #define WHERE STARPU_CUDA
  #endif
  #else
  #define WHERE STARPU_CPU
#endif


// computes the level t of a block. the amount of buffers depends on how many neighbors the block has,
// so they and their modes are set per task, in the order of rtm_handle_pos_* (kernel.h)
struct starpu_codelet rtm_codelet = {
    .cpu_funcs = { rtm_kernel },
#ifdef CUDA_BACKEND
    .cuda_funcs = { rtm_kernel_cuda },
    .cuda_flags = { STARPU_CUDA_ASYNC },
#endif
    .where = WHERE,
    .nbuffers = STARPU_VARIABLE_NBUFFERS,
    .model = &rtm_perf_model
};

// copies the initial level of a block into the slabs it sends, buffers in the order of pack_handle_pos_* (kernel.h)
struct starpu_codelet pack_codelet = {
    .cpu_funcs = { pack_kernel },
#ifdef CUDA_BACKEND
    .cuda_funcs = { pack_kernel_cuda },
    .cuda_flags = { STARPU_CUDA_ASYNC },
#endif
    .where = WHERE,
    .nbuffers = STARPU_VARIABLE_NBUFFERS,
    // no perf model: it runs once per block per execution, too few times to ever be calibrated
};

// sets the buffer at `position` of a task whose codelet has STARPU_VARIABLE_NBUFFERS
static inline void set_buffer(struct starpu_task* task, int position, starpu_data_handle_t handle, enum starpu_data_access_mode mode){
    STARPU_TASK_SET_HANDLE(task, handle, position);
    STARPU_TASK_SET_MODE(task, mode, position);
}

// registers the slab a block sends towards the neighbor at entry table_index of NEIGHBOR_OFFSETS,
// p at field 0 and q at field 1. StarPU allocates it where it is first written (home node -1)
static void register_slab(starpu_data_handle_t* handle, int table_index){
    const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
    const uvec3_t dims = slab_dims(offset);
    starpu_tensor_data_register(handle, -1, 0, dims.x, dims.x * dims.y, uvec3_prod(dims), dims.x, dims.y, dims.z, 2, sizeof(FP));
}

err_t write_wave(int64_t* n_out, starpu_data_handle_t* wave_iter){

    err_t err = 0;
    // salva o primeiro bloco (nulo)
    for(size_t k = 0; k < g.num_blocks.z; k++)
    for(size_t j = 0; j < g.num_blocks.y; j++)
    for(size_t i = 0; i < g.num_blocks.x; i++){

        // call the write task
        dump_block_args_t* dump_args;
        if((err = make_dump_block_args(&dump_args, i, j, k, *n_out)) != 0){
            return err;
        }

        struct starpu_task* dump_block_task = starpu_task_create();
        dump_block_task->name = "write_block";
        dump_block_task->cl = &dump_block_codelet;
        dump_block_task->cl_arg = dump_args;
        dump_block_task->cl_arg_size = sizeof(dump_block_args_t);
        dump_block_task->cl_arg_free = 1;
        //NOTE: mudado para p_wave_iter[0] para se adequar melhor com o trace do fletcher base
        //      no caso não importa pois write(t = 0) é só zeros no original
        dump_block_task->handles[0] = wave_iter[block_idx(i, j, k)];

        if((err = starpu_task_submit(dump_block_task)) != 0){
            return err;
        };
    }

    (*n_out)++;
    return 0;
}


// if x fails (!= 0), goto the end of main and log status;
#define TRY(x, ...)							\
  do {									\
    err_t err_try;							\
    TRYTO(err_try = (x), program_status = EXIT_FAILURE;			\
        printf("[error] Failed at line %d with error %ld: %s\n", __LINE__, \
               err_try, err_name(err_try));				\
        printf("[err-msg] " __VA_ARGS__);, program_end)			\
} while(0)

// if x fails (!= 0) emits a message, but
// keeps the executions as is
#define CHECK(x, ...)                                                          \
  do {                                                                         \
    err_t err_check = x;                                                       \
    if (err_check != 0) {                                                      \
      printf("[warnign] Check at line %d detected error %ld: %s\n", __LINE__,  \
             err_check, err_name(err_check));				\
      printf("[warning-msg] " __VA_ARGS__);			\
	}								\
} while(0)


#ifdef RELEASE
#define DEBUG(...) 
#else
#define DEBUG(...) printf("[debug] " __VA_ARGS__)
#endif

#ifndef RANDOM_SEED
#define RANDOM_SEED 0
#endif

int main(int argc, char **argv){

    // line buffer stdout, so the logs are not lost if the process is killed
    // while the output is redirected to a file
    setvbuf(stdout, NULL, _IOLBF, 0);

    const uint64_t initialization_start_time = get_timestamp_ns();
    srand(RANDOM_SEED);

    //need to be toplevel for the try macro
    int program_status = EXIT_SUCCESS;
    mem_vec_t static_allocs = NULL;
    mem_vec_t medium_allocs = NULL;

#ifdef CUDA_BACKEND
    printf("Using CUDA backend.\n");

    CHECK(get_envvar(&g.cuda_threads.x, "CUDA_THREAD_X"),
            "Not setting CUDA_THREAD_X, using default: %u\n", g.cuda_threads.x);
    CHECK(get_envvar(&g.cuda_threads.y, "CUDA_THREAD_Y"),
            "Not setting CUDA_THREAD_Y, using default: %u\n", g.cuda_threads.y);
    CHECK(get_envvar(&g.cuda_threads.z, "CUDA_THREAD_Z"),
            "Not setting CUDA_THREAD_Z, using default: %u\n", g.cuda_threads.z);

    // the stencil kernel uses many registers per thread, so it can not be launched with blocks as large as other kernels
    const int stencil_max_threads = stencil_max_threads_per_block();
    if(stencil_max_threads > 0){
        TRY(uvec3_prod(g.cuda_threads) <= (size_t) stencil_max_threads ? 0 : ME_IMPROPER_INPUT,
            "CUDA_THREAD_X * CUDA_THREAD_Y * CUDA_THREAD_Z = %ld, but the stencil kernel supports at most %d threads per block.\n",
            uvec3_prod(g.cuda_threads), stencil_max_threads);
    }else{
        printf("[warning-msg] Could not query the CUDA device, the CUDA_THREAD_* values are not checked.\n");
    }
#else
    printf("Not using CUDA.\n");
#endif

#ifdef RELEASE
    printf("Compiled as release.\n");
#else
    printf("Compiled as debug.\n");
#endif
    

    char* output_folder = DEFAULT_OUTPUT_FOLDER;
    CHECK(get_envvar(&output_folder, "OUTPUT_FOLDER"), 
            "Not setting a folder, using default: %s\n", DEFAULT_OUTPUT_FOLDER);

    char* output_filename = DEFAULT_OUTPUT_NAME;
    CHECK(get_envvar(&output_filename, "OUTPUT_FILE"),
            "Not setting a filename, using default: %s\n", DEFAULT_OUTPUT_NAME);

    int64_t enable_io = 1;
    CHECK(get_envvar(&enable_io, "ENABLE_IO"), "By default, using IO.\n");

    enum Form form = 0;
    char* form_str = NULL;
    uint32_t nx, ny, nz, absorb_width;
    FP dx, dy, dz, dt, tmax;

    TRY(read_args(&argc, argv, 14,
        ARG_str, &form_str,
        ARG_u32, &nx, ARG_u32, &ny, ARG_u32, &nz, ARG_u32, &absorb_width,
        FP_ARG, &dx, FP_ARG, &dy, FP_ARG, &dz, FP_ARG, &dt, FP_ARG, &tmax,
        ARG_u32, &g.num_blocks.x, ARG_u32, &g.num_blocks.y, ARG_u32, &g.num_blocks.z,
        FP_ARG, &g.dt_output
    ), "Error in parsing arg at %d.\n"
       "usage: main FORM nx ny nz absorb_width dx dy dz dt tmax blocks_x blocks_y blocks_z dt_output\n", argc);

    TRY(str_to_medium(form_str, &form), "Failed at string to medium conversion");

    //fazendo dessa forma para ficar igual ao fletcher base
    g.volume_dims = (uvec3_t){
        nx + 2 * absorb_width + 2 * BORDER_WIDTH,
        ny + 2 * absorb_width + 2 * BORDER_WIDTH,
        nz + 2 * absorb_width + 2 * BORDER_WIDTH,
    };
    TRY(g.num_blocks.x > 0 && g.num_blocks.y > 0 && g.num_blocks.z > 0 ? 0 : ME_IMPROPER_INPUT,
        "The amount of blocks in each axis must be greater than 0.\n");

    // the number of blocks divides the volume in every axis
    TRY(g.volume_dims.x % g.num_blocks.x == 0 &&
        g.volume_dims.y % g.num_blocks.y == 0 &&
        g.volume_dims.z % g.num_blocks.z == 0 ? 0 : ME_COUNT_DONT_MATCH,
        "A largura do volume + kernel size devem ser divisíveis pela quantidade de blocos em cada eixo.\n");

    g.block_sizes = (uvec3_t){
        g.volume_dims.x / g.num_blocks.x,
        g.volume_dims.y / g.num_blocks.y,
        g.volume_dims.z / g.num_blocks.z,
    };
    g.block_size = uvec3_prod(g.block_sizes);
    g.total_blocks = uvec3_prod(g.num_blocks);

    // the halo is STENCIL_RADIUS wide on the split axes, a slab must fit inside a single block
    TRY(halo_init(g.num_blocks, g.block_sizes, &g.halo_widths),
        "On a split axis the blocks must be at least %d points wide.\n", STENCIL_RADIUS);

    // the block at the center has the most neighbors, its update task has the most buffers
    const int max_neighbor_count = __builtin_popcount(halo_neighbor_mask(g.num_blocks.x / 2, g.num_blocks.y / 2, g.num_blocks.z / 2));
    TRY(rtm_handle_pos_count(max_neighbor_count) <= STARPU_NMAXBUFS ? 0 : ME_COUNT_DONT_MATCH,
        "StarPU was built with %d buffers per task, but an update task needs %d. Raise maxBuffers in star-fletcher.nix.\n",
        STARPU_NMAXBUFS, rtm_handle_pos_count(max_neighbor_count));

    // the stencil uses 32 bit strides: up to 4 planes of the scratch buffer away from the center
    const size_t scratch_plane = (size_t) (g.block_sizes.x + 2 * g.halo_widths.x) * (g.block_sizes.y + 2 * g.halo_widths.y);
    TRY(STENCIL_RADIUS * scratch_plane <= INT32_MAX ? 0 : ME_OVERFLOW_CONVERT,
        "A block plane of %u x %u points is too large for the 32 bit strides of the stencil.\n", g.block_sizes.x, g.block_sizes.y);

#ifdef CUDA_BACKEND
    // the GPU kernels index the blocks, the slabs, the medium and the scratch with 32 bit integers
    const size_t scratch_points = scratch_plane * (g.block_sizes.z + 2 * g.halo_widths.z);
    TRY(MEDIUM_FIELD_COUNT * g.block_size <= UINT32_MAX && 2 * scratch_points <= UINT32_MAX ? 0 : ME_OVERFLOW_CONVERT,
        "Blocks of %u x %u x %u points are too large for the 32 bit indices of the GPU kernels, split the volume in more blocks.\n",
        g.block_sizes.x, g.block_sizes.y, g.block_sizes.z);
#endif

    // amount of steps to reach tmax. a quotient within float rounding of an integer (e.g. 0.0003 / 0.0001 = 3.0000002)
    // is that integer, otherwise it is rounded up so the simulation covers tmax
    const double steps_to_tmax = (double) tmax / (double) dt;
    const int64_t st = fabs(steps_to_tmax - round(steps_to_tmax)) < 1e-4 ? llround(steps_to_tmax) : (int64_t) ceil(steps_to_tmax);

    // the point in the global volume that the source is inserted
    const size_t volume_propagation_idx = volume_idx(g.volume_dims.x / 2, g.volume_dims.y / 2, g.volume_dims.z / 2);
    DEBUG("Index of propagation is %ld.\n", volume_propagation_idx);


    char bin_filename[512] = {'\0'};
    sprintf(bin_filename, "%s/out-%s.rsf@", output_folder, output_filename);

    int ret = starpu_init(NULL);
    STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");

    // the tasks that write the output only have a CPU implementation
    TRY(!enable_io || starpu_cpu_worker_get_count() > 0 ? 0 : ME_IMPROPER_INPUT,
        "Writing the output needs at least one CPU worker: set STARPU_NCPU >= 1 or ENABLE_IO=0.\n"
        "To run every computation on the GPU, build with NO_CPU_KERNEL=1 and run with STARPU_NCPU=1.\n");

    // amount of steps between two saved moments, counted in steps so float rounding can not skip a moment
    int64_t steps_per_output = 1;

    //init the IO here, only if it will be used
    if(enable_io){
        printf("bin filename %s\n", bin_filename);

        // with output time 0, tmax / g.dt_output is inf, and casting it to size_t is UB
        TRY(g.dt_output > FP_LIT(0.0) ? 0 : ME_IMPROPER_INPUT,
            "Output time must be greater than 0 when IO is enabled.\n");

        // the output time is rounded to a whole number of steps, at least one
        steps_per_output = llround((double) g.dt_output / (double) dt);
        if(steps_per_output < 1) steps_per_output = 1;

        // amount of moments the program will save: the level 0, then one every steps_per_output steps
        const size_t total_saved_moments = (size_t) (st / steps_per_output) + 1;
        TRY(io_state_init(bin_filename, 4096, total_saved_moments, uvec3_prod(g.volume_dims)));
    }

    // mem_allocate the buffers that will be used in computing the 
    // intermediary values
    FP *vpz, *vsv, *epsilon, *delta, *phi, *theta;
    const size_t medium_size = sizeof(FP) * uvec3_prod(g.volume_dims);

    TRY(mem_allocate_local(&medium_allocs, (void**) &vpz, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &vsv, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &epsilon, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &delta, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &phi, medium_size));
    TRY(mem_allocate_local(&medium_allocs, (void**) &theta, medium_size));

    // inicialize the buffers above based on the type of medium
    medium_initialize(form, uvec3_prod(g.volume_dims), vpz, vsv, epsilon, delta, phi, theta);
    
    // set the absorption zone for vpz and vsv
    medium_random_velocity_boundary(BORDER_WIDTH, absorb_width, vpz, vsv);

    // run the stability condition for the size
    const FP stability_condition = medium_stability_condition(dx, dy, dz, vpz, epsilon, uvec3_prod(g.volume_dims));
    DEBUG("The stability condition (proper value for dt) for this problem is %lf.\n", stability_condition);


    // one buffer per block with the MEDIUM_FIELD_COUNT fields one after the other, see enum medium_field
    FP **medium;
    TRY(mem_allocate(&static_allocs, (void**) &medium, g.total_blocks * sizeof(FP*)));
    for(size_t i = 0; i < g.total_blocks; i++){
        TRY(mem_allocate(&static_allocs, (void**)(medium + i), MEDIUM_FIELD_COUNT * g.block_size * sizeof(FP)));
    }

    medium_calc_intermediary_values(vpz, vsv, epsilon, delta, phi, theta, medium);

    //at this point the values for the medium will not be used again
    mem_free_local(&medium_allocs);

    // the levels of p and q of every block: [0] is the level t being computed, [1] is t - 1 and [2] is t - 2.
    // every level of a block is a new handle (home node -1), written once by the update task that computes it
    // and unregistered right after the last task that reads it is submitted, so no task ever waits for an older reader
    starpu_data_handle_t* p_wave_iter[3];
    starpu_data_handle_t* q_wave_iter[3];
    for(int level = 0; level < 3; level++){
        TRY(mem_allocate(&static_allocs, (void**) &p_wave_iter[level], g.total_blocks * sizeof(starpu_data_handle_t)));
        TRY(mem_allocate(&static_allocs, (void**) &q_wave_iter[level], g.total_blocks * sizeof(starpu_data_handle_t)));
    }

    // the slab each block sends to each neighbor, at [block * NEIGHBOR_COUNT + table_index].
    // slabs[0] holds the ones of the level t - 1, read by the update of the level t,
    // slabs[1] the ones of the level t, written by it. new handles written once, like the levels
    starpu_data_handle_t* slabs[2];
    TRY(mem_allocate(&static_allocs, (void**) &slabs[0], g.total_blocks * NEIGHBOR_COUNT * sizeof(starpu_data_handle_t)));
    TRY(mem_allocate(&static_allocs, (void**) &slabs[1], g.total_blocks * NEIGHBOR_COUNT * sizeof(starpu_data_handle_t)));

    starpu_data_handle_t *hdl_medium;
    TRY(mem_allocate(&static_allocs, (void**) &hdl_medium, g.total_blocks * sizeof(starpu_data_handle_t)));

    for(size_t idx = 0; idx < g.total_blocks; idx++){
        starpu_tensor_data_register(hdl_medium + idx, STARPU_MAIN_RAM, (uintptr_t) medium[idx],
            g.block_sizes.x, g.block_sizes.x * g.block_sizes.y, g.block_size,
            g.block_sizes.x, g.block_sizes.y, g.block_sizes.z, MEDIUM_FIELD_COUNT, sizeof(FP));
    }

    // the padded working copy of the update task, p at field 0 and q at field 1.
    // with STARPU_SCRATCH each worker gets its own buffer, allocated on first use and never transferred
    const uvec3_t scratch_dims = {
        g.block_sizes.x + 2 * g.halo_widths.x,
        g.block_sizes.y + 2 * g.halo_widths.y,
        g.block_sizes.z + 2 * g.halo_widths.z,
    };
    starpu_data_handle_t scratch_handle;
    starpu_tensor_data_register(&scratch_handle, -1, 0, scratch_dims.x, scratch_dims.x * scratch_dims.y,
        uvec3_prod(scratch_dims), scratch_dims.x, scratch_dims.y, scratch_dims.z, 2, sizeof(FP));

    // the initial levels: the level -1 is zero, the level 0 is zero plus the source value of t = 0 at the source point.
    // the zero buffer is shared by the read only handles, the block of the source has its own p and q buffers
    FP *zero_block, *source_p_block, *source_q_block;
    TRY(mem_allocate(&static_allocs, (void**) &zero_block, g.block_size * sizeof(FP)));
    TRY(mem_allocate(&static_allocs, (void**) &source_p_block, g.block_size * sizeof(FP)));
    TRY(mem_allocate(&static_allocs, (void**) &source_q_block, g.block_size * sizeof(FP)));
    for(size_t point = 0; point < g.block_size; point++){
        zero_block[point] = source_p_block[point] = source_q_block[point] = FP_LIT(0.0);
    }

    const size_t source_block = volume_to_block_idx(volume_propagation_idx);
    const size_t source_local_idx = volume_to_local_idx(volume_propagation_idx);
    const uvec3_t source_point = unindex(source_local_idx, g.block_sizes.x, g.block_sizes.y);
    source_p_block[source_local_idx] = medium_source_value(dt, 0);
    source_q_block[source_local_idx] = medium_source_value(dt, 0);

    #define BLOCK_REGISTER(handle, ptr) starpu_block_data_register((handle), STARPU_MAIN_RAM, (uintptr_t) (ptr), \
        g.block_sizes.x, g.block_sizes.x * g.block_sizes.y, g.block_sizes.x, g.block_sizes.y, g.block_sizes.z, sizeof(FP))

    for(size_t block = 0; block < g.total_blocks; block++){
        BLOCK_REGISTER(&p_wave_iter[2][block], zero_block);
        BLOCK_REGISTER(&q_wave_iter[2][block], zero_block);
        BLOCK_REGISTER(&p_wave_iter[1][block], block == source_block ? source_p_block : zero_block);
        BLOCK_REGISTER(&q_wave_iter[1][block], block == source_block ? source_q_block : zero_block);
    }

    #undef BLOCK_REGISTER

    // copies the level 0 into the slabs read by the update of the level 1
    for(size_t k = 0; k < g.num_blocks.z; k++)
    for(size_t j = 0; j < g.num_blocks.y; j++)
    for(size_t i = 0; i < g.num_blocks.x; i++){
        const size_t block = block_idx(i, j, k);
        const uint32_t neighbor_mask = halo_neighbor_mask(i, j, k);
        if(neighbor_mask == 0){
            continue;
        }

        pack_args_t* pack_args;
        TRY(make_pack_args(&pack_args, neighbor_mask));

        struct starpu_task* task = starpu_task_create();
        task->name = "pack_initial_level";
        task->cl = &pack_codelet;
        task->cl_arg = pack_args;
        task->cl_arg_size = sizeof(pack_args_t);
        task->cl_arg_free = 1;
        task->nbuffers = pack_handle_pos_count(__builtin_popcount(neighbor_mask));
        set_buffer(task, PACK_HANDLE_POS_P, p_wave_iter[1][block], STARPU_R);
        set_buffer(task, PACK_HANDLE_POS_Q, q_wave_iter[1][block], STARPU_R);

        int slab_position = 0;
        for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
            if(!(neighbor_mask & NBIT(table_index))) continue;
            starpu_data_handle_t* slab = &slabs[0][block * NEIGHBOR_COUNT + table_index];
            register_slab(slab, table_index);
            set_buffer(task, pack_handle_pos_sent_slab(slab_position++), *slab, STARPU_W);
        }

        TRY(starpu_task_submit(task));
    }

    const uint64_t initialization_end_time = get_timestamp_ns();

    const double initialization_total_time = elapsed_seconds(initialization_end_time, initialization_start_time);

    printf("Initialization Elapsed time is: %lfs\n", initialization_total_time);


    int64_t n_out = 0;
    // salva o nível 0
    if(enable_io){
        TRY(write_wave(&n_out, p_wave_iter[1]));
    }

    const uint64_t start_time = get_timestamp_ns();

    for(int64_t t = 1; t <= st; t++){
        starpu_iteration_push(t);

        // added to the level t by the update task of the block of the source
        const FP source_value = medium_source_value(dt, t);

        for(size_t k = 0; k < g.num_blocks.z; k++) // z
        for(size_t j = 0; j < g.num_blocks.y; j++) // y
        for(size_t i = 0; i < g.num_blocks.x; i++){// x
            const size_t block = block_idx(i, j, k);
            const uint32_t neighbor_mask = halo_neighbor_mask(i, j, k);
            const int neighbor_count = __builtin_popcount(neighbor_mask);

            //let starpu allocate the data by setting home_node = -1 
            starpu_block_data_register(&p_wave_iter[0][block], -1, 0,
                g.block_sizes.x, g.block_sizes.x * g.block_sizes.y,
                g.block_sizes.x, g.block_sizes.y, g.block_sizes.z, sizeof(FP)
            );
            starpu_block_data_register(&q_wave_iter[0][block], -1, 0,
                g.block_sizes.x, g.block_sizes.x * g.block_sizes.y,
                g.block_sizes.x, g.block_sizes.y, g.block_sizes.z, sizeof(FP)
            );

            // if v == true, return the border size, else 0
            #define AS_BORDER(v) ((v) ? BORDER_WIDTH : 0)

            const rtm_args_t values = {
                .x_start = AS_BORDER(i == 0), .x_end = g.block_sizes.x - AS_BORDER(i == g.num_blocks.x - 1),
                .y_start = AS_BORDER(j == 0), .y_end = g.block_sizes.y - AS_BORDER(j == g.num_blocks.y - 1),
                .z_start = AS_BORDER(k == 0), .z_end = g.block_sizes.z - AS_BORDER(k == g.num_blocks.z - 1),
                .dx = dx, .dy = dy, .dz = dz, .dt = dt,
                .neighbor_mask = neighbor_mask,
                .has_source = block == source_block,
                .source_x = source_point.x, .source_y = source_point.y, .source_z = source_point.z,
                .source_value = source_value,
                .t = t,
            };

            #undef AS_BORDER

            struct rtm_args* rtm_args;
            TRY(make_rtm_args(&rtm_args, &values));

            struct starpu_task* task = starpu_task_create();
            task->name = "wave_propagation";
            task->cl = &rtm_codelet;
            task->cl_arg = rtm_args;
            task->cl_arg_size = sizeof(struct rtm_args);
            task->cl_arg_free = 1; // free the args after use
            task->nbuffers = rtm_handle_pos_count(neighbor_count);

            set_buffer(task, RTM_HANDLE_POS_MEDIUM, hdl_medium[block], STARPU_R);
            set_buffer(task, RTM_HANDLE_POS_P_T, p_wave_iter[0][block], STARPU_W);
            set_buffer(task, RTM_HANDLE_POS_Q_T, q_wave_iter[0][block], STARPU_W);
            set_buffer(task, RTM_HANDLE_POS_P_T1, p_wave_iter[1][block], STARPU_R);
            set_buffer(task, RTM_HANDLE_POS_Q_T1, q_wave_iter[1][block], STARPU_R);
            set_buffer(task, RTM_HANDLE_POS_P_T2, p_wave_iter[2][block], STARPU_R);
            set_buffer(task, RTM_HANDLE_POS_Q_T2, q_wave_iter[2][block], STARPU_R);
            set_buffer(task, RTM_HANDLE_POS_SCRATCH, scratch_handle, STARPU_SCRATCH);

            // the slab received to do some computation
            int slab_position = 0;
            for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
                if(!(neighbor_mask & NBIT(table_index))) continue;
                const ivec3_t offset = NEIGHBOR_OFFSETS[table_index];
                const size_t neighbor = block_idx(i + offset.x, j + offset.y, k + offset.z);
                const starpu_data_handle_t received = slabs[0][neighbor * NEIGHBOR_COUNT + halo_opposite_neighbor(table_index)];
                set_buffer(task, rtm_handle_pos_received_slab(slab_position++), received, STARPU_R);
            }

            // the slab to send after the computation
            slab_position = 0;
            for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
                if(!(neighbor_mask & NBIT(table_index))) continue;
                starpu_data_handle_t* sent = &slabs[1][block * NEIGHBOR_COUNT + table_index];
                register_slab(sent, table_index);
                set_buffer(task, rtm_handle_pos_sent_slab(neighbor_count, slab_position++), *sent, STARPU_W);
            }

            TRY(starpu_task_submit(task));
        }

        // save the level every steps_per_output steps
        if(enable_io && t % steps_per_output == 0){
            TRY(write_wave(&n_out, p_wave_iter[0]));
        }

        // the level t - 2 and the slabs of the level t - 1 were only read by the updates just submitted
        for(size_t k = 0; k < g.num_blocks.z; k++)
        for(size_t j = 0; j < g.num_blocks.y; j++)
        for(size_t i = 0; i < g.num_blocks.x; i++){
            const size_t block = block_idx(i, j, k);
            starpu_data_unregister_submit(p_wave_iter[2][block]);
            starpu_data_unregister_submit(q_wave_iter[2][block]);
            const uint32_t neighbor_mask = halo_neighbor_mask(i, j, k);
            for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
                if(!(neighbor_mask & NBIT(table_index))) continue;
                starpu_data_unregister_submit(slabs[0][block * NEIGHBOR_COUNT + table_index]);
            }
        }

        // the level t becomes the level t - 1 of the next iteration
        rotate_for_next_iter(p_wave_iter);
        rotate_for_next_iter(q_wave_iter);
        starpu_data_handle_t* swap = slabs[0];
        slabs[0] = slabs[1];
        slabs[1] = swap;

        starpu_iteration_pop();
    }
    DEBUG("Submitted all tasks\n");
    //at least after all iterations
    starpu_task_wait_for_all();

    const uint64_t end_time = get_timestamp_ns();

    const double total_elapsed_time = elapsed_seconds(end_time, start_time);

    printf("Computation Elapsed time is: %lfs\n", total_elapsed_time);

    const uint64_t samplesPropagate = (g.volume_dims.x - 2 * BORDER_WIDTH) * (g.volume_dims.y - 2 * BORDER_WIDTH) * (g.volume_dims.z - 2 * BORDER_WIDTH);
    const uint64_t totalSamples = samplesPropagate * (uint64_t) st;

    #define MEGA 1.0e-6
    const double mega_samples = (MEGA * (double) totalSamples) / total_elapsed_time;

    printf("Msamples/s: %lf\n", mega_samples);


    // cleanup the remaining handles: the last two levels and the slabs of the last level, which nobody read
    for(size_t k = 0; k < g.num_blocks.z; k++)
    for(size_t j = 0; j < g.num_blocks.y; j++)
    for(size_t i = 0; i < g.num_blocks.x; i++){
        const size_t block = block_idx(i, j, k);
        starpu_data_unregister(p_wave_iter[1][block]);
        starpu_data_unregister(p_wave_iter[2][block]);
        starpu_data_unregister(q_wave_iter[1][block]);
        starpu_data_unregister(q_wave_iter[2][block]);
        const uint32_t neighbor_mask = halo_neighbor_mask(i, j, k);
        for(int table_index = 0; table_index < NEIGHBOR_COUNT; table_index++){
            if(!(neighbor_mask & NBIT(table_index))) continue;
            starpu_data_unregister(slabs[0][block * NEIGHBOR_COUNT + table_index]);
        }
    }

    starpu_data_unregister(scratch_handle);

    for(size_t idx = 0; idx < g.total_blocks; idx++){
        starpu_data_unregister(hdl_medium[idx]);
    }

    //Write the header file, only if there is a binary file for it to describe
    if(enable_io){
        FILE *rsf_header = NULL;
        char filename[512] = {'\0'};
        sprintf(filename, "%s/out-%s.rsf", output_folder, output_filename);

        TRY((rsf_header = fopen(filename, "w+")) == NULL);

        fprintf(rsf_header, "in=\"%s\"\n", io_state_get_filename());
        fprintf(rsf_header, "data_format=\"native_float\"\n");
        fprintf(rsf_header, "esize=%lu\n", sizeof(float));
        fprintf(rsf_header, "n1=%u\n", g.volume_dims.x /*sx*/);
        fprintf(rsf_header, "n2=%u\n", g.volume_dims.y /*sy*/);
        fprintf(rsf_header, "n3=%u\n", g.volume_dims.z /*sz*/);
        fprintf(rsf_header, "n4=%ld\n", n_out);
        fprintf(rsf_header, "d1=%f\n", dx);
        fprintf(rsf_header, "d2=%f\n", dy);
        fprintf(rsf_header, "d3=%f\n", dz);
        // the time between two saved moments
        fprintf(rsf_header, "d4=%f\n", steps_per_output * dt);
        fprintf(rsf_header, "seg=%u\n", g.num_blocks.x);

        fclose(rsf_header);
    }


    program_end:

    // tasks may still be running if a TRY failed after they were submitted, wait for them before freeing their buffers
    if(starpu_is_initialized()){
        starpu_task_wait_for_all();
    }
    mem_free(&medium_allocs);
    mem_free(&static_allocs);
    starpu_shutdown();
    assert(io_state_finish() == 0);
    return program_status;
}

