#include "io.h"

#define _GNU_SOURCE
#include <errno.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#include <starpu.h>

#include "vector.h"

#define HAS_FLAG(x, f) !(~(x) & f)
#define ROUND_UP(x, a) ((((x) + (a) - 1) / (a)) * (a))
// used when the kernel or the filesystem does not report the O_DIRECT alignment
#define FALLBACK_DIO_ALIGN 4096

#define MAX_PATH 512

// where the block (i, j, k) of moment t is in the file of its worker
typedef struct io_block_entry {
    size_t i, j, k, t;
    size_t offset;
    size_t bytes;
} io_block_entry_t;

// only touched by its own worker, so no locks are needed
typedef struct io_worker_file {
    file_desc_t fd; // -1 for the workers that do not write (non CPU)
    char filename[MAX_PATH]; // relative to the folder of the index
    size_t offset;
    // O_DIRECT requires aligned address, offset and size, which the starpu blocks
    // do not have, so they are copied to this buffer and padded before the write
    void* bounce;
    size_t bounce_size;
    vector(io_block_entry_t) entries;
    err_t err; // first error on write, the kernel can't return it
} io_worker_file_t;

struct io_state {
    bool initialized;
    char filename[MAX_PATH];
    char meta_path[MAX_PATH];
    size_t volume_width;
    size_t align;
    unsigned nworkers;
    io_worker_file_t* workers;
};

io_state_t _state = { .initialized = false };

io_state_t* get_io_state(){
    return &_state;
}
void delete_io_state(){
    io_state_t* state = get_io_state();
    for(unsigned w = 0; state->workers && w < state->nworkers; w++){
        free(state->workers[w].bounce);
        vector_free(state->workers[w].entries);
    }
    free(state->workers);
    *state = (io_state_t){ .initialized = false };
}

bool is_io_state_initialized(){
    return get_io_state()->initialized;
}

// get the alignment the system requires for O_DIRECT in the file
static err_t io_alignment_restrictions(file_desc_t fd, size_t* alignment){
#ifdef STATX_DIOALIGN
    struct statx stx;
    if(statx(fd, "", AT_EMPTY_PATH, STATX_DIOALIGN, &stx) != 0){
        return errno;
    }
    if(!HAS_FLAG(stx.stx_mask, STATX_DIOALIGN) || stx.stx_dio_mem_align == 0){
        return ME_MISSING_ATTR;
    }
    const size_t mem = stx.stx_dio_mem_align, off = stx.stx_dio_offset_align;
    *alignment = mem > off ? mem : off;
    return 0;
#else
    (void) fd; (void) alignment;
    return ME_MISSING_ATTR;
#endif
}

err_t io_state_init(const char* folder, const char* basename, const size_t volume_width){
    if(is_io_state_initialized()){
        return ME_REINITILIZATION;
    }

    io_state_t* state = get_io_state();
    *state = (io_state_t){
        .volume_width = volume_width,
        .align = 0,
        .nworkers = starpu_worker_get_count()
    };
    snprintf(state->filename, MAX_PATH, "%s/out-%s.rsf@", folder, basename);
    snprintf(state->meta_path, MAX_PATH, "%s/out-%s.meta", folder, basename);

    if((state->workers = calloc(state->nworkers, sizeof(io_worker_file_t))) == NULL){
        return errno;
    }

    err_t err = 0;
    for(unsigned w = 0; w < state->nworkers; w++){
        io_worker_file_t* wf = &state->workers[w];
        wf->fd = -1;
        // the dump codelet only has a CPU implementation
        if(starpu_worker_get_type(w) != STARPU_CPU_WORKER){
            continue;
        }
        snprintf(wf->filename, MAX_PATH, "out-%s-w%u.bin", basename, w);

        char path[2 * MAX_PATH];
        snprintf(path, sizeof(path), "%s/%s", folder, wf->filename);
        // O_DIRECT: write direct to disk without the page cache, requires aligned writes
        // O_DSYNC is not used, the files are synced once on finish
        if((wf->fd = open(path, O_WRONLY | O_CREAT | O_TRUNC | O_DIRECT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH)) == -1){
            err = errno;
            goto init_fail;
        }

        if(state->align == 0){
            if(io_alignment_restrictions(wf->fd, &state->align) != 0){
                printf("[warning] Could not get the O_DIRECT alignment of %s, using %d\n", path, FALLBACK_DIO_ALIGN);
                state->align = FALLBACK_DIO_ALIGN;
            }
        }
    }

    state->initialized = true;
    return 0;

init_fail:
    for(unsigned w = 0; w < state->nworkers; w++){
        if(state->workers[w].fd != -1){
            close(state->workers[w].fd);
        }
    }
    delete_io_state();
    return err;
}

void io_state_write_file(
    const size_t i, const size_t j, const size_t k, const size_t t, 
    const size_t nx, const size_t ny, const size_t nz, const FP* block
){
    io_state_t* state = get_io_state();
    io_worker_file_t* wf = &state->workers[starpu_worker_get_id()];

    if(wf->err != 0){
        return;
    }

    const size_t bytes = nx * ny * nz * sizeof(FP);
    const size_t padded = ROUND_UP(bytes, state->align);

    if(wf->bounce_size < padded){
        free(wf->bounce);
        wf->bounce_size = 0;
        if((wf->err = posix_memalign(&wf->bounce, state->align, padded)) != 0){
            wf->bounce = NULL;
            return;
        }
        wf->bounce_size = padded;
    }
    memcpy(wf->bounce, block, bytes);
    memset((char*) wf->bounce + bytes, 0, padded - bytes);

    size_t written = 0;
    while(written < padded){
        const ssize_t res = pwrite(wf->fd, (char*) wf->bounce + written, padded - written, wf->offset + written);
        if(res == -1){
            if(errno == EINTR){
                continue;
            }
            wf->err = errno;
            return;
        }
        written += res;
    }

    vector_push(wf->entries, ((io_block_entry_t){
        .i = i, .j = j, .k = k, .t = t,
        .offset = wf->offset,
        .bytes = bytes
    }));
    wf->offset += padded;
}

err_t io_state_finish(){
    if(!is_io_state_initialized()){
        return 0;
    }

    io_state_t* state = get_io_state();
    err_t err = 0;

    for(unsigned w = 0; w < state->nworkers; w++){
        io_worker_file_t* wf = &state->workers[w];
        if(wf->fd == -1){
            continue;
        }
        if(err == 0 && wf->err != 0){
            err = wf->err;
        }
        if(fsync(wf->fd) != 0 && err == 0){
            err = errno;
        }
        if(close(wf->fd) != 0 && err == 0){
            err = errno;
        }
        wf->fd = -1;
    }

    FILE* meta = NULL;
    if((meta = fopen(state->meta_path, "w")) == NULL){
        if(err == 0){
            err = errno;
        }
        delete_io_state();
        return err;
    }

    // the output is relative to the folder of the index
    const char* output = strrchr(state->filename, '/');
    output = output ? output + 1 : state->filename;

    fprintf(meta, "# star-fletcher io index v1\n");
    fprintf(meta, "volume_width %zu\n", state->volume_width);
    fprintf(meta, "esize %zu\n", sizeof(FP));
    fprintf(meta, "align %zu\n", state->align);
    fprintf(meta, "output %s\n", output);
    fprintf(meta, "# t i j k file offset bytes\n");
    for(unsigned w = 0; w < state->nworkers; w++){
        const io_worker_file_t* wf = &state->workers[w];
        for(long long e = 0; e < vector_count(wf->entries); e++){
            const io_block_entry_t* b = &wf->entries[e];
            fprintf(meta, "%zu %zu %zu %zu %s %zu %zu\n", b->t, b->i, b->j, b->k, wf->filename, b->offset, b->bytes);
        }
    }

    if(fclose(meta) != 0 && err == 0){
        err = errno;
    }

    delete_io_state();
    return err;
}

const char* io_state_get_filename(){
    return get_io_state()->filename;
}

err_t io_available_huge_page_sizes(size_t max_count, size_t* outp_respc, size_t* outp_available_sizes){
    struct dirent* dent;
    DIR* srcdir;
    *outp_respc = 0;

    if ((srcdir = opendir("/sys/kernel/mm/hugepages")) == NULL){
        return errno;
    }

    while((dent = readdir(srcdir)) != NULL){

        const char* folder = dent->d_name;
        const char* prepend = "hugepages-";
        const size_t prepend_size = strlen(prepend);
        if(strlen(folder) < prepend_size){
            continue;
        }
        if(strncmp(folder, prepend, prepend_size) == 0){
            char* out = "\0\0";

            const size_t res = strtoull(folder + prepend_size, &out, 10);
            if(errno){
                return errno;
            }

            if(strncmp(out, "kB", 2) != 0){
                return ME_INCOMPLETE_PARSE;
            }
            if(*outp_respc < max_count){
                outp_available_sizes[*outp_respc]  = res * 1024;
            }
            *outp_respc += 1;
        }
    }
    closedir(srcdir);
    return 0;
}
