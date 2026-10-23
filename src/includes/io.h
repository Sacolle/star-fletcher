#ifndef _IO_GUARD_
#define _IO_GUARD_

#define _GNU_SOURCE
#include <stddef.h>
#include <stdint.h>

#include "err.h"
#include "floatingpoint.h"

typedef int file_desc_t;

typedef struct io_state io_state_t;

// opens one file per CPU worker at `folder`/out-`basename`-w<id>.bin with O_DIRECT,
// so the blocks never pass through the page cache.
// must be called after starpu_init, as it needs the workers.
err_t io_state_init(const char* folder, const char* basename, const size_t volume_width);

// appends the block (i, j, k) of moment t to the file of the calling worker,
// padded to the O_DIRECT alignment, and records where it was written
void io_state_write_file(
    const size_t i, const size_t j, const size_t k, const size_t t,
    const size_t nx, const size_t ny, const size_t nz, const FP* block
);

// closes the worker files and writes the index `folder`/out-`basename`.meta,
// with "t i j k file offset bytes" per block. scripts/congeal.py joins them into the volume
err_t io_state_finish();

// the volume file that scripts/congeal.py builds from the worker files
const char* io_state_get_filename();

err_t io_available_huge_page_sizes(size_t max_count, size_t* outp_respc, size_t* outp_available_sizes);

#endif
