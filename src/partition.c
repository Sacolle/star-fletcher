#include <stdint.h>
#include <starpu.h>

#include "partition.h"

extern const size_t BORDER_WIDTH;

// Canonical child order (index = StarPU child id). See partition.h.
static const cube_faces_t part_idx_to_cube_face[CUBE_FACE_NPARTS] = {
    CFACE_LEFT, CFACE_RIGHT, CFACE_TOP, CFACE_BOTTOM, CFACE_FRONT, CFACE_BACK,
    CFACE_LEFT | CFACE_TOP,     CFACE_LEFT | CFACE_BOTTOM,
    CFACE_RIGHT | CFACE_TOP,     CFACE_RIGHT | CFACE_BOTTOM,
    CFACE_LEFT | CFACE_FRONT,   CFACE_LEFT | CFACE_BACK,
    CFACE_RIGHT | CFACE_FRONT,   CFACE_RIGHT | CFACE_BACK,
    CFACE_TOP  | CFACE_FRONT,   CFACE_TOP  | CFACE_BACK,
    CFACE_BOTTOM | CFACE_FRONT, CFACE_BOTTOM | CFACE_BACK,
};

size_t cube_face_to_part_idx(cube_faces_t face)
{
    for (size_t i = 0; i < CUBE_FACE_NPARTS; i++)
        if (part_idx_to_cube_face[i] == face)
            return i;
    return SIZE_MAX;
}

typedef struct {
    uint32_t x_offset, x_size;
    uint32_t y_offset, y_size;
    uint32_t z_offset, z_size;
} face_extent_t;

static face_extent_t face_to_extent(cube_faces_t face, uint32_t nx, uint32_t ny, uint32_t nz)
{
    const uint32_t border_width = (uint32_t) BORDER_WIDTH;
    face_extent_t e = { 0, nx, 0, ny, 0, nz };

    if (face & CFACE_LEFT)   { e.x_offset = 0;                 e.x_size = border_width; }
    if (face & CFACE_RIGHT)   { e.x_offset = nx - border_width; e.x_size = border_width; }
    if (face & CFACE_TOP)    { e.y_offset = 0;                 e.y_size = border_width; }
    if (face & CFACE_BOTTOM) { e.y_offset = ny - border_width; e.y_size = border_width; }
    if (face & CFACE_FRONT)  { e.z_offset = 0;                 e.z_size = border_width; }
    if (face & CFACE_BACK)  { e.z_offset = nz - border_width; e.z_size = border_width; }
    return e;
}

// Same conventions as StarPU's starpu_block_filter_block: sizes are always set,
// pointer/strides/offset only when the parent is allocated (dev_handle != 0).
static void cube_face_filter_func(void *parent_interface, void *child_interface,
                                  struct starpu_data_filter *f, unsigned id, unsigned nparts)
{
    (void) f;
    struct starpu_block_interface *father = parent_interface;
    struct starpu_block_interface *child = child_interface;

    STARPU_ASSERT_MSG(father->id == STARPU_BLOCK_INTERFACE_ID, "cube_face_filter only applies to block data");
    STARPU_ASSERT_MSG(nparts == CUBE_FACE_NPARTS, "cube_face_filter expects %d parts, got %u", CUBE_FACE_NPARTS, nparts);
    STARPU_ASSERT_MSG(father->nx >= BORDER_WIDTH && father->ny >= BORDER_WIDTH && father->nz >= BORDER_WIDTH,
                      "cube (%u, %u, %u) is thinner than BORDER_WIDTH", father->nx, father->ny, father->nz);

    const face_extent_t e = face_to_extent(part_idx_to_cube_face[id], father->nx, father->ny, father->nz);

    child->id       = father->id;
    child->elemsize = father->elemsize;
    child->nx       = e.x_size;
    child->ny       = e.y_size;
    child->nz       = e.z_size;

    if (father->dev_handle) {
        // sub-block is not contiguous: it keeps the parent's strides
        const size_t off = ((size_t) e.x_offset
                          + (size_t) e.y_offset * father->ldy
                          + (size_t) e.z_offset * father->ldz) * father->elemsize;
        if (father->ptr)
            child->ptr = father->ptr + off;
        child->ldy        = father->ldy;
        child->ldz        = father->ldz;
        child->dev_handle = father->dev_handle;
        child->offset     = father->offset + off;
    }
}

struct starpu_data_filter g_cube_face_filter = {
    .filter_func = cube_face_filter_func,
    .nchildren   = CUBE_FACE_NPARTS,
};
