#ifndef GUARD_MACROS
#define GUARD_MACROS

#include <stddef.h>

// linear idx https://stackoverflow.com/a/34363187
// calculates a + width_a * (b + c * width_b)
static inline size_t idx(size_t a, size_t b, size_t c, size_t width_a, size_t width_b){
    return a + width_a * (b + c * width_b);
}

struct coords3d {
    size_t x;
    size_t y;
    size_t z;
};

// extract the components
// https://stackoverflow.com/a/11712864
// stackoverflow misturou o x com o z no meu contexto, pois o x é o eixo menor
static inline struct coords3d unindex(size_t index, size_t width_x, size_t width_y){
    struct coords3d c;

    c.x = index % width_x;
    c.y = (index / width_x) % width_y;
    c.z = index / (width_x * width_y);

    return c;
}

// the linear index of the point(x, y, z) inside a segmented block
// dependes on the global variables `g_cube_width_{x,y}`;
static inline size_t cube_idx(size_t x, size_t y, size_t z){
    extern size_t g_cube_width_x, g_cube_width_y;
    return idx(x, y, z, g_cube_width_x, g_cube_width_y);
}

// the linear index of the block(i,j,k) in the total volume
// dependes on the global variables `g_width_in_cubes_{x,y}`;
static inline size_t block_idx(size_t i, size_t j, size_t k){
    extern size_t g_width_in_cubes_x, g_width_in_cubes_y;
    return idx(i, j, k, g_width_in_cubes_x, g_width_in_cubes_y);
}

// the linear index of the block(i,j,k) in the padded handle arrays (one halo block per side)
// dependes on the global variables `g_width_in_cubes_{x,y}`;
static inline size_t padded_block_idx(size_t i, size_t j, size_t k){
    extern size_t g_width_in_cubes_x, g_width_in_cubes_y;
    return idx(i, j, k, g_width_in_cubes_x + 2, g_width_in_cubes_y + 2);
}

// the linear index for a point (vx, vy, vz) in the whole volume
// dependes on the global variables `g_volume_width_{x,y}`;
static inline size_t volume_idx(size_t vx, size_t vy, size_t vz){
    extern size_t g_volume_width_x, g_volume_width_y;
    return idx(vx, vy, vz, g_volume_width_x, g_volume_width_y);
}

#define SQUARE(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))

//if exp != 0, run body and goto tag
#define TRYTO(exp, body, tag) if((exp) != 0) { body; goto tag; }

// turn a (x, y, z) point in a (i, j, k) block in a (vx, vy, vz) coordenate for the complete volume
// depends on `g_cube_width_{x,y,z}` and `g_volume_width_{x,y}`
static inline size_t block_cube_to_volume_idx(size_t x, size_t y, size_t z, size_t i, size_t j, size_t k){
    extern size_t g_cube_width_x, g_cube_width_y, g_cube_width_z;
    return volume_idx(x + i * g_cube_width_x, y + j * g_cube_width_y, z + k * g_cube_width_z);
}

// the (vx, vy, vz) coordenates of a global volume linear idx
// depends on `g_volume_width_{x,y}`
static inline struct coords3d volume_unindex(size_t idx){
    extern size_t g_volume_width_x, g_volume_width_y;
    return unindex(idx, g_volume_width_x, g_volume_width_y);
}

// turn a global volume linear idx into a the index of the block that point correspondes
// depends on `g_cube_width_{x,y,z}`, `g_width_in_cubes_{x,y}` and `g_volume_width_{x,y}`
static inline size_t volume_to_block_idx(size_t idx){
    extern size_t g_cube_width_x, g_cube_width_y, g_cube_width_z;

    const struct coords3d c = volume_unindex(idx);
    return block_idx(c.x / g_cube_width_x, c.y / g_cube_width_y, c.z / g_cube_width_z);
}

// turn a global volume linear idx into a the index of the block that point correspondes,
// in the padded handle arrays (shifted by one halo block in every direction)
// depends on `g_cube_width_{x,y,z}`, `g_width_in_cubes_{x,y}` and `g_volume_width_{x,y}`
static inline size_t volume_to_padded_block_idx(size_t idx){
    extern size_t g_cube_width_x, g_cube_width_y, g_cube_width_z;

    const struct coords3d c = volume_unindex(idx);
    return padded_block_idx(c.x / g_cube_width_x + 1, c.y / g_cube_width_y + 1, c.z / g_cube_width_z + 1);
}

// turn a global volume linear index into a the index inside its corresponding block
// depends on `g_cube_width_{x,y,z}` and `g_volume_width_{x,y}`
static inline size_t volume_to_cube_idx(size_t idx){
    extern size_t g_cube_width_x, g_cube_width_y, g_cube_width_z;

    const struct coords3d c = volume_unindex(idx);
    return cube_idx(c.x % g_cube_width_x, c.y % g_cube_width_y, c.z % g_cube_width_z);
}

#endif
