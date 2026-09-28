#ifndef GUARD_PARTITION
#define GUARD_PARTITION

#include <starpu.h>
#include <stdint.h>

// Número de filhos gerados pelo filtro: 6 faces + 12 arestas.
#define CUBE_FACE_NPARTS 18

// Face bitmask: 2 bits por eixo.
// bits 0-1: eixo X (LEFT  = x baixo, RIGHT  = x alto), x cresce da esquerda para a direita
// bits 2-3: eixo Y (TOP   = y baixo, BOTTOM = y alto), y cresce de cima para baixo
// bits 4-5: eixo Z (FRONT = z baixo, BACK   = z alto), z cresce da frente para trás
// Combinações com | de dois eixos representam arestas. Quinas (3 eixos) são inválidas.
typedef uint8_t cube_faces_t;
enum {
    CFACE_LEFT   = 0x01,  // x = [0 .. BORDER_WIDTH)
    CFACE_RIGHT  = 0x02,  // x = [w - BORDER_WIDTH .. w)
    CFACE_TOP    = 0x04,  // y = [0 .. BORDER_WIDTH)
    CFACE_BOTTOM = 0x08,  // y = [w - BORDER_WIDTH .. w)
    CFACE_FRONT  = 0x10,  // z = [0 .. BORDER_WIDTH)
    CFACE_BACK   = 0x20,  // z = [w - BORDER_WIDTH .. w)
};

// Ordem canônica dos 18 filhos StarPU (índice do filho em partition_plan):
//   0:LEFT  1:RIGHT  2:TOP  3:BOTTOM  4:FRONT  5:BACK
//   6:LEFT|TOP     7:LEFT|BOTTOM    8:RIGHT|TOP     9:RIGHT|BOTTOM
//  10:LEFT|FRONT  11:LEFT|BACK     12:RIGHT|FRONT  13:RIGHT|BACK
//  14:TOP|FRONT   15:TOP|BACK      16:BOTTOM|FRONT 17:BOTTOM|BACK

// Filtro StarPU que particiona um handle de bloco em CUBE_FACE_NPARTS sub-handles.
// Os filhos são fatias de espessura BORDER_WIDTH, não contíguas (mantêm ldy/ldz do pai).
extern struct starpu_data_filter g_cube_face_filter;

// Converte uma face/aresta (cube_faces_t) no índice do filho 0..17.
// Retorna SIZE_MAX se a combinação não é uma face ou aresta válida.
size_t cube_face_to_part_idx(cube_faces_t face);

#endif
