from functools import partial

import os

def tab(text, amount):
    return '\n'.join([('\t' * amount) + t for t in text.splitlines()])


# Per axis: coordinate variable, name of the neighbor on the low / high side
# (same names as the cube_faces_t bits in partition.h) and how to get the stride
# of a block view along that axis.
AXES = {
    'x': {'coord': 'x', 'neg': 'left',  'pos': 'right',  'stride': lambda view: '1'},
    'y': {'coord': 'y', 'neg': 'top',   'pos': 'bottom', 'stride': lambda view: f'{view}.ldy'},
    'z': {'coord': 'z', 'neg': 'front', 'pos': 'back',   'stride': lambda view: f'{view}.ldz'},
}

AXIS_PAIRS = [('x', 'y'), ('y', 'z'), ('x', 'z')]

SIDES = ['neg', 'pos', 'center']


def view_name(sides):
    """sides: {axis: 'neg' | 'pos'}; returns e.g. 'left_top', or 'center' if empty."""
    parts = [AXES[a][s] for a, s in sorted(sides.items()) if s != 'center']
    return '_'.join(parts) if parts else 'center'


def neighbor_offset(sides):
    """NEIGHBOR_IDX arguments for the view on the given sides."""
    delta = {'neg': -1, 'pos': 1, 'center': 0}
    return ', '.join(str(delta[sides.get(a, 'center')]) for a in 'xyz')


def linear_idx(view, coords):
    """x + y * view.ldy + z * view.ldz, dropping terms whose coordinate is 0."""
    terms = []
    for axis, stride in (('x', None), ('y', f'{view}.ldy'), ('z', f'{view}.ldz')):
        c = coords[axis]
        if c == '0':
            continue
        if stride is None:
            terms.append(c)
        elif ' ' in c:
            terms.append(f'({c}) * {stride}')
        else:
            terms.append(f'{c} * {stride}')
    return ' + '.join(terms) if terms else '0'


def border_coords(sides):
    """Coordinates inside the neighbor view that touch the central cube."""
    coords = {a: AXES[a]['coord'] for a in 'xyz'}
    for axis, side in sides.items():
        if side == 'neg':
            coords[axis] = 'STENCIL_RADIUS - 1'
        elif side == 'pos':
            coords[axis] = '0'
    return coords


class CrossCase:
    """One (axis pair, side of axis 1, side of axis 2) cross derivative variant."""

    def __init__(self, a1, a2, w1, w2):
        self.a1, self.a2, self.w1, self.w2 = a1, a2, w1, w2
        # view touched by each kind of tap, as named by compute_case()
        self.views = {'block': 'center'}
        if w1 != 'center':
            self.views[f'block_{w1}_center'] = view_name({a1: w1})
        if w2 != 'center':
            self.views[f'block_center_{w2}'] = view_name({a2: w2})
        if w1 != 'center' and w2 != 'center':
            self.views[f'block_{w1}_{w2}'] = view_name({a1: w1, a2: w2})

    @property
    def name(self):
        return f'cross_deriv_{self.a1}{self.a2}_{self.w1}_{self.w2}'

    def view_decls(self):
        lines = []
        sides_of = {'center': {}}
        if self.w1 != 'center':
            sides_of[view_name({self.a1: self.w1})] = {self.a1: self.w1}
        if self.w2 != 'center':
            sides_of[view_name({self.a2: self.w2})] = {self.a2: self.w2}
        if self.w1 != 'center' and self.w2 != 'center':
            sides_of[view_name({self.a1: self.w1, self.a2: self.w2})] = {self.a1: self.w1, self.a2: self.w2}
        for view, sides in sides_of.items():
            lines.append(f'const block_view_t {view} = neighborhood[NEIGHBOR_IDX({neighbor_offset(sides)})];')
        return '\n'.join(lines)

    def depth_calc(self):
        out = [f'const int32_t base_idx = {linear_idx("center", border_coords({}))};']
        for n, (axis, side) in enumerate(((self.a1, self.w1), (self.a2, self.w2)), start=1):
            coord = AXES[axis]['coord']
            if side == 'neg':
                out.append(f'const int32_t depth{n} = {coord};')
            elif side == 'pos':
                out.append(f'const int32_t depth{n} = cube_width - 1 - {coord};')
            if side != 'center':
                view = view_name({axis: side})
                out.append(f'const int32_t border_{n} = {linear_idx(view, border_coords({axis: side}))};')
        if self.w1 != 'center' and self.w2 != 'center':
            sides = {self.a1: self.w1, self.a2: self.w2}
            out.append(f'const int32_t border_3 = {linear_idx(view_name(sides), border_coords(sides))};')
        return '\n'.join(out)

    def tap(self, x, y, d1, d2, lamb):
        """Access to the tap at offset x along axis 1 and y along axis 2."""
        block, idx, n_1, n_2 = lamb(x, y, d1, d2)
        view = self.views[block]
        index = idx
        for n, axis in ((n_2, self.a2), (n_1, self.a1)):
            if n == 0:
                continue
            stride = AXES[axis]['stride'](view)
            step = f'{abs(n)}' if stride == '1' else f'{abs(n)} * {stride}'
            index += f' {symbol(n)} {step}'
        return f'{view}.ptr[{index}]'


def is_case(val, w):
    match w:
        case 'pos':
            return val > 0
        case 'neg':
            return val < 0
        case 'center':
            return False


def proper_dir(val1, endval):
    if val1 < 0:
        return -endval
    else:
        return endval


# por algum motivo essa é a ordem certa, y,x
def compute_case(y, x, depth_x, depth_y, case1, case2):
    dist_x = abs(x) - depth_x
    dist_y = abs(y) - depth_y
    change_x = is_case(x, case1) and dist_x > 0  # nums go from -4 to 4
    change_y = is_case(y, case2) and dist_y > 0
    match (change_x, change_y):
        case (False, False):
            return ('block', 'base_idx', x, y)
        case (True, False):
            return (f'block_{case1}_center', 'border_1', proper_dir(x, dist_x - 1), y)
        case (False, True):
            return (f'block_center_{case2}', 'border_2', x, proper_dir(y, dist_y - 1))
        case (True, True):
            return (f"block_{case1}_{case2}", 'border_3', proper_dir(x, dist_x - 1), proper_dir(y, dist_y - 1))


def function_impl(case, lamb):
    tap = partial(case.tap, lamb=lamb)
    out = ''
    if case.w1 == 'center':
        if case.w2 == 'center':
            out += tab(make_cross_comp(4, 4, tap), 1)
        else:
            out += 'switch (depth2){'
            for d2 in range(0, 4):
                out += f'''
    case {d2}:
{tab(make_cross_comp(4, d2, tap), 2)}'''
            out += '\n\tdefault: UNREACHABLE;\n\t}'
    else:
        if case.w2 == 'center':
            out += '''switch (depth1){'''
            for d1 in range(0, 4):
                out += f'''
    case {d1}:
{tab(make_cross_comp(d1, 4, tap),2)}'''
            out += '\n\tdefault: UNREACHABLE;\n\t}'
        else:
            out += function_impl_full(tap)

    return out


def function_impl_full(tap):
    out = '''switch (depth1 * STENCIL_RADIUS + depth2){'''
    for d1 in range(0, 4):
        for d2 in range(0, 4):
            out += f'''
    case {d1} * STENCIL_RADIUS + {d2}:
{tab(make_cross_comp(d1, d2, tap), 2)}'''
    out += '\n\tdefault: UNREACHABLE; \n\t}'
    return out


symbol = lambda x: '-' if x < 0 else '+'


def debug(p):
    return f'''
printf("Computed operation:\\n %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f\\n    ) +\\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f + \\n        %.9f - %.9f - %.9f + %.9f\\n    ) +\\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f +\\n        %.9f - %.9f - %.9f + %.9f\\n    ) +\\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f +\\n        %.9f - %.9f - %.9f + %.9f\\n    ) +\\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f\\n    ) +\\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f + \\n        %.9f - %.9f - %.9f + %.9f\\n	) +        \\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f + \\n        %.9f - %.9f - %.9f + %.9f\\n	) +      \\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f\\n\\n	) + \\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f + \\n        %.9f - %.9f - %.9f + %.9f\\n	) + \\n    %.9f *  (\\n        %.9f - %.9f - %.9f + %.9f\\n    )) * %.9f\\n",
    L11, {p(+1, +1)}, {p(+1, -1)}, {p(-1, +1)}, {p(-1, -1)} ,
    L12, {p(+1, +2)}, {p(+1, -2)}, {p(-1, +2)}, {p(-1, -2)},
        {p(+2, +1)}, {p(+2, -1)}, {p(-2, +1)}, {p(-2, -1)} ,
    L13, {p(+1, +3)}, {p(+1, -3)}, {p(-1, +3)}, {p(-1, -3)},
        {p(+3, +1)}, {p(+3, -1)}, {p(-3, +1)}, {p(-3, -1)} ,
    L14, {p(+1, +4)}, {p(+1, -4)}, {p(-1, +4)}, {p(-1, -4)},
        {p(+4, +1)}, {p(+4, -1)}, {p(-4, +1)}, {p(-4, -1)} ,
    L22, {p(+2, +2)}, {p(+2, -2)}, {p(-2, +2)}, {p(-2, -2)} ,
    L23, {p(+2, +3)}, {p(+2, -3)}, {p(-2, +3)}, {p(-2, -3)},
        {p(+3, +2)}, {p(+3, -2)}, {p(-3, +2)}, {p(-3, -2)} ,
    L24, {p(+2, +4)}, {p(+2, -4)}, {p(-2, +4)}, {p(-2, -4)},
        {p(+4, +2)}, {p(+4, -2)}, {p(-4, +2)}, {p(-4, -2)} ,
    L33, {p(+3, +3)}, {p(+3, -3)}, {p(-3, +3)}, {p(-3, -3)} ,
    L34, {p(+3, +4)}, {p(+3, -4)}, {p(-3, +4)}, {p(-3, -4)},
        {p(+4, +3)}, {p(+4, -3)}, {p(-4, +3)}, {p(-4, -3)},
    L44, {p(+4, +4)}, {p(+4, -4)}, {p(-4, +4)}, {p(-4, -4)}, dinv
    );
'''

# {debug(p)}
def make_cross_comp(d1, d2, tap):
    p = partial(tap, d1 = d1, d2 = d2)
    return f'''
{debug(p) if os.environ.get('PRINTOUT') is not None else ""}
return ((
    L11 * (
        {p(+1, +1)} - {p(+1, -1)} - {p(-1, +1)} + {p(-1, -1)}
    ) +
    L12 * (
        {p(+1, +2)} - {p(+1, -2)} - {p(-1, +2)} + {p(-1, -2)} +
        {p(+2, +1)} - {p(+2, -1)} - {p(-2, +1)} + {p(-2, -1)}
    ) +
    L13 * (
        {p(+1, +3)} - {p(+1, -3)} - {p(-1, +3)} + {p(-1, -3)} +
        {p(+3, +1)} - {p(+3, -1)} - {p(-3, +1)} + {p(-3, -1)}
    ) +
    L14 * (
        {p(+1, +4)} - {p(+1, -4)} - {p(-1, +4)} + {p(-1, -4)} +
        {p(+4, +1)} - {p(+4, -1)} - {p(-4, +1)} + {p(-4, -1)}
    ) +
    L22 * (
        {p(+2, +2)} - {p(+2, -2)} - {p(-2, +2)} + {p(-2, -2)}
    ) +
    L23 * (
        {p(+2, +3)} - {p(+2, -3)} - {p(-2, +3)} + {p(-2, -3)} +
        {p(+3, +2)} - {p(+3, -2)} - {p(-3, +2)} + {p(-3, -2)}
	) +
    L24 * (
        {p(+2, +4)} - {p(+2, -4)} - {p(-2, +4)} + {p(-2, -4)} +
        {p(+4, +2)} - {p(+4, -2)} - {p(-4, +2)} + {p(-4, -2)}
	) +
    L33 * (
        {p(+3, +3)} - {p(+3, -3)} - {p(-3, +3)} + {p(-3, -3)}

	) +
    L34 * (
        {p(+3, +4)} - {p(+3, -4)} - {p(-3, +4)} + {p(-3, -4)} +
        {p(+4, +3)} - {p(+4, -3)} - {p(-4, +3)} + {p(-4, -3)}
	) +
    L44 * (
        {p(+4, +4)} - {p(+4, -4)} - {p(-4, +4)} + {p(-4, -4)}
    )) * dinv);
'''


SIGNATURE = '''const block_view_t *neighborhood,
    int32_t x, int32_t y, int32_t z,
    int32_t cube_width, FP dinv'''


def dispatcher(a1, a2):
    """cross_deriv_<a1><a2>: picks the variant for the sides the point touches."""
    c1, c2 = AXES[a1]['coord'], AXES[a2]['coord']
    side_const = {'neg': 'SIDE_NEG', 'pos': 'SIDE_POS', 'center': 'SIDE_CENTER'}
    cases = []
    for w1 in SIDES:
        for w2 in SIDES:
            if w1 == 'center' and w2 == 'center':
                continue
            cases.append(f'    case {side_const[w1]} * 3 + {side_const[w2]}: '
                         f'return {CrossCase(a1, a2, w1, w2).name}(neighborhood, x, y, z, cube_width, dinv);')
    cases.append(f'    case SIDE_CENTER * 3 + SIDE_CENTER:\n'
                 f'    default: return {CrossCase(a1, a2, "center", "center").name}(neighborhood, x, y, z, cube_width, dinv);')
    body = '\n'.join(cases)
    return f'''
ATTRIBUTE FP cross_deriv_{a1}{a2}(
    {SIGNATURE}
){{
    switch (stencil_side({c1}, cube_width) * 3 + stencil_side({c2}, cube_width))
    {{
{body}
    }}
}}'''


def main():
    print(f'''/*
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
''')
    for a1, a2 in AXIS_PAIRS:
        for w1 in SIDES:
            for w2 in SIDES:
                case = CrossCase(a1, a2, w1, w2)
                lamb = partial(compute_case, case1=w1, case2=w2)
                print(f'''
ATTRIBUTE static FP {case.name}(
    {SIGNATURE}
){{
    (void) cube_width;
{tab(case.view_decls(), 1)}
{tab(case.depth_calc(), 1)}
    {function_impl(case, lamb)}
}}''')
        print(dispatcher(a1, a2))

main()
