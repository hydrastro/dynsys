# dynsys

A C/OpenGL visualizer for continuous-time 3D dynamical systems.

This version uses the vendored TPCAS parser directly, so equations are written in normal TPCAS/infix mode instead of the old Lisp/Polish notation.

## Build

```sh
make
```

With Nix:

```sh
nix build
nix run
```

The executable is `build/dynsys` when using Make.

## Equation syntax

Use `x`, `y`, and `z` for the current point, and write ordinary infix arithmetic:

```text
Lorenz:
dx/dt = 10 * (y - x)
dy/dt = x * (28 - z) - y
dz/dt = x * y - (8 / 3) * z

Rössler:
dx/dt = 0 - (y + z)
dy/dt = x + (1 / 5) * y
dz/dt = (1 / 5) + z * (x - (57 / 10))
```

Unary minus and decimals are supported, for example `-x + 1.5`. The evaluator also accepts `sin`, `cos`, `tan`, `exp`, `log`, `sqrt`, `abs`, `pow`, `min`, and `max` as function calls.

## Runtime font

The uploaded bitmap/font file is not vendored here. By default, the Make build tries `/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf`. Override it with:

```sh
DYNSYS_FONT_PATH=/path/to/font.ttf ./build/dynsys
```

or at compile time:

```sh
make FONT_PATH=/path/to/font.ttf
```

## Controls

- `w`, `a`, `s`, `d`: rotate camera
- `-`, `=`: zoom
- `z`: toggle cross
- `x`: toggle axes
- `c`: clear/reset
- `spacebar`: pause/resume
- `m`, `,`: decrease speed
- `.`, `/`: increase speed
- `j`, `l`: x translation
- `i`, `k`: y translation
- `p`, `;`: z translation
- `v`, `b`: adjust `dt`
