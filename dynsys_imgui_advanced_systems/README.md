# dynsys

A GLFW/OpenGL visualizer for continuous-time 3D dynamical systems.

This version uses:

- Dear ImGui for the control UI
- TPCAS for equation parsing
- normal TPCAS/infix equations instead of Lisp/Polish notation
- named auxiliary definitions and user-defined helper functions

## Build

Preferred workflow:

```sh
nix develop
make
```

Run it with:

```sh
make run
```

Useful variants:

```sh
make release
make asan
make clean
make print-vars
```

The executable is `build/dynsys`.

The flake intentionally provides a development shell, not a `nix build`/`nix run` application. `make` owns the build.

## Dear ImGui source

`nix develop` exports `IMGUI_DIR` to the pinned Dear ImGui source input declared in `flake.nix`. The Makefile compiles the Dear ImGui core sources plus the GLFW/OpenGL3 backends directly into dynsys.

Outside Nix, provide your own Dear ImGui checkout:

```sh
make IMGUI_DIR=/path/to/imgui
```

The expected source layout is:

```text
/path/to/imgui/imgui.cpp
/path/to/imgui/imgui_draw.cpp
/path/to/imgui/imgui_tables.cpp
/path/to/imgui/imgui_widgets.cpp
/path/to/imgui/backends/imgui_impl_glfw.cpp
/path/to/imgui/backends/imgui_impl_opengl3.cpp
```

## System syntax

The GUI now has one multiline system editor instead of three isolated derivative fields. Each non-empty line is a definition:

```text
name = expression
name(arg1, arg2) = expression
dx = expression
dy = expression
dz = expression
```

`dx`, `dy`, and `dz` are required. `dx/dt`, `dy/dt`, and `dz/dt` are also accepted aliases.

Comments start with `#` or `//`.

Available state variables:

```text
x, y, z, t
```

Built-in constants:

```text
pi, e
```

Built-in functions:

```text
sin, cos, tan, exp, log, sqrt, abs, pow, min, max
```

Example Lorenz system using named parameters:

```text
# Lorenz system
sigma = 10
rho = 28
beta = 8 / 3

dx = sigma * (y - x)
dy = x * (rho - z) - y
dz = x * y - beta * z
```

Example with user-defined helper functions:

```text
# Thomas cyclically symmetric attractor
b = 0.208186
wave(u) = sin(u)

dx = wave(y) - b * x
dy = wave(z) - b * y
dz = wave(x) - b * z
```

Example with explicit time dependence:

```text
# Forced Duffing-like 3D autonomous form
# z is still rendered; t is internal simulation time.
delta = 0.2
alpha = -1
beta = 1
gamma = 0.3
omega = 1.2

drive(tau) = gamma * cos(omega * tau)

dx = y
dy = 0 - delta * y - alpha * x - beta * pow(x, 3) + drive(t)
dz = 1
```

The evaluator rejects cyclic scalar definitions and recursive helper functions. Input is TPCAS infix mode; Lisp/S-expression syntax is not used.

## UI

The control panel supports:

- preset switching between Lorenz, Rössler, Thomas, and Duffing-like systems
- editing the full system live
- named parameters and auxiliary expressions
- helper functions such as `wave(u) = sin(u)`
- parse/runtime error display without crashing the app
- start point, start time, `dt`, point-buffer, and speed controls
- camera rotation, zoom, translation, axes, and crosshair controls

## Controls

The GUI is the primary control surface. Keyboard/mouse shortcuts are still available when ImGui is not capturing input:

- mouse drag: rotate camera
- mouse wheel: zoom
- `w`, `a`, `s`, `d`: rotate camera
- `-`, `=`: zoom
- `z`: toggle crosshair
- `x`: toggle axes
- `c`: clear/reset
- `spacebar`: pause/resume
- `m`, `,`: decrease speed
- `.`, `/`: increase speed
- `j`, `l`: x translation
- `i`, `k`: y translation
- `p`, `;`: z translation
- `v`, `b`: adjust `dt`
- `Esc`: quit

## Updating an existing dirty checkout

If you copied this tree over an older dynsys checkout, remove the legacy FreeType/C entry point before building:

```sh
make prune-legacy
make prune-legacy CLEAN_APPLY=1
make clean
make
```

The old file `src/dynsys.c` is not part of this package. The active entry point is `src/dynsys.cpp`.
