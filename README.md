# dynsys

A GLFW/OpenGL visualizer for continuous-time 3D dynamical systems.

This version uses:

- Dear ImGui for the control UI
- TPCAS for equation parsing
- normal TPCAS/infix equations instead of Lisp/Polish notation

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

## Equation syntax

Equations are edited inside the Dear ImGui control panel. Use `x`, `y`, and `z` for the current point, and write ordinary infix arithmetic:

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

## UI

The control panel supports:

- preset switching between Lorenz and Rössler
- editing `dx/dt`, `dy/dt`, and `dz/dt` live
- parse/runtime error display without crashing the app
- start point, `dt`, point-buffer, and speed controls
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

The old file `src/dynsys.c` is not part of the Dear ImGui build. The active entry point is `src/dynsys.cpp`.

