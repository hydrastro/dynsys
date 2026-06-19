# dynsys

An interactive visualizer and analyzer for dynamical systems — phase planes,
attractors, basins of attraction, and continuation/bifurcation analysis aimed at
MatCont-style parity — built in C++17 on Dear ImGui + GLFW/OpenGL, with an
optional exact-symbolic (CAS) path.

You give it a system of ODEs or a map in a small `.dyn` text format; it renders
the flow, finds and classifies equilibria, traces limit cycles, continues
bifurcations in one and two parameters, and (when the CAS bridge is available)
computes exact eigenvalues and certified derivatives.

## Build

Dependencies are provided **entirely through the Nix flake** — Dear ImGui, GLEW,
GLFW, cglm, the OpenGL/X libraries, the TPCAS parser frontend, and the Lizard +
Sangaku CAS. The repository vendors no third-party source itself.

```
nix develop        # enters the dev shell; exports IMGUI_DIR, TPCAS_DIR, GLEW_DIR, LIZARD, SANGAKU_ROOT
make               # builds build/dynsys
make run           # or ./build/dynsys
make clean         # remove build artifacts
```

If you build without `nix develop`, supply the same variables by hand, e.g.:

```
make IMGUI_DIR=/path/to/imgui TPCAS_DIR=/path/to/tpcas
```

The Makefile checks these are set and errors with a hint if they are not.

### Windows cross-build

```
nix develop .#windows    # MinGW-w64 cross toolchain + static GLFW/GLEW
make windows             # -> build/windows/dynsys.exe
```

## Run

```
./build/dynsys                        # interactive GUI (loads a default system)
./build/dynsys --headless examples/lorenz.dyn --steps 10000   # headless integration
```

In the GUI the plot fills the window; controls are the top toolbar and the left
panel. In 2-D phase view, left-click adds an orbit, left-drag pans, the wheel
(or the zoom buttons / `+`/`-`) zooms, and double-click auto-fits. In 3-D view,
left-drag rotates, right-drag pans, and the wheel zooms.

## The `.dyn` format

A system is a short text file. Example (`examples/lorenz.dyn`):

```
state x, y, z
mode = ode
integrator = rk4

param sigma = 10 [0,50]
param rho   = 28 [0,100]
param beta  = 8 / 3 [0,10]

plot3d = x, y, z
plot2d = x, z
observe r = sqrt(x*x + y*y + z*z)

section = z - 27
section_direction = positive
section_plot = x, y

initial x = 0.1
initial y = 0
initial z = 0

dx = sigma * (y - x)
dy = x * (rho - z) - y
dz = x * y - beta * z
```

`mode = ode` integrates a continuous vector field (`dx = ...` per state);
`mode = map` iterates a discrete map; `mode = ifs` runs an iterated function
system. Parameters take an optional `[lo,hi]` range for the sliders and for
continuation. See `examples/` for maps (Hénon, Thomas), oscillators (Van der
Pol), predator-prey (Lotka-Volterra), and higher-dimensional systems.

## What it computes

- Phase-plane flow, nullclines with flow-direction arrows, and a speed-colored
  vector field; automatic equilibrium finding and linear classification with
  stable/unstable manifolds.
- 3-D attractor view; basins of attraction (parallel, deterministic) and
  fractal/escape-time renders; box-counting dimension; Lyapunov spectrum.
- Equilibrium continuation with fold / Hopf / branch-point detection and branch
  switching; first Lyapunov coefficient; fold and cusp normal forms.
- Two-parameter fold and Hopf curves with codim-2 detection (Bogdanov-Takens,
  cusp, generalized-Hopf, zero-Hopf, Hopf-Hopf) and their normal-form
  coefficients; direct codim-2 point location and continuation of the BT,
  zero-Hopf, cusp, and Hopf-Hopf loci as curves in a third parameter.
- Limit-cycle continuation (collocation + pseudo-arclength, adaptive mesh) with
  Floquet multipliers; fold-of-cycles, period-doubling, and Neimark-Sacker
  curves; branch-point-of-cycles; homoclinic and heteroclinic solving and
  continuation.
- An optional exact-symbolic path (the CAS bridge): when the `lizard` binary and
  Sangaku are available (the flake supplies both), eigenvalues, equilibria, and
  derivatives can be computed symbolically; the bridge degrades gracefully to
  the numeric path when they are absent.

See `docs/` for the design notes, the MatCont comparison, and the roadmap.

## Tests

```
make test          # the full smoke-test suite (numeric core)
make test-cas      # the exact-symbolic (CAS) bridge tests; needs LIZARD + SANGAKU_ROOT
```

## Scope

dynsys is an independent tool, not a drop-in replacement for MatCont or AUTO. It
validates each analysis feature against analytic ground truth in its test suite,
but it does not claim their robustness track record. `docs/MATCONT_COMPARISON.md`
states honestly what is at rough parity and what remains.
