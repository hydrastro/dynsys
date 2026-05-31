# dynsys — master plan

This is the single source of truth for where dynsys is going. We
implement it **little by little, carefully**, one numbered task per
iteration, each ending in: compiles warning-clean + tests pass + you
confirm it on screen.

---

## Where we are (honest assessment)

**Strong:** the numerics. A headless, unit-tested analysis core
(eigensolver, equilibrium classification + scan, pseudo-arclength
continuation with fold/Hopf, forward-mode AD, bifurcation/Lyapunov
sweeps). The 2D phase portrait is genuinely pplane-grade.

**Weak:** UX and code organization.
- `dynsys.cpp` is ~5,100 lines — one giant file.
- The side panel is **12 collapsing headers in one long scroll**.
  Important things (Parameters, Bifurcation diagram) get buried — this is
  literally why "I can't see the bifurcation map" and "I can't see
  parameters" happened. They render; they're just below the fold.
- **1D systems have no view at all.** The phase background returns early
  when there are fewer than 2 state variables, so the tent map and the
  logistic-as-1D draw nothing. These are THE canonical bifurcation
  examples, so this is a real hole.

**Conclusion:** the next phase is not more math — it's making the
existing power visible and the code maintainable.

---

## Root-cause fixes these two bugs imply (Phase A, do first)

1. **A real view system.** Instead of "2D phase OR 3D" with everything
   else in panels, have an explicit set of VIEWS the active system can
   show, chosen in the toolbar, each filling the plot area:
   - **1D line view** (for 1D maps/ODEs): cobweb diagram for maps, flow-
     on-a-line for 1D ODEs. Fixes the tent map.
   - **2D phase** (exists).
   - **3D scene** (exists).
   - **Bifurcation diagram** (exists as a buried panel → promote to a
     full view).
   - **Time series**, **Poincaré** (exist as panels → optionally views).
   The toolbar offers only the views that make sense for the current
   system's dimension/mode. This single change fixes both reported bugs
   and is the backbone for everything later.

2. **Reorganize the side panel** so it's not a 12-item scroll: group into
   a few tabs — **System** (editor, params, dim), **View** (controls for
   the active view), **Analysis** (fixed points, continuation,
   bifurcation, Lyapunov), **Data** (time series, Poincaré, export).
   Only show controls relevant to the active view.

3. **Split `dynsys.cpp`** into a few translation units so it's
   maintainable and faster to compile:
   - `app_state.h` (the AppState struct + enums)
   - `engine.cpp` (compile_system, integrators, stepping, eval)
   - `views.cpp` (all rendering: phase, scene, 1D, bifurcation, series)
   - `ui.cpp` (toolbar, side panel, controls)
   - `main.cpp` (glue, main loop)
   Done incrementally so nothing breaks; the analysis core already lives
   in `analysis.{h,cpp}`.

---

## Phase B — richer numerics (no lizard needed)

**DONE: Full Lyapunov spectrum + Kaplan-Yorke dimension** (test-lyap; validated
against Lorenz, Henon, stable-linear). Wired into the Analysis tab.

- **Full Lyapunov spectrum** (Benettin/Gram-Schmidt on the variational
  flow; we have AD already) → Kaplan-Yorke (Lyapunov) dimension. First
  real fractal number.
- **DONE: Basins of attraction** (test-basin; validated on a bistable ODE and Newton z^3-1)
  ~~Basins of attraction~~ original: (grid of ICs colored by destination) → basin
  fractals.
- **Streamlines / LIC** for dense ODE portraits.
- **DONE: adaptive/multiple integrators** (Euler, Heun, RK2, RK4, RK 3/8, RKF45, Dormand-Prince; test-solver)
  ~~Adaptive integrators~~ (RK45 Dormand-Prince), event detection,
  symplectic option.
- **DONE: 2-parameter scans** (largest-Lyapunov shrimp maps; test-scan)
  ~~2-parameter scans~~ (period / largest-Lyapunov over a (p1,p2) grid)
  → "shrimp"/period maps.
- **Bifurcation diagram polish**: zoom/pan, click-to-set parameter,
  export, multi-observable.

## Phase C — fractals as first-class objects (mostly no lizard)

**DONE: 3D bridge scene** — the Mandelbrot set (flat) with the logistic
bifurcation diagram rising along its real axis; they align via the verified
conjugacy c(r)=r/2-r^2/4 (test-bridge). Always-available 3D view.

**DONE: escape-time fractal view (Mandelbrot/Julia bridge)** — a full-window
view that iterates the system's 2D map over the plane. Parameter space =
Mandelbrot-type, state space = Julia-type, driven by the SAME expression
language as the simulator. New 'Complex quadratic' preset. Escape-time
membership locked by test-fractal.

**DONE: bifurcation diagram polish** — density (log-histogram) rendering,
zoom/pan over both axes, axis tick labels, Lyapunov overlay. Available for
any system with a parameter (not just the 1D maps).

- **Escape-time fractals** (Mandelbrot/Julia) where the iterated map is
  ANY user expression compiled through the existing IR → fractal explorer
  speaks the same language as the simulator.
- **IFS / chaos game** (Barnsley fern, Sierpinski).
- **Newton fractals** (symbolic f, f' once CAS exists).
- **Fractal dimension** (box-counting, correlation dimension).

## Phase D — matcont-grade continuation UI (engine exists)

- **Equilibrium continuation diagram** view (branches, fold/Hopf/branch
  points, click to start).
- **Limit-cycle continuation**, period-doubling cascades.
- **2-parameter codim-1 curves** + codim-2 points (needs CAS for
  robustness).

## Phase E — lizard + CAS integration (the symbolic leap)

- **Symbolic Jacobians/Hessians** replace AD/finite-diff.
- **Exact equilibria** via symbolic solving (find ALL, including unstable).
- **Symbolic nullclines** (analytic, no sampling error).
- **Normal-form coefficients** (Hopf first Lyapunov coefficient, etc.) —
  the matcont killer feature.
- **`cas_bridge`**: translate CAS S-expressions ↔ tpcas AST ↔ `expr_ir`
  VM. One language for modeling, analysis, and rendering; CAS off the hot
  path (cached), IR on it.

---

## Working agreement (how we do this carefully)

- **One numbered task per iteration.** No big-bang rewrites.
- Every iteration: compiles with **zero warnings**, `make test` passes,
  and where math is involved, a headless test locks the behavior.
- **You confirm each on screen** before we move on (the green NEW-UI
  stamp tells you the build is current).
- Keep the analysis core headless and tested; keep growing that test
  suite (it's our safety net).
- Prefer promoting existing-but-hidden features over writing new ones.

---

## Immediate next step (proposed)

**Phase A, Task 1: the 1D view + view system**, because it fixes the two
bugs you just hit (tent map invisible; bifurcation buried) and lays the
backbone. Specifically: add a toolbar view selector that offers only the
views valid for the current system, implement the 1D cobweb/line view,
and promote the bifurcation diagram to a full view.
