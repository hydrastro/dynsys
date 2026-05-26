# Changes

## Optimization pass (in progress, MatCont deferred)

- Added an expression IR (`src/expr_ir.{h,cpp}`) — flat opcode array,
  power-of-two-aligned, with PushConst/PushState/PushParam/PushLocal,
  PushT/PushPi/PushE, Neg/Add/Sub/Mul/Div, CallBuiltin (19 builtins:
  sin/cos/tan/asin/acos/atan/exp/log/log10/sqrt/abs/floor/ceil/sign/
  pow/min/max/mod/clamp), CallDef, and BrIfZero/Jump for lazy
  control flow. The IR module is self-contained — no AppState, no
  ImGui dep — and ships with `test/ir_smoke.cpp`, a 56-case
  standalone test that pins down operator precedence, all builtins,
  0-arity-def memoization, n-arity user-function calls, lazy
  `select` (where the old eager version produced NaNs for
  `select(x, log(x), -1)` when `x = 0`), unknown-name errors,
  arity-mismatch errors, and mutual-recursion cycle detection.

- `compile_system` now lowers every parsed `node_t *body` to an
  `ir::Program` after parsing finishes, in two phases (register all
  def signatures first, then lower bodies — so mutual references
  between defs work). Lowering happens BEFORE the transactional
  AppState swap; if it fails the next-arena is destroyed cleanly
  without touching the live system.

- Hot-path callers route through the IR: `eval_rhs`,
  `step_ode_state`, `step_map_state`, `eval_plot3d`,
  `maybe_record_poincare`, `value_by_name`. Each compile runs a
  self-check that evaluates every program both via the new IR and
  via the original `eval_expr_at` AST walker on a probe state and
  aborts the compile if they disagree by more than 1e-10 relative.

- Pre-allocated integrator scratch states `scratch_k1..k4` and
  `scratch_mid` live in `AppState`, sized to `dim` at compile time.
  RK4 now performs zero `std::vector<double>` allocations per
  integration step (previously eight: four k-states + three mid
  states + the output state).

- New flag `--headless [path.dyn] [--steps N] [--use-ast]
  [--dump]` runs a system without GLFW for differential testing
  and benchmarking. Trajectories are bit-identical between IR and
  AST paths across 200k-step runs on every preset.

- `find_fixed_point` no longer recomputes the Jacobian after
  convergence — keeps the last Newton-iteration matrix.

- `pc_list_new`-style hidden zero-init dependency: no longer
  relevant here, but `param_values` (the flat array the IR reads)
  is explicitly synced from `params[i].value` on every slider
  update, reset, and after `compile_system`.

### Numbers

Release build, x86_64, 200k integration steps:

```
                    AST ns/step    IR ns/step   speedup    final state
Lorenz (3D)              1636.7         348.5      4.70x   identical
Lotka-Volterra (2D)      1330.9         265.2      5.02x   identical
Van der Pol (2D)         1105.1         229.5      4.81x   identical
4D demo                  2282.3         400.3      5.70x   identical
```

ASAN+UBSan clean on 50k-step runs across all presets.

## v4

- Added direct click-to-set initial conditions in the 2D phase plane.
- Added shift-click extra trajectory creation.
- Added phase-plane right-drag panning, cursor-centered wheel zoom, and double-click auto-fit.
- Added multi-trajectory phase-plane rendering plus grid/circle seed tools.
- Added 2D fixed-point classification from finite-difference Jacobian eigenvalues.
- Added fixed-point markers directly in the phase-plane plot.
- Added saddle separatrix tracing for 2D ODEs.
- Added equal-aspect phase-plane mode, tick labels, and normalized/raw vector-field controls.
- Kept bifurcation/continuation deferred; v4 focuses on pplane-quality local analysis.

## v3

- Added a true 2D phase-plane tab separate from the 3D OpenGL renderer.
- Added selectable phase-plane axes over the declared state variables.
- Added phase-plane vector-field arrows.
- Added approximate ODE nullcline drawing for the selected phase-plane components.
- Changed 3D zoom to scale the model instead of pushing the camera through the near plane.
- Added optional orthographic 3D projection and explicit scene-scale/camera-distance controls.
- Removed bifurcation scanning from the active UI; it is now deferred to a future continuation subsystem.
- Added a continuation roadmap panel for the later MatCont-like/Lisp-language direction.

## v2

- Generalized the runtime state from hardcoded `x,y,z` to named state variables via `state ...`.
- Added generic ODE equations for every state variable: `dx`, `dtheta`, `dx/dt`, `dtheta/dt`, etc.
- Added generic map equations for every state variable: `x_next`, `x'`, `next x`, etc.
- Added `plot3d = expr, expr, expr` for rendering arbitrary 3D projections of N-dimensional systems.
- Added `initial var = expr` / `start var = expr` source-level initial conditions.
- Generalized RK/Euler integration over the full state vector.
- Generalized discrete map stepping over the full state vector.
- Generalized Lyapunov distance and renormalization over the full state vector.
- Generalized Newton fixed-point search and finite-difference Jacobian over the full state vector.
- Generalized trajectory CSV export to include all state variables and observables.
- Auto-generates GUI start-value controls for every declared state variable.
- Added N-dimensional example presets and example `.dyn` files.
