# dynsys — fixes round (warnings, zoom, 2D view, eigendirections)

Unzip at repo root, replace all, then `make clean && make run`.

## First: those messages were warnings, not errors
The build was completing. This round silences them anyway:
- `_FORTIFY_SOURCE requires compiling with optimization`: the nix stdenv
  injects `-D_FORTIFY_SOURCE=2`, which needs `-O`. Debug build used
  `-O0`. Changed debug to `-Og` (full debuggability, satisfies fortify).
- `missing field 'state' initializer`: `EvalContext` is now fully
  value-initialized.
- the `unused function` warnings: those helpers are genuinely unused and
  are now marked `[[maybe_unused]]` (matching existing style), including
  `clear_analysis_objects`, unused since its call sites were replaced.

## Fixes for what you reported

1. Zoom now works (phase plane). The wheel-zoom scales the view span
   about the cursor, keeps the cursor point fixed, and writes straight
   to the manual bounds so it sticks. Verified the anchor stays fixed
   and the span changes.

2. 2D systems no longer draw a duplicate orbit in the 3D backdrop. The
   3D view is the GL window background (not a closable panel), so for a
   2D system it now paints a clean "2D system - see the phase plane tab"
   hint instead of the redundant curve. Declare a 3rd state variable and
   the 3D orbit returns. (The launch default is Lorenz, which is 3D;
   switch to Van der Pol / Lotka-Volterra / Damped pendulum for 2D.)

3. Phase-plane auto-fit no longer looks weird: replaced per-frame easing
   with a stable grow-only fit (holds steady, enlarges in chunks only
   when the orbit leaves it, like pplane). Double-click refits.

## New feature

4. Eigendirection rays at 2D fixed points. After "Find fixed point" the
   marker draws the linearized eigenvectors as rays colored by stability
   (green stable, orange unstable) with flow arrowheads. Complex pairs
   draw none, as expected.

## Verification
- dynsys.cpp compiles warning-clean at -O2 and -Og against real Dear
  ImGui / GLFW / cglm / tpcas headers.
- make test: 23 analysis + 11 AD + nullcline geometry pass.
- Zoom and marching-squares checked with standalone numeric tests.

## Next
- Nullcline caching so they're free to leave on.
- Vector-field arrows scaled/colored by speed.
- matcont-style bifurcation diagram UI (gated on lizard symbolic
  Jacobians, per your call).
- Find equilibrium nearest a clicked point; multiple stored equilibria.
