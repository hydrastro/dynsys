# dynsys — bug fixes: NaN bounds, zoom recovery, custom orbits, fixed views

Unzip at repo root, replace all, then **`make clean && make`** (the
clean is important), then `make run`. Confirm the green "dynsys NEW-UI"
label at the top — that means you're on this build.

## Fixes for what you reported

1. **Logistic map "nan/nan coordinates" — fixed at the real cause.**
   I ran the logistic map through the actual IR engine in a test harness:
   it stays bounded (no nan). The nan came from the *view bounds*, not
   the map: a bad zoom or a stray non-finite value collapsed the bounds,
   `screen_to_plot` then produced nan coordinates, and that nan fed back
   into the zoom anchor and poisoned every later frame. Added a single
   `sanitize_bounds` choke point (guarantees finite, correctly-ordered
   bounds with a sane minimum span) applied wherever bounds are set or
   read, plus a non-finite-anchor guard in zoom so it always recovers.

2. **"Sometimes the zoom bugs and shows nothing" — same root cause.**
   Degenerate/collapsed bounds are now impossible, so the view can't go
   blank from zooming. If a zoom ever lands oddly, double-click still
   refits.

3. **Custom orbits — added.** In the Phase plane panel, "Add custom
   orbit" gives a numeric field per state variable, plus "Add orbit at
   these values" and "Use current state". (Clicking the plot still adds
   an orbit at the cursor; orbits accumulate until Clear.)

4. **"2D plot isn't fixed, it enlarges" — presets now open framed.**
   New `view2d = xmin, xmax, ymin, ymax` directive sets a fixed initial
   window. Added to the 2D presets (Van der Pol -4..4 x -6..6,
   Lotka-Volterra 0..4 x 0..4, damped pendulum, Hénon-2D) so they open
   with a good steady view instead of auto-growing. You can add `view2d`
   to your own systems too.

5. **Nullclines + their direction arrows.** The pplane-style arrows along
   nullclines are present and on by default ("nullcline arrows" toggle).
   IMPORTANT: nullclines and auto fixed points are mathematically
   ODE-only — a *map* (logistic, Hénon) has no continuous vector field,
   so there's nothing to draw. If you were testing on maps, that's why
   you saw none. Load **Van der Pol** or **Lotka-Volterra** (ODEs) to see
   nullclines, their flow arrows, and the auto-classified equilibria with
   stable/unstable manifolds. Maps now show an explicit note saying so.

## Carried forward
- Full-window background plot, top toolbar + left control panel.
- Left-click adds orbits (accumulate; Clear to remove), drag pans,
  wheel/+/- zoom, double-click auto-fit; zooming never moves a window.
- Auto-scan & classify all equilibria in view (test-fp), dimension
  auto-detect + Force2D/Force3D (test-dim), marching-squares nullclines
  (test-nullcline), 3D-in-FBO, analysis engine, portable ImTextureID.

## Verification
- All 4 C++ translation units compile with ZERO warnings (-O2) against
  the real Dear ImGui / GLFW / cglm / tpcas headers, ImTextureID set to
  the integer type your nix build uses.
- The logistic iteration was verified bounded via a numeric harness
  driving the real IR.
- make test: analysis(23) + AD(11) + nullcline + dimension + fixed
  points — all pass.
- GUI rendering still needs your `make run`; a screenshot (with the green
  NEW-UI label visible) lets me confirm visuals on the same binary.

## Next
- Speed-colored streamlines for denser ODE portraits.
- Click-to-select the equilibrium nearest the cursor.
- matcont-style bifurcation diagram UI (gated on lizard symbolic math).
