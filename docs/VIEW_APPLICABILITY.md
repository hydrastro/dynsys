# Which views apply to which systems (consistency audit)

You asked the right question: the views are NOT all applicable to every
system, and some restrictions are real and intentional. Here is the honest
matrix, the reasoning, and the one genuine inconsistency (IFS) plus its fix.

## The applicability matrix (source of truth: view_valid + view_requirement)

| View | Applies to | Why restricted |
|---|---|---|
| 1D line / cobweb | 1-D systems only | a single state variable on a line |
| 2D phase | dim >= 2 | needs an x and a y state axis |
| 3D scene | dim == 3 | needs three state axes |
| Bifurcation | any system WITH a parameter | sweeps a parameter; needs one |
| Fractal (escape-time) | 2-D MAPS only | it's z -> f(z) iteration in the complex plane |
| 3D bridge | 1-D MAPS only | the logistic<->Mandelbrot needle is specific to 1-D maps |
| Basins | dim >= 2 (ideally multistable) | a 2-D grid of initial conditions; only interesting with several attractors |
| Param scan | >= 2 parameters | colors a 2-parameter grid |
| Continuation (equilibria) | ODEs WITH a parameter | equilibria & their stability are a flow concept |
| Limit cycle | ODEs WITH a parameter | periodic orbits are a flow concept |
(IFS is NOT a view — it's a system kind; see below)

So most restrictions are correct and principled:
- **Maps vs ODEs matter.** Equilibrium continuation and limit-cycle
  continuation are meaningless for a map (a map has fixed points and
  periodic points, not flows/equilibria in the ODE sense). The escape-time
  fractal and the Mandelbrot bridge are map-only by definition.
- **Dimension matters.** 1-D/2-D/3-D views require the matching number of
  state variables.
- **Parameters matter.** Bifurcation, both continuations, and the parameter
  scan all sweep parameters, so they need 1 (or 2) to exist.

These now show a tooltip explaining the requirement when a view is greyed
out, so the rules are visible rather than mysterious.

## The real inconsistency: IFS

Every other view is a *lens on the loaded system*. IFS is not — it renders a
hardcoded gallery (fern, Sierpinski, dragon, ...) that ignores whatever
system is loaded. So "switch to IFS while the Lorenz system is loaded" is
conceptually incoherent: the view has nothing to do with your model. For
now the view says so explicitly ("a standalone fractal gallery, independent
of the loaded model").

## The fix: make an IFS a first-class model (planned next)

An IFS *is* a discrete dynamical system — a finite set of affine contraction
maps applied stochastically (the chaos game). The clean, consistent design
is a third system mode alongside ODE and Map:

    mode = ifs
    ifs_map = 0.5, 0, 0, 0.5, 0,    0      # a b c d e f  (prob optional)
    ifs_map = 0.5, 0, 0, 0.5, 0.5,  0
    ifs_map = 0.5, 0, 0, 0.5, 0.25, 0.5
    # x' = a*x + b*y + e ;  y' = c*x + d*y + f

With `SystemMode::IFS`:
- **Stepping** = the chaos game (already implemented and validated in
  analysis::chaos_game; box-dimension cross-checked: Sierpinski ~1.585,
  fern ~1.8).
- **Consistency falls out automatically** from view_valid:
  - IFS view -> applies (renders the loaded maps).
  - Box-counting dimension -> applies (an IFS attractor is a fractal).
  - Bifurcation / continuation / limit-cycle -> correctly DO NOT apply
    (no parameter-vs-attractor, no equilibria, no flow).
  - Basins -> does not apply (the chaos game has one attractor by
    construction).
  - 1D/2D/3D phase, vector field -> do not apply (no state-vector ODE/map).
- The current built-in gallery (fern, Sierpinski, dragon, tree, maple leaf)
  becomes a set of **presets** in this mode, exactly like the ODE/map
  presets — so "load Barnsley fern" works like "load Lorenz".
- Optional polish: parameterize a map coefficient (e.g. a rotation angle)
  so even an IFS can have a swept parameter and a small animation — but
  that's a later nicety, not required for consistency.

This turns IFS from an odd detached view into a proper model type, making
the whole tool consistent: load a system of some kind (ODE, map, or IFS),
and the views that apply to that kind light up.

## The cleaner separation (system KIND vs view LENS)

There are two distinct axes, and they should not be mixed:

- **System kind** (what you input): ODE, map, IFS. Chosen by `mode =` in the
  model text. This is the model.
- **View / lens** (how you look at the loaded model): 1D line, 2D phase, 3D,
  bifurcation, basins, continuation, limit cycle, ... These are ways of
  *looking*, shared across systems.

An ODE doesn't get its own "ODE view" button — it's displayed through the
universal lenses (3D, phase, bifurcation, ...). **IFS is the same: it's a
system KIND, not a view.** So:

- An IFS is loaded as a model (`mode = ifs`, or the Barnsley fern / Sierpinski
  presets).
- It is displayed in the **2D phase view** (an IFS attractor is a planar
  object), rendered via the chaos game — exactly as an ODE shows its
  trajectory there.
- There is **no separate "IFS" button** in the view row. (An earlier version
  had one; that conflated the two axes. Removed.)
- In IFS mode, only the 2D phase view applies; bifurcation / continuation /
  limit-cycle / basins / 1D / 3D correctly do not (no vector field,
  parameter, or equilibria).

## Status
- DONE: applicability tooltips on every greyed-out view.
- DONE: IFS is a first-class system KIND (`mode = ifs`, `ifs_map` lines,
  fern/Sierpinski presets).
- DONE: IFS displays through the normal **2D phase view** (not a bespoke
  view button) — system-kind and view-lens are now cleanly separated, so
  IFS is treated like any other dynamical-system input.
