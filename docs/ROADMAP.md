# dynsys roadmap — toward the best dynamical-system analysis tool

This is a living plan. It is organized by capability, not by calendar.
Each phase lists what already exists, what is missing, and the concrete
next move. The north star is a tool that is *better than pplane* for
interactive 2D phase-plane work and *competitive with matcont* for
continuation and bifurcation analysis, built on the dynsys/tpcas/lizard
stack rather than MATLAB.

The single most important architectural fact: every model expression is
lowered to a flat stack-VM program (`dynsys::ir::Program`) once at
compile time. That IR is the leverage point. Numerical analysis targets
it directly, and the eventual symbolic layer (lizard + tpcas) emits into
the same IR. Nothing downstream needs to know whether a Jacobian came
from finite differences, forward-mode AD over the IR, or a symbolic
derivative compiled by lizard.

---

## Phase 0 — Correctness foundation (in progress)

The interactive layer had two structural defects that made analysis
results feel unreliable.

1. **Analysis objects were forgotten on every recompile.** `reset_simulation`
   cleared `phase_trajectories`, `separatrix_curves`, and the located
   fixed point, and it ran after every "Apply system", every preset
   change, and the `C` key. A separate `clear_analysis_objects` existed
   for exactly this purpose but was never the one on the hot path.
   - **Fix:** `reset_simulation` now resets *only* the running orbit and
     its derived buffers. Analysis survival is decided explicitly:
     editing parameters or pressing reset keeps custom orbits; changing
     the *equations* invalidates located equilibria (their coordinates
     no longer satisfy the new RHS) but keeps user-seeded orbits unless
     the state dimension/names changed.

2. **The main orbit lived in two stores with different lifetimes.**
   `history` (deque, `history_limit` long) drove the phase/series views;
   `points[]` (GL ring buffer, `num_points` long) drove the 3D view. The
   two disagreed about where the orbit began and used different colors.
   - **Fix:** a single `history_limit == num_points` invariant plus a
     shared color rule, with `history` remaining the canonical N-D store
     and the GL buffer treated as its projected mirror. (Full
     single-source-of-truth refactor is Phase 5.)

---

## Phase 1 — General stability analysis (this iteration)

**Have:** N-D Newton solver for equilibria; finite-difference Jacobian
(N-D); eigenvalue classification *only for 2D* via closed-form
trace/determinant.

**Missing:** eigenvalues in 3D and higher; therefore no general
node/spiral/saddle/focus classification, and no Hopf detection.

**Move (done here):** a dependency-free real eigensolver for general
N×N matrices — balancing, reduction to upper Hessenberg, then the
double-shift (Francis) QR iteration returning complex eigenvalues. On
top of it, dimension-agnostic classification: count eigenvalues by sign
of real part, flag a complex pair crossing the imaginary axis (candidate
Hopf), flag a real eigenvalue crossing zero (candidate fold). This lifts
every existing 2D feature to arbitrary dimension and is the prerequisite
for continuation event detection.

---

## Phase 2 — Continuation engine (this iteration: equilibrium curves)

**Have:** brute-force parameter sweeps (`run_bifurcation`) that record
attractor samples — a bifurcation *diagram*, not continuation.

**Missing:** the matcont core — tracking a solution branch as a
parameter varies, through folds, with codimension-1 event detection.

**Move (done here):** a pseudo-arclength continuation engine for
equilibria. Given a located fixed point and a chosen continuation
parameter, it predicts along the tangent to the augmented system
`{ f(x, p) = 0, arclength condition }`, corrects with Newton on the
bordered system, and adapts step size. At each accepted point it
computes the Jacobian spectrum (Phase 1) and watches two test
functions — the smallest `|Re λ|` of a complex pair (Hopf) and
`det(f_x)` (fold) — bracketing and bisecting sign changes to locate
bifurcations on the branch. Output is a branch object: a list of
`(p, x, eigenvalues, stability)` plus tagged special points.

**Next after this:** limit-cycle continuation (needs a boundary-value /
collocation solver and a phase condition), branch switching at
detected branch points, and two-parameter continuation of codim-1
curves.

---

## Phase 3 — Symbolic Jacobian via the IR (forward-mode AD: done)

**Have now:** exact Jacobian columns and df/dp via forward-mode
automatic differentiation evaluated directly over the IR. A separate
dual-number executor (`expr_ir_ad.{h,cpp}`) mirrors the plain IR
`run()` opcode for opcode — including `CallDef` recursion, 0-arity
memoization, and the `select()` branch ops — propagating `(value,
derivative)` pairs with per-builtin derivative rules. Seeding the dual
on state variable j yields Jacobian column j; seeding the active
parameter yields df/dp. `build_model` in dynsys.cpp installs these as
`Model::jacobian_x` / `Model::dfdp`, so the continuation engine's
linear algebra is now exact rather than finite-difference-noisy, with
no change to any caller. Verified in `test/ad_smoke.cpp` against
analytic derivatives and central finite differences (`make test-ad`).

**Next after this:** lizard/tpcas *symbolic* differentiation of `.dyn`
expressions, re-lowered to IR, differential-tested against this AD
path (the AD path becomes the oracle); second derivatives for
fold/Hopf normal-form coefficients; using the exact Jacobian inside
the Newton equilibrium solver and a future stiff integrator.

---

## Phase 4 — Integrators

**Have:** fixed-step Euler / RK2 / RK4.

**Missing:** adaptive embedded RK (Dormand–Prince) with error control;
a stiff solver (implicit, using the exact Jacobian from Phase 3);
event-located integration (stop/record exactly at a zero of a user
event function, reusing the bracketing code from Phase 2).

---

## Phase 5 — Model language and persistence

**Have:** a compact `.dyn` line language; presets.

**Missing:** a real language. The intended path is lizard as the
scripting/CAS layer: define systems, scripted parameter studies,
symbolic preprocessing, and saved sessions (orbits, branches, view
state) as s-expressions. Single-source-of-truth refactor of the orbit
representation also lands here.

---

## Phase 6 — Visualization and UX

Cached nullclines (recomputed only on parameter/equation/bounds change,
not every frame), corrected marching-squares for cells with four edge
crossings, branch/eigenvalue overlays in the phase plane, manifold
rendering in 3D, and export of branches and bifurcation data.

---

## Status legend used in code comments

- `PHASE0`/`PHASE1`/... tags mark code introduced for a given phase.
- The new numerical core lives in `src/analysis.{h,cpp}`, deliberately
  free of AppState/ImGui dependencies so it is unit-testable on its own
  (see `test/analysis_smoke.cpp`), mirroring how `expr_ir` is structured.

---

## Phase 6 progress (visible-UX iteration)

- 2D systems no longer render the 3D OpenGL duplicate (`system_is_3d`
  gate in `draw_scene`); the phase plane is the 2D view.
- Nullclines rewritten as proper marching squares (node-sampled,
  16-case, saddle-correct), drawn as solid curves with a naming legend.
  Geometric correctness checked in `test/nullcline_smoke.cpp`.
- Phase-plane auto-bounds eased (grow-fast / shrink-slow) so the view
  no longer jitters or collapses as the orbit evolves.

Still open under Phase 6: nullcline caching (recompute only on
parameter/equation/bounds change), eigendirection overlays at fixed
points, and the matcont-style bifurcation *diagram* UI (gated on
lizard symbolic Jacobians, by request).
