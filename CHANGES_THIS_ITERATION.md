# dynsys — box-counting fractal dimension + vector-field color fix

Unzip at repo root, then `make clean && make && make run`. Green
"dynsys NEW-UI <date>" label = you're on this build.

Now that the parameter-sync fix made every sweep trustworthy, this round
adds a genuinely new analysis capability and fixes a rendering bug I found
while looking around.

## New: box-counting (fractal) dimension

In the Analysis tab, under the Lyapunov spectrum, there's a new
**"Measure box-counting dimension"** button. It estimates the
Minkowski-Bouligand (box-counting) dimension of the set currently in the
2D plane:
- For a MAP it samples the attractor by iterating (200k points after a
  transient).
- For an ODE it measures the live trajectory (let it run a bit first).
It reports D_box and the log-log fit quality R^2.

This is a real piece of nonlinear-dynamics kit that none of pplane / XPPAUT
/ MatCont / AUTO surface, and it complements the Kaplan-Yorke dimension
already there (Kaplan-Yorke comes from the Lyapunov spectrum; box-counting
is measured directly from the geometry — nice to compare the two).

Validated headlessly (`make test-boxdim`) against shapes of known
dimension: filled square -> 1.99, line -> 1.04, circle -> 1.09, Sierpinski
triangle -> 1.60 (theory 1.585), Henon attractor -> 1.31 (literature ~1.26),
all with R^2 >= 0.995.

## Fixed: vector-field arrow coloring

While wiring the above I noticed the phase-plane vector field's speed->color
mapping had a collapsed term (a stray * 0.0) that tinted almost every arrow
the same. It now scales color to the actual max speed in the current view
(two-pass): slow = blue, fast = warm. The flow structure reads much more
clearly. (I also caught myself starting to build a duplicate direction-field
renderer and removed it — the existing one just needed the fix.)

## Roadmap
docs/COMPARISON_AND_ROADMAP.md updated: box-counting dimension marked done
(Phase C). Remaining: IFS/chaos game; the big one is still Phase D step 2
(limit-cycle continuation) for full MatCont parity; Phase E waits on Lizard.

## Verification
- All 4 C++ TUs compile with ZERO warnings (-O2).
- make test: 19 checks/suites pass (added **test-boxdim**).
