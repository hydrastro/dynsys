# dynsys — 3D bridge reset bug fixed, + two MatCont normal-form coefficients

Unzip at repo root, `make clean && make && make run`.

## 1. The 3D bridge (and other caches) now reset on EVERY system (re)compile
You were right: the bridge could keep showing a previous system. The reset was
gated on the variable/parameter NAMES changing — so loading a different system
that happened to reuse the same names (e.g. both use x,y and a parameter a, but
different equations) would NOT clear the bridge, bifurcation diagram, fractal,
basin, or scan caches. Now those stale-data invalidations run unconditionally
on every successful compile (the name-gated block still handles only the
name-structural defaults like the bifurcation parameter range/observable).
Concretely, each recompile now clears: the bifurcation point cloud, the 3D
bridge geometry, and marks the fractal/basin/scan images dirty, plus the
certified-eigenvalue and Hopf/fold classifications.

## 2. MatCont normal-form coefficients at equilibria
Two codim-1 normal-form quantities MatCont reports and most teaching tools
don't — both computed from finite differences of the vector field, so they
work for ANY system, and both verified against textbook normal forms:

- Hopf FIRST LYAPUNOV COEFFICIENT l1 (this had been implemented but never
  shipped). At a Hopf point: l1<0 supercritical (a stable limit cycle is
  born), l1>0 subcritical (hard loss of stability), l1~0 degenerate (Bautin /
  generalized-Hopf codim-2). Button: "Classify Hopf (first Lyapunov coeff)".

- Fold (limit point) NORMAL-FORM COEFFICIENT a (new this round): a = (1/2)
  <p,B(q,q)>/<p,q> from the right/left null vectors at the fold. a != 0 is the
  genuine quadratic-fold condition; a ~ 0 signals a CUSP (codim-2). Button:
  "Classify fold (normal-form coeff a)". It self-checks for a near-zero
  eigenvalue and tells you when the equilibrium isn't actually at a fold.

Both appear in the Analysis tab's fixed-point block once an equilibrium is in
hand. Verification (headless, under ASan/UBSan):
- fold x'=p+x^2 at the fold -> a = 1 exactly; planar saddle-node -> |a| = 1;
  cubic x'=p+x^3 -> a ~ 0 (correctly flags the cusp).
- supercritical Hopf normal form -> l1 < 0; subcritical -> l1 > 0.

These join the existing equilibrium continuation with fold/Hopf detection, so
dynsys now both DETECTS codim-1 points on a branch and CLASSIFIES their
criticality the way MatCont does.

## Verification
- All 4 C++ TUs compile with ZERO warnings (-O2 -Wall -Wextra).
- Brace/tab structure balanced after edits.
- make test: 27 groups (added test-foldnf) + graceful CAS skip.
- New fold coefficient + Hopf l1 verified against known normal forms under
  AddressSanitizer + UBSan.

## Still pending toward fuller MatCont parity (not in this build)
Codim-2 point detection ON a branch (Bogdanov-Takens, cusp, generalized Hopf
located automatically during continuation), branch switching at branch points,
and limit-cycle continuation by collocation. Exact equilibria via the CAS
(solve-poly / Groebner) also remain as the next CAS slice.
