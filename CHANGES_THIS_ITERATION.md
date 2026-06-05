# dynsys — Phase E slice 1b: certified eigenvalues in the equilibrium readout
#           + a Sangaku CAS module that fills real gaps

Unzip at repo root, `make clean && make && make run`. With Nix, `nix develop`
provides the Lizard interpreter + Sangaku so the CAS features light up.

## 1. Certified eigenvalues wired into the fixed-point panel (slice 1b)
After "Find fixed point", when the CAS is reachable, a new button
**"Certified eigenvalues (CAS)"** appears in the Analysis tab. It:
- takes the Jacobian dynsys already computed at the equilibrium
  (app.fixed_jacobian, doubles),
- rationalizes each entry (continued fractions; flags whether every entry is
  exactly rational),
- runs Sangaku via Lizard and shows the EXACT eigenvalues, an exact
  Routh-Hurwitz stability verdict (stable / unstable / marginal), the
  Re<0 / Re>0 / Re=0 counts, and — when a pair sits precisely on the
  imaginary axis — an "exact Hopf point (Re = 0 exactly)" note.
The numeric spectrum above it is untouched and always stands alone; the CAS
section is purely additive and silently absent when LIZARD/SANGAKU_ROOT
aren't set. If the Jacobian isn't exactly rational, the result is clearly
labelled "approximate (nearest rational system)".

New in src/cas_bridge.h: `rationalize()` (double -> exact rational with a
tolerance + denominator bound) and `eigen_report_from_doubles()` (the
UI entry point). Verified live: Lorenz origin -> -8/3 and (-11±sqrt1201)/2,
1 RHP (saddle-type); Van der Pol origin -> (1±i sqrt3)/2, unstable; stable
nodes/focis classified correctly.

## 2. Sangaku CAS module (shipped separately as sangaku_dynsys_contrib.zip)
I audited Sangaku and implemented the biggest dynamical-systems gaps as a new
library module cas/dynsys.lisp (with a golden test + example 392). It adds
multivariate partial derivatives, the Jacobian of a polynomial vector field,
exact eigenvalues at an equilibrium, divergence, gradient, and Hessian — none
of which existed in Sangaku before. Verified on Lorenz (exact Jacobian,
eigenvalues, and the constant divergence -41/3). This means the WHOLE chain —
polynomial field -> exact symbolic Jacobian -> exact eigenvalues — can run
inside the certified CAS, not just dynsys-assembled matrices.

## 3. What's still missing in Sangaku
See docs/SANGAKU_GAPS.md: the remaining candidates (transcendental
multivariate certified derivatives, packaged eigenvectors, an in-CAS
half-plane/Hurwitz decision, a turnkey exact-equilibria wrapper, and the
certified-decimal realizer) with effort/value notes, plus the deliberate
non-gaps (FP linear algebra / numeric ODE — dynsys's job, not Sangaku's).

## Verification
- All 4 C++ TUs compile with ZERO warnings (-O2 -Wall -Wextra).
- make test: 25 groups + graceful CAS skip (26 with the CAS present); the new
  UI path is covered by the cas bridge tests (from-doubles cases).
