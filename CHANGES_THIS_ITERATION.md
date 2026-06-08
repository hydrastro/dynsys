# dynsys — codim-2 depth (cusp cubic + generalized-Hopf l2) and homoclinic seeding

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell). Verified by new regression tests
against analytic results and by building/running the real GUI.

## 1. Cusp normal-form coefficient c (codim-2)
On a fold curve, where the quadratic fold coefficient a -> 0, the cusp is
classified by the CUBIC coefficient c of the reduced dynamics y' = c y^3 + ...:
  c = (1/6) <p, C(q,q,q) - 3 B(q,h2)>,  A h2 = -(B(q,q) - <p,B(q,q)> q),
with A singular (Govaerts bordered solve). New analysis::cusp_normal_form.
  - EXACT on the normal form: x'=2x^3 -> c=2.0000, x'=-0.5x^3 -> c=-0.5000,
    2-D x'=x^3 (stable y) -> c=1.0000.

## 2. Generalized-Hopf (Bautin) second Lyapunov coefficient l2 (codim-2)
On a Hopf curve, where the first Lyapunov coefficient l1 -> 0, the Bautin point
is classified by l2. New analysis::gh_second_lyapunov implements Kuznetsov's
invariant Re<p,...> combination at 5th order on the critical eigenspace, with
new 4th/5th-order directional-derivative stencils, symmetric 4-/5-linear forms
(D, E) by polarization, and a complex multilinear evaluator. The intermediate
vectors h20,h11,h30,h21 are solved from (k i w I - A) systems (the singular
k=+/-1 case by a bordered solve).
  - Validated on the Bautin normal form x'=-y+L x r^4, y'=x+L y r^4 (l1=0): the
    SIGN of l2 is correct for L=+1, -1, +2 and the magnitude scales linearly.
  - HONEST NOTE: l2 leans on 4th/5th finite differences, so the MAGNITUDE
    carries a convention factor and is approximate. The SIGN -- which is what
    classifies the Bautin point and the direction the fold-of-cycles curve opens
    -- is the reliable output, and is reported in the GUI as super/subcritical.

Both coefficients are computed at the corresponding refined codim-2 points
during two-parameter continuation and printed in the Continuation-view readout
alongside the existing Bogdanov-Takens a,b.

## 3. Homoclinic seeding hardening
New analysis::seed_homoclinic_by_integration builds the homoclinic seed by
integrating the unstable manifold: it nudges off the saddle along the dominant
unstable eigenvector, tries BOTH orientations, integrates with RK4, and keeps
the excursion that returns nearest the saddle. This replaces the hand-built
inline seeding in the GUI homoclinic driver and is a far better initial guess.
  - On the unforced Duffing oscillator (genuine homoclinic, peak sqrt(2)) it
    seeds and the BVP converges in 2 Newton steps to residual ~1e-10.
  - HONEST LIMIT (unchanged, now better understood): this does NOT make the
    Bogdanov-Takens homoclinic converge. That is not a seeding bug -- a
    homoclinic is codimension-1, so at an APPROXIMATE BT parameter there is no
    exact connection to shadow, and the manifold trajectory runs to the time
    budget without returning. Closing that needs Lin's method (measure the gap,
    adjust a parameter to close it) -- a substantially larger build, deferred.

## Validation
- New regression tests: codim2_coeffs_smoke (cusp c exact + GH l2 sign) and
  homoclinic_seed_smoke (Duffing manifold seed + solve). Both in `make test`.

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds/runs.
- Full suite passes (incl. codim2_coeffs_smoke, homoclinic_seed_smoke,
  homoclinic_smoke, homoclinic_cont_smoke, validation_smoke, basins_mt_smoke);
  CAS green.
- Continuation + fractal views run headless without crashing.
