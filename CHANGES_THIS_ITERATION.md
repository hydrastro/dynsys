# dynsys — robust limit-cycle & homoclinic continuation; Lin's-method progress

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell). Verified by regression tests and
by building/running the real GUI.

## 1. Limit-cycle continuation: simulation-based self-seeding  [NEW, the headline]
continue_limit_cycle previously REQUIRED the caller to supply a near-converged
cycle: a crude guess (e.g. a circle) failed Newton for far-from-circular cycles
like the van der Pol relaxation oscillation. It now SELF-SEEDS: if the initial
collocation Newton fails, it integrates the ODE to settle onto the attracting
cycle, detects one period via a Poincare-section return, resamples that loop to
the mesh, and retries. Both the arclength and monotone paths use it.
  - VALIDATED (new test lc_selfseed_smoke): van der Pol continued from a crude
    CIRCLE across mu in [1.0, 2.67] (81 cycles); period grows correctly with mu
    (relaxation regime) and the trivial Floquet multiplier stays within 0.04 of
    1 on the stiffening orbit. Previously this guess failed outright.

## 2. Homoclinic continuation: tangent predictor  [IMPROVED]
continue_homoclinic now uses a TANGENT PREDICTOR (extrapolate the primary
parameter p along the curve from the last two accepted points) plus a robust
seed handoff and a one-retry fallback, instead of reusing the last p. This holds
the curve over a wide secondary-parameter range and tolerates gentle folds.
  - VALIDATED: on the conservative family x'=y, y'=q x - x^2 + p y it traces 36
    points across q in [0.25, 2.00] with p-drift 0 from the locus p=0 and
    amplitude error ~0.15%.

## 3. Lin's method for the BT homoclinic  [EXPERIMENTAL, improved but not reliable]
lin_homoclinic's test function was rebuilt: it now relocates the saddle at each
p, integrates the unstable manifold (choosing the bounded-return orientation),
finds the closest return to the saddle AFTER the excursion peak, and measures
the SIGNED distance along the stable left-eigenvector (covector normal to the
stable manifold). A coarse-scan + golden-section minimizer drives |gap| toward
its dip.
  - HONEST STATUS: the signed gap genuinely DIPS near the true Bogdanov-Takens
    homoclinic (visible in a parameter sweep), and the orbits are now real
    returning loops. But the automated locator does NOT robustly pin the stiff
    BT homoclinic: the one-sided manifold return leaves a residual the secant
    can't bracket, and amplitude/sign confounds make a clean minimum elusive.
    The routine now REPORTS "homoclinic not robustly located" and points to
    solve_homoclinic / continue_homoclinic, rather than emitting a false
    positive (a spurious tiny-loop "pass" was caught and removed). It is marked
    EXPERIMENTAL in the header. The genuinely unsolved piece is a two-sided Lin
    gap (forward unstable + backward stable on a shared section) for a true sign
    change -- deferred.

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds and
  runs the LimitCycle and Continuation views headless.
- New lc_selfseed_smoke passes; lc_colloc, lpc_arclength, lpccurve,
  branch_switch, homoclinic_cont all pass; CAS green.
