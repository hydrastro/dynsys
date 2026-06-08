# dynsys — branch point of cycles (BPC) detection

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell).

## Branch point of cycles (BPC)  [NEW — roadmap item B]
continue_limit_cycle now DETECTS branch points of cycles -- parameter values
where two periodic-orbit branches cross (e.g. a symmetry-breaking / transcritical
bifurcation of cycles, the multiplier signature being a SECOND Floquet multiplier
passing through +1 while the branch does not fold).

Method: a signed bordered-determinant test function (cycle_bp_test). The
extended cycle Jacobian dF/d(U,p) (N rows, N+1 unknowns) is bordered with the
branch tangent as the last row; the determinant of that square matrix changes
sign at a branch point. Two numerical points were essential:
  - the raw determinant of a ~200x200 matrix underflows to ~1e-27 even when
    well-conditioned, so we return the SIGNED GEOMETRIC-MEAN determinant
    sign(det)*|det|^(1/N), which stays O(1) and preserves the sign;
  - the bordering tangent's overall sign is a gauge (the two-direction
    continuation negates it), so we fix the gauge (orient the tangent's
    parameter component non-negative) -- otherwise the sign flips spuriously.
Distinguishing BPC from LPC: a fold ALSO makes the cycle determinant vanish, so
a bare sign change is ambiguous. A genuine BPC has the branch passing straight
through (parameter locally monotonic) whereas a fold turns in the parameter; the
detector accepts a BPC bracket only where there is no p-turn, and the fold flag
itself was tightened to require a local extremum in p. Each CycleSample carries
bp_test and, when bracketed, is_bp / bp_p (the refined branch-point parameter).

### Validation (new test bpc_smoke), two cases:
  - van der Pol (a clean branch, NO branch point): bp_test stays one sign across
    the whole branch -> ZERO BPCs reported (no false positives).
  - A constructed transcritical-of-cycles  x'=-y+x(1-r^2), y'=x+y(1-r^2),
    z'=(mu - x^2) z  has its transverse Floquet multiplier pass through +1 at
    mu=1/2 (since <x^2> = 1/2 on the unit cycle). The detector brackets the BPC
    at mu=0.4950, matching the analytic 0.5 to within 0.005.

## Roadmap status (docs/ROADMAP_MATCONT_PARITY.md)
A heteroclinic connections (prior iteration) and B branch-point-of-cycles
(this iteration) are now DONE. Remaining: C codim-2 points on cycle-bifurcation
curves, D a fully robust Lin's method for the stiff Bogdanov-Takens homoclinic.

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds + runs.
- New bpc_smoke passes; lpc_arclength, lpccurve, pd_curve, lc_selfseed,
  heteroclinic all pass; CAS green.
