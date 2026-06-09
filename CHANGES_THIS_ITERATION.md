# dynsys — Bogdanov-Takens point location (final gap-closer) + full checkup

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell).

## Bogdanov-Takens point location by defining-system Newton  [NEW]
Previously dynsys DETECTED Bogdanov-Takens (BT) points along a fold/Hopf curve
(eigenvalue-test bracketing + bisection). It now also LOCATES a BT point
directly, the way MatCont's BT initialization does: Newton on the defining
system
    { f(x,p,q) = 0  (n eqns) ;  mu_prod = 0 ;  mu_sum = 0 }
where mu_prod is the product and mu_sum the sum of the real parts of the two
smallest-magnitude eigenvalues -- both vanish at the double-zero. This is a
square (n+2)x(n+2) Newton in the full (x, p, q) space with Levenberg damping,
pinning the equilibrium AND both parameters to machine precision. On success the
BT normal-form coefficients a,b are computed at the located point.

analysis::locate_bogdanov_takens(model2, x0, p0, q0) -> Codim2Point.
This is the prerequisite for the next rung (continuing codim-2 points as curves
in a third parameter).

### Validation (new test bt_locate_smoke)
On the Bogdanov-Takens normal form x'=y, y'=b1+b2*y+x^2+x*y -- whose BT point is
known EXACTLY at (b1,b2)=(0,0) with equilibrium (0,0) -- starting from a nearby
guess the locator converges in 8 iterations to the exact point: residual
2.7e-15, x=(0,0), b1=0, b2=0, and the BT normal form a=b=1.

## Refreshed MATCONT_COMPARISON.md
Marked BT point location as done; the codim-2 row is now "≈ parity on
coefficients + BT point location". The remaining continuation gap is narrowed to
continuing codim-2 points as curves in a third parameter.

## Final checkup (this iteration)
- All 4 C++ TUs (dynsys, expr_ir, analysis, expr_ir_ad) compile with ZERO
  warnings under -Wall -Wextra; brace check passes.
- The REAL GLFW/OpenGL binary builds (~12 MB) and renders the food-chain preset.
- CAS (exact symbolic Routh-Hurwitz path) green: === PASS ===.
- The full codim-2 suite passes: bt_locate, hopf_hopf_nf, zero_hopf_nf,
  codim2_coeffs, bt_codim2 -- all "all checks pass".
- Full `make test` run: 0 failures (see below).

## State at finish
dynsys now has the complete planar codim-2 normal-form set (BT, cusp, GH, ZH,
HH) with validated coefficients, direct BT point location, two-parameter
fold/Hopf and cycle-bifurcation (LPC/PD/NS) curves with codim-2 detection,
branch-point-of-cycles, full Floquet via the variational monodromy, arclength
limit-cycle continuation, homoclinic + heteroclinic solving/continuation + a
one-call homoclinic locator, plus the exact-symbolic (CAS) differentiator.
Remaining honest gaps: a fully general parameter-as-unknown homoclinic BVP,
continuation of codim-2 points as curves in a 3rd parameter, and 15 years of
robustness/track record.
