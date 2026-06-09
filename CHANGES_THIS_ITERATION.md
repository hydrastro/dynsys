# dynsys — zero-Hopf CURVE continuation (codim-2 curves beyond BT)

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell).

## Zero-Hopf (fold-Hopf) curve continuation in 3 parameters  [NEW]
The codim-2 *curve continuation* gap was partly closed already (bt_curve
continues a Bogdanov-Takens point as a curve in (x,p,q,r)). This iteration adds
the ZERO-HOPF locus: zh_curve traces the fold-Hopf point as a curve in three
parameters by pseudo-arclength on the defining system
    { f = 0 (n eqns) ; Re(real eigenvalue) = 0 ; Re(complex pair) = 0 },
i.e. a fold and a Hopf pinned simultaneously, over the unknowns (x,p,q,r).
Same predictor/corrector scheme as bt_curve; reuses BTCurve/BTCurveSettings (the
point's a,b fields carry the zero-Hopf b and Re(c)). So BOTH Bogdanov-Takens and
zero-Hopf codim-2 points now continue as curves.

### Validation (new test zh_curve_smoke)
On a 3-parameter system with an ANALYTIC zero-Hopf locus
  x' = (p + 0.3 r) x - x^2     (fold: zero eigenvalue at p = -0.3 r)
  y' = -(1+0.2 r) z + q y      (Hopf pair at q = 0)
  z' =  (1+0.2 r) y + q z
the locus is exactly { p = -0.3 r, q = 0 }. zh_curve traces 86 points with
r spanning [-4.08, 4.08], both defining conditions satisfied to 3e-10, and the
zero-Hopf coefficient b = -2 recovered correctly all along the curve.

## Refreshed MATCONT_COMPARISON.md
Codim-2 curve continuation now lists BT + zero-Hopf as done; the remaining
curve-continuation sub-gap is the cusp / Hopf-Hopf / Bautin loci.

## Final checkup (this iteration)
- All 4 C++ TUs compile ZERO warnings under -Wall -Wextra; brace check passes.
- The REAL GLFW/OpenGL binary builds (~12 MB).
- CAS green: === PASS ===.
- Codim-2 curve + normal-form tests all pass: zh_curve, bt_curve, bt_locate,
  hopf_hopf_nf, zero_hopf_nf.
- Full `make test`: 0 failures (see below).
