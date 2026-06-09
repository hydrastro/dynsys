# dynsys — cusp + Hopf-Hopf curve continuation (codim-2 curves COMPLETE)

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell).

## Cusp and Hopf-Hopf curve continuation in 3 parameters  [NEW]
Previous iterations added continuation of Bogdanov-Takens (bt_curve) and
zero-Hopf (zh_curve) codim-2 points as curves. This iteration adds the last two
planar codim-2 loci, so ALL FOUR now continue as curves in three parameters:

- cusp_curve: continues the CUSP locus -- fold (one real eigenvalue = 0) AND the
  fold normal-form coefficient a = 0 -- over the unknowns (x,p,q,r).
- hh_curve: continues the HOPF-HOPF locus -- TWO distinct complex-conjugate
  pairs simultaneously on the imaginary axis (both Re = 0).

Both are built on a shared pseudo-arclength codim-2 curve tracer
(trace_codim2_curve) refactored out of bt_curve/zh_curve, so all four
continuations share one validated predictor/corrector core, parameterized only
by their two defining conditions and a per-point normal-form recorder.

### Validation (new tests cusp_curve_smoke, hh_curve_smoke)
On systems with ANALYTIC codim-2 loci:
- cusp: x' = p + (q + 0.5 r) x - x^3 has its cusp locus exactly { p=0, q=-0.5 r }.
  cusp_curve traces 91 points on it, cusp coefficient c=-1 recovered all along,
  r in [-4.0, 4.0], both defining conditions to 2e-10.
- Hopf-Hopf: two oscillators with dampings p and (q+0.4 r) and frequencies 1 and
  2.3 give the locus { p=0, q=-0.4 r }. hh_curve traces 89 points on it, both
  frequencies (1.0, 2.3) recovered, r in [-4.1, 4.1], conditions to 4e-10.

## Refreshed MATCONT_COMPARISON.md
Codim-2 curve continuation is now marked DONE for all four planar loci
(BT/ZH/cusp/HH). The remaining codim-2 sub-gap is narrowed to codim-2 CYCLE
points (on limit-cycle bifurcation curves) and continuing those as curves.

## Final checkup (this iteration)
- All 4 C++ TUs compile ZERO warnings under -Wall -Wextra; brace check passes.
- The REAL GLFW/OpenGL binary builds (~12.6 MB).
- CAS green: === PASS ===.
- All codim-2 curve + locator tests pass: bt_curve, zh_curve, cusp_curve,
  hh_curve, bt_locate.
- Full `make test`: 0 failures (see below).
