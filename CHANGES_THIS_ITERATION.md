# dynsys — Hopf-Hopf normal-form coefficients + refreshed MatCont comparison

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell).

## Hopf-Hopf (double-Hopf) normal-form coefficients  [NEW — completes the planar codim-2 set]
This was the last codim-2 type dynsys detected but did not quantify. With it,
dynsys now computes normal-form coefficients for ALL of MatCont's planar codim-2
points: BT (a,b), cusp (c), generalized-Hopf (l2), zero-Hopf (b, Re c), and now
Hopf-Hopf.

analysis::hopf_hopf_normal_form, at a point whose Jacobian has two pure-imaginary
pairs {+-i*omega1, +-i*omega2}, computes the four cubic coefficients of the
amplitude normal form
    r1' = r1 (mu1 + p11 r1^2 + p12 r2^2)
    r2' = r2 (mu2 + p21 r1^2 + p22 r2^2):
  p11, p22 = the two self-interaction (first-Lyapunov) coefficients,
  p12, p21 = the cross-coupling coefficients (mode i driven by mode j amplitude).
It uses the existing bilinear/trilinear toolkit (B_cplx, C_real) plus the
quadratic-correction vectors h20, h11, h(w1+w2), h(w1-w2), ... solved from the
resonant linear systems. The signs of p11,p22 and the product p12*p21 classify
the local diagram (whether the two cycles coexist -- the "simple" case -- or
compete, the "difficult" case where a 2-torus / heteroclinic structure can
appear).

### Validation (new test hopf_hopf_nf_smoke)
On two oscillators coupled at the origin with PRESCRIBED self-cubics (A1,A2) and
cross-couplings (G12,G21) and frequencies (W1,W2):
- omega1, omega2 recovered exactly,
- all four coefficients match the prescribed values up to a single common
  convention factor (exactly 2): p11/A1 = p22/A2 = p12/G12 = p21/G21 = 2,
- signs of all four are correct and the cross terms p12 != p21 are correctly
  DISTINGUISHED (the hard part of HH).

### Integration
Wired into the two-parameter fold/Hopf curve detection (HH points get their
coefficients computed on refinement) and the GUI prints them with the
simple/difficult-case interpretation; ZH and HH markers were added to the curve
plot in the previous iteration.

## Refreshed MATCONT_COMPARISON.md
The comparison table was stale (it still listed Floquet, pseudo-arclength cycle
continuation, and the cycle-bifurcation curves as gaps -- all since shipped). It
now honestly marks at "≈ parity": arclength LC continuation, full Floquet,
LPC/PD/NS curves, branch-point-of-cycles, codim-2-on-cycle-curves, the full
planar codim-2 normal-form set (now incl. ZH and HH), and homoclinic +
heteroclinic + find_homoclinic. Remaining gaps narrowed to: a fully general
parameter-as-unknown homoclinic BVP (production HomCont), continuation of codim-2
points as curves, and 15 years of robustness/track record.

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds + runs.
- New hopf_hopf_nf_smoke passes; zero_hopf, bt_codim2, codim2_coeffs all pass
  (the 2-param path is intact); CAS green.
