# dynsys — codim-2 ZH/HH detection, plus a Lin's-method scaffold

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell). Verified by new regression tests
and by building/running the real GUI.

## 1. Zero-Hopf (fold-Hopf) and Hopf-Hopf (double-Hopf) codim-2 points  [NEW]
Two more codim-2 singularities are now detected along a Hopf continuation curve
(3+ D systems), completing the planar codim-2 set alongside BT/GH/cusp:
  - ZERO-HOPF (ZH): a real eigenvalue passes through zero while the
    pure-imaginary Hopf pair persists (eigenvalues {0, +-i w}). Signed test =
    that real eigenvalue; detected by its sign change (or grazing zero where the
    curve terminates).
  - HOPF-HOPF (HH): a SECOND complex pair reaches the imaginary axis (two
    distinct frequencies neutral at once). Signed test = the real part of the
    nearest other pair.
New SpecialPointKind::ZeroHopf and ::HopfHopf; both are bisection-refined like
the other codim-2 points and named in the Continuation-view readout.
  - VALIDATED (new test zhhh_smoke): a 3-D fold-Hopf system (oscillator + real
    mode) detects ZH at the analytic (p,q)=(0,0); a 4-D double-oscillator
    detects HH at (0,0). Both PASS.

## 2. Lin's method for homoclinic location  [EXPERIMENTAL SCAFFOLD]
A new analysis::lin_homoclinic lays the groundwork for locating homoclinics
that integration-seeding cannot reach (the Bogdanov-Takens loop). It relocates
the saddle at each parameter, computes the dominant stable/unstable directions,
integrates each manifold to a Poincare section (trying both orientations,
rejecting the escaping branch), measures an in-section gap, and runs a secant
root-find on the primary parameter.
  - HONEST STATUS: this is PRELIMINARY and not yet reliable for the stiff BT
    loop. The orbits and gap are now real (amplitudes ~1-3, not noise), but
    obtaining a SIGNED gap that brackets the homoclinic -- so the root-find
    converges -- needs careful section placement where BOTH manifolds genuinely
    cross; the backward (stable) manifold does not always reach a section built
    from the forward (unstable) excursion. It compiles cleanly, runs, and is
    documented in the header as experimental. The VALIDATED homoclinic path
    remains solve_homoclinic / continue_homoclinic.

## Everything else in this list was already shipped and is unchanged here:
branch points + branch switching, limit-cycle continuation (collocation +
pseudo-arclength + adaptive mesh), LPC fold-of-cycles curve, BT/GH/cusp codim-2
(BT a,b normal form; cusp cubic c exact; GH l2 sign), Floquet multipliers, and
homoclinic continuation -- all retained with their regression tests.

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds/runs.
- New zhhh_smoke passes; codim2_coeffs_smoke, homoclinic_seed_smoke,
  homoclinic_smoke, homoclinic_cont_smoke, validation_smoke, basins_mt_smoke,
  lpc_arclength_smoke, bt_codim2_smoke all pass; CAS green.
- Continuation view runs headless without crashing.
