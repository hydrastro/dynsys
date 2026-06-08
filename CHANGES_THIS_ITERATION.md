# dynsys — three fronts: 2-parameter homoclinic continuation, parallel fractals, analytic validation suite

Unzip at repo root, `make clean && make && make run`. CAS features need
LIZARD + SANGAKU_ROOT (or the Nix dev shell). Verified by new regression tests
against analytic results, forced-concurrency pixel-identity checks, and by
building/running the real GUI.

## 1. Two-parameter homoclinic CONTINUATION (completes HomCont)
The previous release solved a single homoclinic orbit (truncated BVP +
projection BCs). A homoclinic connection is codim-1, so in a (p,q) plane it
traces a CURVE. `analysis::continue_homoclinic` now follows it: it steps the
secondary parameter q and, at each q, finds the primary parameter p where the
connection closes (BVP residual -> 0), re-using the previous converged orbit as
the seed (natural/warm-started continuation), tracing both directions.
  - In the GUI (Continuation view): after "Find homoclinic orbit", a "Continue
    homoclinic curve" button traces the locus in (cont_param vs twopar_p2) and
    plots it (primary param vs secondary).
  - VALIDATED (new test homoclinic_cont_smoke): on x'=y, y'=q*x-x^2+p*y the
    traced locus matches the conservative line p=0 with peak amplitude 1.5*q at
    every q across 25 points, residuals ~1e-12.
  - Honest limit: warm-started natural continuation following p(q); it relies on
    a good starting homoclinic (the inner solve's seeding). Very stiff / large
    BT-type loops need a better seed than the simple GUI excursion provides.

## 2. Parallel FRACTAL / escape-time rendering (multi-core)
After the basins parallelisation, the escape-time fractal (the other hotspot) is
now multi-threaded too. Each worker thread owns a private ThreadStepper (its own
IR eval scratch AND its own parameter snapshot), so PARAMETER-space fractals --
which vary parameters per pixel -- are now thread-safe (each thread overrides
parameters only in its private copy). Rows are distributed across
hardware_concurrency() threads; the AST-fallback evaluator keeps the serial
path.
  - VERIFIED bit-identical: with the parallel path FORCED on (sandbox is single
    core), the rendered Mandelbrot is pixel-for-pixel identical to the serial
    render (max pixel difference 0 across 1.28M pixels).
  - Speedup is realised on multi-core machines; correctness (identical output)
    is what was checked here.

## 3. Analytic VALIDATION suite (head-to-head ground truth)
A new `validation_smoke` test runs dynsys's analysis core against systems with
KNOWN analytic answers -- the kind of comparison one runs against MatCont:
  - Lorenz origin eigenvalues: 11.8277 and -22.8277 (exact). PASS.
  - Supercritical Hopf normal form: first Lyapunov coefficient l1 < 0 (correctly
    supercritical), frequency omega = 1. PASS.
  - Van der Pol (mu=1): largest Lyapunov exponent ~ 0 (on the limit cycle). PASS.
  - Homoclinic x'=y,y'=x-x^2: peak amplitude 1.4996 (analytic 1.5). PASS.
  6/6 checks pass. (If you supply a specific system + parameters, the same kind
  of check can be run head-to-head against your MatCont output.)

## Verification
- All 4 C++ TUs compile ZERO warnings; the REAL GLFW/OpenGL binary builds/runs.
- Full suite passes (incl. homoclinic_cont_smoke, validation_smoke,
  homoclinic_smoke, basins_mt_smoke, bt_codim2_smoke, branch_switch_smoke,
  lpc_arclength_smoke, bridge_family_smoke); CAS green.
- Parallel fractal proven pixel-identical to serial under forced concurrency;
  fractal + continuation views run headless without crashing.
