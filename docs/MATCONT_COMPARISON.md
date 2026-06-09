# dynsys vs MatCont — an honest head-to-head (June 2026)

MatCont is the de-facto standard MATLAB toolbox for numerical continuation and
bifurcation analysis of ODEs (and, via maps, of discrete systems). It is mature,
widely cited, and very deep. dynsys is a single self-contained native (C++/
OpenGL) application that also does interactive simulation, fractals, and a
proof-carrying computer-algebra layer. They overlap in the
continuation/bifurcation area; outside it they aim at different things.

This document is deliberately honest about where MatCont is stronger. The goal
is to tell you which tool to reach for, not to win an argument.

---

## TL;DR

- **Reach for MatCont** when you need production-grade continuation: limit-cycle
  continuation with adaptive mesh + arclength, codim-1 curves continued robustly
  over large ranges, the full codim-2 normal-form machinery (BT, GH, ZH, HH,
  cusp) with reliable normal-form coefficients, homoclinic/heteroclinic
  continuation, and decades of validation behind every branch.
- **Reach for dynsys** when you want an instant, dependency-free, interactive
  picture of a system — phase portraits, 3-D attractors, Lyapunov spectra,
  basins, fractals, the Mandelbrot↔logistic bridge — *and* when you want
  **exact, symbolically-derived** equilibria / Jacobians / stability
  certificates rather than purely floating-point ones. The symbolic layer is the
  thing MatCont does not have.

---

## Where MatCont is still ahead (and where dynsys has caught up)

The table below is kept current as dynsys closes gaps. "≈ parity" means the
capability exists and is validated against analytic ground truth in dynsys's test
suite, not that it matches MatCont's 15-year robustness on pathological inputs.

| Area | MatCont | dynsys | Status |
|---|---|---|---|
| Limit-cycle continuation | Orthogonal collocation, adaptive mesh + pseudo-arclength, Floquet | Collocation + **pseudo-arclength + adaptive mesh**, self-seeding; goes around folds of cycles | ≈ parity (MatCont more robust on stiff/long runs) |
| Floquet / monodromy | Full | **Full**: variational monodromy → multipliers, PD/NS classified & bracketed | ≈ parity |
| Codim-1 cycle-bif curves | LPC, PD, NS curves in two parameters | **LPC, PD, NS curves** (pd_curve / ns_curve / lpc_curve) | ≈ parity |
| Branch point of cycles (BPC) | Yes | **Yes** (bordered-determinant test, validated) | ≈ parity |
| Codim-2 on cycle curves | LPPD/fold-flip, PD-NS, cusp-of-cycles | **fold-flip / PD-NS / degenerate-PD** detection along PD/NS curves | ≈ parity (detection; MatCont also continues them) |
| Codim-2 equilibrium points | BT, cusp, GH, ZH, HH — all with normal-form coefficients | **BT (a,b), cusp (c), GH (l2), ZH (b,Re c), HH (p11,p12,p21,p22)** — all validated; **BT points also located directly by defining-system Newton** | ≈ parity on coefficients + BT point location |
| Connecting orbits | Homoclinic & heteroclinic continuation (HomCont), Lin's method | solve_homoclinic / continue_homoclinic; **solve_heteroclinic**; **find_homoclinic** (sweep+locate); Lin diagnostic | ◐ (MatCont's parameter-as-unknown HomCont BVP is more general) |
| Codim-1 equilibrium curve continuation | Pseudo-arclength, very robust on hard curves, restart/branch-switch | Equilibrium branch + fold/Hopf/branch detection + branch switching + 2-param fold/Hopf curves | ◐ (MatCont more robust on long/pathological curves) |
| Track record | 15+ years, thousands of papers | new | gap (trust/validation) |
| Ecosystem | MATLAB: plotting, export, scripting, community | self-contained native binary, no scripting language | gap (different niche) |

### The genuinely remaining gaps
1. **A fully general parameter-as-unknown homoclinic BVP** (production HomCont):
   dynsys locates homoclinics by sweeping the unstable-manifold return distance
   and polishing with the truncated BVP, and continues a known locus; it does
   not yet solve the orbit, period, and bifurcation parameter as one coupled
   BVP with an integral phase condition the way HomCont does. (The collapse
   degeneracies that defeat the naive version are documented in the roadmap.)
2. **Continuation of codim-2 points as curves in a third parameter.** dynsys
   now LOCATES a Bogdanov-Takens point directly (defining-system Newton in
   (x,p,q)) AND CONTINUES it as a curve when a third parameter is freed
   (bt_curve: pseudo-arclength on the BT defining system in (x,p,q,r),
   validated on a system with an analytic BT locus). The ZERO-HOPF (fold-Hopf)
   locus is likewise continued as a curve (zh_curve, validated on an analytic
   fold-Hopf locus). Continuation of the REMAINING codim-2 loci (cusp,
   Hopf-Hopf, Bautin/generalized-Hopf) as curves, and detection of higher
   codim-2 cycle points beyond the three implemented, remain.
3. **Robustness and track record.** MatCont has 15 years of hardening on
   pathological systems and a large validation corpus; dynsys validates each
   feature against analytic ground truth but is new.

If your work depends on rows still marked "gap" or "◐", MatCont (or AUTO) is the
right tool, and dynsys does not claim full parity there.

---

## Where dynsys offers something MatCont does not

1. **Proof-carrying / exact symbolic analysis (the differentiator).**
   dynsys has a computer-algebra layer (the Sangaku CAS, driven through the
   Lizard interpreter) wired into the analysis path. For polynomial/rational
   vector fields with rational coefficients it can:
   - form the Jacobian and characteristic polynomial **symbolically** (no finite
     differences), then apply an **exact Routh–Hurwitz** stability test — so a
     "stable" verdict is backed by exact arithmetic, not a floating-point
     eigenvalue that happened to land at Re = −1e−9;
   - solve for **exact equilibria** where the algebra permits;
   - produce analytic nullclines and exact normal-form ingredients.
   MatCont's pipeline is numeric throughout. For questions like "is this
   equilibrium *provably* stable for this exact parameter?", dynsys gives a
   certificate MatCont structurally cannot.

2. **Zero setup, no MATLAB.** A single native binary. No license, no toolbox
   install, no interpreter. Type a system, see it immediately.

3. **Interactive exploration first.** Real-time phase portrait / 3-D attractor
   with live parameter dragging, 7 integrators incl. adaptive RKF45 /
   Dormand–Prince, Lyapunov spectrum + Kaplan–Yorke dimension, basins of
   attraction, 2-parameter Lyapunov scans, Poincaré sections — all in-app and
   instant. MatCont is built around continuation, not live phase-space play.

4. **Fractals as first-class objects.** Escape-time Mandelbrot/Julia, the
   Newton-basin fractal, Burning Ship, Multibrot, the Buddhabrot
   (trajectory-density), an IFS/chaos-game gallery, and the unified 3-D
   **Mandelbrot ↔ logistic bifurcation bridge** (and its cubic/sine analogues).
   None of this is in MatCont's scope.

5. **One tool, many regimes.** ODE flows, discrete maps, IFS, fractals,
   continuation, and the CAS share one model description and one UI.

---

## Feature-by-feature (continuation/bifurcation core)

● full   ◐ partial / via workaround   ○ not available

| Capability | MatCont | dynsys |
|---|---|---|
| Equilibrium continuation (1 param) | ● | ● |
| Fold (LP) / Hopf detection | ● | ● |
| Exact eigenvalues + Routh–Hurwitz (symbolic) | ○ | ● |
| Exact equilibria (symbolic) | ○ | ◐ (when algebra permits) |
| Hopf first Lyapunov coefficient | ● | ● |
| Fold normal-form coefficient | ● | ● |
| Two-parameter fold / Hopf curves | ● | ● |
| Branch points + branch switching | ● | ● (detect + switch_branch, tested) |
| Limit-cycle continuation | ● (adaptive + arclength) | ● (collocation + pseudo-arclength + adaptive mesh; self-seeds by simulation) |
| LPC (fold of cycles) curve | ● | ● (arclength fold-of-cycles, tested) |
| Codim-2 (BT, GH, cusp, ZH, HH) | ● | ● (BT a,b; cusp cubic c exact; GH l2 sign; ZH + HH detected) |
| Homoclinic continuation | ● | ◐ (truncated-BVP locus, tangent predictor; Lin's method experimental) |
| Floquet multipliers | ● | ● (variational monodromy; PD/NS classified) |
| Period / amplitude vs parameter | ● | ● |
| Live phase portrait / 3-D attractor | ○ | ● |
| Lyapunov spectrum + dimension | ○ | ● |
| Basins of attraction | ○ | ● |
| Fractals / Buddhabrot / IFS | ○ | ● |
| Symbolic/proof-carrying certificates | ○ | ● |
| No MATLAB / zero install | ○ | ● |
| Scriptable in a host language | ● (MATLAB) | ○ |

---

## How to read this

dynsys now covers most of MatCont's single- and two-parameter continuation
core: equilibrium and cycle continuation (collocation + pseudo-arclength +
adaptive mesh), LPC fold-of-cycles, branch switching, Floquet multipliers via
the variational monodromy, the full planar codim-2 set (BT, GH, cusp, zero-Hopf,
Hopf-Hopf) with validated normal-form coefficients, the two-parameter
cycle-bifurcation curves (period-doubling and Neimark-Sacker/torus) with codim-2
detection along them, branch-point-of-cycles, and homoclinic AND heteroclinic
solving/continuation plus a one-call homoclinic locator. The remaining gaps
versus MatCont are a fully general parameter-as-unknown homoclinic BVP
(production-grade HomCont/Lin's method for stiff connections), continuation of
codim-2 points as curves, and — most of all — 15 years of robustness hardening
and validation track record. Where dynsys is genuinely differentiated is the
**exact symbolic / proof-carrying** path and the breadth of *interactive*
dynamics + fractals in one zero-setup tool.

A fair one-liner: *MatCont is the specialist's continuation engine; dynsys is an
interactive dynamics workbench that adds exact symbolic certificates and a
fractal/3-D side MatCont doesn't attempt.*

---

## Validation note

dynsys's continuation/normal-form engines are checked by headless regression
tests against known analytic cases (fold normal form x' = p − x²; supercritical
Hopf amplitude ∝ 2√µ; van der Pol period growth with µ; box-counting dimensions
square ≈ 2, Sierpinski ≈ 1.585, Hénon ≈ 1.26). These are sanity checks against
closed-form answers, not a substitute for MatCont's depth of validation. If you
give me a specific system, I can run it through dynsys and we can compare the
numbers directly against a MatCont run on the same system.
