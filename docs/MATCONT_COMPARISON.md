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

## Where MatCont is clearly ahead

| Area | MatCont | dynsys | Honest gap |
|---|---|---|---|
| Limit-cycle continuation | Orthogonal collocation, **adaptive mesh + arclength**, Floquet multipliers | Collocation cycle solver, but **steps the parameter monotonically (no pseudo-arclength)** | dynsys can't go around folds of cycles robustly; MatCont can |
| Codim-1 curve continuation | Pseudo-arclength, handles folds in the curve, very long runs | Equilibrium branch + fold/Hopf detection; two-parameter fold/Hopf curves | MatCont's continuation is far more robust on hard curves |
| Codim-2 points | BT, GH, ZH, HH, cusp with validated normal-form coefficients | BT / cusp / generalized-Hopf detection + first normal-form coefficients | MatCont covers more types and is better validated |
| Homoclinic/heteroclinic | Dedicated continuation (HomCont) | Not available | genuine missing feature |
| Floquet / monodromy | Yes, full | Partial (period & amplitude; not full Floquet) | gap |
| Track record | 15+ years, thousands of papers | new | trust/validation |
| Ecosystem | MATLAB: plotting, export, scripting, community | self-contained, no scripting language | MATLAB integration |

If your work depends on any row above, MatCont (or AUTO) is the right tool, and
dynsys does not claim parity there.

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
| Branch points + branch switching | ● | ◐ |
| Limit-cycle continuation | ● (adaptive + arclength) | ◐ (monotone param, no arclength) |
| LPC (fold of cycles) curve | ● | ◐ |
| Codim-2 (BT, GH, cusp) | ● (+ ZH, HH) | ◐ (BT, GH, cusp; first coeffs) |
| Homoclinic continuation | ● | ○ |
| Floquet multipliers | ● | ◐ |
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

dynsys is **not** a MatCont replacement for serious continuation work — the
limit-cycle solver alone (monotone parameter, no pseudo-arclength) means it
cannot follow cycle folds the way MatCont does, and MatCont has far more codim-2
machinery and validation. Where dynsys is genuinely differentiated is the
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
