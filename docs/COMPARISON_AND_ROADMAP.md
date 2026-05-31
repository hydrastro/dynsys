# dynsys — roadmap status & comparison (May 2026)

*Grounded in the current build: 41 presets, 9 full-window views, 7
integrators, 16 passing headless test suites, dependency-free PNG export.*

---

## Roadmap status at a glance

| Phase | Theme | Status |
|---|---|---|
| **A** | View system, robustness, rendering | DONE |
| **B** | Richer numerics (Lyapunov, basins, scans, solvers) | DONE (minor polish left) |
| **C** | Fractals as first-class objects | MOSTLY DONE |
| **D** | MatCont-grade continuation | STEP 1 DONE (equilibrium branch + fold/Hopf) |
| **E** | Lizard + CAS symbolic layer | BLOCKED ON LIZARD (in progress by you) |

### Phase A — DONE
9-view architecture, robust auto-fit, crash/freeze hardening, debounced +
progressive grid rendering.

### Phase B — DONE (polish remaining)
Full Lyapunov spectrum + Kaplan-Yorke; basins of attraction (with
chaotic/non-convergent handling); 2-parameter Lyapunov scans; full
integrator set incl. adaptive RKF45/Dormand-Prince.
Polish left: streamlines/direction-field arrows; box-counting &
correlation dimension.

### Phase C — MOSTLY DONE
Escape-time Mandelbrot/Julia (progressive render); Mandelbrot<->bifurcation
3D bridge; density-rendered bifurcation diagrams (now with a contrast lift
so the period-doubling tree is clearly visible); PNG export of any view.
DONE: box-counting fractal dimension (validated: square~2, Sierpinski~1.585,
Henon~1.26). Left: IFS / chaos game (Barnsley fern, Sierpinski); optional cellular-automata / life-like module (independent of
Lizard).

### Phase D — STEP 1 DONE
Equilibrium continuation view: traces the equilibrium branch both ways
vs. a parameter, stable (solid green) vs. unstable (dashed red), with Fold
(LP) and Hopf (H) markers. Engine validated (test-continuation) on the
fold normal form x'=p-x^2 and the Hopf normal form.
Left: branch switching at bifurcations; limit-cycle continuation
(needs a BVP collocation solver - the big one); codim-2 (cusp, BT, GH).

### Phase E — BLOCKED ON LIZARD
When Lizard's CAS is ready: cas_bridge (CAS S-expr <-> tpcas AST <-> expr_ir
VM), then exact symbolic Jacobians/Hessians -> exact equilibria -> analytic
nullclines -> normal-form coefficients. This is the differentiator no
competitor has.

---

## Comparison table vs. the established tools

Legend: 
- full
- (partial)
- (none)

| Capability | pplane | XPPAUT | MatCont | AUTO | dynsys |
|---|:--:|:--:|:--:|:--:|:--:|
| MODELLING | | | | | |
| Arbitrary N-D ODEs | no | full | full | full | full |
| Discrete maps | no | full | partial | no | full |
| Symbolic expression input | partial | full | full | partial | full |
| Adaptive integrators | no | full | full | full | full |
| Built-in example library | partial | partial | partial | partial | full (41) |
| PHASE-PLANE / GEOMETRY | | | | | |
| 2D phase portrait | full | full | partial | no | full |
| Nullclines | full | full | no | no | full |
| Equilibria + stability (N-D) | full | partial | full | full | full |
| 1D flow / cobweb | no | partial | no | no | full |
| 3D phase view | no | partial | partial | no | full |
| NONLINEAR DYNAMICS / CHAOS | | | | | |
| Bifurcation diagram (sweep) | no | partial | full | full | full |
| Largest Lyapunov | no | partial | no | no | full |
| Full Lyapunov spectrum | no | no | no | no | full |
| Kaplan-Yorke dimension | no | no | no | no | full |
| Box-counting dimension | no | no | no | no | full |
| Basins of attraction | no | partial | no | no | full |
| 2-parameter scan (shrimps) | no | no | partial | partial | full |
| Escape-time fractals | no | no | no | no | full |
| CONTINUATION | | | | | |
| Equilibrium continuation | no | partial | full | full | partial (view shipped) |
| Fold / Hopf detection | no | partial | full | full | full |
| Limit-cycle continuation | no | partial | full | full | no |
| Codim-2 points | no | no | full | full | no |
| Normal-form coefficients | no | no | full | partial | no (planned Phase E) |
| SYMBOLIC / EXACTNESS | | | | | |
| Jacobian | numeric | numeric | symbolic | numeric | exact (AD) |
| Exact equilibria / nullclines | no | no | no | no | no (planned Phase E) |
| ENGINEERING | | | | | |
| Runtime dependency | MATLAB/Java | none | MATLAB | none | none |
| Modern GUI | no | partial | partial | no | full |
| Image (PNG) export | partial | partial | full | partial | full |
| Headless test suite | no | no | partial | partial | full (16) |
| Open source | partial | full | full | full | full |

---

## Where this leaves dynsys

- Ahead of everyone on the ergodic/chaos side (full Lyapunov spectrum,
  Kaplan-Yorke, basins, 2-parameter scans, escape-time fractals) and on
  packaging (one native binary, no MATLAB, modern GUI, a real test suite).
- Exact AD Jacobians beat the finite-difference Jacobians in XPPAUT/AUTO.
- One subsystem from parity with MatCont: limit-cycle continuation +
  codim-2. The equilibrium half is shipped.
- Phase E is the moat: a live CAS feeding exact symbolic analysis is
  something none of pplane/XPPAUT/MatCont/AUTO offer.

Net: surpassed pplane; at/above XPPAUT; leading on chaos; equilibrium
continuation shipped, limit-cycle continuation is the remaining gap to
MatCont; the symbolic layer (Phase E, via Lizard) is the long-term
differentiator.
