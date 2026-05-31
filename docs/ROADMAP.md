# dynsys roadmap

A plan for growing dynsys from a phase-plane/attractor visualizer into a
pplane + matcont competitor, and then beyond into fractals — with the
lizard Lisp + a Maxima-style CAS feeding symbolic power back in.

Ordering principle: **numeric first (no symbolic dependency), then
symbolic once lizard/CAS land.** Each item is tagged with what it needs.

---

## Done (shipped)

- Full-window background plot (2D phase or 3D), top toolbar + side panel.
- 2D phase portrait: speed-colored vector field, marching-squares
  nullclines WITH flow-direction arrows, auto-scan + classification of
  ALL equilibria in view, stable/unstable eigendirection manifolds,
  separatrices.
- Standard pplane mouse scheme; orbits accumulate; custom orbit by exact
  typed initial conditions; seed grid / seed circle.
- Keyboard + mouse zoom about cursor/center; sanitized bounds (no NaN,
  no blank view); view2d directive + framed presets.
- 3D scene rendered to an FBO; dimension auto-detect + Force2D/Force3D.
- Maps vs ODEs handled distinctly.
- Numerical analysis core (headless, tested): N-D real eigensolver,
  equilibrium classification, pseudo-arclength continuation with
  fold/Hopf detection, forward-mode AD over the expression IR.
- Orbit (bifurcation) diagram with Lyapunov-vs-parameter overlay;
  Tent map; live largest-Lyapunov estimate.

---

## Phase 7 - richer numerics (no lizard needed)

- Bifurcation diagram polish: zoom/pan in the diagram, click a parameter
  value to jump the live system there, export to CSV/PNG, log scale,
  multi-observable overlays.
- Full Lyapunov spectrum (all exponents) via Benettin/Gram-Schmidt on the
  variational flow — we already have AD for exact Jacobian-vector
  products. Gives the Kaplan-Yorke (Lyapunov) dimension: the first real
  "fractal number."
- Streamlines / line-integral-convolution for dense smooth ODE portraits.
- Basins of attraction: integrate a grid of ICs, color each pixel by
  which attractor it reaches -> basin fractals directly.
- Poincare upgrades: first-return maps, interactive section placement,
  return-time coloring.
- Better integrators: adaptive RK45 (Dormand-Prince), symplectic
  integrators, event detection.
- Two-parameter scans: color a (p1,p2) grid by period / largest Lyapunov
  -> the "shrimp"/period map (another fractal).

## Phase 8 - fractals as first-class objects (mostly no lizard)

- Escape-time fractals: Mandelbrot/Julia for z->z^2+c, but the map can be
  ANY user expression compiled through the existing IR — so the fractal
  explorer speaks the same language as the simulator. Deep-zoom with
  perturbation theory later.
- IFS / chaos game: affine iterated function systems (Barnsley fern,
  Sierpinski), chaos-game renderer.
- Fractal dimension: box-counting and correlation dimension on
  attractors/basins, reported next to Lyapunov dimension.
- Newton fractals for f(z)=0 (symbolic f, f' once the CAS exists).

## Phase 9 - matcont-grade continuation UI (engine exists; UI + symbolic Jacobians needed for full power)

- Equilibrium continuation diagram: branches in (parameter, state) space,
  fold/Hopf/branch markers, click a located equilibrium to start
  continuation. (Numeric now; better with symbolic Jacobians.)
- Limit-cycle continuation and period-doubling cascades.
- Two-parameter continuation of fold/Hopf curves; codim-2 points
  (Bogdanov-Takens, cusp, generalized Hopf). Where matcont lives; wants
  symbolic 2nd/3rd derivatives -> needs CAS.
- Normal-form coefficients (e.g. Hopf first Lyapunov coefficient):
  essentially requires the CAS.

---

## Phase 10 - lizard + CAS integration (the symbolic leap)

Once lizard (Lisp) and the Maxima-style CAS on top of it exist, wire them
in. This is the multiplier that makes the above exact and robust.

- Symbolic Jacobians/Hessians: replace finite-difference and AD
  derivatives with exact symbolic ones, for the phase plane and
  continuation. More accurate fold/Hopf/codim-2 detection.
- Exact equilibria: solve f(x)=0 symbolically where possible (polynomial
  systems via resultants/Groebner) -> find ALL equilibria with certainty,
  including unstable ones Newton may miss.
- Symbolic nullclines: where a component can be solved for one variable,
  draw the nullcline analytically (no sampling error) and label branches.
- Normal forms & center-manifold reduction symbolically (the matcont
  killer feature).
- Symbolic Newton fractals: f and f' from a typed expression.
- A shared expression language: feed the IR FROM lizard/CAS expressions,
  so the same symbolic object is analyzed, differentiated, and compiled
  to the fast numeric VM. One language for modeling, analysis, rendering.

### Integration architecture sketch
- Fast path: lizard/CAS simplifies expressions -> lower to expr_ir VM for
  the hot numeric loops.
- Off the hot path: derivatives, solves, normal forms computed once per
  system/parameter change, cached.
- A thin cas_bridge module translating CAS S-expressions <-> tpcas AST
  <-> expr_ir, so the simulator never needs to know the backend.

---

## Cross-cutting / quality
- Save & load sessions (system + view + orbits).
- PNG/SVG export of any view.
- Larger library of canonical systems.
- Multithread the bifurcation/basin/Lyapunov sweeps.
- Keep growing the headless-testable core.
