#pragma once

/* ============================================================
 * dynsys analysis core  (roadmap phases 1-3)
 *
 * Dependency-free numerical analysis for dynamical systems:
 *
 *   - general real eigensolver (balance -> Hessenberg ->
 *     double-shift QR) returning complex eigenvalues for any
 *     N x N matrix;
 *   - dimension-agnostic equilibrium classification with
 *     fold / Hopf candidate flags;
 *   - pseudo-arclength continuation of equilibria with
 *     codimension-1 event (fold / Hopf) detection.
 *
 * Like expr_ir, this module knows nothing about AppState or
 * imgui. It talks to the rest of the program only through the
 * small callback structs below, so it can be unit-tested on its
 * own (see test/analysis_smoke.cpp) and reused headless.
 *
 * All matrices are row-major, size n*n, indexed J[row*n + col].
 * ============================================================ */

#include <complex>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>

namespace dynsys::analysis {

using Complex = std::complex<double>;

/* ---- linear algebra ----------------------------------------- */

/* Solve A x = b for x. A is row-major n*n. Gaussian elimination
 * with partial pivoting. Returns false on a singular pivot. */
bool solve_linear(const std::vector<double> &A, const std::vector<double> &b,
                  std::vector<double> *x);

/* Determinant of a row-major n*n matrix via LU with partial
 * pivoting. Sign-tracked; returns 0.0 on exact singular pivot. */
double determinant(const std::vector<double> &A, std::size_t n);

/* Eigenvalues of a general real n*n matrix (row-major). Real
 * matrices can have complex-conjugate eigenvalue pairs, hence the
 * complex result. Implementation: Parlett-Reinsch balancing, then
 * reduction to upper Hessenberg form, then the Francis double-shift
 * QR algorithm. Returns false only on dimension/|allocation issues;
 * non-convergence of a single eigenvalue is reported by leaving its
 * imaginary part finite but flagging via the return of
 * eigenvalues_converged(). */
bool eigenvalues(const std::vector<double> &A, std::size_t n,
                 std::vector<Complex> *out);

/* ---- equilibrium classification ----------------------------- */

struct Classification {
  std::string label;                 /* human-readable summary       */
  std::vector<Complex> eigenvalues;  /* spectrum of the Jacobian     */
  int n_stable = 0;                  /* eigenvalues with Re < -tol   */
  int n_unstable = 0;                /* eigenvalues with Re > +tol   */
  int n_center = 0;                  /* eigenvalues with |Re| <= tol */
  bool is_saddle = false;            /* both stable and unstable dirs*/
  bool hopf_candidate = false;       /* a complex pair near Re = 0   */
  bool fold_candidate = false;       /* a real eigenvalue near 0     */
};

/* Classify an equilibrium from its Jacobian. Works in any
 * dimension. `tol` is the threshold on Re(lambda) for calling an
 * eigenvalue marginal. */
Classification classify_equilibrium(const std::vector<double> &jacobian,
                                    std::size_t n, double tol = 1e-7);

/* ---- forward-mode AD over a user callback ------------------- */

/* The continuation engine and the AD Jacobian only need two
 * operations on the model, supplied by the caller (dynsys.cpp wires
 * these to eval_rhs and the IR). Keeping them as std::function means
 * analysis.cpp never includes the model headers.
 *
 * The "extended" state passed to these callbacks is the dynamical
 * state augmented with the single active continuation parameter as
 * the last entry: z = [x_0 ... x_{n-1}, p]. The residual has length
 * n (the vector field); the parameter is read from z[n]. */
struct Model {
  std::size_t n = 0;  /* number of state variables */

  /* f: write the n-vector field at (x, p) into `f_out` (length n).
   * x has length n, p is the scalar parameter. Returns false on a
   * numerical/eval error (message in `err`). */
  std::function<bool(const double *x, double p, double *f_out,
                     std::string *err)>
      vector_field;

  /* Optional exact Jacobian d f / d x (row-major n*n) at (x, p).
   * If null, the engine falls back to finite differences built on
   * vector_field. dynsys.cpp supplies a forward-mode-AD version. */
  std::function<bool(const double *x, double p, double *jac_out,
                     std::string *err)>
      jacobian_x;

  /* Optional exact d f / d p (length n) at (x, p). If null, finite
   * differences are used. */
  std::function<bool(const double *x, double p, double *dfdp_out,
                     std::string *err)>
      dfdp;
};

/* Build d f / d x by finite differences using only vector_field.
 * Public so dynsys.cpp can reuse it and so tests can exercise it. */
bool finite_diff_jacobian(const Model &m, const double *x, double p,
                          std::vector<double> *jac_out, std::string *err,
                          double eps = 1e-6);

/* ---- pseudo-arclength continuation of equilibria ------------ */

enum class SpecialPointKind { None, Fold, Hopf, EndOfBranch };

struct BranchPoint {
  double p = 0.0;               /* parameter value                  */
  std::vector<double> x;        /* equilibrium coordinates (len n)  */
  std::vector<Complex> eigenvalues;
  int n_unstable = 0;           /* unstable eigenvalue count        */
  bool stable = false;          /* n_unstable == 0                  */
  SpecialPointKind special = SpecialPointKind::None;
};

struct ContinuationSettings {
  double h0 = 1e-2;             /* initial arclength step           */
  double h_min = 1e-6;
  double h_max = 1.0;
  int max_points = 600;        /* cap on accepted points           */
  int max_corrector_iters = 12;
  double corrector_tol = 1e-9;
  double p_min = -1e9;         /* stop if parameter leaves range   */
  double p_max = 1e9;
  double event_tol = 1e-9;     /* bisection tolerance for events   */
  int direction = +1;          /* +1 or -1: which way along tangent */
  bool detect_fold = true;
  bool detect_hopf = true;
};

struct Branch {
  std::vector<BranchPoint> points;
  std::vector<std::size_t> special_indices; /* into points[] */
  std::string message;                      /* status / why it stopped */
  bool ok = false;
};

/* Continue an equilibrium branch starting from (x0, p0), which
 * should already be an approximate equilibrium (the engine corrects
 * it first). The continuation parameter is the scalar `p`. */
Branch continue_equilibrium(const Model &m, const std::vector<double> &x0,
                            double p0, const ContinuationSettings &settings);

/* ---- 2D fixed-point scanning (pplane-style) ----------------- *
 * Find ALL equilibria of a planar vector field inside a rectangle by
 * launching Newton from a grid of seeds, deduplicating the converged
 * roots, and classifying each. AppState-free: the field is supplied as
 * a callback f(x,y) -> (u,v). */
struct PlanarField {
  /* write (u, v) = f(x, y); return false on eval error. */
  std::function<bool(double x, double y, double *u, double *v)> eval;
};

struct FixedPoint2D {
  double x = 0.0, y = 0.0;
  std::vector<double> jacobian;      /* 2x2 row-major */
  std::vector<Complex> eigenvalues;  /* 2 entries */
  std::string label;                 /* saddle / stable node / ... */
  bool is_saddle = false;
  /* Real eigenvector directions (unit), present when the eigenvalues are
   * real; used to draw stable/unstable manifolds. dir[k] pairs with
   * eigenvalues[k]; empty when complex. */
  std::vector<std::pair<double, double>> directions;
};

/* Scan [xmin,xmax] x [ymin,ymax] with a seeds x seeds grid of Newton
 * starts. `dedup_tol` is in data units (points closer than this are the
 * same root). Returns the distinct classified equilibria found. */
std::vector<FixedPoint2D> scan_fixed_points_2d(const PlanarField &field,
                                               double xmin, double xmax,
                                               double ymin, double ymax,
                                               int seeds = 11,
                                               double dedup_tol_frac = 0.01);

/* ---- Lyapunov spectrum (Benettin / Gram-Schmidt) ----------- *
 * Integrate the system together with n perturbation vectors evolving
 * under the variational (linearized) flow; periodically QR-reorthonormali
 * ze and accumulate the log growth of each direction. The time-averaged
 * logs are the Lyapunov exponents (sorted descending). From them we get
 * the Kaplan-Yorke (Lyapunov) dimension — a real fractal dimension.
 *
 * Works for both ODEs and maps via the Model's vector_field/jacobian_x:
 *   - ODE: x' = f(x);   tangent vectors obey  v' = Df(x) v.
 *   - map: x_{k+1} = f(x_k);  tangent vectors  v_{k+1} = Df(x_k) v_k.
 * AppState-free: the caller supplies the Model and tells us map vs ODE.
 */
struct LyapunovResult {
  bool ok = false;
  std::vector<double> exponents;   /* sorted descending, length n */
  double kaplan_yorke = 0.0;       /* Lyapunov dimension */
  double sum = 0.0;                /* sum of exponents (phase-space contraction) */
  long steps_used = 0;
  std::string message;
};

struct LyapunovOptions {
  bool is_map = false;     /* true: discrete map; false: continuous ODE */
  double dt = 0.01;        /* integration step (ODE only) */
  long transient = 2000;   /* steps to settle onto the attractor first */
  long steps = 20000;      /* steps over which to average */
  long reorth_every = 1;   /* QR reorthonormalize every k steps (>=1) */
};

LyapunovResult lyapunov_spectrum(const Model &m, const std::vector<double> &x0,
                                 double p, const LyapunovOptions &opt);

/* Kaplan-Yorke dimension from an already-sorted-descending spectrum.
 * D_KY = j + (sum_{i<=j} lambda_i) / |lambda_{j+1}|, where j is the
 * largest index with a non-negative partial sum. Returns 0 if the first
 * exponent is negative (point attractor), n if the whole sum is >= 0. */
double kaplan_yorke_dimension(const std::vector<double> &sorted_desc);

/* ---- Basins of attraction ----------------------------------- *
 * For a planar system, integrate from a grid of initial conditions and
 * color each cell by which attractor its trajectory settles onto. The
 * attractors are discovered automatically by clustering the trajectory
 * endpoints (so fixed points and periodic attractors are both handled).
 * The result is the basin-of-attraction picture — a fractal whenever the
 * basin boundaries are fractal (e.g. Newton's method, forced oscillators).
 *
 * A Julia set is the special case where the "attractor" is infinity:
 * cells whose orbit escapes vs. stay bounded.
 *
 * AppState-free and headless-testable: the caller supplies a stepping
 * function advance(x,y) -> (x',y') that performs one integration step (or
 * one map iteration).
 */
struct BasinResult {
  bool ok = false;
  int width = 0, height = 0;
  std::vector<int> cell_attractor;   /* width*height; index into attractors, or -1 = diverged, -2 = did not settle (chaotic) */
  std::vector<float> cell_speed;     /* width*height; 0..1 convergence speed (1 = fast) */
  std::vector<std::pair<double,double>> attractors; /* representative (x,y) of each basin */
  long n_converged = 0, n_diverged = 0, n_nonconvergent = 0;
  std::string message;
};

struct BasinOptions {
  double xmin = -2, xmax = 2, ymin = -2, ymax = 2;
  int width = 200, height = 200;
  long max_steps = 2000;       /* integration steps per cell */
  double settle_tol = 1e-4;    /* movement below this = settled */
  double cluster_tol = 1e-2;   /* endpoints within this = same attractor */
  double diverge_r = 1e6;      /* |state| beyond this = diverged */
  int max_attractors = 16;     /* cap distinct basins */
};

/* advance: one step of the dynamics, (x,y) -> (*nx,*ny); return false on
 * error. For an ODE pass a small fixed-step integrator; for a map, one
 * iteration. */
BasinResult compute_basins(const std::function<bool(double x, double y, double *nx, double *ny)> &advance,
                           const BasinOptions &opt);

/* ---- Box-counting fractal dimension --------------------------- *
 * Estimate the box-counting (Minkowski–Bouligand) dimension of a set of
 * 2D points: cover the bounding box with a grid of boxes of side eps,
 * count how many boxes contain at least one point (N(eps)), and fit the
 * slope of log N(eps) vs log(1/eps) over a range of scales. For a curve
 * D≈1, for an area-filling set D≈2, for a strange attractor a fraction in
 * between (e.g. Hénon ≈1.26). AppState-free and headless-testable. */
struct BoxCountResult {
  bool ok = false;
  double dimension = 0.0;            /* fitted slope */
  double r_squared = 0.0;            /* fit quality 0..1 */
  std::vector<double> log_inv_eps;   /* x of the fit (log 1/eps) */
  std::vector<double> log_count;     /* y of the fit (log N)     */
  std::string message;
};

/* points: interleaved or as pairs? We take separate x and y spans for
 * generality. n_levels grid refinements (each halving eps). The coarsest
 * eps spans the whole bounding box; the finest is box/2^(n_levels-1). */
BoxCountResult box_counting_dimension(const std::vector<double> &xs,
                                      const std::vector<double> &ys,
                                      int n_levels = 10);

}  // namespace dynsys::analysis
