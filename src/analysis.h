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

}  // namespace dynsys::analysis
