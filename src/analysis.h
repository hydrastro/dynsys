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

/* Kinds of special point detected along an equilibrium branch.
 * Codim-1: Fold (limit point), Hopf. Codim-2 (degeneracies of the codim-1
 * conditions): BogdanovTakens (fold & Hopf collide, double-zero eigenvalue),
 * Cusp (fold normal-form coefficient a -> 0), GeneralizedHopf / Bautin (Hopf
 * first Lyapunov coefficient l1 -> 0). BranchPoint = two equilibrium branches
 * cross (the augmented system is singular though the fold test need not be). */
enum class SpecialPointKind { None, Fold, Hopf, BogdanovTakens, Cusp,
                              GeneralizedHopf, ZeroHopf, HopfHopf, BranchPoint, EndOfBranch };

struct BranchPoint {
  double p = 0.0;               /* parameter value                  */
  std::vector<double> x;        /* equilibrium coordinates (len n)  */
  std::vector<Complex> eigenvalues;
  int n_unstable = 0;           /* unstable eigenvalue count        */
  bool stable = false;          /* n_unstable == 0                  */
  SpecialPointKind special = SpecialPointKind::None;
  /* Normal-form data attached at special points (NaN when not computed):
   * at a Hopf, lyapunov1 is the first Lyapunov coefficient; at a Fold,
   * fold_a is the quadratic normal-form coefficient. These let the UI label
   * criticality and flag codim-2 degeneracies (l1~0 -> generalized Hopf,
   * a~0 -> cusp). second_tangent, when non-empty at a BranchPoint, is the
   * other branch's direction in (x,p) space for branch switching. */
  double lyapunov1 = 0.0;
  double fold_a = 0.0;
  bool has_normal_form = false;
  std::vector<double> second_tangent; /* length n+1, only at BranchPoint */
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

/* ---- Hopf first Lyapunov coefficient (normal-form classification) -------- *
 * At a Hopf point the Jacobian has a pure-imaginary pair +/- i*omega. The
 * sign of the first Lyapunov coefficient l1 decides the criticality:
 *   l1 < 0 : SUPERcritical Hopf  (a stable limit cycle is born)
 *   l1 > 0 : SUBcritical  Hopf  (an unstable limit cycle; hard loss of stab.)
 *   l1 ~ 0 : degenerate (Bautin / generalized Hopf — a codim-2 point)
 * This is the quantity MatCont reports and pplane/XPPAUT do not. Computed from
 * the 2nd and 3rd derivatives of the vector field at the equilibrium via the
 * standard projection onto the critical eigenvector (Kuznetsov's formula),
 * using finite differences on m.vector_field (so it works for any system).
 * Returns false if there is no clean imaginary pair at (x,p) or on eval error;
 * on success *l1 holds the coefficient and *omega the frequency. */
bool hopf_first_lyapunov(const Model &m, const std::vector<double> &x, double p,
                         double *l1, double *omega, std::string *err);

/* Fold (limit point) normal-form coefficient a (Kuznetsov):
 *   a = (1/2) <p, B(q,q)> / <p, q>,   A q = 0,  p^T A = 0.
 * a != 0 is the fold non-degeneracy condition; a ~ 0 signals a CUSP (codim-2).
 * This is the limit-point counterpart of the Hopf l1 above and the quantity
 * MatCont reports at an LP. On success *a holds the coefficient and *lambda0
 * the magnitude of the near-zero eigenvalue (small => genuinely at a fold). */
bool fold_normal_form(const Model &m, const std::vector<double> &x, double p,
                      double *a, double *lambda0, std::string *err);

/* Cusp normal-form coefficient c (codim-2 on a fold curve where a -> 0). The
 * center manifold is cubic y' = c y^3 + ...; c != 0 is the non-degeneracy
 * condition. */
bool cusp_normal_form(const Model &m, const std::vector<double> &x, double p,
                      double *c, std::string *err);

/* Generalized-Hopf (Bautin) SECOND Lyapunov coefficient l2 (codim-2 on a Hopf
 * curve where l1 -> 0). Uses multilinear forms up to 5th order by finite
 * differences, so it is approximate -- the SIGN is the reliable, classifying
 * output. Returns false away from a Hopf point. */
bool gh_second_lyapunov(const Model &m, const std::vector<double> &x, double p,
                        double *l2, std::string *err);

/* Zero-Hopf (Gavrilov-Guckenheimer / fold-Hopf) normal-form coefficients.
 *
 * At a zero-Hopf point the Jacobian has eigenvalues {0, +i*omega, -i*omega}.
 * The quadratic normal form on the 3-D center manifold (Kuznetsov, Elements of
 * Applied Bifurcation Theory, sec. 8.5) is governed by two leading coefficients:
 *     b = < p0, B(q0, q0) >                 (real)
 *     c = < p1, B(q0, q1) >                 (complex; Re c is the key part)
 * where q0 is the real null eigenvector (A q0 = 0), q1 the complex eigenvector
 * (A q1 = i*omega q1), and p0, p1 the corresponding adjoint vectors normalized
 * so <p0,q0> = 1 and <p1,q1> = 1. The SIGN of b and of Re(c) (and the product
 * s = b * Re(c)) classify which local bifurcation diagram occurs -- in
 * particular whether an invariant torus and/or secondary connections are born.
 * Returns b, Re(c), omega, and s = b*Re(c). */
bool zero_hopf_normal_form(const Model &m, const std::vector<double> &x, double p,
                           double *b_out, double *recc_out, double *omega_out,
                           double *s_out, std::string *err);

/* Hopf-Hopf (double-Hopf) normal-form coefficients.
 *
 * At a Hopf-Hopf point the Jacobian has TWO pure-imaginary pairs
 * {+-i*omega1, +-i*omega2} (omega1 != omega2, non-resonant). The truncated
 * amplitude normal form on the 4-D center manifold (Kuznetsov, Elements of
 * Applied Bifurcation Theory, sec. 8.6) is
 *     r1' = r1 (mu1 + p11 r1^2 + p12 r2^2)
 *     r2' = r2 (mu2 + p21 r1^2 + p22 r2^2)
 * (real cubic truncation). The four real cubic coefficients are:
 *   p11 = self-interaction of mode 1 (its first Lyapunov coefficient),
 *   p22 = self-interaction of mode 2,
 *   p12 = mode-1 amplitude driven by mode-2 amplitude (cross-coupling),
 *   p21 = mode-2 amplitude driven by mode-1 amplitude (cross-coupling).
 * The signs of p11,p22 and of theta = p12/p22, delta = p21/p11 classify the
 * local diagram (whether the two limit cycles can coexist, and whether an
 * invariant 2-torus / heteroclinic structures appear -- the "simple" vs
 * "difficult" cases). Returns p11,p12,p21,p22 and the two frequencies. */
bool hopf_hopf_normal_form(const Model &m, const std::vector<double> &x, double p,
                           double *p11_out, double *p12_out, double *p21_out,
                           double *p22_out, double *omega1_out, double *omega2_out,
                           std::string *err);

/* Bogdanov-Takens normal-form coefficients (a, b) at a BT point (double-zero
 * eigenvalue, 2x2 Jordan block). The planar BT normal form is w0'=w1,
 * w1'=a w0^2 + b w0 w1; non-degeneracy needs a != 0 and b != 0. */
bool bt_normal_form(const Model &m, const std::vector<double> &x, double p,
                    double *a, double *b, std::string *err);

/* Branch switching at a branch point. Given a BranchPoint that carries a
 * second_tangent (the crossing branch's direction in (x,p) space), take a
 * small step along +/- that direction, re-converge to an equilibrium, and
 * continue from there — tracing the OTHER branch. Returns an empty/!ok Branch
 * if bp has no second tangent or the off-branch point can't be corrected. */
Branch switch_branch(const Model &m, const BranchPoint &bp,
                     const ContinuationSettings &settings);

/* ---- two-parameter continuation of codim-1 curves ----------------------- *
 * A system whose vector field depends on TWO parameters (p, q). Used to trace
 * a curve of fold points or Hopf points in the (p, q) plane: the locus where
 * the codim-1 condition holds as both parameters vary. */
struct Model2 {
  std::size_t n = 0;
  /* f(x, p, q) -> f_out (length n). */
  std::function<bool(const double *x, double p, double q, double *f_out, std::string *err)>
      vector_field;
};

/* One point on a two-parameter curve: the two parameter values and the
 * equilibrium there. kind is Fold or Hopf (which curve this is). */
struct TwoParamPoint {
  double p = 0.0, q = 0.0;
  std::vector<double> x;
  /* codim-2 point detected ON the curve (None for an ordinary point). On a
   * fold curve: Cusp (fold coeff a -> 0), BogdanovTakens (a 2nd eigenvalue
   * reaches zero). On a Hopf curve: GeneralizedHopf (l1 -> 0),
   * BogdanovTakens (frequency omega -> 0). */
  SpecialPointKind special = SpecialPointKind::None;
  /* When special != None, the location is REFINED by bisection along the curve:
   * (p2,q2,x2) is the precise codim-2 point, and the normal-form coefficients
   * below classify its unfolding (NaN when not computed). For BT: bt_a, bt_b
   * (the quadratic normal-form coefficients of the planar BT normal form
   *  w2' = a*w1^2 + b*w1*w2); for Cusp: the fold coeff a is ~0 there and we
   * report the cubic-ish residual; for GeneralizedHopf: l1 ~ 0 there. */
  double p2 = 0.0, q2 = 0.0;
  std::vector<double> x2;
  double bt_a = 0.0, bt_b = 0.0;
  double cusp_c = 0.0;   /* cusp cubic coefficient (Cusp points)            */
  double gh_l2 = 0.0;    /* second Lyapunov coefficient (GeneralizedHopf)   */
  /* zero-Hopf (ZeroHopf) normal-form coefficients: b, Re(c), and their product
   * s = b*Re(c) whose sign classifies the local diagram. */
  double zh_b = 0.0, zh_recc = 0.0, zh_s = 0.0;
  /* Hopf-Hopf (HopfHopf) cubic coefficients of the amplitude normal form
   * r1'=r1(.. + p11 r1^2 + p12 r2^2), r2'=r2(.. + p21 r1^2 + p22 r2^2). */
  double hh_p11 = 0.0, hh_p12 = 0.0, hh_p21 = 0.0, hh_p22 = 0.0;
  bool has_codim2_nf = false;
};
struct TwoParamCurve {
  std::vector<TwoParamPoint> points;
  std::vector<std::size_t> special_indices; /* indices of codim-2 points */
  std::string message;
  bool ok = false;
  SpecialPointKind kind = SpecialPointKind::Fold;
};

enum class TwoParamKind { Fold, Hopf };

struct TwoParamSettings {
  double p_min = -10, p_max = 10;   /* bounds on parameter 1 (x axis)  */
  double q_min = -10, q_max = 10;   /* bounds on parameter 2 (y axis)  */
  double h0 = 0.05;                 /* arclength step in (p,q,x) space */
  int max_points = 800;
  int max_corrector_iters = 20;
  double corrector_tol = 1e-9;
};

/* Trace a fold or Hopf curve in the (p,q) plane, starting from a codim-1 point
 * found at (p0,q0,x0). The extended system is G(x,p,q) = [ f ; g ] = 0 where g
 * is the fold test det(f_x) (Fold) or the real-part-sum proxy (Hopf); this is
 * n+1 equations in n+2 unknowns, whose solution set is the 1-D curve. Returns
 * points sampled along it, both directions from the start. */
TwoParamCurve two_param_curve(const Model2 &m, TwoParamKind kind,
                              const std::vector<double> &x0, double p0, double q0,
                              const TwoParamSettings &settings);

/* ---- CODIM-2 POINT LOCATION (defining-system Newton) -------------------- *
 * A codim-2 point sits at the intersection of two codim-1 conditions and is
 * found, in two free parameters, by solving its DEFINING SYSTEM directly (the
 * way MatCont's BT/ZH "init" does) rather than only being detected along a
 * curve. This pins the point -- equilibrium x AND both parameters (p,q) -- to
 * Newton precision, which is the prerequisite for then continuing it in a third
 * parameter. */
struct Codim2Point {
  std::vector<double> x;   /* equilibrium at the located point */
  double p = 0.0, q = 0.0; /* the two parameter values         */
  double residual = 0.0;   /* defining-system residual norm    */
  int iters = 0;
  bool ok = false;
  SpecialPointKind kind = SpecialPointKind::None;
  /* normal-form coefficients at the located point (filled for BT: a,b) */
  double a = 0.0, b = 0.0;
  std::string message;
};

/* Locate a Bogdanov-Takens point (double-zero eigenvalue) by Newton on the
 * defining system { f(x,p,q) = 0 ; mu1(x,p,q) = 0 ; mu2(x,p,q) = 0 }, where
 * mu1, mu2 are the two bordered "double-zero" test functions (the product and
 * sum of the two smallest eigenvalues vanish at BT). x0,p0,q0 is a guess near
 * the point (e.g. from detection along a fold/Hopf curve). On success the BT
 * normal-form coefficients a,b are also computed at the located point. */
Codim2Point locate_bogdanov_takens(const Model2 &m, const std::vector<double> &x0,
                                   double p0, double q0);

/* ---- CODIM-2 CURVE CONTINUATION (a codim-2 point in 3 parameters) -------- *
 * A Bogdanov-Takens point is codim-2 (two defining conditions). With a THIRD
 * free parameter its solution set is a 1-D CURVE in (p,q,r) space. We continue
 * it by pseudo-arclength on the defining system
 *   { f(x,p,q,r) = 0  (n eqns) ;  mu_prod = 0 ;  mu_sum = 0 }
 * which is n+2 equations in the n+3 unknowns (x, p, q, r). This is MatCont's
 * "BT curve" (a codim-2 locus continued in three parameters). */
struct Model3 {
  std::size_t n = 0;
  /* f(x, p, q, r) -> f_out (length n). */
  std::function<bool(const double *x, double p, double q, double r, double *f_out, std::string *err)>
      vector_field;
};

struct BTCurvePoint {
  std::vector<double> x;          /* equilibrium */
  double p = 0.0, q = 0.0, r = 0.0; /* the three parameters */
  double a = 0.0, b = 0.0;        /* BT normal-form coefficients here */
  double residual = 0.0;
};

struct BTCurve {
  std::vector<BTCurvePoint> points;
  std::string message;
  bool ok = false;
};

struct BTCurveSettings {
  double ds = 0.02;          /* pseudo-arclength step in (x,p,q,r) space */
  int max_points = 300;      /* per direction */
  int max_corrector_iters = 25;
  double corrector_tol = 1e-9;
  /* optional box to stop the curve when it leaves the region of interest */
  double p_min = -1e9, p_max = 1e9, q_min = -1e9, q_max = 1e9, r_min = -1e9, r_max = 1e9;
};

/* Continue a Bogdanov-Takens curve in three parameters, starting from a point
 * (x0,p0,q0,r0) on it (e.g. refined first by locate_bogdanov_takens on the
 * (p,q) slice at the chosen r0). Traces both directions from the start. */
BTCurve bt_curve(const Model3 &m, const std::vector<double> &x0,
                 double p0, double q0, double r0, const BTCurveSettings &settings);

/* Continue a ZERO-HOPF (fold-Hopf) point as a curve in three parameters.
 * The zero-Hopf defining conditions are: one real eigenvalue = 0 AND one
 * complex-conjugate pair on the imaginary axis (Re = 0). Over the unknowns
 * (x,p,q,r) this is n+2 equations -> a 1-D locus, traced by the same
 * pseudo-arclength scheme as bt_curve. Reuses BTCurve/BTCurveSettings (the
 * BTCurvePoint a,b fields carry the zero-Hopf b and Re(c) here). */
BTCurve zh_curve(const Model3 &m, const std::vector<double> &x0,
                 double p0, double q0, double r0, const BTCurveSettings &settings);

/* Continue a CUSP point as a curve in three parameters. A cusp sits on the fold
 * curve where the fold (saddle-node) normal-form coefficient a vanishes. The
 * defining conditions over (x,p,q,r) are: one real eigenvalue = 0 (fold) AND
 * the cusp coefficient a = 0. Traced by the same pseudo-arclength scheme as
 * bt_curve. The BTCurvePoint a,b fields carry the cusp coefficient c (a) and the
 * fold-eigenvalue residual (b). */
BTCurve cusp_curve(const Model3 &m, const std::vector<double> &x0,
                   double p0, double q0, double r0, const BTCurveSettings &settings);

/* Continue a HOPF-HOPF (double-Hopf) point as a curve in three parameters. The
 * defining conditions are: TWO distinct complex-conjugate pairs simultaneously
 * on the imaginary axis (both Re = 0). Over (x,p,q,r) this is n+2 equations ->
 * a 1-D locus, traced by the same scheme. The BTCurvePoint a,b fields carry the
 * two frequencies omega1, omega2 at each point. */
BTCurve hh_curve(const Model3 &m, const std::vector<double> &x0,
                 double p0, double q0, double r0, const BTCurveSettings &settings);

/* ---- HOMOCLINIC ORBITS (connection to a saddle) ------------------------- *
 * A homoclinic orbit leaves a saddle equilibrium along its unstable manifold
 * and returns to the SAME saddle along its stable manifold, so x(t) -> x0 as
 * t -> +/- infinity. We approximate it on a long but finite interval, rescaled
 * to tau in [0,1], with:
 *   - collocation of  x'(tau) = 2T f(x)   (T = the truncation half-time; the
 *     factor 2T maps real time [-T,+T] to tau in [0,1]);
 *   - PROJECTION boundary conditions (Beyn): the deviation x(0)-x0 lies in the
 *     unstable subspace of A=f_x(x0) (its stable-subspace projection is zero),
 *     and x(1)-x0 lies in the stable subspace (its unstable projection is zero);
 *   - a phase condition pinning the translation freedom.
 * This is the rigorous truncated-BVP formulation used by HomCont/AUTO. */
struct HomoclinicSettings {
  int mesh = 150;            /* collocation intervals over [0,1]            */
  double T = 0.0;            /* truncation half-time; 0 => auto from the seed */
  int newton_iters = 40;
  double newton_tol = 1e-9;
  bool free_T = true;        /* let Newton adjust T as an unknown            */
};
struct HomoclinicResult {
  bool ok = false;
  std::string message;
  std::vector<double> saddle;                 /* x0 (length n)               */
  std::vector<std::vector<double>> orbit;     /* mesh+1 points, each length n */
  std::vector<double> tau;                    /* mesh+1 values in [0,1]      */
  double T = 0.0;                             /* converged half-time         */
  double amplitude = 0.0;                     /* max deviation from saddle    */
  int n_unstable = 0;                         /* dim of unstable subspace     */
  double newton_residual = 0.0;
  int newton_steps = 0;
};

/* Solve for a homoclinic orbit to the saddle x0_guess at parameter p. The seed
 * orbit (one long excursion that leaves and returns near the saddle, e.g. from
 * a simulation) is given as seed_orbit (>=4 points, ordered along the orbit);
 * it is resampled onto the collocation mesh. Returns the refined orbit + T. */
HomoclinicResult solve_homoclinic(const Model &m,
                                  const std::vector<double> &x0_guess, double p,
                                  const std::vector<std::vector<double>> &seed_orbit,
                                  const HomoclinicSettings &settings);

/* Heteroclinic connection: an orbit from saddle x0 (leaving along its UNSTABLE
 * manifold) to a DIFFERENT saddle x1 (arriving along x1's STABLE manifold),
 * solved as a truncated BVP with projection boundary conditions split between
 * the two equilibria. Generalizes solve_homoclinic (which is the x1 == x0
 * case). seed_orbit is an ordered guess of the connecting trajectory (e.g. from
 * integrating x0's unstable manifold until it nears x1). */
struct HeteroclinicResult {
  bool ok = false;
  std::string message;
  std::vector<double> saddle0, saddle1;       /* the two equilibria          */
  std::vector<std::vector<double>> orbit;     /* mesh+1 points, each length n */
  std::vector<double> tau;                    /* mesh+1 values in [0,1]      */
  double T = 0.0;                             /* converged half-time         */
  double length = 0.0;                        /* arclength of the connection  */
  int n_unstable0 = 0, n_stable1 = 0;
  double newton_residual = 0.0;
  int newton_steps = 0;
};
HeteroclinicResult solve_heteroclinic(const Model &m,
                                      const std::vector<double> &x0_guess,
                                      const std::vector<double> &x1_guess, double p,
                                      const std::vector<std::vector<double>> &seed_orbit,
                                      const HomoclinicSettings &settings);

/* Build a homoclinic SEED by integrating the unstable manifold: nudge a point
 * off the saddle along its dominant unstable eigenvector and integrate forward
 * with RK4 (trying both orientations) until the trajectory returns near the
 * saddle. For a system with (or near) a homoclinic this shadows the orbit and
 * gives a far better initial guess than a hand-built bump. Returns false if no
 * return is detected within the budget. */
bool seed_homoclinic_by_integration(const Model &m, const std::vector<double> &saddle,
                                     double p, double dt, double max_time,
                                     std::vector<std::vector<double>> *seed_out,
                                     std::string *err);

/* One-call homoclinic FINDER (no need to know the locus parameter up front).
 *
 * Sweeps the parameter p over [p_lo, p_hi]; at each p it integrates the saddle's
 * unstable manifold and measures the CLOSEST RETURN distance back to the saddle
 * after the excursion. That distance is ~0 exactly on the homoclinic and grows
 * to either side (the manifold either spirals in short or escapes), so its
 * minimum over the sweep brackets the connection. p is then refined by golden-
 * section on the return distance, the orbit at that p is seeded by integration,
 * and solve_homoclinic is run to produce the converged BVP orbit. This wraps the
 * validated seeding + BVP solve and works for same-saddle homoclinics where the
 * unstable manifold genuinely returns (it reports honestly when, as in the
 * Bogdanov-Takens regime, the manifold instead connects to another equilibrium).
 *
 * Returns the located p in p_out. */
HomoclinicResult find_homoclinic(const Model2 &m2,
                                 const std::vector<double> &saddle_guess,
                                 double q_fixed, double p_lo, double p_hi,
                                 const HomoclinicSettings &settings, double *p_out);

/* ---- two-parameter continuation of the HOMOCLINIC locus ----------------- *
 * A homoclinic connection is codimension-1, so in a two-parameter (p,q) plane
 * it traces a CURVE. We continue it by stepping the secondary parameter q and,
 * at each q, finding the primary parameter p at which the truncated-BVP
 * homoclinic closes (its boundary/collocation residual is minimised to ~0),
 * re-using the previous converged orbit as the seed (natural continuation).
 * This wraps solve_homoclinic as the inner solve. */
struct HomoclinicCurvePoint {
  double p = 0.0, q = 0.0;
  double T = 0.0, amplitude = 0.0;
  double residual = 0.0;
};
struct HomoclinicCurve {
  std::vector<HomoclinicCurvePoint> points;
  std::vector<std::vector<std::vector<double>>> orbits; /* the orbit at each point (optional, same order) */
  std::string message;
  bool ok = false;
};
struct HomoclinicContSettings {
  int max_steps = 60;        /* number of q-steps in each direction          */
  double dq = 0.02;          /* step in the secondary parameter q            */
  double q_min = -1e9, q_max = 1e9; /* stop if q leaves this window           */
  HomoclinicSettings bvp;    /* inner BVP settings (mesh, T, newton)          */
  bool both_directions = true;
  bool store_orbits = false;
};
/* Trace the homoclinic curve from a starting (p0,q0) where a homoclinic is
 * known (seed_orbit shadows it). m2 is the two-parameter field. */
HomoclinicCurve continue_homoclinic(const Model2 &m2,
                                    const std::vector<double> &x0_guess,
                                    double p0, double q0,
                                    const std::vector<std::vector<double>> &seed_orbit,
                                    const HomoclinicContSettings &settings);

/* ---- Lin's method: locate a homoclinic by closing the Lin gap ----------- *
 * EXPERIMENTAL / PRELIMINARY. The scaffolding is in place -- it relocates the
 * saddle at each parameter, computes the dominant stable/unstable directions,
 * integrates each manifold to a Poincare section, measures an in-section gap,
 * and runs a secant root-find on the primary parameter. However, robustly
 * obtaining a SIGNED gap that brackets the homoclinic (so the root-find
 * converges) requires careful section placement where both manifolds genuinely
 * cross; for stiff Bogdanov-Takens loops this is not yet reliable. Use
 * solve_homoclinic / continue_homoclinic for the validated path. This is kept
 * as a foundation for a future complete Lin's-method implementation.
 *
 * For systems where a good orbit seed cannot be obtained by integration alone
 * (e.g. a Bogdanov-Takens homoclinic, a large slow loop that exists only on an
 * exact parameter curve), Lin's method splits the connection at a Poincare
 * section: integrate the UNSTABLE manifold forward to the section (x+) and the
 * STABLE manifold backward to it (x-), measure the signed in-section gap, and
 * adjust the PRIMARY parameter p (at fixed secondary q) to drive that gap to
 * zero. m2 is f(x,p,q); the saddle is relocated at each p. */
struct LinResult {
  bool ok = false;
  std::string message;
  double p = 0.0;            /* parameter at which the gap closed            */
  double gap = 0.0;          /* residual Lin gap there (small if converged)  */
  std::vector<double> saddle;
  std::vector<std::vector<double>> orbit; /* concatenated u-branch + v-branch */
  double amplitude = 0.0;
  int iterations = 0;
};
struct LinSettings {
  double eps = 1e-4;         /* initial displacement along the eigenvectors  */
  double dt = 0.01;          /* integration step                            */
  double max_time = 200.0;   /* per-branch integration budget               */
  double p_lo = -1e9, p_hi = 1e9; /* bracket for the parameter root-find     */
  int max_iter = 60;
  double tol = 1e-7;         /* gap tolerance                               */
};
LinResult lin_homoclinic(const Model2 &m2, const std::vector<double> &saddle_guess,
                         double p0, double q, const LinSettings &settings);

/* ---- two-parameter continuation of a FOLD-OF-CYCLES (LPC) curve --------- *
 * Declared after the cycle structs below (it needs CycleSettings). */
struct LPCPoint { double p = 0.0, q = 0.0; double period = 0.0; double amplitude = 0.0; };
struct LPCCurve {
  std::vector<LPCPoint> points;
  std::string message;
  bool ok = false;
};

/* Two-parameter loci of cycle codim-1 bifurcations: period-doubling (PD, a real
 * Floquet multiplier = -1) and Neimark-Sacker (NS, a complex multiplier pair on
 * the unit circle). Same point type as LPC. Traced by scanning the secondary
 * parameter q, continuing the cycle in p at each q (with Floquet on), and
 * picking the bracketed PD / NS sample on that branch. */
struct CycleBifPoint {
  double p = 0.0, q = 0.0; double period = 0.0; double amplitude = 0.0;
  /* secondary cycle test-function values AT this curve point, used to flag
   * codim-2 degeneracies along the curve (e.g. an LPC coinciding with the PD). */
  double fold_test = 0.0;     /* cycle fold (LPC) test at this point          */
  double ns_test = 0.0;       /* Neimark-Sacker test at this point            */
  double max_nontrivial_mult = 0.0;
};

/* A codim-2 point on a cycle-bifurcation curve: where a SECOND codim-1 cycle
 * condition is met along the (PD or NS) locus. */
enum class CycleCodim2Kind { None, FoldFlip, PDNS, DegeneratePD, CuspOfCycles };
struct CycleCodim2Point {
  double p = 0.0, q = 0.0;
  CycleCodim2Kind kind = CycleCodim2Kind::None;
  const char *label() const {
    switch (kind) {
      case CycleCodim2Kind::FoldFlip:     return "LPPD (fold-flip: LPC = PD)";
      case CycleCodim2Kind::PDNS:         return "PD-NS interaction";
      case CycleCodim2Kind::DegeneratePD: return "degenerate PD (curve meets Hopf, amp->0)";
      case CycleCodim2Kind::CuspOfCycles: return "cusp of cycles (LPC turning point)";
      default: return "none";
    }
  }
};
struct CycleBifCurve {
  std::vector<CycleBifPoint> points;
  std::vector<CycleCodim2Point> codim2;   /* codim-2 points detected on the curve */
  std::string message;
  bool ok = false;
};

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

/* Parallel variant: `make_advance(tid)` returns an advance function private to
 * worker thread tid (so each thread can own its own evaluation scratch and
 * avoid data races). The expensive per-cell integration runs across hardware
 * threads; the attractor clustering is done serially afterwards so labels stay
 * deterministic. Produces the same result as compute_basins, just faster on
 * multi-core machines. Falls back to serial when only one core is available. */
using AdvanceFn = std::function<bool(double x, double y, double *nx, double *ny)>;
BasinResult compute_basins_mt(const std::function<AdvanceFn(int tid)> &make_advance,
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

/* ---- Iterated Function System (chaos game) -------------------- *
 * An IFS is a set of affine contraction maps with probabilities; its
 * attractor (Barnsley fern, Sierpinski, dragon, ...) is drawn by the
 * "chaos game": start anywhere, repeatedly pick a map at random (weighted
 * by probability) and apply it, plotting each point. AppState-free and
 * testable (the attractor has a known box-counting dimension). */
struct AffineMap {
  double a = 1, b = 0, c = 0, d = 1, e = 0, f = 0; /* x' = a x + b y + e, y' = c x + d y + f */
  double p = 0.0;                                  /* selection probability */
};

struct IFSResult {
  bool ok = false;
  std::vector<float> xs, ys;   /* the plotted points (after burn-in) */
  double xmin = 0, xmax = 0, ymin = 0, ymax = 0;
  std::string message;
};

/* Run the chaos game for `iterations` steps (minus a small burn-in),
 * returning the visited points. `seed` makes it deterministic. */
IFSResult chaos_game(const std::vector<AffineMap> &maps, long iterations,
                     unsigned int seed = 12345u);

/* ---- Limit-cycle period & amplitude (foundation for LC continuation) -- *
 * Given a settled, sampled scalar signal y(t) from an oscillating system,
 * estimate the oscillation PERIOD and AMPLITUDE. The period is found from
 * the dominant peak of the autocorrelation (robust to noise and to the
 * exact section), the amplitude from the settled signal's peak-to-peak
 * range. Returns ok=false for non-oscillatory signals (fixed points). This
 * is the per-parameter measurement that a limit-cycle continuation diagram
 * (period vs parameter, amplitude vs parameter) is built from. */
struct LimitCycleResult {
  bool ok = false;
  double period = 0.0;     /* in units of the sample spacing dt */
  double amplitude = 0.0;  /* peak-to-peak of the settled signal */
  std::string message;
};

LimitCycleResult limit_cycle_period_amplitude(const std::vector<double> &y, double dt);

/* ---- periodic-orbit continuation by collocation ------------------------- *
 * Represents a periodic orbit on a uniform mesh of m points in [0,1) with the
 * BVP  x'(s) = T f(x(s)),  x(0) = x(1),  plus an integral phase condition that
 * pins the orbit's position along itself. Solving this by Newton finds the
 * cycle AND its period T directly, so it follows unstable cycles and folds of
 * cycles (unlike a simulate-and-measure approach). Continuation steps the
 * parameter and re-solves, producing the cycle's period and amplitude vs the
 * parameter. */
struct CycleSample {
  double p = 0.0;          /* parameter value                              */
  double period = 0.0;     /* orbit period T                               */
  double amplitude = 0.0;  /* peak-to-peak of the first state over a cycle */
  double min0 = 0.0, max0 = 0.0; /* min/max of first state along the orbit */
  bool stable = false;     /* nontrivial Floquet multipliers inside unit circle */
  double fold_test = 0.0;  /* fold-of-cycles (LPC) test function value      */
  bool is_fold = false;    /* an LPC (cycle fold) was bracketed at this sample */
  /* Floquet multipliers of the cycle (the monodromy matrix eigenvalues). One
   * multiplier is always ~1 (the trivial one, along the flow); the rest govern
   * stability. Populated when compute_floquet is on. */
  std::vector<double> floquet_re, floquet_im;
  double max_nontrivial_mult = 0.0; /* |largest non-trivial multiplier|     */
  bool is_pd = false;      /* period-doubling (a multiplier passed -1)      */
  bool is_ns = false;      /* Neimark-Sacker / torus (complex pair hit |.|=1)*/
  /* SIGNED test functions whose sign change between samples brackets the
   * bifurcation (refined by bisection, like the fold test):
   *  pd_test = prod over non-trivial multipliers of (mu_i + 1)  -> 0 at PD
   *  ns_test = (max non-trivial |mu| of a COMPLEX pair) - 1     -> 0 at NS  */
  double pd_test = 0.0, ns_test = 0.0;
  double pd_p = 0.0, ns_p = 0.0;  /* refined parameter at the PD / NS point  */
  /* Branch point of cycles (BPC): where two periodic-orbit branches cross. The
   * signed bp_test (a bordered determinant) changes sign at a BPC; is_bp/bp_p
   * are set when a sign change brackets one. */
  double bp_test = 0.0;
  bool is_bp = false;
  double bp_p = 0.0;
};
struct CycleBranch {
  std::vector<CycleSample> samples;
  std::string message;
  bool ok = false;
  bool turned = false;     /* true if arclength continuation went around a fold */
};
struct CycleSettings {
  int mesh = 60;                 /* collocation points around the orbit     */
  double p_min = -1e9, p_max = 1e9;
  double dp = 0.02;              /* parameter step (monotone mode)          */
  int max_steps = 400;
  int newton_iters = 30;
  double newton_tol = 1e-8;
  bool arclength = true;         /* pseudo-arclength (follows folds of cycles) */
  double ds = 0.05;              /* arclength step size                      */
  bool compute_floquet = true;   /* compute Floquet multipliers per sample   */
  bool adaptive_mesh = true;     /* redistribute mesh by arclength (stiff cycles) */
};

/* Trace a branch of periodic orbits starting from an initial guess (a set of
 * `mesh` points around one loop and an initial period guess), continuing in
 * the parameter from p0. The initial guess is typically taken from a simulated
 * cycle; this solver then refines it exactly and follows it. */
CycleBranch continue_limit_cycle(const Model &m,
                                 const std::vector<std::vector<double>> &guess_points,
                                 double period_guess, double p0,
                                 const CycleSettings &settings);

/* Two-parameter fold-of-cycles (LPC) curve: traces, in the (p,q) plane, the
 * locus where a limit cycle folds (a saddle-node of cycles; a nontrivial
 * Floquet multiplier passes +1). For each q it continues the cycle in p and
 * records where the fold-of-cycles test changes sign. Seeded from a cycle
 * guess + period at (p0,q0). */
LPCCurve lpc_curve(const Model2 &m,
                   const std::vector<std::vector<double>> &guess_points,
                   double period_guess, double p0, double q0,
                   const TwoParamSettings &settings, const CycleSettings &cyc);

/* Period-doubling (PD) curve in two parameters. */
CycleBifCurve pd_curve(const Model2 &m,
                       const std::vector<std::vector<double>> &guess_points,
                       double period_guess, double p0, double q0,
                       const TwoParamSettings &settings, const CycleSettings &cyc);

/* Neimark-Sacker (NS / torus) curve in two parameters. */
CycleBifCurve ns_curve(const Model2 &m,
                       const std::vector<std::vector<double>> &guess_points,
                       double period_guess, double p0, double q0,
                       const TwoParamSettings &settings, const CycleSettings &cyc);

}  // namespace dynsys::analysis
