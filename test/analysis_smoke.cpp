/* Standalone smoke test for the dynsys analysis core.
 *
 * Builds with just a C++17 compiler and -lm (no GL / imgui / tpcas),
 * mirroring test/ir_smoke.cpp. Verifies the eigensolver against known
 * spectra, the classifier, and the continuation engine on systems with
 * an analytically known fold and Hopf bifurcation.
 *
 *   make test-analysis
 */

#include "../src/analysis.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace dynsys::analysis;

static int g_failures = 0;
static int g_checks = 0;

static void check(bool cond, const char *what) {
  ++g_checks;
  if (!cond) {
    ++g_failures;
    std::printf("  FAIL: %s\n", what);
  }
}

static bool close(double a, double b, double tol = 1e-6) {
  return std::fabs(a - b) <= tol * (1.0 + std::fabs(a) + std::fabs(b));
}

/* Find an eigenvalue in the list close to a target complex value. */
static bool has_eigenvalue(const std::vector<Complex> &ev, Complex target,
                           double tol = 1e-5) {
  for (const Complex &z : ev)
    if (std::abs(z - target) <= tol * (1.0 + std::abs(target))) return true;
  return false;
}

static void test_eigen_diagonal() {
  std::printf("eigen: diagonal\n");
  std::vector<double> A = {2.0, 0.0, 0.0, 0.0, -3.0, 0.0, 0.0, 0.0, 5.0};
  std::vector<Complex> ev;
  check(eigenvalues(A, 3, &ev), "diag solve ok");
  check(has_eigenvalue(ev, {2, 0}), "has 2");
  check(has_eigenvalue(ev, {-3, 0}), "has -3");
  check(has_eigenvalue(ev, {5, 0}), "has 5");
}

static void test_eigen_rotation() {
  /* [[0,-1],[1,0]] has eigenvalues +/- i. */
  std::printf("eigen: rotation (complex pair)\n");
  std::vector<double> A = {0.0, -1.0, 1.0, 0.0};
  std::vector<Complex> ev;
  check(eigenvalues(A, 2, &ev), "rot solve ok");
  check(has_eigenvalue(ev, {0, 1}), "has +i");
  check(has_eigenvalue(ev, {0, -1}), "has -i");
}

static void test_eigen_known_3x3() {
  /* Companion-style matrix with known characteristic polynomial.
   * Use a triangular matrix whose eigenvalues are the diagonal, but
   * with off-diagonal coupling so it actually exercises QR. */
  std::printf("eigen: 3x3 with complex pair\n");
  /* Block diag: 2x2 rotation-scaled (real part -0.5, imag +/-2) and a
   * real eigenvalue 4.
   *   [[-0.5, -2, 0], [2, -0.5, 0], [0, 0, 4]] */
  std::vector<double> A = {-0.5, -2.0, 0.0, 2.0, -0.5, 0.0, 0.0, 0.0, 4.0};
  std::vector<Complex> ev;
  check(eigenvalues(A, 3, &ev), "3x3 solve ok");
  check(has_eigenvalue(ev, {-0.5, 2.0}), "has -0.5+2i");
  check(has_eigenvalue(ev, {-0.5, -2.0}), "has -0.5-2i");
  check(has_eigenvalue(ev, {4.0, 0.0}), "has 4");
}

static void test_determinant() {
  std::printf("determinant\n");
  std::vector<double> A = {1, 2, 3, 4};  /* det = -2 */
  check(close(determinant(A, 2), -2.0), "det 2x2");
  std::vector<double> B = {6, 1, 1, 4, -2, 5, 2, 8, 7};  /* det = -306 */
  check(close(determinant(B, 3), -306.0, 1e-9), "det 3x3");
}

static void test_classify_saddle() {
  std::printf("classify: 2D saddle\n");
  /* Eigenvalues +1, -1 => saddle. */
  std::vector<double> J = {1.0, 0.0, 0.0, -1.0};
  Classification cl = classify_equilibrium(J, 2);
  check(cl.is_saddle, "is saddle");
  check(cl.n_stable == 1 && cl.n_unstable == 1, "1 stable 1 unstable");
}

static void test_classify_stable_spiral() {
  std::printf("classify: stable spiral\n");
  std::vector<double> J = {-0.5, -2.0, 2.0, -0.5};
  Classification cl = classify_equilibrium(J, 2);
  check(cl.n_unstable == 0, "no unstable dirs");
  check(cl.label.find("stable") != std::string::npos, "labelled stable");
}

/* ---- continuation: fold (saddle-node) -----------------------
 * Scalar system  x' = p + x^2.
 * Equilibria x = +/- sqrt(-p) exist for p < 0 and collide at the fold
 * (x=0, p=0). Continuing from p=-1, x=-1 (one branch) should march the
 * parameter toward 0 and detect a fold near p=0. */
static void test_continuation_fold() {
  std::printf("continuation: scalar fold  x' = p + x^2\n");
  Model m;
  m.n = 1;
  m.vector_field = [](const double *x, double p, double *f,
                      std::string *) -> bool {
    f[0] = p + x[0] * x[0];
    return true;
  };
  ContinuationSettings s;
  s.h0 = 0.05;
  s.max_points = 400;
  s.p_min = -2.0;
  s.p_max = 2.0;
  Branch b = continue_equilibrium(m, {-1.0}, -1.0, s);
  check(b.ok, "branch produced");
  bool found_fold = false;
  double fold_p = 0.0;
  for (std::size_t i : b.special_indices) {
    if (b.points[i].special == SpecialPointKind::Fold) {
      found_fold = true;
      fold_p = b.points[i].p;
    }
  }
  check(found_fold, "fold detected");
  if (found_fold) check(close(fold_p, 0.0, 5e-2), "fold near p=0");
}

/* ---- continuation: Hopf -------------------------------------
 * 2D system (a Hopf normal form with linear parameter):
 *   x' = p*x - y - x*(x^2+y^2)
 *   y' = x + p*y - y*(x^2+y^2)
 * The origin is always an equilibrium; its Jacobian at 0 is
 * [[p,-1],[1,p]] with eigenvalues p +/- i, so a Hopf bifurcation
 * occurs at p=0. Continuing the origin in p should flag a Hopf near
 * p=0. */
static void test_continuation_hopf() {
  std::printf("continuation: Hopf  (origin of normal form)\n");
  Model m;
  m.n = 2;
  m.vector_field = [](const double *x, double p, double *f,
                      std::string *) -> bool {
    const double r2 = x[0] * x[0] + x[1] * x[1];
    f[0] = p * x[0] - x[1] - x[0] * r2;
    f[1] = x[0] + p * x[1] - x[1] * r2;
    return true;
  };
  ContinuationSettings s;
  s.h0 = 0.05;
  s.max_points = 400;
  s.p_min = -1.0;
  s.p_max = 1.0;
  s.detect_fold = false;  /* origin is never a fold here */
  Branch b = continue_equilibrium(m, {0.0, 0.0}, -0.5, s);
  check(b.ok, "branch produced");
  bool found_hopf = false;
  double hopf_p = 0.0;
  for (std::size_t i : b.special_indices) {
    if (b.points[i].special == SpecialPointKind::Hopf) {
      found_hopf = true;
      hopf_p = b.points[i].p;
    }
  }
  check(found_hopf, "Hopf detected");
  if (found_hopf) check(close(hopf_p, 0.0, 5e-2), "Hopf near p=0");
}

int main() {
  std::printf("=== dynsys analysis smoke test ===\n");
  test_eigen_diagonal();
  test_eigen_rotation();
  test_eigen_known_3x3();
  test_determinant();
  test_classify_saddle();
  test_classify_stable_spiral();
  test_continuation_fold();
  test_continuation_hopf();
  std::printf("=== %d/%d checks passed ===\n", g_checks - g_failures, g_checks);
  return g_failures == 0 ? 0 : 1;
}
