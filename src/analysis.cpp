/* ============================================================
 * dynsys analysis core implementation. See analysis.h.
 *
 * The eigensolver is a from-scratch real general-matrix solver:
 * Parlett-Reinsch balancing, Householder reduction to upper
 * Hessenberg, then the Francis double-shift QR algorithm (the
 * classic "hqr" structure from Wilkinson / EISPACK, rewritten in
 * modern C++). It returns complex eigenvalues and needs no LAPACK.
 *
 * The continuation engine solves the augmented system
 *     G(z) = [ f(x, p) ; n . (z - z_prev) - h ] = 0,
 * z = [x; p], by a predictor (tangent) / corrector (bordered
 * Newton) scheme with adaptive arclength step h, watching test
 * functions for fold (det f_x = 0) and Hopf (a complex pair
 * crossing the imaginary axis).
 * ============================================================ */

#include "analysis.h"

#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <thread>
#include <atomic>

namespace dynsys::analysis {

namespace {

/* Run body(i) for i in [0,count) across hardware threads. Falls back to a plain
 * serial loop when only one core is available or count is tiny, so behaviour is
 * identical (just faster) on multi-core machines. body must be thread-safe with
 * respect to disjoint i (each call should touch only its own slice of output).*/
template <typename F>
void parallel_for(int count, F &&body) {
  if (count <= 0) return;
  unsigned hw = std::thread::hardware_concurrency();
  if (hw < 2 || count < 64) { for (int i = 0; i < count; ++i) body(i); return; }
  unsigned nthreads = std::min<unsigned>(hw, (unsigned)count);
  std::atomic<int> next{0};
  auto worker = [&]() {
    for (;;) {
      int i = next.fetch_add(1, std::memory_order_relaxed);
      if (i >= count) break;
      body(i);
    }
  };
  std::vector<std::thread> pool;
  pool.reserve(nthreads);
  for (unsigned t = 0; t < nthreads; ++t) pool.emplace_back(worker);
  for (auto &th : pool) th.join();
}

constexpr double kRadix = 2.0; /* floating-point base for balancing */

bool is_finite_vec(const std::vector<double> &v) {
  for (double x : v)
    if (!std::isfinite(x)) return false;
  return true;
}

/* ---- Parlett-Reinsch balancing (in place) ------------------- */
void balance(std::vector<double> &a, std::size_t n) {
  const double b2 = kRadix * kRadix;
  bool last = false;
  while (!last) {
    last = true;
    for (std::size_t i = 0; i < n; ++i) {
      double r = 0.0, c = 0.0;
      for (std::size_t j = 0; j < n; ++j) {
        if (j != i) {
          c += std::fabs(a[j * n + i]);
          r += std::fabs(a[i * n + j]);
        }
      }
      if (c == 0.0 || r == 0.0) continue;
      double g = r / kRadix;
      double f = 1.0;
      const double s = c + r;
      while (c < g) {
        f *= kRadix;
        c *= b2;
      }
      g = r * kRadix;
      while (c >= g) {
        f /= kRadix;
        c /= b2;
      }
      if ((c + r) / f < 0.95 * s) {
        last = false;
        const double gg = 1.0 / f;
        for (std::size_t j = 0; j < n; ++j) a[i * n + j] *= gg;
        for (std::size_t j = 0; j < n; ++j) a[j * n + i] *= f;
      }
    }
  }
}

/* ---- reduction to upper Hessenberg (in place) --------------- */
void to_hessenberg(std::vector<double> &a, std::size_t n) {
  for (std::size_t m = 1; m + 1 < n; ++m) {
    double x = 0.0;
    std::size_t i = m;
    for (std::size_t j = m; j < n; ++j) {
      if (std::fabs(a[j * n + (m - 1)]) > std::fabs(x)) {
        x = a[j * n + (m - 1)];
        i = j;
      }
    }
    if (i != m) {
      for (std::size_t j = m - 1; j < n; ++j)
        std::swap(a[i * n + j], a[m * n + j]);
      for (std::size_t j = 0; j < n; ++j)
        std::swap(a[j * n + i], a[j * n + m]);
    }
    if (x != 0.0) {
      for (i = m + 1; i < n; ++i) {
        double y = a[i * n + (m - 1)];
        if (y != 0.0) {
          y /= x;
          a[i * n + (m - 1)] = y;
          for (std::size_t j = m; j < n; ++j) a[i * n + j] -= y * a[m * n + j];
          for (std::size_t j = 0; j < n; ++j) a[j * n + m] += y * a[j * n + i];
        }
      }
    }
  }
}

/* ---- Francis double-shift QR on an upper Hessenberg matrix --- *
 * Adapted from the classic EISPACK hqr / Numerical Recipes hqr
 * structure. Eigenvalues are accumulated into `w`. Returns false
 * if any eigenvalue fails to converge within the iteration cap. */
bool hqr(std::vector<double> &a, std::size_t n, std::vector<Complex> *w) {
  w->assign(n, Complex(0.0, 0.0));
  double anorm = 0.0;
  for (std::size_t i = 0; i < n; ++i)
    for (std::size_t j = (i == 0 ? 0 : i - 1); j < n; ++j)
      anorm += std::fabs(a[i * n + j]);

  int nn = static_cast<int>(n) - 1;
  double t = 0.0;
  bool converged = true;

  while (nn >= 0) {
    int its = 0;
    int l;
    do {
      for (l = nn; l >= 1; --l) {
        double s = std::fabs(a[(l - 1) * n + (l - 1)]) +
                   std::fabs(a[l * n + l]);
        if (s == 0.0) s = anorm;
        if (std::fabs(a[l * n + (l - 1)]) + s == s) {
          a[l * n + (l - 1)] = 0.0;
          break;
        }
      }
      double x = a[nn * n + nn];
      if (l == nn) {
        (*w)[nn] = Complex(x + t, 0.0);
        --nn;
      } else {
        double y = a[(nn - 1) * n + (nn - 1)];
        double ww = a[nn * n + (nn - 1)] * a[(nn - 1) * n + nn];
        if (l == nn - 1) {
          double p = 0.5 * (y - x);
          double q = p * p + ww;
          double z = std::sqrt(std::fabs(q));
          x += t;
          if (q >= 0.0) {
            z = p + (p >= 0.0 ? std::fabs(z) : -std::fabs(z));
            (*w)[nn - 1] = Complex(x + z, 0.0);
            (*w)[nn] = Complex(z != 0.0 ? x - ww / z : x + z, 0.0);
          } else {
            (*w)[nn - 1] = Complex(x + p, z);
            (*w)[nn] = Complex(x + p, -z);
          }
          nn -= 2;
        } else {
          if (its == 60) {
            converged = false;
            break;
          }
          double p = 0.0, q = 0.0, r = 0.0;
          if (its == 10 || its == 20) {
            t += x;
            for (int i = 0; i <= nn; ++i) a[i * n + i] -= x;
            double s = std::fabs(a[nn * n + (nn - 1)]) +
                       std::fabs(a[(nn - 1) * n + (nn - 2)]);
            y = x = 0.75 * s;
            ww = -0.4375 * s * s;
          }
          ++its;
          int m;
          for (m = nn - 2; m >= l; --m) {
            double z = a[m * n + m];
            r = x - z;
            double s = y - z;
            p = (r * s - ww) / a[(m + 1) * n + m] + a[m * n + (m + 1)];
            q = a[(m + 1) * n + (m + 1)] - z - r - s;
            r = a[(m + 2) * n + (m + 1)];
            s = std::fabs(p) + std::fabs(q) + std::fabs(r);
            p /= s;
            q /= s;
            r /= s;
            if (m == l) break;
            double u = std::fabs(a[m * n + (m - 1)]) * (std::fabs(q) + std::fabs(r));
            double v = std::fabs(p) *
                       (std::fabs(a[(m - 1) * n + (m - 1)]) + std::fabs(z) +
                        std::fabs(a[(m + 1) * n + (m + 1)]));
            if (u + v == v) break;
          }
          for (int i = m + 2; i <= nn; ++i) {
            a[i * n + (i - 2)] = 0.0;
            if (i != m + 2) a[i * n + (i - 3)] = 0.0;
          }
          for (int k = m; k <= nn - 1; ++k) {
            if (k != m) {
              p = a[k * n + (k - 1)];
              q = a[(k + 1) * n + (k - 1)];
              r = 0.0;
              if (k != nn - 1) r = a[(k + 2) * n + (k - 1)];
              x = std::fabs(p) + std::fabs(q) + std::fabs(r);
              if (x != 0.0) {
                p /= x;
                q /= x;
                r /= x;
              }
            }
            double s = std::sqrt(p * p + q * q + r * r);
            if (p < 0.0) s = -s;
            if (s != 0.0) {
              if (k == m) {
                if (l != m) a[k * n + (k - 1)] = -a[k * n + (k - 1)];
              } else {
                a[k * n + (k - 1)] = -s * x;
              }
              p += s;
              double xx = p / s;
              double yy = q / s;
              double zz = r / s;
              q /= p;
              r /= p;
              for (int j = k; j <= nn; ++j) {
                p = a[k * n + j] + q * a[(k + 1) * n + j];
                if (k != nn - 1) {
                  p += r * a[(k + 2) * n + j];
                  a[(k + 2) * n + j] -= p * zz;
                }
                a[(k + 1) * n + j] -= p * yy;
                a[k * n + j] -= p * xx;
              }
              int mmin = nn < k + 3 ? nn : k + 3;
              for (int i = l; i <= mmin; ++i) {
                p = xx * a[i * n + k] + yy * a[i * n + (k + 1)];
                if (k != nn - 1) {
                  p += zz * a[i * n + (k + 2)];
                  a[i * n + (k + 2)] -= p * r;
                }
                a[i * n + (k + 1)] -= p * q;
                a[i * n + k] -= p;
              }
            }
          }
        }
      }
    } while (l < nn - 1 && nn >= 0 && converged);
    if (!converged) break;
  }
  return converged;
}

}  // namespace

/* ---- public linear algebra ---------------------------------- */

bool solve_linear(const std::vector<double> &A, const std::vector<double> &b,
                  std::vector<double> *x) {
  const std::size_t n = b.size();
  if (A.size() != n * n || n == 0) return false;
  std::vector<double> m(n * (n + 1), 0.0);
  for (std::size_t r = 0; r < n; ++r) {
    for (std::size_t c = 0; c < n; ++c) m[r * (n + 1) + c] = A[r * n + c];
    m[r * (n + 1) + n] = b[r];
  }
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    for (std::size_t r = col + 1; r < n; ++r)
      if (std::fabs(m[r * (n + 1) + col]) > std::fabs(m[pivot * (n + 1) + col]))
        pivot = r;
    if (std::fabs(m[pivot * (n + 1) + col]) < 1e-300) return false;
    if (pivot != col)
      for (std::size_t c = col; c <= n; ++c)
        std::swap(m[pivot * (n + 1) + c], m[col * (n + 1) + c]);
    const double div = m[col * (n + 1) + col];
    for (std::size_t c = col; c <= n; ++c) m[col * (n + 1) + c] /= div;
    for (std::size_t r = 0; r < n; ++r)
      if (r != col) {
        const double factor = m[r * (n + 1) + col];
        for (std::size_t c = col; c <= n; ++c)
          m[r * (n + 1) + c] -= factor * m[col * (n + 1) + c];
      }
  }
  x->assign(n, 0.0);
  for (std::size_t r = 0; r < n; ++r) (*x)[r] = m[r * (n + 1) + n];
  return true;
}

double determinant(const std::vector<double> &A, std::size_t n) {
  if (A.size() != n * n || n == 0) return 0.0;
  std::vector<double> m(A);
  double det = 1.0;
  for (std::size_t col = 0; col < n; ++col) {
    std::size_t pivot = col;
    for (std::size_t r = col + 1; r < n; ++r)
      if (std::fabs(m[r * n + col]) > std::fabs(m[pivot * n + col])) pivot = r;
    if (std::fabs(m[pivot * n + col]) < 1e-300) return 0.0;
    if (pivot != col) {
      for (std::size_t c = 0; c < n; ++c)
        std::swap(m[pivot * n + c], m[col * n + c]);
      det = -det;
    }
    det *= m[col * n + col];
    const double inv = 1.0 / m[col * n + col];
    for (std::size_t r = col + 1; r < n; ++r) {
      const double factor = m[r * n + col] * inv;
      if (factor != 0.0)
        for (std::size_t c = col; c < n; ++c)
          m[r * n + c] -= factor * m[col * n + c];
    }
  }
  return det;
}

bool eigenvalues(const std::vector<double> &A, std::size_t n,
                 std::vector<Complex> *out) {
  if (A.size() != n * n || n == 0) return false;
  if (n == 1) {
    out->assign(1, Complex(A[0], 0.0));
    return true;
  }
  std::vector<double> a(A);
  balance(a, n);
  to_hessenberg(a, n);
  return hqr(a, n, out);
}

/* ---- classification ----------------------------------------- */

Classification classify_equilibrium(const std::vector<double> &jacobian,
                                    std::size_t n, double tol) {
  Classification cl;
  if (!eigenvalues(jacobian, n, &cl.eigenvalues)) {
    cl.label = "eigenvalue computation failed";
    return cl;
  }
  bool has_complex_marginal = false;
  bool has_real_marginal = false;
  for (const Complex &z : cl.eigenvalues) {
    const double re = z.real();
    if (re < -tol)
      ++cl.n_stable;
    else if (re > tol)
      ++cl.n_unstable;
    else {
      ++cl.n_center;
      if (std::fabs(z.imag()) > tol)
        has_complex_marginal = true;
      else
        has_real_marginal = true;
    }
  }
  cl.is_saddle = cl.n_stable > 0 && cl.n_unstable > 0;
  cl.hopf_candidate = has_complex_marginal;
  cl.fold_candidate = has_real_marginal;

  /* Build a label. In 2D we keep the familiar node/spiral words; in
   * higher dimensions we report the stable/unstable split, which is
   * the honest dimension-agnostic description. */
  bool any_complex = false;
  for (const Complex &z : cl.eigenvalues)
    if (std::fabs(z.imag()) > tol) any_complex = true;

  if (cl.n_center > 0 && cl.n_stable == 0 && cl.n_unstable == 0) {
    cl.label = has_complex_marginal ? "center (marginal)" : "non-hyperbolic";
  } else if (cl.is_saddle) {
    cl.label = "saddle (" + std::to_string(cl.n_stable) + " stable, " +
               std::to_string(cl.n_unstable) + " unstable)";
  } else if (cl.n_unstable == 0) {
    cl.label = any_complex ? "stable spiral/focus" : "stable node";
  } else {
    cl.label = any_complex ? "unstable spiral/focus" : "unstable node";
  }
  if (cl.hopf_candidate) cl.label += " [Hopf candidate]";
  if (cl.fold_candidate) cl.label += " [fold candidate]";
  return cl;
}

/* ---- finite-difference Jacobian ----------------------------- */

bool finite_diff_jacobian(const Model &m, const double *x, double p,
                          std::vector<double> *jac_out, std::string *err,
                          double eps) {
  const std::size_t n = m.n;
  jac_out->assign(n * n, 0.0);
  std::vector<double> xp(x, x + n), fp(n), fm(n);
  for (std::size_t col = 0; col < n; ++col) {
    const double h = eps * (std::fabs(x[col]) + eps);
    xp[col] = x[col] + h;
    if (!m.vector_field(xp.data(), p, fp.data(), err)) return false;
    xp[col] = x[col] - h;
    if (!m.vector_field(xp.data(), p, fm.data(), err)) return false;
    xp[col] = x[col];
    const double inv = 1.0 / (2.0 * h);
    for (std::size_t row = 0; row < n; ++row)
      (*jac_out)[row * n + col] = (fp[row] - fm[row]) * inv;
  }
  return true;
}

namespace {

bool jacobian_x(const Model &m, const double *x, double p,
                std::vector<double> *jac, std::string *err) {
  if (m.jacobian_x) {
    /* The caller may pass an empty vector; the exact callback writes into
     * jac->data(), so it MUST be sized n*n first (otherwise data() is null
     * or undersized -> out-of-bounds write / crash). The finite-difference
     * fallback resizes internally, which masked this on the FD path. */
    jac->assign(m.n * m.n, 0.0);
    return m.jacobian_x(x, p, jac->data(), err);
  }
  return finite_diff_jacobian(m, x, p, jac, err);
}

bool dfdp(const Model &m, const double *x, double p, std::vector<double> *out,
          std::string *err) {
  const std::size_t n = m.n;
  out->assign(n, 0.0);
  if (m.dfdp) return m.dfdp(x, p, out->data(), err);
  const double h = 1e-6 * (std::fabs(p) + 1e-6);
  std::vector<double> fp(n), fm(n);
  if (!m.vector_field(x, p + h, fp.data(), err)) return false;
  if (!m.vector_field(x, p - h, fm.data(), err)) return false;
  const double inv = 1.0 / (2.0 * h);
  for (std::size_t i = 0; i < n; ++i) (*out)[i] = (fp[i] - fm[i]) * inv;
  return true;
}

/* Newton-correct an equilibrium at fixed p (used to clean up the
 * starting point before continuation). */
bool correct_equilibrium(const Model &m, std::vector<double> *x, double p,
                         int max_iters, double tol, std::string *err) {
  const std::size_t n = m.n;
  std::vector<double> f(n), J, delta;
  for (int it = 0; it < max_iters; ++it) {
    if (!m.vector_field(x->data(), p, f.data(), err)) return false;
    double res = 0.0;
    for (double v : f) res += v * v;
    if (std::sqrt(res) < tol) return true;
    if (!jacobian_x(m, x->data(), p, &J, err)) return false;
    std::vector<double> rhs(n);
    for (std::size_t i = 0; i < n; ++i) rhs[i] = -f[i];
    if (!solve_linear(J, rhs, &delta)) {
      if (err) *err = "singular Jacobian during equilibrium correction";
      return false;
    }
    for (std::size_t i = 0; i < n; ++i) (*x)[i] += delta[i];
  }
  /* One last residual check. */
  if (!m.vector_field(x->data(), p, f.data(), err)) return false;
  double res = 0.0;
  for (double v : f) res += v * v;
  return std::sqrt(res) < tol * 1e3;
}

/* Fold test function: det(f_x). Sign change => a real eigenvalue
 * crossed zero. */
double fold_test(const Model &m, const double *x, double p) {
  std::vector<double> J;
  std::string err;
  if (!jacobian_x(m, x, p, &J, &err)) return std::nan("");
  return determinant(J, m.n);
}

/* Hopf test function: product over complex-conjugate pairs of the
 * sum of their real parts. We use a robust proxy: the minimum
 * |Re(lambda)| among complex (nonzero-imag) eigenvalues, signed by
 * the sign of that Re. A sign change with imag != 0 brackets a Hopf.
 * Returns NaN when there is no complex pair. */
double hopf_test(const Model &m, const double *x, double p, double tol) {
  std::vector<double> J;
  std::string err;
  if (!jacobian_x(m, x, p, &J, &err)) return std::nan("");
  std::vector<Complex> ev;
  if (!eigenvalues(J, m.n, &ev)) return std::nan("");
  bool found = false;
  double best_abs = 0.0, signed_val = 0.0;
  for (const Complex &z : ev) {
    if (std::fabs(z.imag()) > tol) {
      const double a = std::fabs(z.real());
      if (!found || a < best_abs) {
        best_abs = a;
        signed_val = z.real();
        found = true;
      }
    }
  }
  return found ? signed_val : std::nan("");
}

/* Branch-point test function. A branch point (where two equilibrium branches
 * cross) is a point on the curve where the bordered (n+1)x(n+1) matrix
 * [ f_x  f_p ; v^T ] (v = curve tangent) becomes singular, WHILE the
 * equilibrium itself persists. At a regular point this determinant is nonzero;
 * at a fold it stays nonzero (a fold is a turning point, not a crossing). So a
 * sign change of this determinant that is NOT accompanied by a fold-test sign
 * change brackets a branch point. Returns NaN on evaluation failure. */
double branch_point_test(const Model &m, const double *x, double p,
                         const std::vector<double> &tangent, std::size_t n) {
  std::vector<double> J, fp;
  std::string err;
  if (!jacobian_x(m, x, p, &J, &err)) return std::nan("");
  if (!dfdp(m, x, p, &fp, &err)) return std::nan("");
  const std::size_t N = n + 1;
  std::vector<double> A(N * N, 0.0);
  for (std::size_t r = 0; r < n; ++r) {
    for (std::size_t c = 0; c < n; ++c) A[r * N + c] = J[r * n + c];
    A[r * N + n] = fp[r];
  }
  /* last row: the curve tangent (length N). If unavailable, fall back to e_last. */
  if (tangent.size() == N)
    for (std::size_t c = 0; c < N; ++c) A[n * N + c] = tangent[c];
  else
    A[n * N + n] = 1.0;
  return determinant(A, N);
}

}  // namespace

/* ---- pseudo-arclength continuation -------------------------- */

Branch continue_equilibrium(const Model &m, const std::vector<double> &x0,
                            double p0, const ContinuationSettings &settings) {
  Branch branch;
  const std::size_t n = m.n;
  if (x0.size() != n || n == 0) {
    branch.message = "starting point dimension mismatch";
    return branch;
  }

  std::string err;

  /* Clean up the starting equilibrium. */
  std::vector<double> x = x0;
  if (!correct_equilibrium(m, &x, p0, settings.max_corrector_iters,
                           settings.corrector_tol, &err)) {
    branch.message =
        "no equilibrium found from this start — Newton did not converge (" +
        (err.empty() ? std::string("residual too large") : err) +
        "). Try 'Find fixed point' first, or note that some systems "
        "(e.g. conservative ones like Nose-Hoover) have no isolated equilibrium.";
    return branch;
  }

  /* z = [x; p], length n+1. */
  const std::size_t N = n + 1;
  std::vector<double> z(N);
  for (std::size_t i = 0; i < n; ++i) z[i] = x[i];
  z[n] = p0;

  auto record_point = [&](const std::vector<double> &zz,
                          SpecialPointKind kind) {
    BranchPoint bp;
    bp.p = zz[n];
    bp.x.assign(zz.begin(), zz.begin() + n);
    std::vector<double> J;
    if (jacobian_x(m, bp.x.data(), bp.p, &J, &err)) {
      Classification cl = classify_equilibrium(J, n, settings.event_tol * 1e2);
      bp.eigenvalues = cl.eigenvalues;
      bp.n_unstable = cl.n_unstable;
      bp.stable = cl.n_unstable == 0 && cl.n_center == 0;
    }
    bp.special = kind;
    /* Classify criticality + upgrade to codim-2 at special points. */
    if (kind == SpecialPointKind::Fold) {
      double a = 0.0, l0 = 0.0; std::string e2;
      if (fold_normal_form(m, bp.x, bp.p, &a, &l0, &e2)) {
        bp.fold_a = a; bp.has_normal_form = true;
        /* Bogdanov-Takens: a SECOND eigenvalue also near zero (a double zero
         * eigenvalue) turns the fold into a BT point. Cusp: the fold's own
         * quadratic coefficient a vanishes. BT takes precedence. */
        int near_zero = 0;
        for (const Complex &z : bp.eigenvalues)
          if (std::abs(z) < 1e-2) ++near_zero;
        if (near_zero >= 2) bp.special = SpecialPointKind::BogdanovTakens;
        else if (std::fabs(a) < 1e-2) bp.special = SpecialPointKind::Cusp;
      }
    } else if (kind == SpecialPointKind::Hopf) {
      double l1 = 0.0, w = 0.0; std::string e2;
      if (hopf_first_lyapunov(m, bp.x, bp.p, &l1, &w, &e2)) {
        bp.lyapunov1 = l1; bp.has_normal_form = true;
        /* Generalized Hopf (Bautin): the first Lyapunov coefficient vanishes. */
        if (std::fabs(l1) < 1e-3) bp.special = SpecialPointKind::GeneralizedHopf;
      }
    }
    if (bp.special != SpecialPointKind::None)
      branch.special_indices.push_back(branch.points.size());
    branch.points.push_back(std::move(bp));
  };

  /* Initial tangent: solve the bordered system for the null vector of
   * [ f_x | f_p ]. We build [f_x f_p; e_last^T] tangent = e_last and
   * normalize. The last component sign sets the initial direction. */
  auto compute_tangent = [&](const std::vector<double> &zz,
                             const std::vector<double> &prev_tan,
                             std::vector<double> *tan) -> bool {
    std::vector<double> J, fp;
    std::vector<double> xx(zz.begin(), zz.begin() + n);
    if (!jacobian_x(m, xx.data(), zz[n], &J, &err)) return false;
    if (!dfdp(m, xx.data(), zz[n], &fp, &err)) return false;
    /* Augmented (n+1)x(n+1): rows 0..n-1 are [J | fp];
     * last row is the previous tangent (or e_last on the first call)
     * to fix the scale and pick a consistent direction. */
    std::vector<double> A(N * N, 0.0), b(N, 0.0);
    for (std::size_t r = 0; r < n; ++r) {
      for (std::size_t c = 0; c < n; ++c) A[r * N + c] = J[r * n + c];
      A[r * N + n] = fp[r];
    }
    if (prev_tan.empty()) {
      A[n * N + n] = 1.0;  /* e_last */
      b[n] = 1.0;
    } else {
      for (std::size_t c = 0; c < N; ++c) A[n * N + c] = prev_tan[c];
      b[n] = 1.0;
    }
    if (!solve_linear(A, b, tan)) return false;
    double norm = 0.0;
    for (double v : *tan) norm += v * v;
    norm = std::sqrt(norm);
    if (norm < 1e-300) return false;
    for (double &v : *tan) v /= norm;
    /* Keep direction consistent with the previous tangent. */
    if (!prev_tan.empty()) {
      double dot = 0.0;
      for (std::size_t i = 0; i < N; ++i) dot += (*tan)[i] * prev_tan[i];
      if (dot < 0.0)
        for (double &v : *tan) v = -v;
    } else {
      if (((*tan)[n] < 0.0 && settings.direction > 0) ||
          ((*tan)[n] > 0.0 && settings.direction < 0))
        for (double &v : *tan) v = -v;
    }
    return true;
  };

  std::vector<double> tangent;
  std::vector<double> empty_prev;
  if (!compute_tangent(z, empty_prev, &tangent)) {
    branch.message = "could not compute initial tangent (singular system)";
    return branch;
  }

  record_point(z, SpecialPointKind::None);
  double prev_fold = fold_test(m, x.data(), p0);
  double prev_hopf = hopf_test(m, x.data(), p0, settings.event_tol * 1e2);
  double prev_bp = branch_point_test(m, x.data(), p0, tangent, n);

  double h = settings.h0;

  /* Corrector: Newton on G(z) = [f(x,p); tangent.(z - z_pred) - 0]
   * where the predicted point is z_pred = z + h*tangent and the
   * pseudo-arclength condition pins the step. */
  auto corrector = [&](std::vector<double> *zc,
                       const std::vector<double> &z_pred,
                       const std::vector<double> &tan) -> bool {
    std::vector<double> J, fp, f(n);
    for (int it = 0; it < settings.max_corrector_iters; ++it) {
      std::vector<double> xx(zc->begin(), zc->begin() + n);
      const double pp = (*zc)[n];
      if (!m.vector_field(xx.data(), pp, f.data(), &err)) return false;
      double arc = 0.0;
      for (std::size_t i = 0; i < N; ++i) arc += tan[i] * ((*zc)[i] - z_pred[i]);
      double res = arc * arc;
      for (double v : f) res += v * v;
      if (std::sqrt(res) < settings.corrector_tol) return true;
      if (!jacobian_x(m, xx.data(), pp, &J, &err)) return false;
      if (!dfdp(m, xx.data(), pp, &fp, &err)) return false;
      std::vector<double> A(N * N, 0.0), b(N, 0.0);
      for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t c = 0; c < n; ++c) A[r * N + c] = J[r * n + c];
        A[r * N + n] = fp[r];
        b[r] = -f[r];
      }
      for (std::size_t c = 0; c < N; ++c) A[n * N + c] = tan[c];
      b[n] = -arc;
      std::vector<double> delta;
      if (!solve_linear(A, b, &delta)) return false;
      for (std::size_t i = 0; i < N; ++i) (*zc)[i] += delta[i];
    }
    return false;
  };

  /* Refine an event bracketed between z_lo and z_hi (both corrected
   * branch points) for a scalar test function `g`. Bisects the secant
   * parameter, re-correcting each candidate onto the branch with the
   * current tangent so the located event sits exactly on the solution
   * curve, not on the chord. Writes the located point into z_out. */
  auto refine_event =
      [&](const std::vector<double> &z_lo, const std::vector<double> &z_hi,
          const std::function<double(const double *, double)> &g,
          const std::vector<double> &tan,
          std::vector<double> *z_out) -> bool {
    std::vector<double> a = z_lo, b = z_hi;
    auto eval_g = [&](const std::vector<double> &zz) {
      std::vector<double> xx(zz.begin(), zz.begin() + n);
      return g(xx.data(), zz[n]);
    };
    double ga = eval_g(a);
    double gb = eval_g(b);
    if (!std::isfinite(ga) || !std::isfinite(gb) || ga * gb > 0.0) {
      *z_out = z_hi;
      return false;
    }
    for (int it = 0; it < 60; ++it) {
      std::vector<double> mid(N);
      for (std::size_t i = 0; i < N; ++i) mid[i] = 0.5 * (a[i] + b[i]);
      std::vector<double> z_pred = mid;
      /* Correct the midpoint back onto the branch using a fresh
       * arclength anchor at the midpoint itself (h = 0). */
      corrector(&mid, z_pred, tan);
      const double gm = eval_g(mid);
      if (!std::isfinite(gm)) {
        *z_out = mid;
        return true;
      }
      if (std::fabs(gm) <= settings.event_tol ||
          std::fabs(b[n] - a[n]) <= settings.event_tol) {
        *z_out = mid;
        return true;
      }
      if (ga * gm < 0.0) {
        b = mid;
        gb = gm;
      } else {
        a = mid;
        ga = gm;
      }
    }
    for (std::size_t i = 0; i < N; ++i) (*z_out)[i] = 0.5 * (a[i] + b[i]);
    return true;
  };

  int produced = 1;
  while (produced < settings.max_points) {
    std::vector<double> z_pred(N);
    for (std::size_t i = 0; i < N; ++i) z_pred[i] = z[i] + h * tangent[i];
    std::vector<double> z_new = z_pred;
    if (!corrector(&z_new, z_pred, tangent) || !is_finite_vec(z_new)) {
      /* Shrink step and retry. */
      h *= 0.5;
      if (h < settings.h_min) {
        branch.message = "corrector failed; step underflow";
        break;
      }
      continue;
    }

    const double p_new = z_new[n];
    if (p_new < settings.p_min || p_new > settings.p_max) {
      branch.message = "reached parameter bound";
      record_point(z_new, SpecialPointKind::EndOfBranch);
      break;
    }

    /* Event detection between z (old) and z_new. When a test
     * function changes sign we bisect to locate the event precisely
     * and record it as its own tagged point, then record the regular
     * continuation point. */
    std::vector<double> x_new(z_new.begin(), z_new.begin() + n);

    /* Compute the branch-point test up front so the fold block can tell a real
     * fold from a branch point (at a pitchfork the fold test also vanishes, but
     * the point is a BP, not a fold/cusp). */
    const double bp_new_pre = branch_point_test(m, x_new.data(), p_new, tangent, n);
    const bool bp_coincides = std::isfinite(prev_bp) && std::isfinite(bp_new_pre) &&
                              prev_bp * bp_new_pre < 0.0;

    if (settings.detect_fold) {
      const double f_new = fold_test(m, x_new.data(), p_new);
      if (std::isfinite(prev_fold) && std::isfinite(f_new) &&
          prev_fold * f_new < 0.0 && !bp_coincides) {
        std::vector<double> z_ev = z_new;
        auto g = [&](const double *xx, double pp) {
          return fold_test(m, xx, pp);
        };
        refine_event(z, z_new, g, tangent, &z_ev);
        record_point(z_ev, SpecialPointKind::Fold);
      }
      prev_fold = f_new;
    }
    if (settings.detect_hopf) {
      const double hp_new = hopf_test(m, x_new.data(), p_new,
                                      settings.event_tol * 1e2);
      if (std::isfinite(prev_hopf) && std::isfinite(hp_new) &&
          prev_hopf * hp_new < 0.0) {
        std::vector<double> z_ev = z_new;
        const double htol = settings.event_tol * 1e2;
        auto g = [&, htol](const double *xx, double pp) {
          return hopf_test(m, xx, pp, htol);
        };
        refine_event(z, z_new, g, tangent, &z_ev);
        record_point(z_ev, SpecialPointKind::Hopf);
      }
      prev_hopf = hp_new;
    }

    /* Branch-point detection: the augmented-system determinant changes sign
     * WITHOUT a simultaneous fold (which would also flip it). When found,
     * locate it, compute the SECOND tangent (the other branch's direction) for
     * branch switching, and tag the point. */
    {
      const double bp_new = bp_new_pre;
      if (std::isfinite(prev_bp) && std::isfinite(bp_new) && prev_bp * bp_new < 0.0) {
        std::vector<double> z_ev = z_new;
        auto g = [&](const double *xx, double pp) {
          return branch_point_test(m, xx, pp, tangent, n);
        };
        refine_event(z, z_new, g, tangent, &z_ev);
        /* second tangent: a direction in the null space of [f_x f_p] that is
         * independent of the current curve tangent. Solve the bordered system
         * with the RHS picking out a complementary direction, then orthogonal-
         * ize against the primary tangent. */
        std::vector<double> x_ev(z_ev.begin(), z_ev.begin() + n);
        std::vector<double> Jb, fpb;
        std::vector<double> t2;
        if (jacobian_x(m, x_ev.data(), z_ev[n], &Jb, &err) &&
            dfdp(m, x_ev.data(), z_ev[n], &fpb, &err)) {
          /* M = [f_x | f_p] is n x (n+1). At a branch point it has a 2-D null
           * space spanned by the curve tangent and the crossing direction. We
           * find the crossing direction as the smallest eigenvector of the
           * (n+1)x(n+1) Gram matrix M^T M, with the primary tangent projected
           * out (so we get the transverse null vector, not the tangent). */
          const std::size_t Nn = n + 1;
          std::vector<double> M(n * Nn, 0.0);
          for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t c = 0; c < n; ++c) M[r * Nn + c] = Jb[r * n + c];
            M[r * Nn + n] = fpb[r];
          }
          std::vector<double> G(Nn * Nn, 0.0);
          for (std::size_t i = 0; i < Nn; ++i)
            for (std::size_t j = 0; j < Nn; ++j) {
              double s = 0.0;
              for (std::size_t r = 0; r < n; ++r) s += M[r * Nn + i] * M[r * Nn + j];
              G[i * Nn + j] = s;
            }
          /* deflate the primary tangent: G' = G + C * (t t^T) with large C so
           * the tangent direction becomes the LARGEST, leaving the crossing
           * direction as the smallest eigenvector. */
          double gmax = 0.0; for (double v : G) gmax = std::max(gmax, std::fabs(v));
          const double C = (gmax + 1.0) * 1e3;
          for (std::size_t i = 0; i < Nn; ++i)
            for (std::size_t j = 0; j < Nn; ++j)
              G[i * Nn + j] += C * tangent[i] * tangent[j];
          /* smallest eigenvector of G' via inverse power iteration on (G'+epsI) */
          std::vector<double> Gr = G;
          for (std::size_t i = 0; i < Nn; ++i) Gr[i * Nn + i] += 1e-9;
          std::vector<double> v(Nn, 0.0); v[ (n>0?n:0) ] = 1.0; /* seed toward the p-direction */
          bool ok_it = true;
          for (int it = 0; it < 100; ++it) {
            std::vector<double> nv;
            if (!solve_linear(Gr, v, &nv)) { ok_it = false; break; }
            double nrm = 0.0; for (double q : nv) nrm += q * q; nrm = std::sqrt(nrm);
            if (nrm < 1e-300) { ok_it = false; break; }
            for (double &q : nv) q /= nrm;
            v = nv;
          }
          if (ok_it) {
            /* remove any residual primary-tangent component and normalize */
            double dot = 0.0; for (std::size_t i = 0; i < Nn; ++i) dot += v[i] * tangent[i];
            for (std::size_t i = 0; i < Nn; ++i) v[i] -= dot * tangent[i];
            double nrm = 0.0; for (double q : v) nrm += q * q; nrm = std::sqrt(nrm);
            if (nrm > 1e-6) { for (double &q : v) q /= nrm; t2 = v; }
          }
        }
        const std::size_t idx_before = branch.points.size();
        record_point(z_ev, SpecialPointKind::BranchPoint);
        /* attach the second tangent to the just-recorded branch point */
        if (!t2.empty() && branch.points.size() > idx_before)
          branch.points[idx_before].second_tangent = t2;
      }
      prev_bp = bp_new;
    }

    record_point(z_new, SpecialPointKind::None);
    ++produced;

    /* Advance: new tangent from the corrected point, then adapt step.
     * Grow the step when the corrector converged easily, otherwise
     * keep it. (We approximate "easily" by always nudging up; the
     * retry path above handles hard cases by shrinking.) */
    std::vector<double> new_tan;
    if (!compute_tangent(z_new, tangent, &new_tan)) {
      branch.message = "tangent went singular (possible bifurcation); stopping";
      break;
    }
    z = z_new;
    tangent = new_tan;
    h = std::min(settings.h_max, h * 1.2);
  }

  if (branch.message.empty())
    branch.message = "completed (" + std::to_string(produced) + " points)";
  branch.ok = branch.points.size() > 1;
  return branch;
}

/* ---- branch switching --------------------------------------------------- */

Branch switch_branch(const Model &m, const BranchPoint &bp,
                     const ContinuationSettings &settings) {
  Branch out;
  const std::size_t n = m.n;
  if (bp.x.size() != n) { out.message = "invalid branch point"; return out; }
  const double base_h = (settings.h0 > 0 ? settings.h0 : 0.05);

  /* Candidate kick directions in (x,p) space: the precomputed second_tangent
   * first (if any), then each state-coordinate axis and the parameter axis.
   * For each, step off the branch point, re-converge to an equilibrium at the
   * new parameter, and accept the first that lands on a DISTINCT branch and
   * continues. Trying coordinate kicks (not only the null tangent) makes this
   * robust even at degenerate branch points where the null space is hard to
   * pin down numerically (e.g. a pitchfork, where both f_x and f_p vanish). */
  std::vector<std::vector<double>> dirs;
  if (bp.second_tangent.size() == n + 1) dirs.push_back(bp.second_tangent);
  for (std::size_t i = 0; i < n; ++i) {
    std::vector<double> d(n + 1, 0.0); d[i] = 1.0; dirs.push_back(d);
  }
  { std::vector<double> d(n + 1, 0.0); d[n] = 1.0; dirs.push_back(d); }

  for (const auto &dir : dirs) {
    for (double sgn : {+1.0, -1.0}) {
      for (double scale : {4.0, 8.0, 16.0, 2.0}) {
        const double h = sgn * scale * base_h;
        std::vector<double> x_new(n);
        for (std::size_t i = 0; i < n; ++i) x_new[i] = bp.x[i] + h * dir[i];
        double p_new = bp.p + h * dir[n];
        /* if this direction doesn't move p, nudge it so we leave the singular
         * point instead of re-finding the same equilibrium. */
        if (std::fabs(dir[n]) < 1e-12) p_new = bp.p + sgn * scale * base_h * 0.5;
        std::string err;
        if (!correct_equilibrium(m, &x_new, p_new, settings.max_corrector_iters,
                                 settings.corrector_tol, &err))
          continue;
        double d2 = 0.0;
        for (std::size_t i = 0; i < n; ++i) { const double dd = x_new[i] - bp.x[i]; d2 += dd * dd; }
        if (std::sqrt(d2) < 1e-5) continue; /* landed back on the BP / same branch */
        Branch b = continue_equilibrium(m, x_new, p_new, settings);
        if (b.ok && b.points.size() > 3) {
          b.message = "switched branch at p=" + std::to_string(bp.p) + "; " + b.message;
          return b;
        }
      }
    }
  }
  out.message = "branch switch failed to converge onto a distinct branch";
  return out;
}

/* ---- two-parameter continuation of fold/Hopf curves --------------------- */
namespace {

/* the codim-1 defining function g(x,p,q): det(f_x) for a fold; for Hopf, the
 * minimum |Re| over complex pairs signed by Re (zero when a pair is on the
 * imaginary axis). Built on finite-difference Jacobians of the 2-param field. */
double g_fold2(const Model2 &m, const std::vector<double> &x, double p, double q) {
  const std::size_t n = m.n;
  std::vector<double> J(n * n);
  std::vector<double> f0(n), fp(n);
  std::string err;
  if (!m.vector_field(x.data(), p, q, f0.data(), &err)) return std::nan("");
  const double h = 1e-7;
  std::vector<double> xt = x;
  for (std::size_t j = 0; j < n; ++j) {
    const double save = xt[j]; const double dh = h * (std::fabs(save) + 1.0);
    xt[j] = save + dh;
    if (!m.vector_field(xt.data(), p, q, fp.data(), &err)) return std::nan("");
    for (std::size_t i = 0; i < n; ++i) J[i * n + j] = (fp[i] - f0[i]) / dh;
    xt[j] = save;
  }
  return determinant(J, n);
}

double g_hopf2(const Model2 &m, const std::vector<double> &x, double p, double q) {
  const std::size_t n = m.n;
  std::vector<double> J(n * n), f0(n), fp(n);
  std::string err;
  if (!m.vector_field(x.data(), p, q, f0.data(), &err)) return std::nan("");
  const double h = 1e-7;
  std::vector<double> xt = x;
  for (std::size_t j = 0; j < n; ++j) {
    const double save = xt[j]; const double dh = h * (std::fabs(save) + 1.0);
    xt[j] = save + dh;
    if (!m.vector_field(xt.data(), p, q, fp.data(), &err)) return std::nan("");
    for (std::size_t i = 0; i < n; ++i) J[i * n + j] = (fp[i] - f0[i]) / dh;
    xt[j] = save;
  }
  std::vector<std::complex<double>> ev;
  if (!eigenvalues(J, n, &ev)) return std::nan("");
  bool found = false; double best = 0, sgn = 0;
  for (const auto &z : ev)
    if (std::fabs(z.imag()) > 1e-7) {
      const double a = std::fabs(z.real());
      if (!found || a < best) { best = a; sgn = z.real(); found = true; }
    }
  return found ? sgn : std::nan("");
}

/* the extended residual G(u) where u = [x(0..n-1), p, q], length n+2.
 * G has n+1 components: f (n of them) and the codim-1 condition g. */
void G_residual(const Model2 &m, TwoParamKind kind, const std::vector<double> &u,
                std::vector<double> *out) {
  const std::size_t n = m.n;
  std::vector<double> x(u.begin(), u.begin() + n);
  const double p = u[n], q = u[n + 1];
  out->assign(n + 1, 0.0);
  std::string err;
  std::vector<double> f(n);
  if (m.vector_field(x.data(), p, q, f.data(), &err))
    for (std::size_t i = 0; i < n; ++i) (*out)[i] = f[i];
  (*out)[n] = (kind == TwoParamKind::Fold) ? g_fold2(m, x, p, q) : g_hopf2(m, x, p, q);
}

/* finite-difference Jacobian of G wrt u: (n+1) x (n+2). */
void G_jacobian(const Model2 &m, TwoParamKind kind, const std::vector<double> &u,
                std::vector<double> *Jout) {
  const std::size_t n = m.n, R = n + 1, C = n + 2;
  Jout->assign(R * C, 0.0);
  std::vector<double> g0, gp;
  G_residual(m, kind, u, &g0);
  std::vector<double> ut = u;
  for (std::size_t j = 0; j < C; ++j) {
    const double save = ut[j]; const double dh = 1e-6 * (std::fabs(save) + 1.0);
    ut[j] = save + dh;
    G_residual(m, kind, ut, &gp);
    for (std::size_t i = 0; i < R; ++i) (*Jout)[i * C + j] = (gp[i] - g0[i]) / dh;
    ut[j] = save;
  }
}

}  // namespace

TwoParamCurve two_param_curve(const Model2 &m, TwoParamKind kind,
                              const std::vector<double> &x0, double p0, double q0,
                              const TwoParamSettings &settings) {
  TwoParamCurve curve;
  curve.kind = (kind == TwoParamKind::Fold) ? SpecialPointKind::Fold : SpecialPointKind::Hopf;
  const std::size_t n = m.n;
  if (x0.size() != n || n == 0) { curve.message = "bad starting point"; return curve; }
  const std::size_t C = n + 2;        /* unknowns: x, p, q */
  const std::size_t Rr = n + 1;       /* equations */

  /* Newton corrector onto G=0 with an arclength constraint (project the step
   * orthogonal to the current tangent). We use a least-squares solve of the
   * (n+1) x (n+2) system by appending the tangent row to make it square. */
  auto correct = [&](std::vector<double> *u, const std::vector<double> &tan) -> bool {
    for (int it = 0; it < settings.max_corrector_iters; ++it) {
      std::vector<double> G; G_residual(m, kind, *u, &G);
      double res = 0; for (double v : G) res += v * v; res = std::sqrt(res);
      if (res < settings.corrector_tol) return true;
      std::vector<double> J; G_jacobian(m, kind, *u, &J);
      /* square (n+2)x(n+2): G-rows + tangent row (=0 in residual) */
      std::vector<double> A(C * C, 0.0), b(C, 0.0);
      for (std::size_t i = 0; i < Rr; ++i) {
        for (std::size_t j = 0; j < C; ++j) A[i * C + j] = J[i * C + j];
        b[i] = -G[i];
      }
      for (std::size_t j = 0; j < C; ++j) A[Rr * C + j] = tan[j];
      b[Rr] = 0.0;
      std::vector<double> dz;
      if (!solve_linear(A, b, &dz)) return false;
      for (std::size_t j = 0; j < C; ++j) (*u)[j] += dz[j];
    }
    std::vector<double> G; G_residual(m, kind, *u, &G);
    double res = 0; for (double v : G) res += v * v;
    return std::sqrt(res) < 1e-6;
  };

  /* tangent: null vector of the (n+1)x(n+2) Jacobian (the curve direction). */
  auto tangent_of = [&](const std::vector<double> &u, const std::vector<double> &prev,
                        std::vector<double> *t) -> bool {
    std::vector<double> J; G_jacobian(m, kind, u, &J);
    std::vector<double> A(C * C, 0.0), b(C, 0.0);
    for (std::size_t i = 0; i < Rr; ++i)
      for (std::size_t j = 0; j < C; ++j) A[i * C + j] = J[i * C + j];
    if (prev.empty()) { A[Rr * C + (n + 1)] = 1.0; b[Rr] = 1.0; }  /* seed: move in q */
    else { for (std::size_t j = 0; j < C; ++j) A[Rr * C + j] = prev[j]; b[Rr] = 1.0; }
    if (!solve_linear(A, b, t)) return false;
    double nrm = 0; for (double v : *t) nrm += v * v; nrm = std::sqrt(nrm);
    if (nrm < 1e-300) return false;
    for (double &v : *t) v /= nrm;
    if (!prev.empty()) { double d = 0; for (std::size_t j = 0; j < C; ++j) d += (*t)[j]*prev[j];
      if (d < 0) for (double &v : *t) v = -v; }
    return true;
  };

  std::vector<double> u(C);
  for (std::size_t i = 0; i < n; ++i) u[i] = x0[i];
  u[n] = p0; u[n + 1] = q0;
  std::vector<double> empty_prev;
  if (!correct(&u, std::vector<double>(C, 0.0))) {
    /* still record the start even if it wasn't a perfect codim-1 point */
  }

  auto push = [&](const std::vector<double> &uu) {
    TwoParamPoint pt; pt.p = uu[n]; pt.q = uu[n + 1];
    pt.x.assign(uu.begin(), uu.begin() + n);
    curve.points.push_back(std::move(pt));
  };

  /* trace both directions */
  for (int dir = 0; dir < 2; ++dir) {
    std::vector<double> uu = u, tan;
    if (!tangent_of(uu, empty_prev, &tan)) continue;
    if (dir == 1) for (double &v : tan) v = -v;
    int steps = settings.max_points / 2;
    for (int s = 0; s < steps; ++s) {
      std::vector<double> un = uu;
      for (std::size_t j = 0; j < C; ++j) un[j] += settings.h0 * tan[j];
      std::vector<double> ntan;
      if (!tangent_of(uu, tan, &ntan)) break;  /* tangent at current point */
      if (!correct(&un, ntan)) break;
      const double pp = un[n], qq = un[n + 1];
      if (pp < settings.p_min || pp > settings.p_max ||
          qq < settings.q_min || qq > settings.q_max) { push(un); break; }
      push(un);
      /* advance */
      std::vector<double> t2;
      if (!tangent_of(un, ntan, &t2)) break;
      uu = un; tan = t2;
    }
  }

  /* ---- detect codim-2 points ALONG the curve ----------------------------
   * Two scalar test functions are evaluated at each curve point; a sign change
   * between consecutive points brackets a codim-2 point. For a FOLD curve:
   *   - t_cusp = fold normal-form coefficient a  (a -> 0 is a CUSP)
   *   - t_bt   = 2nd-smallest |eigenvalue|        (-> 0 is Bogdanov-Takens)
   * For a HOPF curve:
   *   - t_gh   = first Lyapunov coefficient l1    (-> 0 is generalized Hopf)
   *   - t_bt   = Hopf frequency omega             (-> 0 is Bogdanov-Takens)
   * The single-parameter normal-form routines are evaluated on a temporary
   * Model that pins q to the point's value and varies p. */
  if (curve.points.size() > 2) {
    auto model_at_q = [&](double qfix) {
      Model mm; mm.n = n;
      mm.vector_field = [&m, qfix](const double *xx, double pp, double *fo, std::string *er) {
        return m.vector_field(xx, pp, qfix, fo, er);
      };
      return mm;
    };
    auto test_funcs = [&](const TwoParamPoint &pt, double *t_cusp_or_gh, double *t_bt) -> bool {
      Model mm = model_at_q(pt.q);
      if (kind == TwoParamKind::Fold) {
        double a = 0, l0 = 0; std::string e;
        if (!fold_normal_form(mm, pt.x, pt.p, &a, &l0, &e)) return false;
        *t_cusp_or_gh = a;
        /* BT indicator on a fold curve: at a fold one eigenvalue is ~0. A BT is
         * where a SECOND eigenvalue also hits zero. The sum of the non-zero
         * eigenvalues = trace(J) - (the near-zero eigenvalue) is a SIGNED test
         * that crosses zero at the BT (unlike |eigenvalue|, which only touches
         * zero and so can't be bracketed by a sign change). For a 2-D system
         * this is just the non-trivial real eigenvalue. */
        std::vector<double> J; if (!finite_diff_jacobian(mm, pt.x.data(), pt.p, &J, &e, 1e-7)) return false;
        std::vector<std::complex<double>> ev; if (!eigenvalues(J, n, &ev)) return false;
        /* find the eigenvalue nearest zero (the fold's), sum the real parts of
         * the rest -> signed BT test */
        std::size_t iz = 0; double best = 1e300;
        for (std::size_t k = 0; k < ev.size(); ++k) { double mg = std::abs(ev[k]); if (mg < best) { best = mg; iz = k; } }
        double sum_rest = 0.0;
        for (std::size_t k = 0; k < ev.size(); ++k) if (k != iz) sum_rest += ev[k].real();
        *t_bt = sum_rest;
        return true;
      } else {
        double l1 = 0, w = 0; std::string e;
        if (!hopf_first_lyapunov(mm, pt.x, pt.p, &l1, &w, &e)) return false;
        *t_cusp_or_gh = l1;
        *t_bt = w;     /* omega -> 0 is Bogdanov-Takens */
        return true;
      }
    };
    /* Extra codim-2 tests on a HOPF curve for 3+ D systems:
     *  ZH (zero-Hopf / fold-Hopf): a REAL eigenvalue passes through 0 while the
     *    pure-imaginary pair persists  -> signed test = that real eigenvalue.
     *  HH (Hopf-Hopf / double-Hopf): a SECOND pair reaches the imaginary axis
     *    -> signed test = the real part of the nearest off-axis complex pair.
     * Both return NaN when not applicable (Fold curve, or n<3, or no clear
     * Hopf pair), so they simply never trigger there. */
    auto zhhh_tests = [&](const TwoParamPoint &pt, double *t_zh, double *t_hh) -> bool {
      *t_zh = std::nan(""); *t_hh = std::nan("");
      if (kind != TwoParamKind::Hopf || n < 3) return false;
      Model mm = model_at_q(pt.q);
      std::vector<double> J; std::string e;
      if (!finite_diff_jacobian(mm, pt.x.data(), pt.p, &J, &e, 1e-7)) return false;
      std::vector<std::complex<double>> ev; if (!eigenvalues(J, n, &ev)) return false;
      /* identify the active Hopf pair: the complex pair with the smallest
       * |Re|. Exclude it, then look at what remains. */
      int hopf_i = -1; double hopf_re = 1e300;
      for (std::size_t k = 0; k < ev.size(); ++k)
        if (std::fabs(ev[k].imag()) > 1e-6 && std::fabs(ev[k].real()) < hopf_re) { hopf_re = std::fabs(ev[k].real()); hopf_i = (int)k; }
      if (hopf_i < 0) return false;
      const double hw = std::fabs(ev[hopf_i].imag());
      std::vector<bool> used(ev.size(), false);
      used[hopf_i] = true;
      for (std::size_t k = 0; k < ev.size(); ++k)   /* mark the conjugate */
        if (!used[k] && std::fabs(ev[k].imag() + ev[hopf_i].imag()) < 1e-6 && std::fabs(ev[k].real() - ev[hopf_i].real()) < 1e-6) { used[k] = true; break; }
      /* ZH: smallest-magnitude REAL eigenvalue among the rest (signed) */
      double zh = std::nan(""); double zh_mag = 1e300;
      for (std::size_t k = 0; k < ev.size(); ++k)
        if (!used[k] && std::fabs(ev[k].imag()) < 1e-6 && std::fabs(ev[k].real()) < zh_mag) { zh_mag = std::fabs(ev[k].real()); zh = ev[k].real(); }
      /* HH: real part of the nearest OTHER complex pair, but only if its
       * frequency differs from the active pair (a genuinely different pair) */
      double hh = std::nan(""); double hh_mag = 1e300;
      for (std::size_t k = 0; k < ev.size(); ++k)
        if (!used[k] && std::fabs(ev[k].imag()) > 1e-6 && std::fabs(std::fabs(ev[k].imag()) - hw) > 1e-3) {
          if (std::fabs(ev[k].real()) < hh_mag) { hh_mag = std::fabs(ev[k].real()); hh = ev[k].real(); }
        }
      *t_zh = zh; *t_hh = hh;
      return true;
    };
    std::vector<double> A(curve.points.size(), std::nan("")), B(curve.points.size(), std::nan(""));
    std::vector<double> Czh(curve.points.size(), std::nan("")), Dhh(curve.points.size(), std::nan(""));
    for (std::size_t i = 0; i < curve.points.size(); ++i) {
      double a = 0, b = 0;
      if (test_funcs(curve.points[i], &a, &b)) { A[i] = a; B[i] = b; }
      double cz = 0, dh = 0;
      if (zhhh_tests(curve.points[i], &cz, &dh)) { Czh[i] = cz; Dhh[i] = dh; }
    }
    for (std::size_t i = 1; i < curve.points.size(); ++i) {
      /* skip the join between the two trace directions (large parameter jump) */
      const double jump = std::fabs(curve.points[i].p - curve.points[i-1].p) +
                          std::fabs(curve.points[i].q - curve.points[i-1].q);
      if (jump > 1.0) continue;
      SpecialPointKind k = SpecialPointKind::None;
      bool use_bt_test = false;
      if (std::isfinite(B[i-1]) && std::isfinite(B[i]) && B[i-1] * B[i] < 0.0) {
        k = SpecialPointKind::BogdanovTakens; use_bt_test = true;   /* BT takes precedence */
      } else if (std::isfinite(A[i-1]) && std::isfinite(A[i]) && A[i-1] * A[i] < 0.0) {
        /* A sign change of the cusp/GH test (fold coeff a, or l1) only counts if
         * it crosses ZERO, not infinity: the fold coefficient diverges (and
         * flips sign) where the cycle/center-manifold reduction is singular,
         * which is NOT a cusp. Require the magnitude on BOTH sides be modest
         * (a genuine zero crossing) -- a pole shows |value| blowing up. */
        const double m0 = std::fabs(A[i-1]), m1 = std::fabs(A[i]);
        const double near = std::min(m0, m1), far = std::max(m0, m1);
        const bool looks_like_zero = (near < 0.5) && (far < 50.0);
        if (looks_like_zero)
          k = (kind == TwoParamKind::Fold) ? SpecialPointKind::Cusp
                                           : SpecialPointKind::GeneralizedHopf;
      } else if (std::isfinite(Czh[i-1]) && std::isfinite(Czh[i]) && Czh[i-1] * Czh[i] < 0.0) {
        /* a real eigenvalue crossed zero with the Hopf pair present -> ZH */
        const double near = std::min(std::fabs(Czh[i-1]), std::fabs(Czh[i]));
        if (near < 0.5) k = SpecialPointKind::ZeroHopf;
      } else if (std::isfinite(Dhh[i-1]) && std::isfinite(Dhh[i]) && Dhh[i-1] * Dhh[i] < 0.0) {
        /* a second complex pair reached the imaginary axis -> HH */
        const double near = std::min(std::fabs(Dhh[i-1]), std::fabs(Dhh[i]));
        if (near < 0.5) k = SpecialPointKind::HopfHopf;
      } else if (std::isfinite(Dhh[i]) && std::fabs(Dhh[i]) < 0.06 &&
                 (i + 1 >= curve.points.size() || !std::isfinite(Dhh[i+1])) &&
                 std::isfinite(Dhh[i-1]) && std::fabs(Dhh[i-1]) > std::fabs(Dhh[i])) {
        /* HH where the Hopf curve TERMINATES: the second pair's real part is
         * decreasing toward zero and the curve ends here (the two neutral pairs
         * make the corrector stall exactly at the double-Hopf). Tag it. */
        k = SpecialPointKind::HopfHopf;
      } else if (std::isfinite(Czh[i]) && std::fabs(Czh[i]) < 0.06 &&
                 (i + 1 >= curve.points.size() || !std::isfinite(Czh[i+1])) &&
                 std::isfinite(Czh[i-1]) && std::fabs(Czh[i-1]) > std::fabs(Czh[i])) {
        /* ZH where the curve terminates (real eigenvalue grazing zero). */
        k = SpecialPointKind::ZeroHopf;
      }
      if (k != SpecialPointKind::None) {
        /* BISECTION refine the codim-2 location along the curve segment
         * [i-1, i]. Parametrize by t in [0,1] linearly interpolating (p,q,x);
         * at each t re-solve the equilibrium in p at the interpolated q, then
         * evaluate the relevant test function (BT: t_bt, else t_cusp_or_gh) and
         * bisect on its sign. This pins the point far more precisely than just
         * tagging the nearer curve node. */
        const TwoParamPoint &P0 = curve.points[i-1], &P1 = curve.points[i];
        auto eval_t = [&](double t, TwoParamPoint *refp, double *tf) -> bool {
          const double q = P0.q + t * (P1.q - P0.q);
          double p = P0.p + t * (P1.p - P0.p);
          std::vector<double> x(n);
          for (std::size_t kk = 0; kk < n; ++kk) x[kk] = P0.x[kk] + t * (P1.x[kk] - P0.x[kk]);
          Model mm = model_at_q(q);
          std::string e;
          if (!correct_equilibrium(mm, &x, p, settings.max_corrector_iters, settings.corrector_tol, &e))
            return false;
          double a = 0, b = 0;
          /* recompute both test functions at the corrected point */
          TwoParamPoint tp; tp.p = p; tp.q = q; tp.x = x;
          if (!test_funcs(tp, &a, &b)) return false;
          *tf = use_bt_test ? b : a;
          if (refp) { *refp = tp; }
          return true;
        };
        double lo = 0.0, hi = 1.0, flo = use_bt_test ? B[i-1] : A[i-1];
        TwoParamPoint best; bool havebest = false;
        for (int it = 0; it < 40; ++it) {
          const double mid = 0.5 * (lo + hi);
          double fm; TwoParamPoint rp;
          if (!eval_t(mid, &rp, &fm)) break;
          havebest = true; best = rp;
          if (flo * fm <= 0.0) { hi = mid; } else { lo = mid; flo = fm; }
          if (hi - lo < 1e-10) break;
        }
        /* tag the nearer node (for the diagram), and attach the refined point */
        std::size_t at = (std::fabs(use_bt_test ? B[i-1] : A[i-1]) <
                          std::fabs(use_bt_test ? B[i] : A[i])) ? i-1 : i;
        if (curve.points[at].special == SpecialPointKind::None) {
          curve.points[at].special = k;
          if (havebest) {
            curve.points[at].p2 = best.p; curve.points[at].q2 = best.q;
            curve.points[at].x2 = best.x;
            /* codim-2 normal-form coefficients at the refined point */
            if (k == SpecialPointKind::BogdanovTakens) {
              Model mm = model_at_q(best.q);
              double a = 0, b = 0; std::string e;
              if (bt_normal_form(mm, best.x, best.p, &a, &b, &e)) {
                curve.points[at].bt_a = a; curve.points[at].bt_b = b;
                curve.points[at].has_codim2_nf = true;
              }
            } else if (k == SpecialPointKind::Cusp) {
              Model mm = model_at_q(best.q);
              double cc = 0; std::string e;
              if (cusp_normal_form(mm, best.x, best.p, &cc, &e)) {
                curve.points[at].cusp_c = cc;
                curve.points[at].has_codim2_nf = true;
              }
            } else if (k == SpecialPointKind::GeneralizedHopf) {
              Model mm = model_at_q(best.q);
              double ll2 = 0; std::string e;
              if (gh_second_lyapunov(mm, best.x, best.p, &ll2, &e)) {
                curve.points[at].gh_l2 = ll2;
                curve.points[at].has_codim2_nf = true;
              }
            }
          } else {
            curve.points[at].p2 = curve.points[at].p;
            curve.points[at].q2 = curve.points[at].q;
            curve.points[at].x2 = curve.points[at].x;
          }
          curve.special_indices.push_back(at);
        }
      }
    }
  }

  curve.ok = curve.points.size() > 2;
  curve.message = curve.ok ? ("traced " + std::to_string(curve.points.size()) + " points" +
                              (curve.special_indices.empty() ? "" :
                               ", " + std::to_string(curve.special_indices.size()) + " codim-2 point(s)"))
                           : "could not trace a curve from this start";
  return curve;
}

/* ---- periodic-orbit continuation by collocation ------------------------- */
namespace {

/* The unknown vector U packs the m mesh points (each length n) followed by the
 * period T:  U = [X0_0..X0_{n-1}, X1_..., ..., X{m-1}_..., T], length m*n+1.
 * Residual F(U) (same length): trapezoidal collocation on each interval
 *   X[i+1] - X[i] - (T/m) * 0.5*(f(X[i]) + f(X[i+1])) = 0   (i = 0..m-1, wrap)
 * which is m*n equations, plus one phase condition fixing the orbit's phase:
 *   <X[0] - Xprev[0], f(Xprev[0])> = 0  (orthogonality to the reference field;
 *   for the first solve we instead pin the first coordinate of X[0] to its
 *   guess, which is robust). */
struct CycleCtx {
  const Model *m;
  double p;
  std::size_t n, mesh;
  std::vector<double> phase_ref;   /* reference point X0 for the phase cond. */
  std::vector<double> phase_dir;   /* f(phase_ref): direction for orthogonality */
  bool pin_mode = true;            /* true: pin X0[0]; false: integral phase   */
  double pin_val = 0.0;            /* value to pin X0[0] to (pin_mode)         */
  /* Adaptive mesh: fraction of the period spanned by each interval i (length
   * `mesh`, summing to 1). Empty => uniform 1/mesh (the original behaviour, so
   * existing callers/tests are unchanged). Remeshing redistributes these to put
   * more points where the orbit moves fast / curves sharply (stiff relaxation
   * cycles), which a uniform-in-time mesh resolves poorly. */
  std::vector<double> frac;
};

bool cyc_field(const CycleCtx &c, const double *x, std::vector<double> *f) {
  f->assign(c.n, 0.0);
  std::string err;
  return c.m->vector_field(x, c.p, f->data(), &err);
}

/* residual of the collocation system at U */
bool cyc_residual(const CycleCtx &c, const std::vector<double> &U, std::vector<double> *Fout) {
  const std::size_t n = c.n, M = c.mesh;
  const double T = U[M * n];
  Fout->assign(M * n + 1, 0.0);
  std::vector<double> fi(n), fj(n), xi(n), xj(n);
  for (std::size_t i = 0; i < M; ++i) {
    const std::size_t j = (i + 1) % M;
    for (std::size_t k = 0; k < n; ++k) { xi[k] = U[i*n+k]; xj[k] = U[j*n+k]; }
    if (!cyc_field(c, xi.data(), &fi)) return false;
    if (!cyc_field(c, xj.data(), &fj)) return false;
    /* interval duration: adaptive fraction if provided, else uniform T/M */
    const double dti = (c.frac.size() == M) ? (T * c.frac[i]) : (T / (double)M);
    for (std::size_t k = 0; k < n; ++k)
      (*Fout)[i*n+k] = xj[k] - xi[k] - dti * 0.5 * (fi[k] + fj[k]);
  }
  /* phase condition */
  if (c.pin_mode || c.phase_ref.size() < n || c.phase_dir.size() < n) {
    /* pin the first coordinate (also the safe fallback when no integral phase
     * reference has been set yet) */
    (*Fout)[M*n] = U[0] - (c.pin_mode ? c.pin_val : U[0]);
  } else {
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k) s += (U[k] - c.phase_ref[k]) * c.phase_dir[k];
    (*Fout)[M*n] = s;
  }
  return true;
}

/* Adaptive remeshing. Given a converged orbit U (M points + period) and a
 * CycleCtx, redistribute the mesh points so they equidistribute a monitor
 * function (here arclength + a curvature weight), and write the resulting
 * per-interval TIME fractions into c.frac. Returns a new orbit Unew sampled at
 * the redistributed points (same M, same period). This concentrates points on
 * the fast/sharp parts of relaxation cycles, which a uniform-in-time mesh
 * resolves poorly. No-op-safe: on any failure it leaves things unchanged. */
bool cyc_remesh(CycleCtx &c, const std::vector<double> &U, std::vector<double> *Unew) {
  const std::size_t n = c.n, M = c.mesh;
  if (M < 4) return false;
  const double T = U[M * n];
  if (!(T > 0)) return false;

  /* current interval time fractions (uniform or existing adaptive) */
  std::vector<double> fr(M);
  for (std::size_t i = 0; i < M; ++i) fr[i] = (c.frac.size() == M) ? c.frac[i] : 1.0 / (double)M;

  /* monitor density per interval: arclength of the segment plus a small floor,
   * raised to a power<1 so adaptation is gentle (avoids starving slow arcs). */
  std::vector<double> w(M, 0.0);
  double wsum = 0.0, len_total = 0.0;
  for (std::size_t i = 0; i < M; ++i) {
    const std::size_t j = (i + 1) % M;
    double seg = 0.0;
    for (std::size_t k = 0; k < n; ++k) { const double d = U[j*n+k] - U[i*n+k]; seg += d * d; }
    seg = std::sqrt(seg);
    len_total += seg;
    w[i] = seg;
  }
  if (!(len_total > 0)) return false;
  /* monitor = (arclength density)^alpha, smoothed, with a floor */
  const double alpha = 0.5;
  for (std::size_t i = 0; i < M; ++i) {
    double dens = w[i] / (len_total / (double)M); /* relative speed */
    dens = std::pow(dens + 0.15, alpha);
    w[i] = dens; wsum += dens;
  }
  if (!(wsum > 0)) return false;

  /* cumulative monitor at node boundaries (length M+1, periodic) */
  std::vector<double> cum(M + 1, 0.0);
  for (std::size_t i = 0; i < M; ++i) cum[i + 1] = cum[i] + w[i] / wsum; /* 0..1 */

  /* also need cumulative TIME at each original node for resampling the orbit */
  std::vector<double> tcum(M + 1, 0.0);
  for (std::size_t i = 0; i < M; ++i) tcum[i + 1] = tcum[i] + fr[i];

  /* new nodes equidistributed in the monitor variable xi in [0,1). For each new
   * node find where cum crosses k/M, interpolate the orbit (linear) and the
   * time there. */
  Unew->assign(M * n + 1, 0.0);
  std::vector<double> new_tnode(M + 1, 0.0);
  for (std::size_t inew = 0; inew <= M; ++inew) {
    const double target = (double)inew / (double)M; /* in [0,1] */
    /* locate interval [cum[s], cum[s+1]] containing target */
    std::size_t s = 0;
    while (s < M && cum[s + 1] < target) ++s;
    if (s >= M) s = M - 1;
    const double seg = cum[s + 1] - cum[s];
    const double local = seg > 1e-15 ? (target - cum[s]) / seg : 0.0;
    if (inew < M) {
      const std::size_t a = s % M, b = (s + 1) % M;
      for (std::size_t k = 0; k < n; ++k)
        (*Unew)[inew * n + k] = (1.0 - local) * U[a * n + k] + local * U[b * n + k];
    }
    new_tnode[inew] = tcum[s] + local * (tcum[s + 1] - tcum[s]); /* fractional time in [0,1] */
  }
  (*Unew)[M * n] = T;

  /* new per-interval time fractions from the new node times */
  std::vector<double> newfrac(M, 0.0);
  double fsum = 0.0;
  for (std::size_t i = 0; i < M; ++i) {
    double d = new_tnode[i + 1] - new_tnode[i];
    if (d <= 1e-9) d = 1e-9;
    newfrac[i] = d; fsum += d;
  }
  if (!(fsum > 0)) return false;
  for (std::size_t i = 0; i < M; ++i) newfrac[i] /= fsum;
  c.frac = newfrac;
  return true;
}

/* Newton solve of the collocation system; returns refined U. Uses a finite-
 * difference Jacobian (dense (M*n+1)^2) — fine for the modest mesh sizes here.*/
bool cyc_newton(const CycleCtx &c, std::vector<double> *U, int iters, double tol) {
  const std::size_t N = c.mesh * c.n + 1;
  std::vector<double> F0, Fp;
  for (int it = 0; it < iters; ++it) {
    if (!cyc_residual(c, *U, &F0)) return false;
    double res = 0; for (double v : F0) res += v*v; res = std::sqrt(res);
    if (res < tol) return true;
    std::vector<double> J(N * N, 0.0);
    std::vector<double> Ut = *U;
    for (std::size_t j = 0; j < N; ++j) {
      const double save = Ut[j]; const double dh = 1e-7 * (std::fabs(save) + 1.0);
      Ut[j] = save + dh;
      if (!cyc_residual(c, Ut, &Fp)) return false;
      for (std::size_t i = 0; i < N; ++i) J[i*N+j] = (Fp[i] - F0[i]) / dh;
      Ut[j] = save;
    }
    std::vector<double> b(N); for (std::size_t i = 0; i < N; ++i) b[i] = -F0[i];
    std::vector<double> du;
    if (!solve_linear(J, b, &du)) return false;
    /* damped step for robustness */
    double nrm = 0; for (double v : du) nrm += v*v; nrm = std::sqrt(nrm);
    double damp = (nrm > 1.0) ? 1.0 / nrm : 1.0;
    for (std::size_t j = 0; j < N; ++j) (*U)[j] += damp * du[j];
  }
  if (!cyc_residual(c, *U, &F0)) return false;
  double res = 0; for (double v : F0) res += v*v;
  return std::sqrt(res) < tol * 100;
}

/* amplitude (peak-to-peak) and min/max of the first coordinate over the mesh */
void cyc_amp(const std::vector<double> &U, std::size_t n, std::size_t M,
             double *amp, double *mn, double *mx) {
  double lo = 1e300, hi = -1e300;
  for (std::size_t i = 0; i < M; ++i) { const double v = U[i*n]; lo = std::min(lo, v); hi = std::max(hi, v); }
  *mn = lo; *mx = hi; *amp = hi - lo;
}

/* Fold-of-cycles (LPC) test function: the determinant of the collocation
 * Jacobian at a converged cycle U. Away from a cycle fold the periodic BVP is
 * regular and this is bounded away from zero; at an LPC (a nontrivial Floquet
 * multiplier crossing +1, where the cycle branch turns) the system becomes
 * singular and the determinant changes sign. We rescale rows to keep the
 * determinant in a sane numeric range (only its SIGN is used downstream). */
double cycle_fold_test(const CycleCtx &c, const std::vector<double> &U) {
  const std::size_t N = c.mesh * c.n + 1;
  std::vector<double> F0, Fp;
  if (!cyc_residual(c, U, &F0)) return std::nan("");
  std::vector<double> J(N * N, 0.0), Ut = U;
  for (std::size_t j = 0; j < N; ++j) {
    const double save = Ut[j]; const double dh = 1e-7 * (std::fabs(save) + 1.0);
    Ut[j] = save + dh;
    if (!cyc_residual(c, Ut, &Fp)) return std::nan("");
    for (std::size_t i = 0; i < N; ++i) J[i*N+j] = (Fp[i] - F0[i]) / dh;
    Ut[j] = save;
  }
  /* row-normalize so det doesn't under/overflow for large meshes */
  for (std::size_t i = 0; i < N; ++i) {
    double s = 0; for (std::size_t j = 0; j < N; ++j) s += J[i*N+j]*J[i*N+j];
    s = std::sqrt(s); if (s < 1e-300) s = 1.0;
    for (std::size_t j = 0; j < N; ++j) J[i*N+j] /= s;
  }
  return determinant(J, N);
}

/* ---- Floquet multipliers of a converged cycle -------------------------------
 * Integrate the variational equation  Phi' = T * Df(x(s)) * Phi,  Phi(0)=I,
 * over one period (s: 0->1) using the orbit samples in U and the model Jacobian.
 * The monodromy matrix Phi(1) has eigenvalues = the Floquet multipliers. One is
 * always ~1 (along the flow); the others govern stability. RK4 on the matrix
 * ODE between mesh points (interpolating x along the orbit). */
void cycle_floquet(const CycleCtx &c, const std::vector<double> &U,
                   std::vector<Complex> *mult_out) {
  const std::size_t n = c.n, M = c.mesh;
  const double T = U[M * n];
  mult_out->clear();
  if (n == 0 || M < 2 || !(T > 0)) return;

  auto node = [&](std::size_t i, std::vector<double> &x) {
    i %= M; for (std::size_t k = 0; k < n; ++k) x[k] = U[i*n+k];
  };
  /* cumulative node times in [0,1]: uniform (i/M) unless an adaptive mesh frac
   * is present, in which case node i sits at the cumulative sum of fractions. */
  std::vector<double> tnode(M + 1, 0.0);
  for (std::size_t i = 0; i < M; ++i)
    tnode[i + 1] = tnode[i] + ((c.frac.size() == M) ? c.frac[i] : 1.0 / (double)M);
  std::vector<double> xa(n), xb(n);
  /* interpolate the orbit at normalized period-time s in [0,1) using tnode */
  auto orbit_at_s = [&](double s, std::vector<double> &x) {
    s -= std::floor(s);
    std::size_t i = 0;
    while (i < M && tnode[i + 1] < s) ++i;
    if (i >= M) i = M - 1;
    const double seg = tnode[i + 1] - tnode[i];
    const double t = seg > 1e-15 ? (s - tnode[i]) / seg : 0.0;
    node(i, xa); node(i + 1, xb);
    for (std::size_t k = 0; k < n; ++k) x[k] = (1.0 - t) * xa[k] + t * xb[k];
  };
  std::vector<double> Jf(n * n);
  auto jac_at = [&](const std::vector<double> &x) -> bool {
    std::string err;
    if (c.m->jacobian_x && c.m->jacobian_x(x.data(), c.p, Jf.data(), &err)) return true;
    std::vector<double> f0(n), f1(n), xx = x;
    if (!c.m->vector_field(x.data(), c.p, f0.data(), &err)) return false;
    for (std::size_t j = 0; j < n; ++j) {
      const double save = xx[j], dh = 1e-7 * (std::fabs(save) + 1.0);
      xx[j] = save + dh;
      if (!c.m->vector_field(xx.data(), c.p, f1.data(), &err)) return false;
      for (std::size_t i = 0; i < n; ++i) Jf[i*n+j] = (f1[i] - f0[i]) / dh;
      xx[j] = save;
    }
    return true;
  };
  std::vector<double> Phi(n * n, 0.0);
  for (std::size_t k = 0; k < n; ++k) Phi[k*n+k] = 1.0;
  auto matmul = [&](const std::vector<double> &A, const std::vector<double> &B, std::vector<double> &out) {
    out.assign(n * n, 0.0);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t kk = 0; kk < n; ++kk) {
        const double a = A[i*n+kk]; if (a == 0.0) continue;
        for (std::size_t j = 0; j < n; ++j) out[i*n+j] += a * B[kk*n+j];
      }
  };
  std::vector<double> xcur(n), TJ(n * n), k1, k2, k3, k4, tmp(n * n);
  const int steps = (int)M * 8;   /* 8 RK4 substeps per mesh interval for stiff cycles */
  const double h = 1.0 / (double)steps;
  auto deriv = [&](double s, const std::vector<double> &Ph, std::vector<double> &dPh) -> bool {
    orbit_at_s(s, xcur);
    if (!jac_at(xcur)) return false;
    for (std::size_t i = 0; i < n*n; ++i) TJ[i] = T * Jf[i];
    matmul(TJ, Ph, dPh);
    return true;
  };
  for (int st = 0; st < steps; ++st) {
    const double s = (double)st * h;
    if (!deriv(s, Phi, k1)) { mult_out->clear(); return; }
    for (std::size_t i = 0; i < n*n; ++i) tmp[i] = Phi[i] + 0.5*h*k1[i];
    if (!deriv(s + 0.5*h, tmp, k2)) { mult_out->clear(); return; }
    for (std::size_t i = 0; i < n*n; ++i) tmp[i] = Phi[i] + 0.5*h*k2[i];
    if (!deriv(s + 0.5*h, tmp, k3)) { mult_out->clear(); return; }
    for (std::size_t i = 0; i < n*n; ++i) tmp[i] = Phi[i] + h*k3[i];
    if (!deriv(s + h, tmp, k4)) { mult_out->clear(); return; }
    for (std::size_t i = 0; i < n*n; ++i)
      Phi[i] += (h/6.0) * (k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
  }
  eigenvalues(Phi, n, mult_out);
}

/* Bordered residual for pseudo-arclength: unknowns V = (U, p), p in last slot. */
bool cyc_residual_p(CycleCtx c, const std::vector<double> &V,
                    std::size_t Ulen, std::vector<double> *Fout) {
  c.p = V[Ulen];
  std::vector<double> U(V.begin(), V.begin() + Ulen);
  return cyc_residual(c, U, Fout);
}

}  // namespace

/* Classify a sample's Floquet multipliers into stability + bifurcation flags.
 * Drops the trivial multiplier (~1) before judging the rest. */
static void classify_floquet(const std::vector<Complex> &mult, CycleSample *s) {
  s->floquet_re.clear(); s->floquet_im.clear();
  if (mult.empty()) return;
  /* find and drop the one multiplier closest to +1 (the trivial one) */
  std::size_t triv = 0; double best = 1e300;
  for (std::size_t i = 0; i < mult.size(); ++i) {
    const double d = std::hypot(mult[i].real() - 1.0, mult[i].imag());
    if (d < best) { best = d; triv = i; }
  }
  double maxmag = 0.0; bool pd = false, ns = false;
  for (std::size_t i = 0; i < mult.size(); ++i) {
    s->floquet_re.push_back(mult[i].real());
    s->floquet_im.push_back(mult[i].imag());
    if (i == triv) continue;
    const double mag = std::hypot(mult[i].real(), mult[i].imag());
    if (mag > maxmag) maxmag = mag;
    if (mult[i].real() < -1.0 + 0.02 && std::fabs(mult[i].imag()) < 0.02) pd = true;
    if (std::fabs(mult[i].imag()) > 0.02 && std::fabs(mag - 1.0) < 0.03) ns = true;
  }
  s->max_nontrivial_mult = maxmag;
  s->stable = (maxmag < 1.0 + 1e-3);
  s->is_pd = pd;
  s->is_ns = ns;
}

/* Refine a rough periodic-orbit guess by SIMULATION: integrate the ODE from a
 * point of the guess to settle onto the attracting cycle, detect one full
 * period by a Poincare-section return (a coordinate crossing in a consistent
 * direction), and resample that loop to `mesh` points. This lets cycle
 * continuation start from a crude guess (e.g. a circle) even when the true
 * cycle is far from circular (van der Pol relaxation), instead of demanding a
 * near-converged orbit up front. Returns false if no clean period is found
 * (e.g. the cycle is unstable, so forward integration won't converge to it). */
static bool refine_cycle_seed_by_simulation(
    const Model &m, double p, const std::vector<std::vector<double>> &guess,
    std::size_t mesh, double period_guess,
    std::vector<std::vector<double>> *out_points, double *out_period) {
  const std::size_t n = m.n;
  if (n < 2 || guess.empty() || mesh < 4) return false;
  double dt = period_guess > 0 ? period_guess / 2000.0 : 0.005;
  if (!(dt > 0) || !std::isfinite(dt)) dt = 0.005;
  dt = std::min(dt, 0.02);
  auto field = [&](const std::vector<double> &x, std::vector<double> &f) -> bool {
    f.assign(n, 0.0); std::string e; return m.vector_field(x.data(), p, f.data(), &e);
  };
  std::vector<double> k1(n), k2(n), k3(n), k4(n), tmp(n);
  auto rk4 = [&](std::vector<double> &x) -> bool {
    if (!field(x, k1)) return false;
    for (std::size_t i=0;i<n;i++) tmp[i]=x[i]+0.5*dt*k1[i];
    if(!field(tmp,k2)) return false;
    for (std::size_t i=0;i<n;i++) tmp[i]=x[i]+0.5*dt*k2[i];
    if(!field(tmp,k3)) return false;
    for (std::size_t i=0;i<n;i++) tmp[i]=x[i]+dt*k3[i];
    if(!field(tmp,k4)) return false;
    for (std::size_t i=0;i<n;i++) x[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0;
    double s=0; for(double v:x) s+=v*v; return std::isfinite(s) && s < 1e18;
  };
  std::vector<double> x = guess[guess.size()/2];
  const double Tg = period_guess > 0 ? period_guess : 2*M_PI;
  long settle = (long)(40.0 * Tg / dt);
  settle = std::min<long>(settle, 4000000);
  for (long s=0;s<settle;s++) if (!rk4(x)) return false;
  std::vector<double> lo(n, 1e300), hi(n, -1e300), sum(n, 0.0);
  std::vector<double> xt = x; long probe = std::min<long>((long)(2.0*Tg/dt), 400000);
  if (probe < 10) probe = 10;
  for (long s=0;s<probe;s++){ for(std::size_t i=0;i<n;i++){lo[i]=std::min(lo[i],xt[i]);hi[i]=std::max(hi[i],xt[i]);sum[i]+=xt[i];} if(!rk4(xt)) return false; }
  std::size_t sc=0; double sw=-1; for(std::size_t i=0;i<n;i++){ double r=hi[i]-lo[i]; if(r>sw){sw=r;sc=i;} }
  if (!(sw>1e-9)) return false;
  const double level = sum[sc]/(double)probe;
  std::vector<std::vector<double>> loop;
  long maxsteps = (long)(20.0*Tg/dt);
  double prev = x[sc] - level; bool armed=false; int crossings=0; double Tlen=0;
  loop.push_back(x);
  for (long s=0;s<maxsteps;s++){
    if (!rk4(x)) return false;
    Tlen += dt;
    double cur = x[sc]-level;
    if (prev < 0.0 && cur >= 0.0) {
      crossings++;
      if (crossings==1) { loop.clear(); Tlen=0; loop.push_back(x); armed=true; }
      else if (crossings==2 && armed) { loop.push_back(x); break; }
    } else if (armed) {
      loop.push_back(x);
    }
    prev = cur;
  }
  if (crossings < 2 || loop.size() < 8 || Tlen <= 0) return false;
  out_points->assign(mesh, std::vector<double>(n, 0.0));
  for (std::size_t i=0;i<mesh;i++){
    double f = (double)i/(double)mesh * (double)(loop.size()-1);
    std::size_t j = (std::size_t)f; double a = f - (double)j;
    std::size_t j1 = std::min(j+1, loop.size()-1);
    for (std::size_t k=0;k<n;k++) (*out_points)[i][k] = (1-a)*loop[j][k] + a*loop[j1][k];
  }
  *out_period = Tlen;
  return true;
}

/* Pseudo-arclength continuation of a periodic orbit in the extended space
 * V = (U, p). Predictor: step along the branch tangent (null vector of the
 * bordered Jacobian). Corrector: Newton on [ cycle BVP ; arclength constraint ]
 * which keeps the solution on a plane a distance ds along the tangent. Because
 * p is an unknown (not marched), this turns around folds of cycles (LPC) where
 * monotone-in-p continuation stalls. */
static CycleBranch continue_limit_cycle_arclength(
    const Model &m, const std::vector<std::vector<double>> &guess_points,
    double period_guess, double p0, const CycleSettings &settings) {
  CycleBranch out;
  const std::size_t n = m.n, M = guess_points.size();
  if (n == 0 || M < 4 || period_guess <= 0) { out.message = "need a periodic-orbit guess (>=4 points) and T>0"; return out; }
  for (const auto &g : guess_points) if (g.size() != n) { out.message = "guess point dimension mismatch"; return out; }
  const std::size_t Ulen = M * n + 1;       /* orbit + period */
  const std::size_t NV = Ulen + 1;          /* + parameter    */

  CycleCtx c; c.m = &m; c.n = n; c.mesh = M; c.p = p0;
  c.pin_mode = false;

  /* pack V0 = (U, p) and converge the initial cycle at p0 with a pinned phase */
  std::vector<double> V(NV, 0.0);
  for (std::size_t i = 0; i < M; ++i) for (std::size_t k = 0; k < n; ++k) V[i*n+k] = guess_points[i][k];
  V[M*n] = period_guess; V[Ulen] = p0;
  {
    CycleCtx c0 = c; c0.pin_mode = true; c0.pin_val = V[0];
    std::vector<double> U0(V.begin(), V.begin() + Ulen);
    if (!cyc_newton(c0, &U0, settings.newton_iters, settings.newton_tol)) {
      /* The supplied guess didn't converge. Try to obtain a better seed by
       * SIMULATION (settle onto the attracting cycle, detect one period) and
       * retry -- this rescues crude guesses for far-from-circular cycles. */
      std::vector<std::vector<double>> sg; double sT = 0;
      bool reseeded = refine_cycle_seed_by_simulation(m, p0, guess_points, M, period_guess, &sg, &sT);
      if (reseeded) {
        for (std::size_t i = 0; i < M; ++i) for (std::size_t k = 0; k < n; ++k) V[i*n+k] = sg[i][k];
        V[M*n] = sT; V[Ulen] = p0;
        c0.pin_val = V[0];
        U0.assign(V.begin(), V.begin() + Ulen);
      }
      if (!reseeded || !cyc_newton(c0, &U0, settings.newton_iters, settings.newton_tol)) {
        out.message = reseeded
          ? "Newton did not converge on the initial cycle even after a simulation reseed (try a larger mesh)"
          : "Newton did not converge on the initial cycle and simulation reseed found no stable cycle (try a better guess / larger mesh)";
        return out;
      }
    }
    for (std::size_t i = 0; i < Ulen; ++i) V[i] = U0[i];
  }

  /* phase condition: orthogonality to the field at a fixed reference X0 */
  auto set_phase = [&](const std::vector<double> &Vref) {
    c.pin_mode = false;
    c.phase_ref.assign(Vref.begin(), Vref.begin() + n);
    c.p = Vref[Ulen];
    std::vector<double> fr;
    if (cyc_field(c, c.phase_ref.data(), &fr)) c.phase_dir = fr; else c.phase_dir.assign(n, 0.0);
  };

  /* build the (Ulen) x (NV) Jacobian of the cycle residual wrt V by finite diff */
  auto resid_jac = [&](const std::vector<double> &Vc, std::vector<double> *F, std::vector<double> *J) -> bool {
    if (!cyc_residual_p(c, Vc, Ulen, F)) return false; /* length Ulen */
    J->assign(Ulen * NV, 0.0);
    std::vector<double> Vt = Vc, Fp;
    for (std::size_t j = 0; j < NV; ++j) {
      const double save = Vt[j], dh = 1e-7 * (std::fabs(save) + 1.0);
      Vt[j] = save + dh;
      if (!cyc_residual_p(c, Vt, Ulen, &Fp)) return false;
      for (std::size_t i = 0; i < Ulen; ++i) (*J)[i*NV+j] = (Fp[i] - (*F)[i]) / dh;
      Vt[j] = save;
    }
    return true;
  };

  /* tangent: null vector of the Ulen x NV Jacobian (one-dim kernel generically).
   * Solve [J; e_k^T] t = e_last for some pivot, then normalize; pick orientation
   * consistent with prev. */
  std::vector<double> tangent(NV, 0.0); tangent[Ulen] = 1.0; /* initial guess: increase p */
  auto compute_tangent = [&](const std::vector<double> &Vc, std::vector<double> &tan_io) -> bool {
    std::vector<double> F, J;
    set_phase(Vc);
    if (!resid_jac(Vc, &F, &J)) return false;
    /* augment with the previous tangent as the last row to fix the kernel scale:
     * [ J ] t = [ 0 ]
     * [ t_prev^T ]   [ 1 ]  */
    std::vector<double> A(NV * NV, 0.0), b(NV, 0.0);
    for (std::size_t i = 0; i < Ulen; ++i)
      for (std::size_t j = 0; j < NV; ++j) A[i*NV+j] = J[i*NV+j];
    for (std::size_t j = 0; j < NV; ++j) A[Ulen*NV+j] = tan_io[j];
    b[Ulen] = 1.0;
    std::vector<double> t;
    if (!solve_linear(A, b, &t)) return false;
    double nrm = 0; for (double v : t) nrm += v*v; nrm = std::sqrt(nrm);
    if (!(nrm > 0) || !std::isfinite(nrm)) return false;
    for (double &v : t) v /= nrm;
    /* keep orientation continuous */
    double dot = 0; for (std::size_t j = 0; j < NV; ++j) dot += t[j]*tan_io[j];
    if (dot < 0) for (double &v : t) v = -v;
    tan_io = t;
    return true;
  };

  auto record = [&](const std::vector<double> &Vc) {
    CycleSample s; s.p = Vc[Ulen]; s.period = Vc[M*n];
    std::vector<double> U(Vc.begin(), Vc.begin() + Ulen);
    cyc_amp(U, n, M, &s.amplitude, &s.min0, &s.max0);
    CycleCtx cf = c; cf.p = Vc[Ulen];
    s.fold_test = cycle_fold_test(cf, U);
    if (settings.compute_floquet) {
      std::vector<Complex> mult; cycle_floquet(cf, U, &mult); classify_floquet(mult, &s);
    }
    out.samples.push_back(s);
  };

  set_phase(V);   /* establish the integral phase reference before recording */
  const std::vector<double> V_seed = V;
  record(V);
  /* initial tangent from the converged point */
  if (!compute_tangent(V, tangent)) { out.message = "could not form branch tangent"; out.ok = out.samples.size() > 1; return out; }
  const std::vector<double> tangent_seed = tangent;

  const double ds = settings.ds > 0 ? settings.ds : 0.05;
  bool turned = false;
  const int steps_per_dir = std::max(1, settings.max_steps / 2);

  for (int dir = 0; dir < 2; ++dir) {
    V = V_seed;
    tangent = tangent_seed;
    if (dir == 1) for (double &v : tangent) v = -v; /* opposite branch direction */
    double prev_dp_sign = (tangent[Ulen] >= 0 ? 1.0 : -1.0);

    for (int step = 0; step < steps_per_dir; ++step) {
      /* predictor */
      std::vector<double> Vp = V;
      for (std::size_t j = 0; j < NV; ++j) Vp[j] += ds * tangent[j];

      /* corrector: Newton on [ cycle_resid(V) ; tangent . (V - Vp) ] = 0 */
      std::vector<double> Vc = Vp;
      bool ok = false;
      for (int it = 0; it < settings.newton_iters; ++it) {
        std::vector<double> F, J;
        set_phase(Vc);
        if (!resid_jac(Vc, &F, &J)) break;
        double g = -ds;
        for (std::size_t j = 0; j < NV; ++j) g += tangent[j] * (Vc[j] - V[j]);
        double res = 0; for (double v : F) res += v*v; res = std::sqrt(res + g*g);
        if (res < settings.newton_tol) { ok = true; break; }
        std::vector<double> A(NV * NV, 0.0), rhs(NV, 0.0);
        for (std::size_t i = 0; i < Ulen; ++i) {
          for (std::size_t j = 0; j < NV; ++j) A[i*NV+j] = J[i*NV+j];
          rhs[i] = -F[i];
        }
        for (std::size_t j = 0; j < NV; ++j) A[Ulen*NV+j] = tangent[j];
        rhs[Ulen] = -g;
        std::vector<double> dV;
        if (!solve_linear(A, rhs, &dV)) break;
        double nrm = 0; for (double v : dV) nrm += v*v; nrm = std::sqrt(nrm);
        double damp = (nrm > 1.0) ? 1.0 / nrm : 1.0;
        for (std::size_t j = 0; j < NV; ++j) Vc[j] += damp * dV[j];
      }
      if (!ok) break;
      if (!(Vc[M*n] > 0) || !std::isfinite(Vc[M*n])) break;
      if (Vc[Ulen] < settings.p_min - 1e-9 || Vc[Ulen] > settings.p_max + 1e-9) {
        record(Vc); break;
      }
      V = Vc;
      /* adaptive remesh: redistribute the mesh points by arclength so stiff /
       * relaxation cycles (sharp corners) stay well-resolved as the branch
       * moves. Updates c.frac and the orbit points; period & parameter keep. */
      if (settings.adaptive_mesh) {
        std::vector<double> U(V.begin(), V.begin() + Ulen), Unew;
        CycleCtx cm = c; cm.p = V[Ulen];
        if (cyc_remesh(cm, U, &Unew)) {
          c.frac = cm.frac;
          for (std::size_t i = 0; i < Ulen; ++i) V[i] = Unew[i];
          /* re-converge on the new mesh so V is an exact solution there */
          set_phase(V);
          std::vector<double> Ufix(V.begin(), V.begin() + Ulen);
          CycleCtx cf = c; cf.p = V[Ulen]; cf.pin_mode = false;
          cf.phase_ref = c.phase_ref; cf.phase_dir = c.phase_dir;
          if (cyc_newton(cf, &Ufix, settings.newton_iters, settings.newton_tol))
            for (std::size_t i = 0; i < Ulen; ++i) V[i] = Ufix[i];
        }
      }
      if (!compute_tangent(V, tangent)) break;
      const double dp_sign = (tangent[Ulen] >= 0 ? 1.0 : -1.0);
      record(V);
      if (dp_sign * prev_dp_sign < 0) { out.samples.back().is_fold = true; turned = true; }
      prev_dp_sign = dp_sign;
    }
  }

  out.turned = turned;
  std::sort(out.samples.begin(), out.samples.end(),
            [](const CycleSample &a, const CycleSample &b){ return a.p < b.p; });
  /* also flag LPC by fold-test sign change (belt and suspenders) */
  for (std::size_t i = 1; i < out.samples.size(); ++i) {
    const double a = out.samples[i-1].fold_test, b = out.samples[i].fold_test;
    if (std::isfinite(a) && std::isfinite(b) && a * b < 0.0)
      out.samples[std::fabs(a) < std::fabs(b) ? i-1 : i].is_fold = true;
  }
  out.ok = out.samples.size() > 1;
  out.message = out.ok ? ("traced " + std::to_string(out.samples.size()) + " cycles" +
                          (turned ? " (turned around a fold of cycles)" : ""))
                       : "could not continue the cycle";
  return out;
}

CycleBranch continue_limit_cycle(const Model &m,
                                 const std::vector<std::vector<double>> &guess_points,
                                 double period_guess, double p0,
                                 const CycleSettings &settings) {
  if (settings.arclength)
    return continue_limit_cycle_arclength(m, guess_points, period_guess, p0, settings);
  CycleBranch out;
  const std::size_t n = m.n;
  const std::size_t M = guess_points.size();
  if (n == 0 || M < 4 || period_guess <= 0) { out.message = "need a periodic-orbit guess (>=4 points) and T>0"; return out; }
  for (const auto &g : guess_points) if (g.size() != n) { out.message = "guess point dimension mismatch"; return out; }

  /* pack initial U */
  std::vector<double> U(M * n + 1, 0.0);
  for (std::size_t i = 0; i < M; ++i) for (std::size_t k = 0; k < n; ++k) U[i*n+k] = guess_points[i][k];
  U[M*n] = period_guess;

  CycleCtx c; c.m = &m; c.n = n; c.mesh = M; c.p = p0;
  c.pin_mode = true; c.pin_val = U[0]; /* pin first coord of X0 for the first solve */

  if (!cyc_newton(c, &U, settings.newton_iters, settings.newton_tol)) {
    /* rescue a crude guess via a simulation-based seed, then retry */
    std::vector<std::vector<double>> sg; double sT = 0;
    bool reseeded = refine_cycle_seed_by_simulation(m, p0, guess_points, M, period_guess, &sg, &sT);
    if (reseeded) {
      for (std::size_t i = 0; i < M; ++i) for (std::size_t k = 0; k < n; ++k) U[i*n+k] = sg[i][k];
      U[M*n] = sT; c.pin_val = U[0];
    }
    if (!reseeded || !cyc_newton(c, &U, settings.newton_iters, settings.newton_tol)) {
      out.message = reseeded
        ? "Newton did not converge on the initial cycle even after a simulation reseed (try a larger mesh)"
        : "Newton did not converge on the initial cycle and simulation reseed found no stable cycle (try a better guess / larger mesh)";
      return out;
    }
  }

  auto record = [&](double p) {
    CycleSample s; s.p = p; s.period = U[M*n];
    cyc_amp(U, n, M, &s.amplitude, &s.min0, &s.max0);
    s.stable = false; /* Floquet stability left as future work */
    s.fold_test = cycle_fold_test(c, U);
    out.samples.push_back(s);
  };

  /* switch to an integral-style phase condition for continuation robustness:
   * orthogonality of X0's increment to the field at the previous X0. */
  auto set_phase = [&](const std::vector<double> &Uref) {
    c.pin_mode = false;
    c.phase_ref.assign(Uref.begin(), Uref.begin() + n);
    std::vector<double> fr;
    if (cyc_field(c, c.phase_ref.data(), &fr)) c.phase_dir = fr;
    else c.phase_dir.assign(n, 0.0);
  };

  /* trace both directions in the parameter */
  for (int dir = 0; dir < 2; ++dir) {
    std::vector<double> Uc = U;
    double p = p0;
    const double dp = (dir == 0 ? settings.dp : -settings.dp);
    if (dir == 0) record(p);
    int steps = settings.max_steps / 2;
    for (int s = 0; s < steps; ++s) {
      const double pnext = p + dp;
      if (pnext < settings.p_min || pnext > settings.p_max) break;
      CycleCtx cc = c; cc.p = pnext;
      set_phase(Uc); cc.pin_mode = c.pin_mode; cc.phase_ref = c.phase_ref; cc.phase_dir = c.phase_dir;
      std::vector<double> Un = Uc;
      if (!cyc_newton(cc, &Un, settings.newton_iters, settings.newton_tol)) break;
      /* sanity: period must stay positive and finite */
      if (!(Un[M*n] > 0) || !std::isfinite(Un[M*n])) break;
      Uc = Un; p = pnext; c = cc;
      CycleSample smp; smp.p = p; smp.period = Uc[M*n];
      cyc_amp(Uc, n, M, &smp.amplitude, &smp.min0, &smp.max0);
      smp.fold_test = cycle_fold_test(c, Uc);
      out.samples.push_back(smp);
    }
  }
  /* sort samples by parameter for a clean plot */
  std::sort(out.samples.begin(), out.samples.end(),
            [](const CycleSample &a, const CycleSample &b){ return a.p < b.p; });
  /* mark fold-of-cycles (LPC). Two signatures: (a) a sign change of the fold
   * test between adjacent samples, or (b) the branch TERMINATES at a fold —
   * the cycle solver can't step past a saddle-node of cycles, so the test
   * collapses toward zero at the last reachable sample. We flag an endpoint
   * whose |fold_test| is far below the branch's typical magnitude. */
  for (std::size_t i = 1; i < out.samples.size(); ++i) {
    const double a = out.samples[i-1].fold_test, b = out.samples[i].fold_test;
    if (std::isfinite(a) && std::isfinite(b) && a * b < 0.0)
      out.samples[std::fabs(a) < std::fabs(b) ? i-1 : i].is_fold = true;
  }
  if (out.samples.size() >= 4) {
    /* median |fold_test| as the branch's scale */
    std::vector<double> mags;
    for (const auto &s : out.samples) if (std::isfinite(s.fold_test)) mags.push_back(std::fabs(s.fold_test));
    if (!mags.empty()) {
      std::sort(mags.begin(), mags.end());
      const double med = mags[mags.size()/2];
      /* check both endpoints (the branch stops at a fold) */
      for (std::size_t e : {std::size_t(0), out.samples.size()-1}) {
        const double f = std::fabs(out.samples[e].fold_test);
        if (std::isfinite(f) && med > 0 && f < med * 1e-3)
          out.samples[e].is_fold = true;
      }
    }
  }
  out.ok = out.samples.size() > 1;
  out.message = out.ok ? ("traced " + std::to_string(out.samples.size()) + " cycles")
                       : "could not continue the cycle";
  return out;
}

/* Two-parameter fold-of-cycles (LPC) curve. For each q across [q_min,q_max] we
 * build the single-parameter cycle model at that q, continue the cycle in p,
 * and look for a fold-of-cycles (a sign change of the cycle fold test). The
 * (p,q) where it occurs is one point of the LPC curve. The seed cycle from the
 * previous q is reused as the next guess, so the locus is followed smoothly. */
LPCCurve lpc_curve(const Model2 &m,
                   const std::vector<std::vector<double>> &guess_points,
                   double period_guess, double p0, double q0,
                   const TwoParamSettings &settings, const CycleSettings &cyc) {
  LPCCurve out;
  const std::size_t n = m.n;
  if (n < 2 || guess_points.size() < 4) { out.message = "need a 2+ D system and a cycle guess"; return out; }

  const int nq = std::max(8, std::min(settings.max_points, 80));
  std::vector<std::vector<double>> guess = guess_points;
  double per = period_guess;
  (void)q0; /* the seed cycle was simulated at q0; the scan covers [q_min,q_max] */

  auto model_at_q = [&](double qfix) {
    Model mm; mm.n = n;
    mm.vector_field = [&m, qfix](const double *xx, double pp, double *fo, std::string *er) {
      return m.vector_field(xx, pp, qfix, fo, er);
    };
    return mm;
  };

  for (int iq = 0; iq <= nq; ++iq) {
    const double q = settings.q_min + (settings.q_max - settings.q_min) * (double)iq / nq;
    Model mm = model_at_q(q);
    CycleSettings cs = cyc;
    cs.p_min = settings.p_min; cs.p_max = settings.p_max;
    CycleBranch br = continue_limit_cycle(mm, guess, per, p0, cs);
    if (!br.ok || br.samples.empty()) continue;
    /* find the LPC sample on this branch */
    for (const auto &smp : br.samples)
      if (smp.is_fold) {
        LPCPoint pt; pt.p = smp.p; pt.q = q; pt.period = smp.period; pt.amplitude = smp.amplitude;
        out.points.push_back(pt);
        break;
      }
    /* reseed the next q from the middle sample of this branch (keeps the guess
     * close as the locus moves) and update period/p0 toward the LPC if found. */
    const CycleSample &mid = br.samples[br.samples.size()/2];
    per = mid.period; p0 = mid.p;
  }
  out.ok = out.points.size() >= 2;
  out.message = out.ok ? ("traced " + std::to_string(out.points.size()) + " LPC points")
                       : "no fold-of-cycles found in this (p,q) window (the cycle may not fold here)";
  return out;
}

/* ---- Hopf first Lyapunov coefficient ------------------------------------ */
namespace {

using Cplx = std::complex<double>;

/* Solve the complex n x n system M z = b by writing it as a real 2n x 2n
 * system [[Re,-Im],[Im,Re]] [zr; zi] = [br; bi] and reusing solve_linear. */
bool solve_complex(const std::vector<Cplx> &M, const std::vector<Cplx> &b,
                   std::size_t n, std::vector<Cplx> *z) {
  std::vector<double> A(4 * n * n, 0.0), rhs(2 * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const double re = M[i * n + j].real(), im = M[i * n + j].imag();
      A[(i) * (2 * n) + (j)] = re;          /* top-left  Re */
      A[(i) * (2 * n) + (j + n)] = -im;     /* top-right -Im */
      A[(i + n) * (2 * n) + (j)] = im;      /* bot-left  Im */
      A[(i + n) * (2 * n) + (j + n)] = re;  /* bot-right Re */
    }
    rhs[i] = b[i].real();
    rhs[i + n] = b[i].imag();
  }
  std::vector<double> sol;
  if (!solve_linear(A, rhs, &sol)) return false;
  z->assign(n, Cplx(0, 0));
  for (std::size_t i = 0; i < n; ++i) (*z)[i] = Cplx(sol[i], sol[i + n]);
  return true;
}

/* B(u,v): the symmetric bilinear form of second derivatives of f at x,
 *   B_i(u,v) = sum_{j,k} d^2 f_i/dx_j dx_k * u_j v_k,
 * by central finite differences of the directional second derivative
 *   B(u,v) = (1/2)[ D^2 f(x; u+v) - D^2 f(x; u) - D^2 f(x; v) ]  (polarization)
 * where D^2 f(x; w) = (f(x+hw) - 2 f(x) + f(x-hw)) / h^2 (real directions).
 * u,v are complex; we expand bilinearly from real/imag parts. */
struct DerivCtx {
  const Model *m; const double *x0; double p; std::size_t n; double h;
  std::vector<double> f_p, f_m, f_0, xt; std::string *err;
};

bool dir_second_real(DerivCtx &c, const std::vector<double> &w, std::vector<double> *out) {
  /* (f(x+hw) - 2 f0 + f(x-hw)) / h^2 */
  const std::size_t n = c.n; const double h = c.h;
  c.xt.assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) c.xt[i] = c.x0[i] + h * w[i];
  if (!c.m->vector_field(c.xt.data(), c.p, c.f_p.data(), c.err)) return false;
  for (std::size_t i = 0; i < n; ++i) c.xt[i] = c.x0[i] - h * w[i];
  if (!c.m->vector_field(c.xt.data(), c.p, c.f_m.data(), c.err)) return false;
  out->assign(n, 0.0);
  for (std::size_t i = 0; i < n; ++i) (*out)[i] = (c.f_p[i] - 2.0 * c.f_0[i] + c.f_m[i]) / (h * h);
  return true;
}

/* real bilinear B(a,b) for real vectors via polarization */
bool B_real(DerivCtx &c, const std::vector<double> &a, const std::vector<double> &b,
            std::vector<double> *out) {
  std::vector<double> spp, sa, sb, wpp(c.n);
  for (std::size_t i = 0; i < c.n; ++i) wpp[i] = a[i] + b[i];
  if (!dir_second_real(c, wpp, &spp)) return false;
  if (!dir_second_real(c, a, &sa)) return false;
  if (!dir_second_real(c, b, &sb)) return false;
  out->assign(c.n, 0.0);
  for (std::size_t i = 0; i < c.n; ++i) (*out)[i] = 0.5 * (spp[i] - sa[i] - sb[i]);
  return true;
}

/* complex B(u,v) expanded from real/imag: with u=ur+i ui, v=vr+i vi,
 * B(u,v) = [B(ur,vr)-B(ui,vi)] + i[B(ur,vi)+B(ui,vr)]. */
bool B_cplx(DerivCtx &c, const std::vector<Cplx> &u, const std::vector<Cplx> &v,
            std::vector<Cplx> *out) {
  std::vector<double> ur(c.n), ui(c.n), vr(c.n), vi(c.n);
  for (std::size_t i = 0; i < c.n; ++i) { ur[i]=u[i].real(); ui[i]=u[i].imag(); vr[i]=v[i].real(); vi[i]=v[i].imag(); }
  std::vector<double> rr, ii, ri, ir;
  if (!B_real(c, ur, vr, &rr)) return false;
  if (!B_real(c, ui, vi, &ii)) return false;
  if (!B_real(c, ur, vi, &ri)) return false;
  if (!B_real(c, ui, vr, &ir)) return false;
  out->assign(c.n, Cplx(0,0));
  for (std::size_t i = 0; i < c.n; ++i) (*out)[i] = Cplx(rr[i]-ii[i], ri[i]+ir[i]);
  return true;
}

/* real trilinear C(a,b,d) via mixed differences using the cubic directional
 * derivative and polarization. We use C(w,w,w) = D^3 f(x;w) and the identity
 *   6 C(a,b,d) = D3(a+b+d) - D3(a+b) - D3(a+d) - D3(b+d) + D3(a)+D3(b)+D3(d)
 * with D3 f(x;w) = (f(x+2hw) -2 f(x+hw) +2 f(x-hw) - f(x-2hw)) / (2 h^3)
 * (the standard 3rd-derivative stencil, all in the real direction w). */
bool dir_third_real(DerivCtx &c, const std::vector<double> &w, std::vector<double> *out) {
  const std::size_t n = c.n; const double h = c.h;
  std::vector<double> fp2(n), fp1(n), fm1(n), fm2(n);
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+2*h*w[i];
  if (!c.m->vector_field(c.xt.data(), c.p, fp2.data(), c.err)) return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+h*w[i];
  if (!c.m->vector_field(c.xt.data(), c.p, fp1.data(), c.err)) return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-h*w[i];
  if (!c.m->vector_field(c.xt.data(), c.p, fm1.data(), c.err)) return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-2*h*w[i];
  if (!c.m->vector_field(c.xt.data(), c.p, fm2.data(), c.err)) return false;
  out->assign(n,0.0);
  for (std::size_t i=0;i<n;i++) (*out)[i] = (fp2[i]-2*fp1[i]+2*fm1[i]-fm2[i])/(2*h*h*h);
  return true;
}
/* 4th and 5th directional derivatives D^k f(x;w) along a single real direction
 * w, central stencils. Used (with polarization) to build the symmetric
 * multilinear forms D(.,.,.,.) and E(.,.,.,.,.) that enter the second Hopf
 * Lyapunov coefficient. High-order finite differences are noisy (h^k in the
 * denominator), so a larger step is used for these and l2's SIGN is the
 * reliable output. */
bool dir_fourth_real(DerivCtx &c, const std::vector<double> &w, std::vector<double> *out) {
  const std::size_t n=c.n; const double h=c.h;
  std::vector<double> fp2(n),fp1(n),f0(n),fm1(n),fm2(n);
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+2*h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fp2.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fp1.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i];
  if(!c.m->vector_field(c.xt.data(),c.p,f0.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fm1.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-2*h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fm2.data(),c.err))return false;
  out->assign(n,0.0);
  for (std::size_t i=0;i<n;i++) (*out)[i]=(fp2[i]-4*fp1[i]+6*f0[i]-4*fm1[i]+fm2[i])/(h*h*h*h);
  return true;
}
bool dir_fifth_real(DerivCtx &c, const std::vector<double> &w, std::vector<double> *out) {
  const std::size_t n=c.n; const double h=c.h;
  std::vector<double> f3(n),f2(n),f1(n),fm1(n),fm2(n),fm3(n);
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+3*h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,f3.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+2*h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,f2.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]+h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,f1.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fm1.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-2*h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fm2.data(),c.err))return false;
  for (std::size_t i=0;i<n;i++) c.xt[i]=c.x0[i]-3*h*w[i];
  if(!c.m->vector_field(c.xt.data(),c.p,fm3.data(),c.err))return false;
  out->assign(n,0.0);
  for (std::size_t i=0;i<n;i++) (*out)[i]=(f3[i]-4*f2[i]+5*f1[i]-5*fm1[i]+4*fm2[i]-fm3[i])/(2*h*h*h*h*h);
  return true;
}

bool C_real(DerivCtx &c, const std::vector<double> &a, const std::vector<double> &b,
            const std::vector<double> &d, std::vector<double> *out) {
  auto add=[&](const std::vector<double>&u,const std::vector<double>&v,const std::vector<double>&w){
    std::vector<double> s(c.n); for(std::size_t i=0;i<c.n;i++) s[i]=u[i]+v[i]+w[i]; return s; };
  std::vector<double> zero(c.n,0.0);
  std::vector<double> d_abd,d_ab,d_ad,d_bd,d_a,d_b,d_d;
  if(!dir_third_real(c, add(a,b,d), &d_abd)) return false;
  if(!dir_third_real(c, add(a,b,zero), &d_ab)) return false;
  if(!dir_third_real(c, add(a,d,zero), &d_ad)) return false;
  if(!dir_third_real(c, add(b,d,zero), &d_bd)) return false;
  if(!dir_third_real(c, a, &d_a)) return false;
  if(!dir_third_real(c, b, &d_b)) return false;
  if(!dir_third_real(c, d, &d_d)) return false;
  out->assign(c.n,0.0);
  for(std::size_t i=0;i<c.n;i++)
    (*out)[i] = (d_abd[i]-d_ab[i]-d_ad[i]-d_bd[i]+d_a[i]+d_b[i]+d_d[i])/6.0;
  return true;
}
/* Symmetric 4-linear form D(a,b,d,e) by polarization of the 4th directional
 * derivative; 24 = 4!. */
bool D_real(DerivCtx &c, const std::vector<double> &a, const std::vector<double> &b,
            const std::vector<double> &d, const std::vector<double> &e, std::vector<double> *out) {
  const std::vector<double>* v[4] = {&a,&b,&d,&e};
  out->assign(c.n, 0.0);
  std::vector<double> s(c.n), term;
  for (int mask=1; mask<16; ++mask) {
    int cnt=0; for (std::size_t i=0;i<c.n;i++) s[i]=0.0;
    for (int k=0;k<4;k++) if (mask&(1<<k)) { ++cnt; for (std::size_t i=0;i<c.n;i++) s[i]+=(*v[k])[i]; }
    if (!dir_fourth_real(c, s, &term)) return false;
    double sign = ((4-cnt)%2==0)? 1.0 : -1.0;
    for (std::size_t i=0;i<c.n;i++) (*out)[i] += sign*term[i];
  }
  for (std::size_t i=0;i<c.n;i++) (*out)[i] /= 24.0;
  return true;
}
/* Symmetric 5-linear form E by polarization of the 5th directional derivative;
 * 120 = 5!. */
bool E_real(DerivCtx &c, const std::vector<double> &a, const std::vector<double> &b,
            const std::vector<double> &d, const std::vector<double> &e, const std::vector<double> &g,
            std::vector<double> *out) {
  const std::vector<double>* v[5] = {&a,&b,&d,&e,&g};
  out->assign(c.n, 0.0);
  std::vector<double> s(c.n), term;
  for (int mask=1; mask<32; ++mask) {
    int cnt=0; for (std::size_t i=0;i<c.n;i++) s[i]=0.0;
    for (int k=0;k<5;k++) if (mask&(1<<k)) { ++cnt; for (std::size_t i=0;i<c.n;i++) s[i]+=(*v[k])[i]; }
    if (!dir_fifth_real(c, s, &term)) return false;
    double sign = ((5-cnt)%2==0)? 1.0 : -1.0;
    for (std::size_t i=0;i<c.n;i++) (*out)[i] += sign*term[i];
  }
  for (std::size_t i=0;i<c.n;i++) (*out)[i] /= 120.0;
  return true;
}
/* Generic complex multilinear evaluator: expand mform over the 2^k real/imag
 * choices, each carrying i^(#imag). Used to assemble D/E on combinations of
 * q, qbar, and real intermediate vectors for l2. */
template <typename F>
bool cmultilinear(DerivCtx &c, const std::vector<std::vector<Cplx>> &vs, F mform, std::vector<Cplx> *out) {
  const std::size_t k = vs.size(), n = c.n;
  out->assign(n, Cplx(0,0));
  std::vector<std::vector<double>> re(k, std::vector<double>(n)), im(k, std::vector<double>(n));
  for (std::size_t j=0;j<k;j++) for (std::size_t i=0;i<n;i++){ re[j][i]=vs[j][i].real(); im[j][i]=vs[j][i].imag(); }
  std::vector<const std::vector<double>*> args(k);
  for (int mask=0; mask<(1<<k); ++mask) {
    int nim=0;
    for (std::size_t j=0;j<k;j++){ if (mask&(1<<j)){ args[j]=&im[j]; ++nim; } else args[j]=&re[j]; }
    Cplx coef(1,0); for (int t=0;t<nim;t++) coef*=Cplx(0,1);
    std::vector<double> term;
    if (!mform(args, &term)) return false;
    for (std::size_t i=0;i<n;i++) (*out)[i] += coef*term[i];
  }
  return true;
}

/* complex C(q,q,qbar): expand from parts. With all args built from qr,qi this
 * is most simply assembled by summing real C over the binary expansion. */
bool C_qqqbar(DerivCtx &c, const std::vector<Cplx> &q, std::vector<Cplx> *out) {
  std::vector<double> qr(c.n), qi(c.n);
  for (std::size_t i=0;i<c.n;i++){ qr[i]=q[i].real(); qi[i]=q[i].imag(); }
  /* C(q,q,qbar) with q=qr+i qi, qbar=qr-i qi. Expand trilinearly:
   * real part  = C(qr,qr,qr) + C(qr,qi,qi) + C(qi,qr,qi) - C(qi,qi,qr)
   *            (collect terms with even # of imag factors, sign from i^2)
   * Rather than track by hand, sum over the 8 sign combinations. */
  out->assign(c.n, Cplx(0,0));
  /* q index sign: +1 -> qr, +i -> qi. For each of the three slots, the factor
   * is (qr + s*i*qi) with s=+1 for q and s=-1 for qbar (third slot). */
  const int s3[3] = {+1,+1,-1};
  for (int b0=0;b0<2;b0++) for (int b1=0;b1<2;b1++) for (int b2=0;b2<2;b2++){
    /* b=0 picks real part qr (factor 1), b=1 picks imag part qi (factor i*s) */
    const std::vector<double> &a0 = b0? qi: qr;
    const std::vector<double> &a1 = b1? qi: qr;
    const std::vector<double> &a2 = b2? qi: qr;
    Cplx coef(1,0);
    if (b0) coef *= Cplx(0, s3[0]);
    if (b1) coef *= Cplx(0, s3[1]);
    if (b2) coef *= Cplx(0, s3[2]);
    std::vector<double> cr;
    if (!C_real(c, a0, a1, a2, &cr)) return false;
    for (std::size_t i=0;i<c.n;i++) (*out)[i] += coef * cr[i];
  }
  return true;
}

Cplx cdot(const std::vector<Cplx> &a, const std::vector<Cplx> &b, std::size_t n) {
  Cplx s(0,0); for (std::size_t i=0;i<n;i++) s += std::conj(a[i]) * b[i]; return s; /* <a,b> = conj(a).b */
}

} // namespace

bool hopf_first_lyapunov(const Model &m, const std::vector<double> &x, double p,
                         double *l1, double *omega, std::string *err) {
  const std::size_t n = m.n;
  if (n < 2 || x.size() < n) { if (err) *err = "need a 2+ dim system"; return false; }
  /* Jacobian at the point */
  std::vector<double> J;
  if (!finite_diff_jacobian(m, x.data(), p, &J, err, 1e-7)) return false;
  /* eigenvalues -> find the imaginary pair */
  std::vector<Cplx> ev;
  if (!eigenvalues(J, n, &ev)) { if (err) *err="eigenvalue failure"; return false; }
  double w = 0.0; bool found = false;
  for (const auto &lam : ev) {
    if (lam.imag() > 1e-6 && std::fabs(lam.real()) < 1e-3 * (1.0 + std::fabs(lam.imag()))) {
      w = lam.imag(); found = true; break;
    }
  }
  if (!found) { if (err) *err = "no pure-imaginary pair at this point (not a Hopf point)"; return false; }

  /* critical eigenvector q: (J - i w I) q = 0  -> nullspace via inverse
   * iteration on the complex matrix. Adjoint p: (J^T + i w I) pp = 0. */
  std::vector<Cplx> Jc(n * n), JTc(n * n);
  for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++){
    Jc[i*n+j]  = Cplx(J[i*n+j], 0) - ((i==j)? Cplx(0,w):Cplx(0,0));
    JTc[i*n+j] = Cplx(J[j*n+i], 0) + ((i==j)? Cplx(0,w):Cplx(0,0));
  }
  /* inverse iteration: solve (M + eps) q_{k+1} = q_k, normalize */
  auto inv_iter = [&](std::vector<Cplx> M, std::vector<Cplx> *vec)->bool{
    for (std::size_t i=0;i<n;i++) M[i*n+i] += Cplx(1e-8,0); /* regularize singular M */
    std::vector<Cplx> v(n, Cplx(1.0/std::sqrt((double)n),0)); v[0]=Cplx(1,0);
    for (int it=0; it<60; ++it) {
      std::vector<Cplx> nv;
      if (!solve_complex(M, v, n, &nv)) return false;
      double nrm=0; for (auto &z:nv) nrm += std::norm(z); nrm=std::sqrt(nrm);
      if (nrm < 1e-300) return false;
      for (auto &z:nv) z/=nrm;
      v = nv;
    }
    *vec = v; return true;
  };
  std::vector<Cplx> q, pp;
  if (!inv_iter(Jc, &q)) { if (err) *err="eigenvector solve failed"; return false; }
  if (!inv_iter(JTc, &pp)) { if (err) *err="adjoint solve failed"; return false; }
  /* normalize so <pp,q> = 1 (complex) */
  Cplx pq = cdot(pp, q, n);
  if (std::abs(pq) < 1e-300) { if (err) *err="degenerate eigenvectors"; return false; }
  for (auto &z : pp) z = std::conj(Cplx(1,0)/pq) * z; /* scale adjoint: <pp,q>=1 */
  /* recheck and rescale robustly */
  pq = cdot(pp, q, n);
  for (auto &z : pp) z /= pq;

  /* derivative context */
  DerivCtx c; c.m=&m; c.x0=x.data(); c.p=p; c.n=n; c.h=1e-3;
  c.f_p.assign(n,0); c.f_m.assign(n,0); c.f_0.assign(n,0); c.xt.assign(n,0); c.err=err;
  if (!m.vector_field(x.data(), p, c.f_0.data(), err)) return false;

  std::vector<Cplx> qbar(n); for (std::size_t i=0;i<n;i++) qbar[i]=std::conj(q[i]);

  /* B(q,qbar) -> solve J a = -B(q,qbar) ; (real linear system since result real) */
  std::vector<Cplx> Bqqb, Bqq, Cqqqb;
  if (!B_cplx(c, q, qbar, &Bqqb)) return false;
  if (!B_cplx(c, q, q, &Bqq)) return false;
  if (!C_qqqbar(c, q, &Cqqqb)) return false;

  /* a = -J^{-1} B(q,qbar)  (J real, B(q,qbar) real up to rounding) */
  std::vector<Cplx> Jreal(n*n);
  for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++) Jreal[i*n+j]=Cplx(J[i*n+j],0);
  std::vector<Cplx> negBqqb(n); for (std::size_t i=0;i<n;i++) negBqqb[i]=-Bqqb[i];
  std::vector<Cplx> avec;
  if (!solve_complex(Jreal, negBqqb, n, &avec)) { if(err)*err="singular J in l1"; return false; }

  /* b = (2 i w I - J)^{-1} B(q,q) */
  std::vector<Cplx> M2(n*n);
  for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++)
    M2[i*n+j] = ((i==j)?Cplx(0,2*w):Cplx(0,0)) - Cplx(J[i*n+j],0);
  std::vector<Cplx> bvec;
  if (!solve_complex(M2, Bqq, n, &bvec)) { if(err)*err="singular (2iwI-J)"; return false; }

  /* h1 = C(q,q,qbar) - 2 B(q,a) + B(qbar,b) ; l1 = (1/2w) Re <pp, h1> */
  std::vector<Cplx> Bqa, Bqbb;
  if (!B_cplx(c, q, avec, &Bqa)) return false;
  if (!B_cplx(c, qbar, bvec, &Bqbb)) return false;
  std::vector<Cplx> h1(n);
  for (std::size_t i=0;i<n;i++) h1[i] = Cqqqb[i] - 2.0*Bqa[i] + Bqbb[i];
  Cplx g = cdot(pp, h1, n);
  *l1 = g.real() / (2.0 * w);
  *omega = w;
  return true;
}

/* Bogdanov-Takens normal-form coefficients (a, b).
 *
 * At a BT point the Jacobian A has a DOUBLE-ZERO eigenvalue with a 2x2 Jordan
 * block: right generalized eigenvectors q0, q1 with A q0 = 0, A q1 = q0; left
 * p0, p1 with A^T p1 = 0, A^T p0 = p1, normalized <p1,q1>=<p0,q0>=1 (and the
 * cross terms ~0). The planar BT normal form is
 *      w0' = w1,
 *      w1' = a w0^2 + b w0 w1 + ...,
 * with (Kuznetsov, Elements of Applied Bifurcation Theory, 3rd ed., eq. 8.61):
 *      a = (1/2) <p1, B(q0,q0)>,
 *      b = <p0, B(q0,q0)> + <p1, B(q0,q1)>,
 * where B is the bilinear form of second derivatives. The non-degeneracy of the
 * BT bifurcation requires a != 0 and b != 0. Computed at a point that is assumed
 * to already be (numerically) a BT point. */
bool bt_normal_form(const Model &m, const std::vector<double> &x, double p,
                    double *a_out, double *b_out, std::string *err) {
  const std::size_t n = m.n;
  if (n < 2 || x.size() < n) { if (err) *err = "need a 2+ dim system"; return false; }
  std::vector<double> J;
  if (!finite_diff_jacobian(m, x.data(), p, &J, err, 1e-7)) return false;
  std::vector<double> JT(n * n);
  for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) JT[i*n+j] = J[j*n+i];

  /* right null vector q0 (A q0 = 0) and left null vector p1 (A^T p1 = 0) via
   * inverse iteration on the regularized matrices. */
  auto inv_null = [&](const std::vector<double> &Min, std::vector<double> *vec) -> bool {
    std::vector<double> M = Min;
    for (std::size_t i = 0; i < n; ++i) M[i*n+i] += 1e-9;
    std::vector<double> v(n, 1.0/std::sqrt((double)n));
    for (int it = 0; it < 200; ++it) {
      std::vector<double> nv;
      if (!solve_linear(M, v, &nv)) return false;
      double nrm = 0; for (double z : nv) nrm += z*z; nrm = std::sqrt(nrm);
      if (nrm < 1e-300) return false;
      for (double &z : nv) z /= nrm;
      v = nv;
    }
    *vec = v; return true;
  };
  std::vector<double> q0, p1;
  if (!inv_null(J, &q0)) { if (err) *err = "BT: q0 solve failed"; return false; }
  if (!inv_null(JT, &p1)) { if (err) *err = "BT: p1 solve failed"; return false; }

  auto dot = [&](const std::vector<double> &u, const std::vector<double> &v) {
    double s = 0; for (std::size_t i = 0; i < n; ++i) s += u[i]*v[i]; return s;
  };
  { double qn = std::sqrt(dot(q0,q0)); if (qn < 1e-300) { if(err)*err="BT: q0 zero"; return false; } for (double &z:q0) z/=qn; }
  { double pn = std::sqrt(dot(p1,p1)); if (pn < 1e-300) { if(err)*err="BT: p1 zero"; return false; } for (double &z:p1) z/=pn; }
  /* fix sign conventions for reproducibility: largest-magnitude component of q0
   * positive, and <p1,...> oriented so the chain scales consistently. */
  { std::size_t im = 0; for (std::size_t i = 1; i < n; ++i) if (std::fabs(q0[i]) > std::fabs(q0[im])) im = i;
    if (q0[im] < 0) for (double &z : q0) z = -z; }
  { std::size_t im = 0; for (std::size_t i = 1; i < n; ++i) if (std::fabs(p1[i]) > std::fabs(p1[im])) im = i;
    if (p1[im] < 0) for (double &z : p1) z = -z; }

  /* Generalized eigenvectors via BORDERED systems. To make the bordered matrix
   * nonsingular for a defective A, border the A-block with the LEFT null vector
   * in the extra column and the RIGHT null vector in the extra row (Govaerts,
   * Numerical Methods for Bifurcations..., the standard bordering):
   *   [ A    p1 ] [q1]   [q0]        [ A^T  q0 ] [p0]   [p1]
   *   [ q0^T 0  ] [s ] = [0 ]   and  [ p1^T 0  ] [s ] = [0 ].
   * The q0^T row fixes the gauge <q0,q1>=0; p1 in the column spans coker A. */
  auto bordered_solve = [&](const std::vector<double> &A, const std::vector<double> &bord_col,
                            const std::vector<double> &bord_row, const std::vector<double> &rhs_top,
                            std::vector<double> *sol) -> bool {
    const std::size_t N = n + 1;
    std::vector<double> M(N * N, 0.0), r(N, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) M[i*N+j] = A[i*n+j];
      M[i*N+n] = bord_col[i];
      r[i] = rhs_top[i];
    }
    for (std::size_t j = 0; j < n; ++j) M[n*N+j] = bord_row[j];
    M[n*N+n] = 0.0; r[n] = 0.0;
    std::vector<double> z;
    if (!solve_linear(M, r, &z)) return false;
    sol->assign(z.begin(), z.begin() + n);
    return true;
  };

  std::vector<double> q1, p0;
  /* q1: A q1 = q0, border col = p1 (coker A), row = q0 (gauge <q0,q1>=0) */
  if (!bordered_solve(J, p1, q0, q0, &q1)) { if (err) *err = "BT: q1 bordered solve failed"; return false; }
  /* p0: A^T p0 = p1, border col = q0 (coker A^T), row = p1 (gauge <p1,p0>=0) */
  if (!bordered_solve(JT, q0, p1, p1, &p0)) { if (err) *err = "BT: p0 bordered solve failed"; return false; }

  /* Jordan-chain normalization (Kuznetsov 8.3.2, eq. 8.61): scale the chains so
   * that <p1,q1> = 1 and <p0,q0> = 1. */
  double s_q = dot(p1, q1);
  if (std::fabs(s_q) < 1e-300) { if (err) *err = "BT: degenerate Jordan chain (<p1,q1>=0)"; return false; }
  for (double &z : q1) z /= s_q;     /* now <p1,q1> = 1 */
  double s_p = dot(p0, q0);
  if (std::fabs(s_p) < 1e-300) { if (err) *err = "BT: degenerate (<p0,q0>=0)"; return false; }
  for (double &z : p0) z /= s_p;     /* now <p0,q0> = 1 */

  /* bilinear form B(.,.) of second derivatives */
  DerivCtx c; c.m=&m; c.x0=x.data(); c.p=p; c.n=n; c.h=1e-3;
  c.f_p.assign(n,0); c.f_m.assign(n,0); c.f_0.assign(n,0); c.xt.assign(n,0); c.err=err;
  if (!m.vector_field(x.data(), p, c.f_0.data(), err)) return false;
  std::vector<double> Bq0q0, Bq0q1;
  if (!B_real(c, q0, q0, &Bq0q0)) return false;
  if (!B_real(c, q0, q1, &Bq0q1)) return false;

  const double a = 0.5 * dot(p1, Bq0q0);
  const double b = dot(p0, Bq0q0) + dot(p1, Bq0q1);
  if (a_out) *a_out = a;
  if (b_out) *b_out = b;
  return true;
}

/* Fold (limit point) normal-form coefficient a.
 *
 * At a fold/saddle-node the Jacobian A has a simple zero eigenvalue with right
 * null vector q (A q = 0) and left null vector p (p^T A = 0). The restriction
 * of the dynamics to the center manifold is, to leading order,
 *   y' = a y^2 + ...,     a = (1/2) <p, B(q,q)> / <p, q>,
 * where B is the bilinear form of second derivatives (Kuznetsov, Elements of
 * Applied Bifurcation Theory). a != 0 is the non-degeneracy condition for a
 * quadratic fold; a ~ 0 signals a CUSP (codim-2). This is the limit-point
 * counterpart of the Hopf first Lyapunov coefficient and the quantity MatCont
 * reports at an LP. Reuses the same finite-difference bilinear B as the Hopf
 * routine. Returns the coefficient in *a and an estimate of the zero
 * eigenvalue's magnitude in *lambda0 (small => genuinely at the fold). */
bool fold_normal_form(const Model &m, const std::vector<double> &x, double p,
                      double *a, double *lambda0, std::string *err) {
  const std::size_t n = m.n;
  if (n < 1 || x.size() < n) { if (err) *err = "need a 1+ dim system"; return false; }
  std::vector<double> J;
  if (!finite_diff_jacobian(m, x.data(), p, &J, err, 1e-7)) return false;

  /* find the smallest-magnitude eigenvalue (should be ~0 at a fold) */
  std::vector<Cplx> ev;
  if (eigenvalues(J, n, &ev)) {
    double best = 1e300;
    for (const auto &lam : ev) best = std::min(best, std::abs(lam));
    if (lambda0) *lambda0 = best;
  } else if (lambda0) *lambda0 = -1.0;

  /* real null vector via inverse iteration on (M + eps I) */
  auto inv_iter_real = [&](const std::vector<double> &Min, std::vector<double> *vec)->bool{
    std::vector<double> M = Min;
    for (std::size_t i=0;i<n;i++) M[i*n+i] += 1e-9;            /* regularize */
    std::vector<double> v(n, 1.0/std::sqrt((double)n));
    for (int it=0; it<80; ++it) {
      std::vector<double> nv;
      if (!solve_linear(M, v, &nv)) return false;
      double nrm=0; for (double z:nv) nrm+=z*z; nrm=std::sqrt(nrm);
      if (nrm < 1e-300) return false;
      for (double &z:nv) z/=nrm;
      v = nv;
    }
    *vec = v; return true;
  };
  std::vector<double> JT(n*n);
  for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++) JT[i*n+j]=J[j*n+i];
  std::vector<double> q, pvec;
  if (!inv_iter_real(J, &q))  { if (err) *err="null-vector solve failed"; return false; }
  if (!inv_iter_real(JT, &pvec)) { if (err) *err="adjoint null-vector solve failed"; return false; }

  /* normalize <q,q>=1 and scale p so <p,q>=1 */
  double qn=0; for (double z:q) qn+=z*z; qn=std::sqrt(qn);
  if (qn < 1e-300) { if (err) *err="degenerate q"; return false; }
  for (double &z:q) z/=qn;
  double pq=0; for (std::size_t i=0;i<n;i++) pq+=pvec[i]*q[i];
  if (std::abs(pq) < 1e-300) { if (err) *err="p and q orthogonal (not a simple fold)"; return false; }
  for (double &z:pvec) z/=pq;

  /* B(q,q) via the shared finite-difference bilinear form */
  DerivCtx c; c.m=&m; c.x0=x.data(); c.p=p; c.n=n; c.h=1e-3;
  c.f_p.assign(n,0); c.f_m.assign(n,0); c.f_0.assign(n,0); c.xt.assign(n,0); c.err=err;
  if (!m.vector_field(x.data(), p, c.f_0.data(), err)) return false;
  std::vector<double> Bqq;
  if (!B_real(c, q, q, &Bqq)) return false;

  double val=0; for (std::size_t i=0;i<n;i++) val += pvec[i]*Bqq[i];
  if (a) *a = 0.5 * val;
  return true;
}

/* Cusp normal-form coefficient c (codim-2, on a fold curve where a -> 0). The
 * center-manifold reduction is cubic y' = c y^3 + ...,
 *   c = (1/6) <p, C(q,q,q) - 3 B(q,h2)>,  A h2 = -(B(q,q) - <p,B(q,q)> q),
 * with gauge <p,h2>=0 (Govaerts bordered solve, since A is singular).
 * c != 0 is the cusp non-degeneracy condition. */
bool cusp_normal_form(const Model &m, const std::vector<double> &x, double p,
                      double *c_out, std::string *err) {
  const std::size_t n = m.n;
  if (n < 1 || x.size() < n) { if (err) *err = "need a 1+ dim system"; return false; }
  std::vector<double> J;
  if (!finite_diff_jacobian(m, x.data(), p, &J, err, 1e-7)) return false;
  std::vector<double> JT(n*n);
  for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++) JT[i*n+j]=J[j*n+i];
  auto inv_iter_real = [&](const std::vector<double> &Min, std::vector<double> *vec)->bool{
    std::vector<double> M = Min;
    for (std::size_t i=0;i<n;i++) M[i*n+i] += 1e-9;
    std::vector<double> v(n, 1.0/std::sqrt((double)n));
    for (int it=0; it<120; ++it) {
      std::vector<double> nv;
      if (!solve_linear(M, v, &nv)) return false;
      double nrm=0; for (double z:nv) nrm+=z*z; nrm=std::sqrt(nrm);
      if (nrm < 1e-300) return false;
      for (double &z:nv) z/=nrm;
      v = nv;
    }
    *vec = v; return true;
  };
  std::vector<double> q, pvec;
  if (!inv_iter_real(J, &q))  { if (err) *err="cusp: null-vector solve failed"; return false; }
  if (!inv_iter_real(JT, &pvec)) { if (err) *err="cusp: adjoint null-vector solve failed"; return false; }
  double qn=0; for (double z:q) qn+=z*z; qn=std::sqrt(qn);
  if (qn<1e-300){ if(err)*err="cusp: degenerate q"; return false; }
  for (double &z:q) z/=qn;
  { std::size_t im=0; for (std::size_t i=1;i<n;i++) if (std::fabs(q[i])>std::fabs(q[im])) im=i; if (q[im]<0) for(double&z:q)z=-z; }
  double pq=0; for (std::size_t i=0;i<n;i++) pq+=pvec[i]*q[i];
  if (std::fabs(pq)<1e-300){ if(err)*err="cusp: p,q orthogonal (not a simple fold)"; return false; }
  for (double &z:pvec) z/=pq;

  DerivCtx c; c.m=&m; c.x0=x.data(); c.p=p; c.n=n; c.h=1e-3;
  c.f_p.assign(n,0); c.f_m.assign(n,0); c.f_0.assign(n,0); c.xt.assign(n,0); c.err=err;
  if (!m.vector_field(x.data(), p, c.f_0.data(), err)) return false;
  std::vector<double> Bqq, Cqqq;
  if (!B_real(c, q, q, &Bqq)) return false;
  if (!C_real(c, q, q, q, &Cqqq)) return false;
  double pBqq=0; for (std::size_t i=0;i<n;i++) pBqq += pvec[i]*Bqq[i];
  std::vector<double> rhs(n);
  for (std::size_t i=0;i<n;i++) rhs[i] = -(Bqq[i] - pBqq*q[i]);
  const std::size_t N=n+1;
  std::vector<double> Mb(N*N,0.0), rb(N,0.0);
  for (std::size_t i=0;i<n;i++){ for (std::size_t j=0;j<n;j++) Mb[i*N+j]=J[i*n+j];
    Mb[i*N+n]=q[i]; rb[i]=rhs[i]; }
  for (std::size_t j=0;j<n;j++) Mb[n*N+j]=pvec[j];
  Mb[n*N+n]=0.0; rb[n]=0.0;
  std::vector<double> z;
  if (!solve_linear(Mb, rb, &z)) { if(err)*err="cusp: bordered solve failed"; return false; }
  std::vector<double> h2(z.begin(), z.begin()+n);
  std::vector<double> Bqh2;
  if (!B_real(c, q, h2, &Bqh2)) return false;
  double val=0; for (std::size_t i=0;i<n;i++) val += pvec[i]*(Cqqq[i] - 3.0*Bqh2[i]);
  if (c_out) *c_out = val/6.0;
  return true;
}

/* Generalized-Hopf (Bautin) SECOND Lyapunov coefficient l2 (codim-2 on a Hopf
 * curve, where l1 -> 0). Kuznetsov's invariant Re<p,...> combination at 5th
 * order on the critical eigenspace. Uses 4th/5th finite differences, so the
 * MAGNITUDE is approximate; the SIGN (which classifies the Bautin point and the
 * direction the fold-of-cycles curve opens) is the reliable output. Returns
 * false away from a Hopf point. */
bool gh_second_lyapunov(const Model &m, const std::vector<double> &x, double p,
                        double *l2, std::string *err) {
  const std::size_t n = m.n;
  if (n < 2 || x.size() < n) { if (err) *err = "need a 2+ dim system"; return false; }
  std::vector<double> J;
  if (!finite_diff_jacobian(m, x.data(), p, &J, err, 1e-7)) return false;
  std::vector<Cplx> ev;
  if (!eigenvalues(J, n, &ev)) { if (err) *err="eigenvalue failure"; return false; }
  double w=0; bool found=false;
  for (const auto &lam : ev)
    if (lam.imag() > 1e-6 && std::fabs(lam.real()) < 1e-2*(1.0+std::fabs(lam.imag()))) { w=lam.imag(); found=true; break; }
  if (!found) { if (err) *err="no pure-imaginary pair (not a Hopf point)"; return false; }

  std::vector<Cplx> Jc(n*n), JTc(n*n);
  for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++){
    Jc[i*n+j]  = Cplx(J[i*n+j],0) - ((i==j)?Cplx(0,w):Cplx(0,0));
    JTc[i*n+j] = Cplx(J[j*n+i],0) + ((i==j)?Cplx(0,w):Cplx(0,0));
  }
  auto inv_iter=[&](std::vector<Cplx> M, std::vector<Cplx>*vec)->bool{
    for (std::size_t i=0;i<n;i++) M[i*n+i]+=Cplx(1e-8,0);
    std::vector<Cplx> v(n,Cplx(1.0/std::sqrt((double)n),0)); v[0]=Cplx(1,0);
    for (int it=0;it<80;it++){ std::vector<Cplx> nv; if(!solve_complex(M,v,n,&nv))return false;
      double nr=0; for(auto&z:nv)nr+=std::norm(z); nr=std::sqrt(nr); if(nr<1e-300)return false;
      for(auto&z:nv){z/=nr;} v=nv; }
    *vec=v; return true; };
  std::vector<Cplx> q,pp;
  if(!inv_iter(Jc,&q)){ if(err)*err="eigvec failed"; return false; }
  if(!inv_iter(JTc,&pp)){ if(err)*err="adjoint failed"; return false; }
  { double qn=0; for(auto&z:q)qn+=std::norm(z); qn=std::sqrt(qn); for(auto&z:q)z/=qn; }
  Cplx pq=cdot(pp,q,n); if(std::abs(pq)<1e-300){ if(err)*err="degenerate eigvecs"; return false; }
  for(auto&z:pp)z/=pq;
  std::vector<Cplx> qbar(n); for(std::size_t i=0;i<n;i++) qbar[i]=std::conj(q[i]);

  DerivCtx c; c.m=&m; c.x0=x.data(); c.p=p; c.n=n; c.h=2e-2;
  c.f_p.assign(n,0); c.f_m.assign(n,0); c.f_0.assign(n,0); c.xt.assign(n,0); c.err=err;
  if (!m.vector_field(x.data(), p, c.f_0.data(), err)) return false;

  auto solveK=[&](int s, const std::vector<Cplx>&rhs, std::vector<Cplx>*X)->bool{
    std::vector<Cplx> M(n*n);
    for(std::size_t i=0;i<n;i++)for(std::size_t j=0;j<n;j++)
      M[i*n+j]=((i==j)?Cplx(0,(double)s*w):Cplx(0,0))-Cplx(J[i*n+j],0);
    if (s!=1 && s!=-1) { return solve_complex(M,rhs,n,X); }
    const std::vector<Cplx>& col = (s==1)? q : qbar;
    std::vector<Cplx> rowv(n);
    if (s==1) for(std::size_t i=0;i<n;i++) rowv[i]=std::conj(pp[i]);
    else      for(std::size_t i=0;i<n;i++) rowv[i]=pp[i];
    const std::size_t N=n+1;
    std::vector<Cplx> Mb(N*N,Cplx(0,0)), rb(N,Cplx(0,0));
    for(std::size_t i=0;i<n;i++){ for(std::size_t j=0;j<n;j++) Mb[i*N+j]=M[i*n+j];
      Mb[i*N+n]=col[i]; rb[i]=rhs[i]; }
    for(std::size_t j=0;j<n;j++) Mb[n*N+j]=rowv[j];
    Mb[n*N+n]=Cplx(0,0); rb[n]=Cplx(0,0);
    std::vector<Cplx> z;
    if(!solve_complex(Mb,rb,N,&z)) return false;
    X->assign(z.begin(), z.begin()+n); return true;
  };

  std::vector<Cplx> Bqq, Bqqb;
  if(!B_cplx(c,q,q,&Bqq)) return false;
  if(!B_cplx(c,q,qbar,&Bqqb)) return false;
  std::vector<Cplx> h20,h11;
  if(!solveK(2,Bqq,&h20)){ if(err)*err="solve (2iw-A) failed"; return false; }
  { std::vector<Cplx> nb(n); for(std::size_t i=0;i<n;i++) nb[i]=Bqqb[i];
    if(!solveK(0,nb,&h11)){ if(err)*err="solve A failed"; return false; } }

  std::vector<Cplx> Cqqq, Cqqqb, Bqh20, Bqbh20, Bqh11;
  { std::vector<std::vector<Cplx>> a={q,q,q};
    if(!cmultilinear(c,a,[&](const std::vector<const std::vector<double>*>&ar,std::vector<double>*o){return C_real(c,*ar[0],*ar[1],*ar[2],o);},&Cqqq)) return false; }
  if(!C_qqqbar(c,q,&Cqqqb)) return false;
  if(!B_cplx(c,q,h20,&Bqh20)) return false;
  if(!B_cplx(c,qbar,h20,&Bqbh20)) return false;
  if(!B_cplx(c,q,h11,&Bqh11)) return false;
  std::vector<Cplx> r30(n); for(std::size_t i=0;i<n;i++) r30[i]=Cqqq[i]+3.0*Bqh20[i];
  std::vector<Cplx> h30; if(!solveK(3,r30,&h30)){ if(err)*err="solve (3iw-A) failed"; return false; }
  std::vector<Cplx> r21(n);
  for(std::size_t i=0;i<n;i++) r21[i]=Cqqqb[i]+Bqbh20[i]+2.0*Bqh11[i];
  Cplx c1=cdot(pp,r21,n);
  std::vector<Cplx> r21p(n); for(std::size_t i=0;i<n;i++) r21p[i]=r21[i]-c1*q[i];
  std::vector<Cplx> h21; if(!solveK(1,r21p,&h21)){ if(err)*err="solve (iw-A) bordered failed"; return false; }

  std::vector<Cplx> Eq, Dqqqh, Dqqqbh, Cterm1, Cterm2, Bterm1, Bterm2;
  { std::vector<std::vector<Cplx>> a={q,q,q,qbar,qbar};
    if(!cmultilinear(c,a,[&](const std::vector<const std::vector<double>*>&ar,std::vector<double>*o){return E_real(c,*ar[0],*ar[1],*ar[2],*ar[3],*ar[4],o);},&Eq)) return false; }
  { std::vector<std::vector<Cplx>> a={q,q,q,h11};
    if(!cmultilinear(c,a,[&](const std::vector<const std::vector<double>*>&ar,std::vector<double>*o){return D_real(c,*ar[0],*ar[1],*ar[2],*ar[3],o);},&Dqqqh)) return false; }
  { std::vector<std::vector<Cplx>> a={q,q,qbar,h20};
    if(!cmultilinear(c,a,[&](const std::vector<const std::vector<double>*>&ar,std::vector<double>*o){return D_real(c,*ar[0],*ar[1],*ar[2],*ar[3],o);},&Dqqqbh)) return false; }
  { std::vector<std::vector<Cplx>> a={q,qbar,h21};
    if(!cmultilinear(c,a,[&](const std::vector<const std::vector<double>*>&ar,std::vector<double>*o){return C_real(c,*ar[0],*ar[1],*ar[2],o);},&Cterm1)) return false; }
  { std::vector<std::vector<Cplx>> a={qbar,h20,h11};
    if(!cmultilinear(c,a,[&](const std::vector<const std::vector<double>*>&ar,std::vector<double>*o){return C_real(c,*ar[0],*ar[1],*ar[2],o);},&Cterm2)) return false; }
  if(!B_cplx(c,qbar,h30,&Bterm1)) return false;
  if(!B_cplx(c,h11,h21,&Bterm2)) return false;

  std::vector<Cplx> hsum(n);
  for(std::size_t i=0;i<n;i++)
    hsum[i] = Eq[i] + 5.0*Dqqqh[i] + 5.0*Dqqqbh[i]
            + 6.0*Cterm1[i] + 6.0*Cterm2[i] + 3.0*Bterm1[i] + 6.0*Bterm2[i];
  Cplx g = cdot(pp, hsum, n);
  if (l2) *l2 = g.real() / (12.0 * w);
  return true;
}

/* ---- 2D fixed-point scanning -------------------------------- */

std::vector<FixedPoint2D> scan_fixed_points_2d(const PlanarField &field,
                                               double xmin, double xmax,
                                               double ymin, double ymax,
                                               int seeds,
                                               double dedup_tol_frac) {
  std::vector<FixedPoint2D> out;
  if (xmax <= xmin || ymax <= ymin) return out;
  seeds = std::max(3, std::min(seeds, 41));
  const double w = xmax - xmin, h = ymax - ymin;
  const double dedup = dedup_tol_frac * std::sqrt(w * w + h * h);

  /* finite-difference 2x2 Jacobian of the planar field */
  auto jac = [&](double x, double y, double J[4]) -> bool {
    const double hx = 1e-6 * (std::fabs(x) + 1e-3);
    const double hy = 1e-6 * (std::fabs(y) + 1e-3);
    double up, vp, um, vm;
    if (!field.eval(x + hx, y, &up, &vp) || !field.eval(x - hx, y, &um, &vm))
      return false;
    J[0] = (up - um) / (2 * hx);  /* du/dx */
    J[2] = (vp - vm) / (2 * hx);  /* dv/dx */
    if (!field.eval(x, y + hy, &up, &vp) || !field.eval(x, y - hy, &um, &vm))
      return false;
    J[1] = (up - um) / (2 * hy);  /* du/dy */
    J[3] = (vp - vm) / (2 * hy);  /* dv/dy */
    return true;
  };

  for (int gj = 0; gj < seeds; ++gj) {
    for (int gi = 0; gi < seeds; ++gi) {
      double x = xmin + (gi + 0.5) * w / seeds;
      double y = ymin + (gj + 0.5) * h / seeds;
      bool ok = true;
      /* Newton, up to 40 iters */
      for (int it = 0; it < 40; ++it) {
        double u, v;
        if (!field.eval(x, y, &u, &v)) { ok = false; break; }
        if (std::hypot(u, v) < 1e-11) break;
        double J[4];
        if (!jac(x, y, J)) { ok = false; break; }
        const double det = J[0] * J[3] - J[1] * J[2];
        if (std::fabs(det) < 1e-14) { ok = false; break; }
        /* solve J [dx;dy] = -[u;v] */
        const double dx = (-u * J[3] + v * J[1]) / det;
        const double dy = (-v * J[0] + u * J[2]) / det;
        x += dx;
        y += dy;
        if (!std::isfinite(x) || !std::isfinite(y)) { ok = false; break; }
        if (std::hypot(dx, dy) < 1e-12) break;
      }
      if (!ok || !std::isfinite(x) || !std::isfinite(y)) continue;
      /* must be inside the region (allow small margin) and an actual root */
      double u, v;
      if (!field.eval(x, y, &u, &v) || std::hypot(u, v) > 1e-6) continue;
      if (x < xmin - 0.05 * w || x > xmax + 0.05 * w ||
          y < ymin - 0.05 * h || y > ymax + 0.05 * h)
        continue;
      /* dedup */
      bool dup = false;
      for (const auto &p : out)
        if (std::hypot(p.x - x, p.y - y) < dedup) { dup = true; break; }
      if (dup) continue;

      FixedPoint2D fp;
      fp.x = x;
      fp.y = y;
      double J[4];
      if (!jac(x, y, J)) continue;
      fp.jacobian = {J[0], J[1], J[2], J[3]};
      Classification cl = classify_equilibrium(fp.jacobian, 2);
      fp.eigenvalues = cl.eigenvalues;
      fp.label = cl.label;
      fp.is_saddle = cl.is_saddle;
      /* real eigendirections for manifold drawing */
      for (const Complex &lam : fp.eigenvalues) {
        if (std::fabs(lam.imag()) > 1e-9) {
          fp.directions.clear();
          break;
        }
        /* (J - lambda I) v = 0  -> v = (J01, lambda-J00) or (lambda-J11, J10) */
        double vx = J[1];
        double vy = lam.real() - J[0];
        if (std::fabs(vx) + std::fabs(vy) < 1e-12) {
          vx = lam.real() - J[3];
          vy = J[2];
        }
        const double n = std::hypot(vx, vy);
        if (n > 1e-12) fp.directions.push_back({vx / n, vy / n});
      }
      out.push_back(std::move(fp));
    }
  }
  return out;
}

/* ---- Lyapunov spectrum --------------------------------------- */

double kaplan_yorke_dimension(const std::vector<double> &l) {
  const std::size_t n = l.size();
  if (n == 0) return 0.0;
  if (l[0] < 0.0) return 0.0; /* attracting fixed point */
  double partial = 0.0;
  std::size_t j = 0;
  /* find largest j such that sum_{i=0}^{j} lambda_i >= 0 */
  for (std::size_t i = 0; i < n; ++i) {
    const double next = partial + l[i];
    if (next < 0.0) break;
    partial = next;
    j = i + 1; /* count of exponents included */
  }
  if (j >= n) return static_cast<double>(n); /* whole sum >= 0 */
  /* D = j + (sum of first j) / |lambda_j|  (0-indexed: l[j] is next) */
  const double denom = std::fabs(l[j]);
  if (denom < 1e-300) return static_cast<double>(j);
  return static_cast<double>(j) + partial / denom;
}

namespace {
/* y += a * x  over n entries */
inline void axpy(std::size_t n, double a, const double *x, double *y) {
  for (std::size_t i = 0; i < n; ++i) y[i] += a * x[i];
}
/* dot product */
inline double dot(std::size_t n, const double *a, const double *b) {
  double s = 0.0;
  for (std::size_t i = 0; i < n; ++i) s += a[i] * b[i];
  return s;
}
/* J*v -> out (J row-major n x n) */
inline void matvec(std::size_t n, const double *J, const double *v, double *out) {
  for (std::size_t i = 0; i < n; ++i) {
    double s = 0.0;
    for (std::size_t k = 0; k < n; ++k) s += J[i * n + k] * v[k];
    out[i] = s;
  }
}
/* Modified Gram-Schmidt on the columns of Q (n rows, m cols, column-major
 * storage: col c occupies Q[c*n .. c*n+n-1]). Writes the diagonal R
 * values (the norms after orthogonalization) into r_diag[c]. */
void mgs(std::size_t n, std::size_t m, std::vector<double> &Q,
         std::vector<double> &r_diag) {
  for (std::size_t c = 0; c < m; ++c) {
    double *vc = &Q[c * n];
    for (std::size_t k = 0; k < c; ++k) {
      const double *vk = &Q[k * n];
      const double proj = dot(n, vk, vc);
      axpy(n, -proj, vk, vc);
    }
    const double nrm = std::sqrt(std::fmax(0.0, dot(n, vc, vc)));
    r_diag[c] = nrm;
    if (nrm > 1e-300) {
      const double inv = 1.0 / nrm;
      for (std::size_t i = 0; i < n; ++i) vc[i] *= inv;
    }
  }
}
}  // namespace

LyapunovResult lyapunov_spectrum(const Model &m, const std::vector<double> &x0,
                                 double p, const LyapunovOptions &opt) {
  LyapunovResult R;
  const std::size_t n = m.n;
  if (n == 0 || x0.size() < n || !m.vector_field) {
    R.message = "lyapunov: invalid model/state";
    return R;
  }
  const long reorth = std::max<long>(1, opt.reorth_every);
  std::string err;

  std::vector<double> x(x0.begin(), x0.begin() + n);
  std::vector<double> Jbuf(n * n, 0.0);
  auto jac = [&](const double *xx, double *Jout) -> bool {
    if (m.jacobian_x) return m.jacobian_x(xx, p, Jout, &err);
    if (!finite_diff_jacobian(m, xx, p, &Jbuf, &err, 1e-7)) return false;
    std::copy(Jbuf.begin(), Jbuf.end(), Jout);
    return true;
  };

  /* one step of the base state; returns false on error */
  auto step_state = [&](std::vector<double> &xx) -> bool {
    if (opt.is_map) {
      std::vector<double> f(n);
      if (!m.vector_field(xx.data(), p, f.data(), &err)) return false;
      xx = f; /* map: x <- f(x) */
      return true;
    }
    /* RK4 for the ODE */
    std::vector<double> k1(n), k2(n), k3(n), k4(n), tmp(n);
    const double h = opt.dt;
    if (!m.vector_field(xx.data(), p, k1.data(), &err)) return false;
    for (std::size_t i = 0; i < n; ++i) tmp[i] = xx[i] + 0.5 * h * k1[i];
    if (!m.vector_field(tmp.data(), p, k2.data(), &err)) return false;
    for (std::size_t i = 0; i < n; ++i) tmp[i] = xx[i] + 0.5 * h * k2[i];
    if (!m.vector_field(tmp.data(), p, k3.data(), &err)) return false;
    for (std::size_t i = 0; i < n; ++i) tmp[i] = xx[i] + h * k3[i];
    if (!m.vector_field(tmp.data(), p, k4.data(), &err)) return false;
    for (std::size_t i = 0; i < n; ++i)
      xx[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return true;
  };

  /* settle onto the attractor */
  for (long t = 0; t < opt.transient; ++t) {
    if (!step_state(x)) { R.message = "lyapunov: " + err; return R; }
    bool fin = true;
    for (double v : x) fin = fin && std::isfinite(v);
    if (!fin) { R.message = "lyapunov: trajectory diverged in transient"; return R; }
  }

  /* tangent frame Q: n columns of length n, init to identity */
  std::vector<double> Q(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) Q[i * n + i] = 1.0;
  std::vector<double> accum(n, 0.0);
  std::vector<double> r_diag(n, 0.0);

  /* propagate one tangent vector by the linearization for `reorth` base
   * steps. For a map: v <- J(x_k) v each step. For an ODE: integrate
   * v' = J(x(t)) v with RK4 using the Jacobian at the current x (frozen
   * across the substep — adequate for small dt). */
  std::vector<double> J(n * n);

  long used = 0;
  for (long t = 0; t < opt.steps; ++t) {
    if (opt.is_map) {
      if (!jac(x.data(), J.data())) { R.message = "lyapunov(jac): " + err; return R; }
      /* advance each column: v <- J v */
      std::vector<double> tmp(n);
      for (std::size_t c = 0; c < n; ++c) {
        matvec(n, J.data(), &Q[c * n], tmp.data());
        for (std::size_t i = 0; i < n; ++i) Q[c * n + i] = tmp[i];
      }
      if (!step_state(x)) { R.message = "lyapunov: " + err; return R; }
    } else {
      /* ODE: RK4 on the linear variational system, Jacobian frozen at x
       * for this step, then advance the base state. */
      if (!jac(x.data(), J.data())) { R.message = "lyapunov(jac): " + err; return R; }
      const double h = opt.dt;
      std::vector<double> k1(n), k2(n), k3(n), k4(n), tmp(n);
      for (std::size_t c = 0; c < n; ++c) {
        double *v = &Q[c * n];
        matvec(n, J.data(), v, k1.data());
        for (std::size_t i = 0; i < n; ++i) tmp[i] = v[i] + 0.5 * h * k1[i];
        matvec(n, J.data(), tmp.data(), k2.data());
        for (std::size_t i = 0; i < n; ++i) tmp[i] = v[i] + 0.5 * h * k2[i];
        matvec(n, J.data(), tmp.data(), k3.data());
        for (std::size_t i = 0; i < n; ++i) tmp[i] = v[i] + h * k3[i];
        matvec(n, J.data(), tmp.data(), k4.data());
        for (std::size_t i = 0; i < n; ++i)
          v[i] += (h / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
      }
      if (!step_state(x)) { R.message = "lyapunov: " + err; return R; }
    }

    bool fin = true;
    for (double v : x) fin = fin && std::isfinite(v);
    if (!fin) { R.message = "lyapunov: trajectory diverged"; return R; }

    if (((t + 1) % reorth) == 0 || t + 1 == opt.steps) {
      mgs(n, n, Q, r_diag);
      for (std::size_t c = 0; c < n; ++c)
        if (r_diag[c] > 1e-300) accum[c] += std::log(r_diag[c]);
      ++used;
    }
  }

  if (used == 0) { R.message = "lyapunov: no samples"; return R; }
  /* time spanned per accumulation block */
  const double block_time = opt.is_map ? static_cast<double>(reorth)
                                       : static_cast<double>(reorth) * opt.dt;
  const double total_time = block_time * static_cast<double>(used);

  R.exponents.resize(n);
  for (std::size_t c = 0; c < n; ++c) R.exponents[c] = accum[c] / total_time;
  std::sort(R.exponents.begin(), R.exponents.end(), std::greater<double>());
  R.sum = 0.0;
  for (double l : R.exponents) R.sum += l;
  R.kaplan_yorke = kaplan_yorke_dimension(R.exponents);
  R.steps_used = used;
  R.ok = true;
  R.message = "ok";
  return R;
}

/* ---- Homoclinic orbits: truncated BVP with projection BCs --------------- */
namespace {

/* Left eigenvectors of A split by the sign of the real part of the eigenvalue.
 * The rows of `stable_left` span the stable subspace of A^T (= orthogonal
 * complement of A's UNSTABLE subspace); rows of `unstable_left` span the
 * unstable subspace of A^T (= orthogonal complement of A's STABLE subspace).
 * Real vectors: a complex pair contributes its real and imaginary parts. Built
 * by inverse iteration on (A^T - lambda I) using the already-computed spectrum.
 * Returns false if the equilibrium is not hyperbolic (an eigenvalue on the
 * imaginary axis) since the projection BCs are then undefined. */
bool saddle_left_subspaces(const std::vector<double> &A, std::size_t n,
                           std::vector<std::vector<double>> *stable_left,
                           std::vector<std::vector<double>> *unstable_left,
                           int *n_unstable, std::string *err) {
  std::vector<Cplx> ev;
  if (!eigenvalues(A, n, &ev)) { if (err) *err = "eigenvalue failure at saddle"; return false; }
  /* A^T */
  std::vector<double> AT(n * n);
  for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) AT[i*n+j] = A[j*n+i];

  stable_left->clear(); unstable_left->clear();
  int nu = 0;
  std::vector<bool> done(ev.size(), false);
  for (std::size_t k = 0; k < ev.size(); ++k) {
    if (done[k]) continue;
    const double re = ev[k].real(), im = ev[k].imag();
    if (std::fabs(re) < 1e-7) { if (err) *err = "saddle is not hyperbolic (eigenvalue on imaginary axis)"; return false; }

    /* inverse iteration on (A^T - lambda I) for the left eigenvector */
    std::vector<Cplx> M(n * n);
    for (std::size_t i = 0; i < n; ++i)
      for (std::size_t j = 0; j < n; ++j)
        M[i*n+j] = Cplx(AT[i*n+j], 0.0) - ((i==j) ? Cplx(re, im) : Cplx(0,0));
    for (std::size_t i = 0; i < n; ++i) M[i*n+i] += Cplx(1e-9, 0); /* regularize */
    std::vector<Cplx> v(n, Cplx(1.0/std::sqrt((double)n), 0)); v[0] = Cplx(1,0);
    for (int it = 0; it < 120; ++it) {
      std::vector<Cplx> nv;
      if (!solve_complex(M, v, n, &nv)) { if (err) *err = "left eigenvector solve failed"; return false; }
      double nrm = 0; for (auto &z : nv) nrm += std::norm(z); nrm = std::sqrt(nrm);
      if (nrm < 1e-300) break;
      for (auto &z : nv) z /= nrm;
      v = nv;
    }
    auto &dst = (re > 0) ? *unstable_left : *stable_left;
    if (std::fabs(im) < 1e-9) {
      std::vector<double> rv(n); for (std::size_t i = 0; i < n; ++i) rv[i] = v[i].real();
      dst.push_back(rv);
      if (re > 0) nu += 1;
      done[k] = true;
    } else {
      /* complex pair -> real & imaginary parts span the 2-D real subspace */
      std::vector<double> rr(n), ri(n);
      for (std::size_t i = 0; i < n; ++i) { rr[i] = v[i].real(); ri[i] = v[i].imag(); }
      dst.push_back(rr); dst.push_back(ri);
      if (re > 0) nu += 2;
      done[k] = true;
      /* mark the conjugate as done */
      for (std::size_t j = k+1; j < ev.size(); ++j)
        if (!done[j] && std::fabs(ev[j].real() - re) < 1e-7 && std::fabs(ev[j].imag() + im) < 1e-7) { done[j] = true; break; }
    }
  }
  if (n_unstable) *n_unstable = nu;
  return true;
}

} /* anonymous namespace */

bool seed_homoclinic_by_integration(const Model &m, const std::vector<double> &saddle,
                                    double p, double dt, double max_time,
                                    std::vector<std::vector<double>> *seed_out,
                                    std::string *err) {
  const std::size_t n = m.n;
  if (saddle.size() < n || n < 2) { if (err) *err = "bad saddle for seeding"; return false; }
  if (dt <= 0) dt = 0.01;
  std::vector<double> A;
  if (!finite_diff_jacobian(m, saddle.data(), p, &A, err, 1e-7)) return false;
  std::vector<Cplx> ev; eigenvalues(A, n, &ev);
  int npos = 0; for (auto &z : ev) if (z.real() > 1e-7) ++npos;
  if (npos == 0 || npos == (int)n) { if (err) *err = "equilibrium is not a saddle"; return false; }
  std::vector<double> d(n, 0.0); d[0] = 1.0;
  for (int it = 0; it < 300; ++it) {
    std::vector<double> nd(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) for (std::size_t j = 0; j < n; ++j) nd[i] += A[i*n+j]*d[j];
    double nr = 0; for (double v : nd) nr += v*v; nr = std::sqrt(nr);
    if (nr < 1e-300) break;
    for (double &v : nd) v /= nr;
    d = nd;
  }
  const double eps = 1e-5;
  auto fld = [&](const std::vector<double> &xx, std::vector<double> &ff) {
    ff.assign(n, 0.0); std::string e; return m.vector_field(xx.data(), p, ff.data(), &e);
  };
  std::vector<std::vector<double>> best; double best_return = 1e300;
  for (int orient = 0; orient < 2; ++orient) {
    const double sgn = orient ? -1.0 : 1.0;
    std::vector<double> x(n); for (std::size_t i = 0; i < n; ++i) x[i] = saddle[i] + sgn*eps*d[i];
    std::vector<std::vector<double>> traj; traj.reserve(40000);
    std::vector<double> k1(n),k2(n),k3(n),k4(n),tmp(n);
    bool left = false; double closest_after_leaving = 1e300;
    const long maxsteps = (long)(max_time/dt);
    for (long s = 0; s < maxsteps; ++s) {
      traj.push_back(x);
      double dev = 0; for (std::size_t i = 0; i < n; ++i){ double dd=x[i]-saddle[i]; dev+=dd*dd; } dev = std::sqrt(dev);
      if (dev > 100.0*eps) left = true;
      if (left) closest_after_leaving = std::min(closest_after_leaving, dev);
      if (left && dev < 50.0*eps && traj.size() > 50) break;
      if (dev > 1e6 || !std::isfinite(dev)) break;
      if (!fld(x,k1)) break;
      for (std::size_t i=0;i<n;i++){ tmp[i]=x[i]+0.5*dt*k1[i]; }
      if(!fld(tmp,k2)) break;
      for (std::size_t i=0;i<n;i++){ tmp[i]=x[i]+0.5*dt*k2[i]; }
      if(!fld(tmp,k3)) break;
      for (std::size_t i=0;i<n;i++){ tmp[i]=x[i]+dt*k3[i]; }
      if(!fld(tmp,k4)) break;
      for (std::size_t i=0;i<n;i++){ x[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0; }
    }
    if (traj.size() > 20 && closest_after_leaving < best_return) { best_return = closest_after_leaving; best = traj; }
  }
  if (best.size() < 20) { if (err) *err = "no near-homoclinic excursion off the saddle"; return false; }
  *seed_out = best;
  return true;
}

HomoclinicResult solve_homoclinic(const Model &m,
                                  const std::vector<double> &x0_guess, double p,
                                  const std::vector<std::vector<double>> &seed_orbit,
                                  const HomoclinicSettings &settings) {
  HomoclinicResult R;
  const std::size_t n = m.n;
  const int M = std::max(8, settings.mesh);
  if (x0_guess.size() != n || seed_orbit.size() < 4) { R.message = "bad input to solve_homoclinic"; return R; }

  /* 1. Refine the saddle equilibrium and get its Jacobian + left subspaces. */
  std::vector<double> x0 = x0_guess; std::string err;
  if (!correct_equilibrium(m, &x0, p, 60, 1e-11, &err)) { R.message = "saddle did not converge: " + err; return R; }
  std::vector<double> A;
  if (!finite_diff_jacobian(m, x0.data(), p, &A, &err, 1e-7)) { R.message = "Jacobian failed at saddle"; return R; }
  std::vector<std::vector<double>> Ls, Lu; int nu = 0;
  if (!saddle_left_subspaces(A, n, &Ls, &Lu, &nu, &err)) { R.message = err; return R; }
  R.saddle = x0; R.n_unstable = nu;
  if (nu == 0 || nu == (int)n) { R.message = "equilibrium is a sink or source, not a saddle (no homoclinic)"; return R; }

  /* 2. Initial guess: resample the seed orbit onto mesh+1 points. The seed is
   *    assumed to be ordered along the orbit and sampled ~uniformly in time
   *    (e.g. from a fixed-step simulation), so we resample by INDEX fraction to
   *    match the collocation's uniform-in-tau grid. (Arclength resampling would
   *    misalign the phase condition's midpoint with the orbit's midpoint.) */
  const int Np = (int)seed_orbit.size();
  std::vector<double> U((size_t)(M+1)*n + 1, 0.0);   /* unknowns: orbit (M+1)*n, then T */
  auto resample = [&](double s, std::vector<double> &xo) {  /* s in [0,1] index fraction */
    double f = s * (Np - 1); int i = (int)std::floor(f); if (i > Np-2) i = Np-2; if (i < 0) i = 0;
    double t = f - i;
    xo.resize(n); for (std::size_t k = 0; k < n; ++k) xo[k] = (1-t)*seed_orbit[i][k] + t*seed_orbit[i+1][k];
  };
  std::vector<double> xtmp(n);
  for (int j = 0; j <= M; ++j) { resample((double)j/M, xtmp); for (std::size_t k = 0; k < n; ++k) U[(size_t)j*n+k] = xtmp[k]; }
  /* arclength of the seed, for a default time scale only */
  double L = 0.0; for (int i = 1; i < Np; ++i) { double d=0; for (std::size_t k=0;k<n;k++){double dd=seed_orbit[i][k]-seed_orbit[i-1][k]; d+=dd*dd;} L += std::sqrt(d); }
  double Tguess = settings.T;
  if (Tguess <= 0) Tguess = std::max(1.0, 0.5 * L);
  U[(size_t)(M+1)*n] = Tguess;

  /* phase condition (Poincare section): the orbit midpoint must lie on the
   * hyperplane through the seed midpoint normal to the seed velocity there,
   *   < x(tau_mid) - x_seed_mid , f(x_seed_mid) > = 0.
   * This pins the time-translation freedom without forcing a coordinate value
   * and converges fast for the single-orbit solve. (NOTE: a pure section also
   * admits the trivial collapse onto the saddle, which is why the two-parameter
   * continuation below is still experimental.) */
  std::vector<double> xmidseed(n); resample(0.5, xmidseed);
  std::vector<double> phase_normal(n);
  { std::string e; m.vector_field(xmidseed.data(), p, phase_normal.data(), &e); }
  std::vector<double> phase_ref = xmidseed;

  /* 3. Residual of the truncated BVP.
   *    F has M*n collocation eqs + n_s + n_u boundary eqs + 1 phase eq. */
  auto field = [&](const double *x, std::vector<double> &f) -> bool {
    f.resize(n); std::string e; return m.vector_field(x, p, f.data(), &e);
  };
  const int n_s = (int)Ls.size(), n_u = (int)Lu.size();
  const int NF = M*(int)n + n_s + n_u + 1;
  const int NV = (M+1)*(int)n + 1;
  auto residual = [&](const std::vector<double> &Uv, std::vector<double> *Fo) -> bool {
    Fo->assign(NF, 0.0);
    const double T = Uv[(size_t)(M+1)*n];
    std::vector<double> fi(n), fj(n), xi(n), xj(n);
    /* trapezoidal collocation of x' = 2T f(x) on tau in [0,1], M intervals */
    for (int i = 0; i < M; ++i) {
      for (std::size_t k = 0; k < n; ++k) { xi[k]=Uv[(size_t)i*n+k]; xj[k]=Uv[(size_t)(i+1)*n+k]; }
      if (!field(xi.data(), fi)) return false;
      if (!field(xj.data(), fj)) return false;
      const double h = 1.0/(double)M;
      for (std::size_t k = 0; k < n; ++k)
        (*Fo)[(size_t)i*n+k] = xj[k]-xi[k] - h*T*(fi[k]+fj[k]); /* 2T*0.5*(fi+fj) = T*(fi+fj) */
    }
    int row = M*(int)n;
    /* left BC: (x(0)-x0) orthogonal to stable left-eigvecs => in unstable subsp */
    for (int r = 0; r < n_s; ++r) {
      double s = 0; for (std::size_t k = 0; k < n; ++k) s += (Uv[k]-x0[k])*Ls[r][k];
      (*Fo)[row++] = s;
    }
    /* right BC: (x(1)-x0) orthogonal to unstable left-eigvecs => in stable subsp */
    for (int r = 0; r < n_u; ++r) {
      double s = 0; for (std::size_t k = 0; k < n; ++k) s += (Uv[(size_t)M*n+k]-x0[k])*Lu[r][k];
      (*Fo)[row++] = s;
    }
    /* phase: orbit midpoint on the seed's Poincare section */
    { double s = 0; for (std::size_t k = 0; k < n; ++k) s += (Uv[(size_t)(M/2)*n+k]-phase_ref[k])*phase_normal[k];
      (*Fo)[row++] = s; }
    return true;
  };

  /* 4. Gauss-Newton with finite-difference Jacobian (NF x NV, NF<=NV). Solve
   *    the normal equations (J^T J) dU = -J^T F with a little regularization
   *    (the system is rectangular: one fewer equation than unknowns, the extra
   *    freedom is the orbit's overall position which the phase pins softly). */
  std::vector<double> F(NF), Uv = U;
  double resn = 0;
  int step = 0;
  for (; step < settings.newton_iters; ++step) {
    if (!residual(Uv, &F)) { R.message = "field eval failed during Newton"; return R; }
    resn = 0; for (double f : F) resn += f*f; resn = std::sqrt(resn);
    if (resn < settings.newton_tol) break;

    /* finite-difference Jacobian J (NF x NV) */
    std::vector<double> J((size_t)NF*NV, 0.0), F2(NF);
    std::vector<double> Up = Uv;
    for (int c = 0; c < NV; ++c) {
      const double save = Up[c];
      const double dh = 1e-7 * (std::fabs(save) + 1e-3);
      Up[c] = save + dh;
      if (!residual(Up, &F2)) { Up[c] = save; continue; }
      Up[c] = save;
      for (int r = 0; r < NF; ++r) J[(size_t)r*NV + c] = (F2[r]-F[r])/dh;
    }
    if (!settings.free_T) {
      /* freeze T: zero its whole column so the normal equations give dT ~ 0 */
      for (int r = 0; r < NF; ++r) J[(size_t)r*NV + (NV-1)] = 0.0;
    }
    /* normal equations: (J^T J + mu I) dU = -J^T F */
    std::vector<double> JTJ((size_t)NV*NV, 0.0), JTF(NV, 0.0);
    for (int a = 0; a < NV; ++a) {
      double g = 0; for (int r = 0; r < NF; ++r) g += J[(size_t)r*NV+a]*F[r];
      JTF[a] = g;
      for (int b = 0; b < NV; ++b) {
        double s = 0; for (int r = 0; r < NF; ++r) s += J[(size_t)r*NV+a]*J[(size_t)r*NV+b];
        JTJ[(size_t)a*NV+b] = s;
      }
    }
    const double mu = 1e-9;
    for (int a = 0; a < NV; ++a) JTJ[(size_t)a*NV+a] += mu;
    std::vector<double> dU, rhs(NV);
    for (int a = 0; a < NV; ++a) rhs[a] = -JTF[a];
    if (!solve_linear(JTJ, rhs, &dU)) { R.message = "Newton linear solve failed"; break; }
    double dn = 0; for (double v : dU) dn += v*v; dn = std::sqrt(dn);
    double damp = dn > 0.5 ? 0.5/dn : 1.0;   /* damp big steps */
    for (int a = 0; a < NV; ++a) Uv[a] += damp*dU[a];
    if (Uv[(size_t)(M+1)*n] < 1e-3) Uv[(size_t)(M+1)*n] = 1e-3; /* keep T positive */
  }

  /* 5. Pack the result. */
  R.newton_residual = resn; R.newton_steps = step;
  R.T = Uv[(size_t)(M+1)*n];
  R.orbit.assign(M+1, std::vector<double>(n));
  R.tau.assign(M+1, 0.0);
  double amp = 0;
  for (int j = 0; j <= M; ++j) {
    R.tau[j] = (double)j/M;
    double dev = 0;
    for (std::size_t k = 0; k < n; ++k) { R.orbit[j][k] = Uv[(size_t)j*n+k]; double d = Uv[(size_t)j*n+k]-x0[k]; dev += d*d; }
    amp = std::max(amp, std::sqrt(dev));
  }
  R.amplitude = amp;
  R.ok = (resn < 1e-5);
  R.message = R.ok ? ("converged homoclinic orbit (residual " + std::to_string(resn) + ", T=" + std::to_string(R.T) + ")")
                   : ("did not fully converge (residual " + std::to_string(resn) + ")");
  return R;
}

HomoclinicCurve continue_homoclinic(const Model2 &m2,
                                    const std::vector<double> &x0_guess,
                                    double p0, double q0,
                                    const std::vector<std::vector<double>> &seed_orbit,
                                    const HomoclinicContSettings &settings) {
  HomoclinicCurve C;
  const std::size_t n = m2.n;
  if (x0_guess.size() != n || seed_orbit.size() < 4) { C.message = "bad input to continue_homoclinic"; return C; }

  /* wrap the two-parameter field as a one-parameter Model at a frozen q */
  auto model_at_q = [&](double q) {
    Model m; m.n = n;
    m.vector_field = [&m2, q](const double *x, double p, double *f, std::string *e) {
      return m2.vector_field(x, p, q, f, e);
    };
    return m;
  };

  /* At fixed q, find p where the BVP homoclinic closes (residual ~0). Because
   * the residual is a smooth, roughly-V-shaped function of p with its minimum
   * (~0) at the homoclinic, we use a short SECANT iteration on p driving the
   * residual down, re-using the current seed orbit. A handful of BVP solves per
   * q-step (vs a full scan) keeps continuation affordable. */
  auto solve_at_q = [&](double q, double p_guess,
                        const std::vector<std::vector<double>> &seed,
                        HomoclinicResult *out_best, double *out_p) -> bool {
    Model m = model_at_q(q);
    auto resid_at_p = [&](double p, HomoclinicResult *R) -> double {
      HomoclinicSettings bs = settings.bvp;
      *R = solve_homoclinic(m, x0_guess, p, seed, bs);
      if (R->orbit.empty()) return 1e9;
      return R->newton_residual;
    };
    /* Probe p_guess and two neighbours to estimate the local slope of the
     * residual, then take damped secant steps toward the minimum. */
    const double hp = std::max(1e-4, 0.25 * settings.dq);
    double bestp = p_guess, bestr = 1e18; HomoclinicResult bestR;
    double pa = p_guess;
    HomoclinicResult Ra; double ra = resid_at_p(pa, &Ra);
    if (ra < bestr) { bestr = ra; bestp = pa; bestR = Ra; }
    for (int it = 0; it < 8 && bestr > 1e-7; ++it) {
      /* numerical derivative of residual wrt p around the current best */
      HomoclinicResult Rp, Rm; double rp = resid_at_p(bestp + hp, &Rp), rm = resid_at_p(bestp - hp, &Rm);
      if (rp < bestr) { bestr = rp; bestp = bestp + hp; bestR = Rp; }
      if (rm < bestr) { bestr = rm; bestp = bestp - hp; bestR = Rm; }
      double g = (rp - rm) / (2 * hp);          /* d resid / d p   */
      double gg = (rp - 2 * ra + rm) / (hp * hp); /* curvature      */
      if (std::fabs(gg) < 1e-12) break;
      double dp = -g / gg;                       /* Newton step on the residual min */
      if (dp > 5 * settings.dq) dp = 5 * settings.dq;
      if (dp < -5 * settings.dq) dp = -5 * settings.dq;
      double pn = bestp + dp;
      HomoclinicResult Rn; double rn = resid_at_p(pn, &Rn);
      ra = rn; pa = pn;
      if (rn < bestr) { bestr = rn; bestp = pn; bestR = Rn; }
      else break; /* no improvement: stop */
    }
    *out_best = bestR; *out_p = bestp;
    /* Accept only a GENUINE loop: small residual AND a non-trivial amplitude
     * (the projection BCs + section also admit the orbit collapsing onto the
     * saddle, which must be rejected). */
    return bestr < 5e-4 && bestR.amplitude > 1e-3;
  };

  auto orbit_to_seed = [&](const HomoclinicResult &R) {
    return R.orbit;   /* the converged orbit reseeds the next step */
  };

  /* record the starting point first */
  {
    HomoclinicResult R0; double pfit = p0;
    if (solve_at_q(q0, p0, seed_orbit, &R0, &pfit)) {
      HomoclinicCurvePoint pt; pt.p = pfit; pt.q = q0; pt.T = R0.T; pt.amplitude = R0.amplitude; pt.residual = R0.newton_residual;
      C.points.push_back(pt);
      if (settings.store_orbits) C.orbits.push_back(R0.orbit);
    } else {
      C.message = "could not establish a homoclinic at the starting (p,q)";
      return C;
    }
  }

  /* trace in +/- q directions with a TANGENT PREDICTOR: after two accepted
   * points we linearly extrapolate p along the curve (p grows/curves with q),
   * giving the corrector a far better starting guess than reusing the last p.
   * The converged orbit reseeds the next step. This makes continuation hold
   * over a wide q range and through gentle folds of the locus. */
  for (int dir = 0; dir < (settings.both_directions ? 2 : 1); ++dir) {
    const double sgn = (dir == 0) ? 1.0 : -1.0;
    double q = q0, p = C.points.front().p;
    std::vector<std::vector<double>> seed = seed_orbit;
    HomoclinicResult Rprev;
    { double pf; if (solve_at_q(q0, p, seed, &Rprev, &pf)) { p = pf; seed = orbit_to_seed(Rprev); } }
    double prev_q = q, prev_p = p;           /* for the tangent predictor      */
    bool have_prev = false;
    int consec_fail = 0;
    for (int s = 0; s < settings.max_steps; ++s) {
      const double qn = q + sgn * settings.dq;
      if (qn < settings.q_min || qn > settings.q_max) break;
      /* predictor: extrapolate p along the curve from the last two points */
      double p_pred = p;
      if (have_prev && std::fabs(q - prev_q) > 1e-12)
        p_pred = p + (p - prev_p) * (sgn * settings.dq) / (q - prev_q);
      HomoclinicResult R; double pfit = p_pred;
      if (!solve_at_q(qn, p_pred, seed, &R, &pfit)) {
        /* one retry from the un-extrapolated p before giving up */
        if (!solve_at_q(qn, p, seed, &R, &pfit)) {
          if (++consec_fail >= 2) break;     /* lost the curve */
          q = qn; continue;
        }
      }
      consec_fail = 0;
      HomoclinicCurvePoint pt; pt.p = pfit; pt.q = qn; pt.T = R.T; pt.amplitude = R.amplitude; pt.residual = R.newton_residual;
      C.points.push_back(pt);
      if (settings.store_orbits) C.orbits.push_back(R.orbit);
      prev_q = q; prev_p = p; have_prev = true;
      q = qn; p = pfit; seed = orbit_to_seed(R);
    }
  }
  C.ok = C.points.size() > 1;
  C.message = C.ok ? ("traced " + std::to_string(C.points.size()) + " homoclinic points")
                   : "could not trace the homoclinic curve";
  return C;
}

LinResult lin_homoclinic(const Model2 &m2, const std::vector<double> &saddle_guess,
                         double p0, double q, const LinSettings &settings) {
  LinResult R;
  const std::size_t n = m2.n;
  if (saddle_guess.size() < n || n < 2) { R.message = "bad input to lin_homoclinic"; return R; }
  const double dt = settings.dt > 0 ? settings.dt : 0.01;

  /* one-parameter field at fixed q */
  auto field = [&](const double *x, double p, std::vector<double> &f) -> bool {
    f.assign(n, 0.0); std::string e; return m2.vector_field(x, p, q, f.data(), &e);
  };
  auto fdjac = [&](const std::vector<double> &x, double p, std::vector<double> &J) -> bool {
    J.assign(n*n, 0.0); std::vector<double> fp(n), fm(n); std::vector<double> xt = x;
    for (std::size_t j = 0; j < n; ++j) {
      double h = 1e-7*(std::fabs(x[j])+1e-3);
      xt[j]=x[j]+h; if(!field(xt.data(),p,fp)) return false;
      xt[j]=x[j]-h; if(!field(xt.data(),p,fm)) return false;
      xt[j]=x[j];
      for (std::size_t i=0;i<n;i++) J[i*n+j]=(fp[i]-fm[i])/(2*h);
    }
    return true;
  };
  /* Newton-correct the saddle at parameter p, starting from x0 */
  auto correct = [&](std::vector<double> x0, double p, std::vector<double> *out) -> bool {
    std::vector<double> x = x0, f(n), J, d;
    for (int it=0; it<80; ++it) {
      if (!field(x.data(), p, f)) return false;
      double fn=0; for(double v:f) fn+=v*v; fn=std::sqrt(fn);
      if (fn < 1e-12) { *out=x; return true; }
      if (!fdjac(x, p, J)) return false;
      std::vector<double> rhs(n); for(std::size_t i=0;i<n;i++) rhs[i]=-f[i];
      if (!solve_linear(J, rhs, &d)) return false;
      double dn=0; for(double v:d) dn+=v*v; dn=std::sqrt(dn);
      double damp = dn>1.0?1.0/dn:1.0;
      for(std::size_t i=0;i<n;i++) x[i]+=damp*d[i];
    }
    std::vector<double> f2; if(!field(x.data(),p,f2)) return false;
    double fn=0; for(double v:f2) fn+=v*v; fn=std::sqrt(fn);
    if (fn<1e-8){ *out=x; return true; }
    return false;
  };
  /* dominant real eigenvector of M with eigenvalue of a given sign, via power
   * iteration on (M - shift I) suitably; here we use: unstable = power iterate
   * M; stable = power iterate -M (largest of -M = most negative of M). */
  auto dominant_dir = [&](const std::vector<double> &M, bool unstable, std::vector<double> *d) -> bool {
    std::vector<double> v(n, 0.0); v[0]=1.0;
    const double s = unstable ? 1.0 : -1.0;
    for (int it=0; it<400; ++it) {
      std::vector<double> nv(n,0.0);
      for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++) nv[i]+=s*M[i*n+j]*v[j];
      /* shift to keep dominant positive: add large multiple of v */
      for (std::size_t i=0;i<n;i++) nv[i]+=5.0*v[i];
      double nr=0; for(double z:nv) nr+=z*z; nr=std::sqrt(nr);
      if (nr<1e-300) return false;
      for(double &z:nv) z/=nr;
      v=nv;
    }
    *d=v; return true;
  };

  /* Build a Poincare section from a reference excursion: integrate the unstable
   * manifold forward a while to find a "far" point, use the field direction
   * there as the section normal and that point as the section anchor. */
  auto integrate = [&](std::vector<double> x, double p, double sgn_time, double tmax,
                       const std::vector<double> *sect_anchor, const std::vector<double> *sect_normal,
                       std::vector<double> *hit, std::vector<std::vector<double>> *traj) -> bool {
    std::vector<double> k1(n),k2(n),k3(n),k4(n),tmp(n);
    const double h = sgn_time*dt;
    double prev_sd = 0; bool have_prev=false;
    std::vector<double> xprev = x;
    const std::vector<double> start = x;
    long steps=(long)(tmax/dt);
    for (long s=0;s<steps;s++) {
      if (traj) traj->push_back(x);
      if (sect_anchor && sect_normal) {
        double sd=0; for(std::size_t i=0;i<n;i++) sd+=(x[i]-(*sect_anchor)[i])*(*sect_normal)[i];
        /* require the orbit to have travelled away from its start before a
         * crossing counts (otherwise the very first steps, still near the
         * saddle, can graze the section). */
        double moved=0; for(std::size_t i=0;i<n;i++){double dd=x[i]-start[i];moved+=dd*dd;} moved=std::sqrt(moved);
        if (have_prev && prev_sd*sd<0.0 && moved>1e-2) {  /* genuine crossing */
          double t = prev_sd/(prev_sd-sd);
          hit->assign(n,0.0);
          for(std::size_t i=0;i<n;i++) (*hit)[i]= xprev[i] + t*(x[i]-xprev[i]); /* interpolate */
          return true;
        }
        prev_sd=sd; have_prev=true;
      }
      xprev = x;
      if(!field(x.data(),p,k1)) return false;
      for(std::size_t i=0;i<n;i++){tmp[i]=x[i]+0.5*h*k1[i];} if(!field(tmp.data(),p,k2))return false;
      for(std::size_t i=0;i<n;i++){tmp[i]=x[i]+0.5*h*k2[i];} if(!field(tmp.data(),p,k3))return false;
      for(std::size_t i=0;i<n;i++){tmp[i]=x[i]+h*k3[i];}     if(!field(tmp.data(),p,k4))return false;
      double dev=0; for(std::size_t i=0;i<n;i++){ x[i]+=h*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6.0; dev+=x[i]*x[i]; }
      if (!std::isfinite(dev) || dev>1e12) return false;
    }
    return false; /* never hit the section */
  };

  /* The Lin gap as a function of p: relocate saddle, get eigvecs, build the
   * section from the unstable excursion's far point, then measure (x+ - x-)
   * projected onto an in-section direction. */
  /* The homoclinic TEST FUNCTION g(p): integrate the unstable manifold forward;
   * the orbit makes a large excursion and returns toward the saddle. At the
   * point of CLOSEST return, measure the SIGNED distance along the saddle's
   * STABLE left-eigenvector (the covector normal to the stable manifold). On a
   * true homoclinic the return lands exactly on the stable manifold, so g=0;
   * off it, g has a definite sign that flips as p crosses the homoclinic value.
   * This is far more robust than a two-sided section gap because it only needs
   * the (well-behaved) forward unstable manifold to come back near the saddle.
   *
   * Orientation of the unstable direction is chosen as the one that makes a
   * BOUNDED returning excursion (the other side may escape to infinity). */
  auto homoclinic_test = [&](double p, double *gval, std::vector<std::vector<double>> *orbit_out) -> bool {
    std::vector<double> x0;
    if (!correct(saddle_guess, p, &x0)) return false;
    std::vector<double> J; if (!fdjac(x0, p, J)) return false;
    std::vector<Cplx> ev; eigenvalues(J, n, &ev);
    int npos=0,nneg=0; for(auto&z:ev){ if(z.real()>1e-7)npos++; if(z.real()<-1e-7)nneg++; }
    if (npos==0||nneg==0) return false;          /* not a saddle */

    /* unstable right eigvec of J; stable LEFT eigvec = unstable right eigvec of
     * J^T for the negative eigenvalue. dominant_dir(M,true) -> most positive;
     * (M,false) -> most negative. For the stable-manifold normal we want the
     * left eigenvector for the stable eigenvalue: dominant_dir(J^T, false). */
    std::vector<double> du, wst;
    if (!dominant_dir(J, true, &du)) return false;
    std::vector<double> JT(n*n);
    for (std::size_t i=0;i<n;i++) for (std::size_t j=0;j<n;j++) JT[i*n+j]=J[j*n+i];
    if (!dominant_dir(JT, false, &wst)) return false;   /* stable left-eigvec (covector) */

    /* pick the unstable orientation that returns (bounded) and get its orbit */
    std::vector<std::vector<double>> best_traj; double best_far=-1; double best_sign=1;
    for (int orient=0; orient<2; ++orient) {
      const double so = orient? -1.0 : 1.0;
      std::vector<double> xs(n); for(std::size_t i=0;i<n;i++) xs[i]=x0[i]+so*settings.eps*du[i];
      std::vector<std::vector<double>> tj; std::vector<double> dummy;
      integrate(xs, p, +1.0, settings.max_time, nullptr, nullptr, &dummy, &tj);
      if (tj.size()<20) continue;
      double fd=0; for(auto&pt:tj){double d=0;for(std::size_t i=0;i<n;i++){double dd=pt[i]-x0[i];d+=dd*dd;}fd=std::max(fd,std::sqrt(d));}
      if (fd>1e4 || !std::isfinite(fd)) continue;       /* escaped */
      if (fd>best_far){ best_far=fd; best_traj=tj; best_sign=so; }
    }
    if (best_traj.size()<20) return false;

    /* locate the far point (max distance) first; the homoclinic return is the
     * closest approach to the saddle AFTER that excursion peak. */
    std::size_t fi=0; double far_d=0;
    for (std::size_t s=0;s<best_traj.size();++s){ double d=0; for(std::size_t i=0;i<n;i++){double dd=best_traj[s][i]-x0[i];d+=dd*dd;} d=std::sqrt(d); if(d>far_d){far_d=d;fi=s;} }
    double closest=1e300; std::size_t ci=fi;
    for (std::size_t s=fi;s<best_traj.size();++s){ double d=0; for(std::size_t i=0;i<n;i++){double dd=best_traj[s][i]-x0[i];d+=dd*dd;} d=std::sqrt(d); if(d<closest){closest=d;ci=s;} }
    /* The homoclinic test value is the CLOSEST-RETURN DISTANCE to the saddle,
     * normalized by the loop size: on a true homoclinic the unstable manifold
     * comes back to the saddle, so this ratio -> 0. This is sign-free and
     * directly measures reconnection, avoiding the sign/amplitude confounds of
     * a projected gap. We still also report the signed projection in `gval` via
     * its sign for diagnostics, but the magnitude used for locating is
     * closest/far_d. */
    double signed_proj=0; for(std::size_t i=0;i<n;i++) signed_proj += (best_traj[ci][i]-x0[i])*wst[i];
    *gval = signed_proj;   /* signed return distance along the stable covector */
    if (orbit_out) { *orbit_out = best_traj; orbit_out->resize(ci+1); }
    (void)best_sign; (void)far_d;
    return true;
  };

  /* 1-D root-find on p (secant with bracketing) to drive the test function to
   * zero -- that is the homoclinic. */
  auto lin_gap = [&](double p, double *gap, std::vector<std::vector<double>> *orbit_out) -> bool {
    return homoclinic_test(p, gap, orbit_out);
  };
  double pa = p0, ga;
  if (!lin_gap(pa, &ga, nullptr)) { R.message = "homoclinic test could not be evaluated at p0 (no saddle / no return)"; return R; }
  double pb = p0 + 0.05*(std::fabs(p0)+1.0);
  double gb;
  if (!lin_gap(pb, &gb, nullptr)) { pb = p0 - 0.05*(std::fabs(p0)+1.0); if(!lin_gap(pb,&gb,nullptr)){ R.message="test eval failed near p0"; return R; } }

  /* The signed return distance |g(p)| dips toward 0 near the homoclinic (the
   * unstable manifold returns closest to the stable manifold). We locate that
   * dip by minimizing |g| over the bracket. NOTE: for stiff Bogdanov-Takens
   * loops the dip does not reach machine zero (the leading-order seed and the
   * one-sided manifold return leave a residual), so this is reported as a
   * located-approximately result with the residual, not a converged root. */
  (void)ga; (void)gb;
  double plo = std::min(pa, pb), phi = std::max(pa, pb);
  { double w = phi - plo; plo -= 1.5*w; phi += 1.5*w; }
  if (plo < settings.p_lo) plo = settings.p_lo;
  if (phi > settings.p_hi) phi = settings.p_hi;

  auto absg = [&](double p)->double{ double g; if(!lin_gap(p,&g,nullptr)) return 1e9; return std::fabs(g); };

  const int NS = 40;
  double best_p = plo; double best_r = 1e300;
  for (int i=0;i<=NS;i++){ double p = plo + (phi-plo)*i/NS; double r = absg(p); if (r<best_r){best_r=r;best_p=p;} }

  double step = (phi-plo)/NS;
  double a = best_p - step, b = best_p + step;
  if (a<settings.p_lo) a=settings.p_lo;
  if (b>settings.p_hi) b=settings.p_hi;
  const double gr = 0.6180339887;
  double c = b - gr*(b-a), d = a + gr*(b-a);
  double fc=absg(c), fd=absg(d);
  int it=0;
  for (; it<settings.max_iter; ++it) {
    if (std::fabs(b-a) < 1e-6) break;
    if (fc < fd) { b=d; d=c; fd=fc; c=b-gr*(b-a); fc=absg(c); }
    else         { a=c; c=d; fc=fd; d=a+gr*(b-a); fd=absg(d); }
  }
  double p = (fc<fd)? c : d;
  double g; lin_gap(p, &g, nullptr);

  std::vector<std::vector<double>> orbit;
  std::vector<double> x0fin;
  lin_gap(p, &g, &orbit);
  correct(saddle_guess, p, &x0fin);
  R.p = p; R.gap = g; R.iterations = it;
  R.saddle = x0fin; R.orbit = orbit;
  double amp=0; for (auto &pt : orbit){ double d=0; for(std::size_t i=0;i<n;i++){double dd=pt[i]-(x0fin.empty()?0:x0fin[i]);d+=dd*dd;} amp=std::max(amp,std::sqrt(d)); }
  R.amplitude = amp;
  /* `gap` is the signed return distance along the stable covector at the |gap|
   * minimum. A genuine homoclinic needs BOTH a substantial loop AND a small
   * return relative to that loop. Tiny-amplitude minima (near a p where no real
   * excursion forms) are rejected as spurious -- they make |gap| small only
   * because the orbit barely moves. */
  R.ok = (amp > 0.2) && (std::fabs(g) < 0.03 * amp);
  R.message = R.ok ? ("homoclinic located approximately (p=" + std::to_string(p) + ", |gap|=" + std::to_string(std::fabs(g)) + ", amp=" + std::to_string(amp) + ")")
                   : ("homoclinic not robustly located (best |gap|=" + std::to_string(std::fabs(g)) + ", amp=" + std::to_string(amp) + ") -- use solve_homoclinic/continue_homoclinic");
  return R;
}

/* ---- Basins of attraction ------------------------------------ */
BasinResult compute_basins(const std::function<bool(double, double, double *, double *)> &advance,
                           const BasinOptions &opt) {
  BasinResult R;
  const int W = std::max(2, opt.width), H = std::max(2, opt.height);
  R.width = W; R.height = H;
  R.cell_attractor.assign((size_t)W * H, -1);
  R.cell_speed.assign((size_t)W * H, 0.0f);
  const double R2 = opt.diverge_r * opt.diverge_r;
  const double clus2 = opt.cluster_tol * opt.cluster_tol;

  for (int j = 0; j < H; ++j) {
    const double y0 = opt.ymin + (opt.ymax - opt.ymin) * (double)j / (H - 1);
    for (int i = 0; i < W; ++i) {
      const double x0 = opt.xmin + (opt.xmax - opt.xmin) * (double)i / (W - 1);

      double x = x0, y = y0;
      long steps = 0;
      bool diverged = false, settled = false;
      double px = x, py = y;
      /* check "settled" on a short stride so periodic attractors (which
       * don't stop moving) still register once the orbit stops drifting
       * its average — we use slow movement over a stride as the signal. */
      const long stride = 8;
      double acc_move = 0.0;
      for (; steps < opt.max_steps; ++steps) {
        double nx = x, ny = y;
        if (!advance(x, y, &nx, &ny)) { diverged = true; break; }
        if (!std::isfinite(nx) || !std::isfinite(ny) || nx * nx + ny * ny > R2) { diverged = true; break; }
        acc_move += std::fabs(nx - x) + std::fabs(ny - y);
        x = nx; y = ny;
        if ((steps % stride) == (stride - 1)) {
          const double drift = std::fabs(x - px) + std::fabs(y - py);
          if (drift < opt.settle_tol) { settled = true; }
          px = x; py = y;
          if (settled) break;
        }
      }

      const size_t idx = (size_t)j * W + i;
      if (diverged) { R.cell_attractor[idx] = -1; R.cell_speed[idx] = 0.0f; ++R.n_diverged; continue; }
      if (!settled) {
        /* Orbit never stopped drifting within max_steps. For a chaotic
         * attractor (e.g. Lorenz) the endpoint is a different point of the
         * attractor each time, so clustering it would manufacture spurious
         * basins. Mark it as non-convergent (-2) instead — an honest,
         * coherent region rather than noise. */
        R.cell_attractor[idx] = -2; R.cell_speed[idx] = 0.0f; ++R.n_nonconvergent; continue;
      }

      /* match endpoint to an existing attractor cluster, else add one */
      int label = -1;
      for (size_t a = 0; a < R.attractors.size(); ++a) {
        const double dx = x - R.attractors[a].first, dy = y - R.attractors[a].second;
        if (dx * dx + dy * dy < clus2) { label = (int)a; break; }
      }
      if (label < 0) {
        if ((int)R.attractors.size() < opt.max_attractors) {
          R.attractors.push_back({x, y});
          label = (int)R.attractors.size() - 1;
        } else {
          /* over the cap: attach to nearest existing */
          double best = 1e300; int bi = -1;
          for (size_t a = 0; a < R.attractors.size(); ++a) {
            const double dx = x - R.attractors[a].first, dy = y - R.attractors[a].second;
            const double d = dx * dx + dy * dy;
            if (d < best) { best = d; bi = (int)a; }
          }
          label = bi;
        }
      }
      R.cell_attractor[idx] = label;
      ++R.n_converged;
      /* convergence speed: faster settle => brighter (1 at 0 steps, ->0 at max) */
      R.cell_speed[idx] = (float)(1.0 - (double)steps / (double)opt.max_steps);
    }
  }
  R.ok = true;
  R.message = "ok";
  return R;
}

BasinResult compute_basins_mt(const std::function<AdvanceFn(int tid)> &make_advance,
                              const BasinOptions &opt) {
  BasinResult R;
  const int W = std::max(2, opt.width), H = std::max(2, opt.height);
  R.width = W; R.height = H;
  R.cell_attractor.assign((size_t)W * H, -1);
  R.cell_speed.assign((size_t)W * H, 0.0f);
  const double R2 = opt.diverge_r * opt.diverge_r;
  const double clus2 = opt.cluster_tol * opt.cluster_tol;

  /* Phase 1 (parallel): integrate each cell independently, recording endpoint,
   * status and step count. Each worker thread owns a private advance (private
   * eval scratch) so there is no shared mutable state in the hot loop. */
  struct CellOut { double ex = 0, ey = 0; long steps = 0; int state = 1; }; /* 0 settled,1 diverged,2 nonconv */
  std::vector<CellOut> cells((size_t)W * H);

  unsigned hw = std::thread::hardware_concurrency();
  unsigned nthreads = (hw < 2 || H < 8) ? 1u : std::min<unsigned>(hw, (unsigned)H);

  auto do_rows = [&](int tid, int j0, int j1) {
    AdvanceFn advance = make_advance(tid);
    for (int j = j0; j < j1; ++j) {
      const double y0 = opt.ymin + (opt.ymax - opt.ymin) * (double)j / (H - 1);
      for (int i = 0; i < W; ++i) {
        const double x0 = opt.xmin + (opt.xmax - opt.xmin) * (double)i / (W - 1);
        double x = x0, y = y0;
        long steps = 0; bool diverged = false, settled = false;
        double px = x, py = y; const long stride = 8;
        for (; steps < opt.max_steps; ++steps) {
          double nx = x, ny = y;
          if (!advance(x, y, &nx, &ny)) { diverged = true; break; }
          if (!std::isfinite(nx) || !std::isfinite(ny) || nx * nx + ny * ny > R2) { diverged = true; break; }
          x = nx; y = ny;
          if ((steps % stride) == (stride - 1)) {
            const double drift = std::fabs(x - px) + std::fabs(y - py);
            if (drift < opt.settle_tol) settled = true;
            px = x; py = y;
            if (settled) break;
          }
        }
        CellOut &c = cells[(size_t)j * W + i];
        c.ex = x; c.ey = y; c.steps = steps;
        c.state = diverged ? 1 : (settled ? 0 : 2);
      }
    }
  };

  if (nthreads <= 1) {
    do_rows(0, 0, H);
  } else {
    std::vector<std::thread> pool;
    pool.reserve(nthreads);
    int per = (H + (int)nthreads - 1) / (int)nthreads;
    for (unsigned t = 0; t < nthreads; ++t) {
      int j0 = (int)t * per, j1 = std::min(H, j0 + per);
      if (j0 >= j1) break;
      pool.emplace_back(do_rows, (int)t, j0, j1);
    }
    for (auto &th : pool) th.join();
  }

  /* Phase 2 (serial): cluster endpoints into attractors with deterministic
   * labels (identical to the serial compute_basins). */
  for (int j = 0; j < H; ++j) {
    for (int i = 0; i < W; ++i) {
      const size_t idx = (size_t)j * W + i;
      const CellOut &c = cells[idx];
      if (c.state == 1) { R.cell_attractor[idx] = -1; ++R.n_diverged; continue; }
      if (c.state == 2) { R.cell_attractor[idx] = -2; ++R.n_nonconvergent; continue; }
      const double x = c.ex, y = c.ey;
      int label = -1;
      for (size_t a = 0; a < R.attractors.size(); ++a) {
        const double dx = x - R.attractors[a].first, dy = y - R.attractors[a].second;
        if (dx * dx + dy * dy < clus2) { label = (int)a; break; }
      }
      if (label < 0) {
        if ((int)R.attractors.size() < opt.max_attractors) {
          R.attractors.push_back({x, y});
          label = (int)R.attractors.size() - 1;
        } else {
          double best = 1e300; int bi = -1;
          for (size_t a = 0; a < R.attractors.size(); ++a) {
            const double dx = x - R.attractors[a].first, dy = y - R.attractors[a].second;
            const double d = dx * dx + dy * dy;
            if (d < best) { best = d; bi = (int)a; }
          }
          label = bi;
        }
      }
      R.cell_attractor[idx] = label;
      ++R.n_converged;
      R.cell_speed[idx] = (float)(1.0 - (double)c.steps / (double)opt.max_steps);
    }
  }
  R.ok = true; R.message = "ok";
  return R;
}
BoxCountResult box_counting_dimension(const std::vector<double> &xs,
                                      const std::vector<double> &ys,
                                      int n_levels) {
  BoxCountResult R;
  const size_t N = std::min(xs.size(), ys.size());
  if (N < 8) { R.message = "too few points"; return R; }
  if (n_levels < 3) n_levels = 3;
  if (n_levels > 24) n_levels = 24;

  double xmin = xs[0], xmax = xs[0], ymin = ys[0], ymax = ys[0];
  for (size_t i = 1; i < N; ++i) {
    if (std::isfinite(xs[i])) { xmin = std::min(xmin, xs[i]); xmax = std::max(xmax, xs[i]); }
    if (std::isfinite(ys[i])) { ymin = std::min(ymin, ys[i]); ymax = std::max(ymax, ys[i]); }
  }
  double span = std::max(xmax - xmin, ymax - ymin);
  if (!(span > 0.0)) { R.message = "degenerate point set (zero extent)"; return R; }
  /* pad slightly so boundary points fall inside the coarse box */
  span *= 1.0000001;

  for (int level = 0; level < n_levels; ++level) {
    const int grid = 1 << level;          /* 1, 2, 4, ... boxes per side */
    const double eps = span / grid;
    /* count occupied cells via a hash set of packed (ix,iy) */
    std::unordered_set<long long> occupied;
    occupied.reserve(N * 2);
    for (size_t i = 0; i < N; ++i) {
      if (!std::isfinite(xs[i]) || !std::isfinite(ys[i])) continue;
      int ix = (int)((xs[i] - xmin) / eps);
      int iy = (int)((ys[i] - ymin) / eps);
      if (ix < 0) ix = 0;
      if (ix >= grid) ix = grid - 1;
      if (iy < 0) iy = 0;
      if (iy >= grid) iy = grid - 1;
      occupied.insert(((long long)ix << 32) | (unsigned int)iy);
    }
    const double count = (double)occupied.size();
    if (count > 0) {
      R.log_inv_eps.push_back(std::log(1.0 / eps));
      R.log_count.push_back(std::log(count));
    }
  }

  /* Least-squares slope of log N vs log(1/eps). Use the middle of the
   * range: the coarsest levels (few boxes) and the finest (each point in
   * its own box -> N saturates at #points) bias the slope, so trim a
   * couple from each end when we have enough levels. */
  size_t lo = 0, hi = R.log_count.size();
  if (hi >= 7) { lo = 1; hi = R.log_count.size() - 2; }
  else if (hi >= 5) { lo = 1; hi = R.log_count.size() - 1; }
  const size_t m = (hi > lo) ? (hi - lo) : 0;
  if (m < 2) { R.message = "insufficient scales for a fit"; return R; }

  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (size_t i = lo; i < hi; ++i) {
    const double X = R.log_inv_eps[i], Y = R.log_count[i];
    sx += X; sy += Y; sxx += X * X; sxy += X * Y;
  }
  const double denom = m * sxx - sx * sx;
  if (std::fabs(denom) < 1e-300) { R.message = "degenerate fit"; return R; }
  const double slope = (m * sxy - sx * sy) / denom;
  const double intercept = (sy - slope * sx) / m;

  /* R^2 */
  const double ymean = sy / m;
  double ss_tot = 0, ss_res = 0;
  for (size_t i = lo; i < hi; ++i) {
    const double X = R.log_inv_eps[i], Y = R.log_count[i];
    const double pred = slope * X + intercept;
    ss_res += (Y - pred) * (Y - pred);
    ss_tot += (Y - ymean) * (Y - ymean);
  }
  R.dimension = slope;
  R.r_squared = (ss_tot > 0) ? (1.0 - ss_res / ss_tot) : 1.0;
  R.ok = true;
  R.message = "ok";
  return R;
}

/* ---- Iterated Function System (chaos game) -------------------- */
IFSResult chaos_game(const std::vector<AffineMap> &maps, long iterations,
                     unsigned int seed) {
  IFSResult R;
  if (maps.empty()) { R.message = "no maps"; return R; }
  if (iterations < 100) iterations = 100;

  /* build a cumulative probability table; if probabilities are all zero or
   * don't sum to ~1, fall back to uniform selection. */
  std::vector<double> cum(maps.size());
  double sum = 0.0;
  for (size_t i = 0; i < maps.size(); ++i) { sum += std::max(0.0, maps[i].p); cum[i] = sum; }
  const bool uniform = !(sum > 1e-12);

  /* simple deterministic LCG so results are reproducible and dependency-free */
  unsigned int state = seed ? seed : 1u;
  auto rnd = [&]() -> double {
    state = state * 1664525u + 1013904223u;
    return (double)(state >> 8) / (double)(1u << 24);
  };

  double x = 0.0, y = 0.0;
  const long burn = 50;
  R.xs.reserve((size_t)iterations);
  R.ys.reserve((size_t)iterations);
  double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
  for (long it = 0; it < iterations; ++it) {
    size_t k = 0;
    if (uniform) {
      k = (size_t)(rnd() * maps.size());
      if (k >= maps.size()) k = maps.size() - 1;
    } else {
      const double r = rnd() * sum;
      k = 0;
      while (k + 1 < maps.size() && r > cum[k]) ++k;
    }
    const AffineMap &m = maps[k];
    const double nx = m.a * x + m.b * y + m.e;
    const double ny = m.c * x + m.d * y + m.f;
    x = nx; y = ny;
    if (!std::isfinite(x) || !std::isfinite(y)) { x = 0; y = 0; continue; }
    if (it >= burn) {
      R.xs.push_back((float)x); R.ys.push_back((float)y);
      xmin = std::min(xmin, x); xmax = std::max(xmax, x);
      ymin = std::min(ymin, y); ymax = std::max(ymax, y);
    }
  }
  if (R.xs.empty()) { R.message = "no points"; return R; }
  R.xmin = xmin; R.xmax = xmax; R.ymin = ymin; R.ymax = ymax;
  R.ok = true; R.message = "ok";
  return R;
}

/* ---- Limit-cycle period & amplitude --------------------------- */
LimitCycleResult limit_cycle_period_amplitude(const std::vector<double> &y, double dt) {
  LimitCycleResult R;
  const size_t n = y.size();
  if (n < 16 || !(dt > 0)) { R.message = "signal too short"; return R; }

  /* mean and peak-to-peak amplitude over the (assumed settled) signal */
  double ymin = y[0], ymax = y[0], mean = 0.0;
  for (size_t i = 0; i < n; ++i) {
    if (!std::isfinite(y[i])) { R.message = "non-finite signal"; return R; }
    ymin = std::min(ymin, y[i]); ymax = std::max(ymax, y[i]); mean += y[i];
  }
  mean /= (double)n;
  const double ptp = ymax - ymin;
  R.amplitude = ptp;

  /* If the signal barely varies it's a fixed point, not a cycle. Scale the
   * threshold to the signal magnitude. */
  const double scale = std::max(1.0, std::fabs(ymax) + std::fabs(ymin));
  if (ptp < 1e-6 * scale) { R.message = "no oscillation (fixed point)"; return R; }

  /* Autocorrelation of the mean-removed signal; the first strong peak after
   * the zero-lag descent gives the period. */
  std::vector<double> c(y.begin(), y.end());
  for (double &v : c) v -= mean;
  const size_t maxlag = n / 2;
  double c0 = 0.0;
  for (size_t i = 0; i < n; ++i) c0 += c[i] * c[i];
  if (c0 <= 0) { R.message = "degenerate signal"; return R; }

  auto acf = [&](size_t lag) {
    double s = 0.0;
    for (size_t i = 0; i + lag < n; ++i) s += c[i] * c[i + lag];
    return s / c0;
  };

  /* find first lag where ACF rises back to a local max above a threshold,
   * after it has first dropped below zero (one full oscillation). */
  bool dropped = false;
  size_t best_lag = 0; double best_val = 0.0;
  for (size_t lag = 1; lag < maxlag; ++lag) {
    const double a = acf(lag);
    if (!dropped && a < 0.0) dropped = true;
    if (dropped) {
      const double am = acf(lag - 1), ap = (lag + 1 < maxlag) ? acf(lag + 1) : a;
      if (a > am && a >= ap && a > 0.2) { best_lag = lag; best_val = a; break; }
    }
    (void)best_val;
  }
  if (best_lag == 0) { R.message = "no clear period (aperiodic or too few cycles)"; return R; }

  /* refine the peak with a parabolic fit around best_lag */
  const double am = acf(best_lag - 1), a0 = acf(best_lag), ap = acf(best_lag + 1);
  double offset = 0.0;
  const double denom = (am - 2 * a0 + ap);
  if (std::fabs(denom) > 1e-12) offset = 0.5 * (am - ap) / denom;
  R.period = ((double)best_lag + offset) * dt;
  R.ok = true;
  R.message = "ok";
  return R;
}

}  // namespace dynsys::analysis