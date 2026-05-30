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

#include <algorithm>
#include <cmath>

namespace dynsys::analysis {

namespace {

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
  if (m.jacobian_x) return m.jacobian_x(x, p, jac->data(), err);
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
    branch.message = "could not converge starting equilibrium: " +
                     (err.empty() ? std::string("residual too large") : err);
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
    branch.points.push_back(std::move(bp));
    if (kind != SpecialPointKind::None)
      branch.special_indices.push_back(branch.points.size() - 1);
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

    if (settings.detect_fold) {
      const double f_new = fold_test(m, x_new.data(), p_new);
      if (std::isfinite(prev_fold) && std::isfinite(f_new) &&
          prev_fold * f_new < 0.0) {
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

}  // namespace dynsys::analysis
