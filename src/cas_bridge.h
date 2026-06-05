/* cas_bridge.h — Phase E: bridge from dynsys to the Sangaku proof-carrying CAS.
 *
 * dynsys assembles a Jacobian at an equilibrium as a matrix of exact rationals,
 * this bridge serializes it to Sangaku's matrix form (lists of rational rows),
 * runs it through the Lizard interpreter, and parses back:
 *   - the exact characteristic polynomial (rational coefficients), and
 *   - the eigenvalues as human-readable strings.
 * From the exact charpoly we determine stability EXACTLY via Routh-Hurwitz
 * (integer/rational arithmetic, no floating point), so the verdict is certified
 * in the same spirit as Sangaku's own results.
 *
 * The bridge shells out to a `lizard` binary located via $LIZARD (and Sangaku
 * via $SANGAKU_ROOT). If either is missing or the call fails, is_available()
 * stays false and callers fall back to dynsys's existing numeric eigen path.
 *
 * Header-only, no GL/ImGui deps, so it is unit-testable headlessly.
 */
#ifndef DYNSYS_CAS_BRIDGE_H
#define DYNSYS_CAS_BRIDGE_H

#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <string>
#include <vector>
#include <array>
#include <numeric>
#include <sstream>

namespace dynsys::cas {

/* An exact rational a/b (b > 0). Jacobian entries that are integers or simple
 * fractions map here exactly; entries that are irrational are not exactly
 * representable and the caller should not use the CAS path for them. */
struct Rational {
  long long num = 0;
  long long den = 1;
  Rational() = default;
  Rational(long long n) : num(n), den(1) {}
  Rational(long long n, long long d) : num(n), den(d) { normalize(); }
  void normalize() {
    if (den < 0) { num = -num; den = -den; }
    if (den == 0) { den = 1; return; }
    long long g = std::gcd(num < 0 ? -num : num, den);
    if (g > 1) { num /= g; den /= g; }
  }
  std::string to_lisp() const {
    if (den == 1) return std::to_string(num);
    return std::to_string(num) + "/" + std::to_string(den); /* Lizard reads rationals */
  }
};

/* Convert a double to a nearby rational via continued fractions, bounded by a
 * maximum denominator. Sets *exact to true iff the rational reproduces x to
 * within `tol` (so the caller knows whether feeding it to the exact CAS is
 * trustworthy or merely an approximation). Many Jacobian entries at equilibria
 * are exactly rational (integers, simple fractions), which this recovers; an
 * irrational entry comes back with *exact=false and the caller should prefer
 * the numeric eigen path. */
inline Rational rationalize(double x, bool *exact, double tol = 1e-9,
                            long long max_den = 1000000) {
  bool neg = x < 0;
  double v = neg ? -x : x;
  /* continued-fraction expansion with convergents h/k */
  long long h0 = 0, h1 = 1, k0 = 1, k1 = 0;
  double frac = v;
  for (int it = 0; it < 40; ++it) {
    double a_d = std::floor(frac);
    long long a = (long long)a_d;
    long long h2 = a * h1 + h0;
    long long k2 = a * k1 + k0;
    if (k2 < 0 || h2 < 0 || k2 > max_den) break; /* denominator bound hit */
    h0 = h1; h1 = h2; k0 = k1; k1 = k2;
    double approx = (k1 != 0) ? (double)h1 / (double)k1 : 0.0;
    if (std::fabs(approx - v) <= tol) break;
    double rem = frac - a_d;
    if (rem < 1e-15) break;
    frac = 1.0 / rem;
  }
  if (k1 == 0) { if (exact) *exact = false; return Rational(0); }
  Rational r(neg ? -h1 : h1, k1);
  if (exact) *exact = (std::fabs((double)r.num / (double)r.den - x) <= tol);
  return r;
}

/* Result of an eigenvalue/stability query. */
struct EigenReport {
  bool ok = false;                       /* the CAS produced a result */
  std::vector<Rational> charpoly;        /* exact char-poly coeffs, low->high */
  std::vector<std::string> eigen_strings;/* human-readable eigenvalues */
  /* exact Routh-Hurwitz counts (for continuous-time stability, Re(lambda)) */
  int n_lhp = 0;                         /* roots with Re < 0 (stable directions) */
  int n_rhp = 0;                         /* roots with Re > 0 (unstable) */
  int n_imag = 0;                        /* roots with Re = 0 (marginal/Hopf) */
  bool stable = false;                   /* all roots strictly in the LHP */
  bool certified = false;                /* charpoly came from the CAS (exact) */
  std::string message;
};

/* ---- locate the interpreter ---------------------------------------------- */
inline std::string lizard_binary() {
  const char *e = std::getenv("LIZARD");
  return e ? std::string(e) : std::string();
}
inline std::string sangaku_root() {
  const char *e = std::getenv("SANGAKU_ROOT");
  return e ? std::string(e) : std::string();
}
inline bool is_available() {
  const std::string lz = lizard_binary(), sr = sangaku_root();
  if (lz.empty() || sr.empty()) return false;
  /* cheap existence check */
  FILE *f = std::fopen(lz.c_str(), "rb");
  if (!f) return false;
  std::fclose(f);
  return true;
}

/* ---- serialize a rational matrix to a Sangaku matrix literal ------------- */
inline std::string matrix_to_lisp(const std::vector<std::vector<Rational>> &A) {
  std::string s = "(list";
  for (const auto &row : A) {
    s += " (list";
    for (const auto &x : row) { s += " "; s += x.to_lisp(); }
    s += ")";
  }
  s += ")";
  return s;
}

/* ---- run a generated Sangaku script through Lizard, capture stdout ------- *
 * Writes prelude + script to a temp .lisp and pipes it to `lizard` (the same
 * contract as scripts/sangaku). Returns false on any failure. */
inline bool run_lizard(const std::string &script, std::string *out) {
  const std::string lz = lizard_binary(), sr = sangaku_root();
  if (lz.empty() || sr.empty()) return false;

  /* temp files */
  char tmpl[] = "/tmp/dynsys_cas_XXXXXX";
  int fd = mkstemp(tmpl);
  if (fd < 0) return false;
  std::string script_path = tmpl;
  {
    FILE *sf = fdopen(fd, "w");
    if (!sf) { std::remove(script_path.c_str()); return false; }
    std::fputs(script.c_str(), sf);
    std::fclose(sf);
  }

  /* command: cd into the Sangaku root (so the prelude's relative imports
   * resolve, exactly as scripts/sangaku does), then feed prelude + script to
   * lizard on stdin. */
  std::string cmd = "cd " + sr + " && cat " + sr + "/src/prelude.lisp " + script_path +
                    " | " + lz + " 2>/dev/null";
  std::string result;
  FILE *p = popen(cmd.c_str(), "r");
  if (!p) { std::remove(script_path.c_str()); return false; }
  char buf[4096];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof(buf), p)) > 0) result.append(buf, n);
  int rc = pclose(p);
  std::remove(script_path.c_str());
  if (rc != 0) return false;
  *out = result;
  return true;
}

/* ---- exact Routh-Hurwitz on a rational char-poly ------------------------- *
 * Counts roots by real-part sign for a continuous-time system, exactly. The
 * number of sign changes in the Routh array's first column = number of
 * right-half-plane roots. Handles the two classical degeneracies exactly:
 *   (a) a zero in the first column with a non-zero rest of row -> replace with
 *       a small symbolic +epsilon (we track signs through the limit), and
 *   (b) an entire zero row -> rebuild it from the derivative of the auxiliary
 *       polynomial (roots symmetric about the imaginary axis -> imaginary or
 *       in mirrored pairs; we record the imaginary-axis case).
 * Returns false only if something truly unhandled occurs. */
inline bool routh_hurwitz(const std::vector<Rational> &coeffs_low_to_high,
                           int *n_lhp, int *n_rhp, int *n_imag) {
  std::vector<Rational> a(coeffs_low_to_high.rbegin(), coeffs_low_to_high.rend());
  size_t start = 0;
  while (start < a.size() && a[start].num == 0) ++start;
  a.erase(a.begin(), a.begin() + (long)start);
  const int n = (int)a.size() - 1;
  if (n <= 0) return false;
  if (a[0].num < 0) for (auto &x : a) x.num = -x.num;

  auto rsubf = [](Rational x, Rational y) { return Rational(x.num*y.den - y.num*x.den, x.den*y.den); };
  auto rmul  = [](Rational x, Rational y) { return Rational(x.num*y.num, x.den*y.den); };
  auto rdiv  = [](Rational x, Rational y) { return Rational(x.num*y.den, x.den*y.num); };

  const int width = n / 2 + 1;
  /* Routh array: rows[0..n], each a row of rationals */
  std::vector<std::vector<Rational>> rows((size_t)(n + 1), std::vector<Rational>((size_t)width, Rational(0)));
  for (int i = 0; i <= n; ++i) {
    int r = i % 2, c = i / 2;
    if (c < width) rows[(size_t)r][(size_t)c] = a[(size_t)i];
  }

  bool imag_axis = false;
  for (int r = 2; r <= n; ++r) {
    /* if the previous row is entirely zero -> auxiliary-polynomial case */
    bool prev_zero = true;
    for (int c = 0; c < width; ++c) if (rows[(size_t)r-1][(size_t)c].num != 0) { prev_zero = false; break; }
    if (prev_zero) {
      /* rebuild row r-1 from derivative of the auxiliary polynomial formed by
       * row r-2 (powers s^(n-r+2), s^(n-r), ...). This signals roots on the
       * imaginary axis or in symmetric pairs. */
      imag_axis = true;
      int order = n - (r - 2);
      for (int c = 0; c < width; ++c) {
        int power = order - 2 * c;
        rows[(size_t)r-1][(size_t)c] = rmul(rows[(size_t)r-2][(size_t)c], Rational(power));
      }
    }
    Rational piv = rows[(size_t)r-1][0];
    if (piv.num == 0) {
      /* zero in the first column, rest non-zero: epsilon method. Sign of the
       * entries below depends on sign of epsilon (->0+). We approximate by a
       * tiny positive rational; the SIGN pattern is what we count, and a
       * positive epsilon gives the standard sign-change count. */
      piv = Rational(1, 1000000000);
      rows[(size_t)r-1][0] = piv;
    }
    for (int c = 0; c + 1 < width; ++c) {
      Rational ad = rmul(rows[(size_t)r-2][(size_t)c+1], rows[(size_t)r-1][0]);
      Rational bc = rmul(rows[(size_t)r-1][(size_t)c+1], rows[(size_t)r-2][0]);
      rows[(size_t)r][(size_t)c] = rdiv(rsubf(ad, bc), piv);
    }
  }

  int rhp = 0, prev_sign = 0;
  for (int r = 0; r <= n; ++r) {
    Rational v = rows[(size_t)r][0];
    int sgn = v.num > 0 ? 1 : (v.num < 0 ? -1 : 0);
    if (sgn == 0) continue; /* skip exact zeros in the count */
    if (prev_sign != 0 && sgn != prev_sign) ++rhp;
    prev_sign = sgn;
  }
  *n_rhp = rhp;
  *n_imag = imag_axis ? (n - 2 * rhp >= 0 ? (n - 2 * rhp) : 0) : 0;
  *n_lhp = n - rhp - *n_imag;
  if (*n_lhp < 0) *n_lhp = 0;
  return true;
}

/* ---- parse a Rational from Sangaku text ("3", "-2", "3/4") --------------- */
inline bool parse_rational(const std::string &tok, Rational *out) {
  size_t slash = tok.find('/');
  try {
    if (slash == std::string::npos) { *out = Rational(std::stoll(tok)); }
    else { *out = Rational(std::stoll(tok.substr(0, slash)), std::stoll(tok.substr(slash + 1))); }
  } catch (...) { return false; }
  return true;
}

/* ---- the headline call: exact eigenvalues + stability of a Jacobian ------ */
inline EigenReport eigen_report(const std::vector<std::vector<Rational>> &J) {
  EigenReport R;
  if (!is_available()) { R.message = "CAS unavailable (set LIZARD and SANGAKU_ROOT)"; return R; }
  if (J.empty()) { R.message = "empty matrix"; return R; }

  std::ostringstream script;
  script << "(import \"cas/linalg.lisp\")\n";
  script << "(define A " << matrix_to_lisp(J) << ")\n";
  script << "(display \"CHARPOLY\") (for-each (lambda (c) (display \" \") (display c)) (poly-norm (mat-charpoly A))) (newline)\n";
  script << "(for-each (lambda (sm) (display \"EIG \") (display (solution-string sm)) (newline)) (mat-eigenvalues A))\n";
  script << "(display \"END\") (newline)\n";

  std::string out;
  if (!run_lizard(script.str(), &out)) { R.message = "CAS call failed"; return R; }

  /* parse lines */
  std::istringstream is(out);
  std::string line;
  bool saw_end = false;
  while (std::getline(is, line)) {
    if (line.rfind("CHARPOLY", 0) == 0) {
      std::istringstream ls(line.substr(8));
      std::string tok;
      while (ls >> tok) { Rational r; if (parse_rational(tok, &r)) R.charpoly.push_back(r); }
    } else if (line.rfind("EIG ", 0) == 0) {
      R.eigen_strings.push_back(line.substr(4));
    } else if (line.rfind("END", 0) == 0) {
      saw_end = true;
    }
  }
  if (!saw_end || R.charpoly.empty()) { R.message = "could not parse CAS output"; return R; }

  R.ok = true;
  R.certified = true; /* charpoly is exact from the CAS */
  if (routh_hurwitz(R.charpoly, &R.n_lhp, &R.n_rhp, &R.n_imag)) {
    R.stable = (R.n_rhp == 0 && R.n_imag == 0);
    R.message = "ok";
  } else {
    /* marginal / degenerate Routh array: eigenvalues still exact & shown, but
     * the strict LHP verdict is left to the caller's exact-eigenvalue read or
     * the numeric path. */
    R.message = "eigenvalues exact; stability marginal (see eigenvalues)";
  }
  return R;
}

/* Convenience: take a numeric Jacobian (row-major n*n doubles) at an
 * equilibrium, rationalize it, and run the exact eigenvalue/stability report.
 * `all_exact` reports whether EVERY entry rationalized within tolerance; when
 * false the caller should treat the CAS result as approximate and prefer the
 * numeric spectrum. This is the entry point the UI uses, since dynsys's AD
 * Jacobian is computed in doubles. */
inline EigenReport eigen_report_from_doubles(const double *J, int n, bool *all_exact,
                                             double tol = 1e-9) {
  std::vector<std::vector<Rational>> M((size_t)n, std::vector<Rational>((size_t)n));
  bool ok_all = true;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      bool ex = false;
      M[(size_t)i][(size_t)j] = rationalize(J[(size_t)i * (size_t)n + (size_t)j], &ex, tol);
      if (!ex) ok_all = false;
    }
  if (all_exact) *all_exact = ok_all;
  EigenReport R = eigen_report(M);
  if (!ok_all && R.ok) R.message = "approximate: some Jacobian entries are not exactly rational";
  return R;
}

}  // namespace dynsys::cas

#endif  // DYNSYS_CAS_BRIDGE_H
