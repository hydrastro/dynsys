/* Regression guard for the codim-2 normal-form coefficients: the CUSP cubic
 * coefficient c and the generalized-Hopf SECOND Lyapunov coefficient l2.
 *
 * CUSP: on the 1-D normal form x' = c x^3 the cusp coefficient is exactly c.
 * GH (Bautin): on x' = -y + L x (x^2+y^2)^2, y' = x + L y (x^2+y^2)^2 (first
 *   Lyapunov coefficient = 0), the SECOND Lyapunov coefficient has the SIGN of
 *   L (magnitude carries a convention factor, so only the sign is asserted --
 *   the sign is what classifies the Bautin point).
 */
#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

static double Lg;

int main() {
  int fails = 0;

  {
    Model m; m.n = 1;
    m.vector_field = [](const double* X, double, double* o, std::string*) { o[0] = 2.0*X[0]*X[0]*X[0]; return true; };
    double cc = 0; std::string e;
    bool ok = cusp_normal_form(m, {0.0}, 0.0, &cc, &e);
    bool t = ok && std::fabs(cc - 2.0) < 1e-2;
    printf("  cusp c (x'=2x^3): c=%.4f (expect 2) %s\n", cc, t ? "PASS" : "FAIL");
    if (!t) fails++;
  }
  {
    Model m; m.n = 1;
    m.vector_field = [](const double* X, double, double* o, std::string*) { o[0] = -0.5*X[0]*X[0]*X[0]; return true; };
    double cc = 0; std::string e;
    bool ok = cusp_normal_form(m, {0.0}, 0.0, &cc, &e);
    bool t = ok && std::fabs(cc + 0.5) < 1e-2;
    printf("  cusp c (x'=-0.5x^3): c=%.4f (expect -0.5) %s\n", cc, t ? "PASS" : "FAIL");
    if (!t) fails++;
  }
  for (double L : {1.0, -1.0}) {
    Lg = L;
    Model m; m.n = 2;
    m.vector_field = [](const double* X, double, double* o, std::string*) {
      double x = X[0], y = X[1], r2 = x*x + y*y;
      o[0] = -y + Lg*x*r2*r2; o[1] = x + Lg*y*r2*r2; return true;
    };
    double l2 = 0; std::string e;
    bool ok = gh_second_lyapunov(m, {0.0, 0.0}, 0.0, &l2, &e);
    int sgot = (l2 > 0) - (l2 < 0), sexp = (L > 0) - (L < 0);
    bool t = ok && sgot == sexp && std::fabs(l2) > 1e-3;
    printf("  GH l2 (L=%+.0f): l2=%.4f sign=%+d (expect %+d) %s\n", L, l2, sgot, sexp, t ? "PASS" : "FAIL");
    if (!t) fails++;
  }

  if (fails) { printf("codim2_coeffs_smoke: %d FAIL\n", fails); return 1; }
  printf("codim2_coeffs_smoke: all checks pass\n");
  return 0;
}
