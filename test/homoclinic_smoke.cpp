/* Regression guard for the homoclinic-orbit solver (truncated BVP + projection
 * boundary conditions onto the saddle's stable/unstable subspaces).
 *
 * Two systems with KNOWN homoclinic orbits to a saddle at the origin:
 *   T1: x' = y, y' = x - x^2.   Exact homoclinic x(t) = (3/2) sech^2(t/2),
 *       peak x = 1.5 at t=0. Seeded from a CRUDE Gaussian bump (wrong shape and
 *       height) to test convergence, not just residual on the exact orbit.
 *   T2: x' = y, y' = 2x - 3x^2.  Peak at x* = 3a/(2b) = 1.0. Different
 *       eigenvalues (+-sqrt(2)), so this checks the projection BCs aren't
 *       specialised to the symmetric case.
 * Each must converge (small residual), recover the analytic peak, place the
 * endpoints near the saddle, and report one unstable direction.
 */
#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

static bool f1(const double* X, double, double* o, std::string*) { o[0]=X[1]; o[1]=X[0]-X[0]*X[0]; return true; }
static bool f2(const double* X, double, double* o, std::string*) { o[0]=X[1]; o[1]=2.0*X[0]-3.0*X[0]*X[0]; return true; }

int main() {
  int fails = 0;

  /* T1 */
  {
    Model m; m.n = 2; m.vector_field = f1;
    double Tt = 8.0; int Np = 300;
    std::vector<std::vector<double>> seed(Np, std::vector<double>(2));
    for (int i = 0; i < Np; ++i) {
      double t = -Tt + 2*Tt*i/(Np-1);
      double bump = 1.3*std::exp(-0.35*t*t);            /* deliberately wrong */
      seed[i][0] = bump; seed[i][1] = 1.3*(-0.7*t)*std::exp(-0.35*t*t);
    }
    HomoclinicSettings s; s.mesh = 200; s.T = Tt; s.newton_iters = 200; s.free_T = false;
    std::vector<double> x0 = {0,0};
    HomoclinicResult R = solve_homoclinic(m, x0, 0.0, seed, s);
    double peak = 0; for (auto &pt : R.orbit) peak = std::max(peak, pt[0]);
    double ed = R.orbit.empty() ? 1e9 : std::hypot(R.orbit.front()[0], R.orbit.front()[1]);
    bool ok = R.ok && std::fabs(peak - 1.5) < 0.05 && R.n_unstable == 1 && ed < 0.02;
    printf("  T1 sech^2 homoclinic: ok=%d steps=%d res=%.2e peak=%.4f (1.5) end=%.4f nu=%d %s\n",
           R.ok, R.newton_steps, R.newton_residual, peak, ed, R.n_unstable, ok ? "PASS" : "FAIL");
    if (!ok) fails++;
  }

  /* T2 */
  {
    Model m; m.n = 2; m.vector_field = f2;
    double lam = std::sqrt(2.0), xstar = 1.0, Tt = 10.0/lam; int Np = 300;
    std::vector<std::vector<double>> seed(Np, std::vector<double>(2));
    for (int i = 0; i < Np; ++i) {
      double t = -Tt + 2*Tt*i/(Np-1);
      double bump = xstar*std::exp(-0.5*lam*lam*t*t);
      seed[i][0] = bump; seed[i][1] = -lam*lam*t*bump;
    }
    HomoclinicSettings s; s.mesh = 200; s.T = Tt; s.newton_iters = 200; s.free_T = false;
    std::vector<double> x0 = {0,0};
    HomoclinicResult R = solve_homoclinic(m, x0, 0.0, seed, s);
    double peak = 0; for (auto &pt : R.orbit) peak = std::max(peak, pt[0]);
    double ed = R.orbit.empty() ? 1e9 : std::hypot(R.orbit.front()[0], R.orbit.front()[1]);
    bool ok = R.ok && std::fabs(peak - xstar) < 0.05 && R.n_unstable == 1 && ed < 0.02;
    printf("  T2 asymmetric homoclinic: ok=%d steps=%d res=%.2e peak=%.4f (1.0) end=%.4f nu=%d %s\n",
           R.ok, R.newton_steps, R.newton_residual, peak, ed, R.n_unstable, ok ? "PASS" : "FAIL");
    if (!ok) fails++;
  }

  if (fails) { printf("homoclinic_smoke: %d FAIL\n", fails); return 1; }
  printf("homoclinic_smoke: all checks pass\n");
  return 0;
}
