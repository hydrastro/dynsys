/* Regression guard for the manifold-based homoclinic seeding helper
 * (seed_homoclinic_by_integration) followed by the BVP solve.
 *
 * Unforced Duffing x'=y, y'=x-x^3 has a genuine homoclinic to the saddle at the
 * origin with peak amplitude sqrt(2). The helper nudges off the unstable
 * manifold, integrates until the orbit returns near the saddle, and the BVP
 * then converges to the connection. */
#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

static bool fld(const double* X, double, double* o, std::string*) { o[0] = X[1]; o[1] = X[0] - X[0]*X[0]*X[0]; return true; }

int main() {
  Model m; m.n = 2; m.vector_field = fld;
  std::vector<double> sad = {0, 0};
  std::vector<std::vector<double>> seed; std::string e;
  bool sok = seed_homoclinic_by_integration(m, sad, 0.0, 0.01, 100.0, &seed, &e);
  printf("  Duffing manifold seed: ok=%d points=%zu\n", sok, seed.size());
  if (!sok) { printf("homoclinic_seed_smoke: FAIL\n"); return 1; }
  HomoclinicSettings hs; hs.mesh = 200; hs.T = 0.5*0.01*(double)seed.size(); hs.free_T = false; hs.newton_iters = 80;
  HomoclinicResult R = solve_homoclinic(m, sad, 0.0, seed, hs);
  double peak = 0; for (auto& pt : R.orbit) peak = std::max(peak, std::fabs(pt[0]));
  double ed = R.orbit.empty() ? 1e9 : std::hypot(R.orbit.front()[0], R.orbit.front()[1]);
  printf("  solve: ok=%d steps=%d res=%.2e peak=%.4f (sqrt2=%.4f) end=%.2e\n",
         R.ok, R.newton_steps, R.newton_residual, peak, std::sqrt(2.0), ed);
  bool ok = R.ok && std::fabs(peak - std::sqrt(2.0)) < 0.05 && ed < 0.02;
  printf("%s\n", ok ? "homoclinic_seed_smoke: all checks pass" : "homoclinic_seed_smoke: FAIL");
  return ok ? 0 : 1;
}
