/* Regression guard for equilibrium branch-point detection + branch switching.
 *
 * Pitchfork: x' = p*x - x^3. The trivial branch x=0 exists for all p; at p=0 a
 * BRANCH POINT occurs where the nontrivial branch x = +/- sqrt(p) bifurcates.
 * We continue the trivial branch, require a branch point is detected near p=0,
 * then switch_branch must land on a DISTINCT branch following |x| ~ sqrt(p).
 */
#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

static bool f_pitch(const double* x, double p, double* o, std::string*) {
  o[0] = p * x[0] - x[0]*x[0]*x[0];
  return true;
}

int main() {
  int fails = 0;
  Model m; m.n = 1; m.vector_field = f_pitch;
  ContinuationSettings s; s.p_min = -1.0; s.p_max = 1.0; s.h0 = 0.02; s.max_points = 400;
  std::vector<double> x0 = {0.0};
  Branch b = continue_equilibrium(m, x0, -0.8, s);

  const BranchPoint* bp = nullptr;
  for (std::size_t i : b.special_indices)
    if (b.points[i].special == SpecialPointKind::BranchPoint) bp = &b.points[i];

  bool t1 = b.ok && bp && std::fabs(bp->p) < 0.05 && std::fabs(bp->x[0]) < 0.05;
  printf("  T1 detect branch point: ok=%d found=%d p=%.4f x=%.4f %s\n",
         b.ok, bp != nullptr, bp ? bp->p : 0.0, bp ? bp->x[0] : 0.0, t1 ? "PASS" : "FAIL");
  if (!t1) { fails++; }

  if (bp) {
    Branch sw = switch_branch(m, *bp, s);
    bool good_shape = false;
    if (sw.ok && sw.points.size() > 3) {
      const BranchPoint& mid = sw.points[sw.points.size()/2];
      const double expect = mid.p > 0 ? std::sqrt(mid.p) : 0.0;
      good_shape = (expect > 0) && std::fabs(std::fabs(mid.x[0]) - expect) < 0.15;
    }
    bool t2 = sw.ok && sw.points.size() > 3 && good_shape;
    printf("  T2 switch branch: ok=%d points=%zu shape~sqrt(p)=%d %s\n",
           sw.ok, sw.points.size(), good_shape, t2 ? "PASS" : "FAIL");
    if (!t2) { fails++; }
  } else {
    printf("  T2 switch branch: SKIPPED (no branch point)\n"); fails++;
  }

  if (fails) { printf("branch_switch_smoke: %d FAIL\n", fails); return 1; }
  printf("branch_switch_smoke: all checks pass\n");
  return 0;
}
