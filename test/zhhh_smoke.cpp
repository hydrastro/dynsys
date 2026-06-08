/* Regression guard for the zero-Hopf (ZH / fold-Hopf) and Hopf-Hopf (HH /
 * double-Hopf) codim-2 detections on a Hopf curve.
 *
 * ZH: 3-D system, oscillator (Re=p, freq 1) plus a real mode z'=q z - z^2.
 *     The Hopf curve is p=0; along it the real eigenvalue q crosses zero at
 *     q=0 -> zero-Hopf at (0,0).
 * HH: 4-D system, two oscillators (Re=p freq 1; Re=q freq 2). The Hopf curve
 *     (active pair Re=p) is p=0; the second pair reaches the axis at q=0 ->
 *     Hopf-Hopf at (0,0).
 */
#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

static bool f3(const double* X, double p, double q, double* o, std::string*) {
  double x = X[0], y = X[1], z = X[2], r2 = x*x + y*y;
  o[0] = p*x - y - x*r2; o[1] = x + p*y - y*r2; o[2] = q*z - z*z; return true;
}
static bool f4(const double* X, double p, double q, double* o, std::string*) {
  double x1 = X[0], y1 = X[1], x2 = X[2], y2 = X[3], r1 = x1*x1 + y1*y1, r2 = x2*x2 + y2*y2;
  o[0] = p*x1 - y1 - x1*r1; o[1] = x1 + p*y1 - y1*r1;
  o[2] = q*x2 - 2*y2 - x2*r2; o[3] = 2*x2 + q*y2 - y2*r2; return true;
}

int main() {
  int fails = 0;

  {
    Model2 m2; m2.n = 3; m2.vector_field = f3;
    TwoParamSettings st; st.h0 = 0.05; st.max_points = 200;
    TwoParamCurve C = two_param_curve(m2, TwoParamKind::Hopf, {0,0,0}, 0.0, 0.5, st);
    int zh = 0;
    for (size_t idx : C.special_indices)
      if (C.points[idx].special == SpecialPointKind::ZeroHopf && std::fabs(C.points[idx].q2) < 0.1) zh++;
    printf("  ZH (3-D fold-Hopf): curve points=%zu specials=%zu zh_found=%d %s\n",
           C.points.size(), C.special_indices.size(), zh, zh >= 1 ? "PASS" : "FAIL");
    if (zh < 1) fails++;
  }
  {
    Model2 m2; m2.n = 4; m2.vector_field = f4;
    TwoParamSettings st; st.h0 = 0.05; st.max_points = 200;
    TwoParamCurve C = two_param_curve(m2, TwoParamKind::Hopf, {0,0,0,0}, 0.0, 0.5, st);
    int hh = 0;
    for (size_t idx : C.special_indices)
      if (C.points[idx].special == SpecialPointKind::HopfHopf && std::fabs(C.points[idx].q2) < 0.1) hh++;
    printf("  HH (4-D double-Hopf): curve points=%zu specials=%zu hh_found=%d %s\n",
           C.points.size(), C.special_indices.size(), hh, hh >= 1 ? "PASS" : "FAIL");
    if (hh < 1) fails++;
  }

  if (fails) { printf("zhhh_smoke: %d FAIL\n", fails); return 1; }
  printf("zhhh_smoke: all checks pass\n");
  return 0;
}
