#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

/* Part 1: DIRECT BT normal-form coefficient check.
 * Use x'=y, y' = p + x^2 + x*y (single param p, BT at p=0, origin).
 * J at origin = [[0,1],[0,0]] (double-zero, Jordan block). Analytic coeffs:
 *   B(q0,q0) with q0=(1,0): only 2nd-order term in y' is x^2 -> Hess=2 -> (0,2)
 *   a = 1/2 <p1,(0,2)> with p1=(0,1) (after <p1,q1>=1) = 1
 *   b = <p0,B(q0,q0)> + <p1,B(q0,q1)>; B(q0,q1)=B((1,0),(0,1)) from x*y -> (0,1)
 *       -> b = 1 (the x*y coefficient). So expect a~1, b~1. */
static bool f_dir(const double*X,double p,double*o,std::string*){
  double x=X[0],y=X[1]; o[0]=y; o[1]=p + x*x + x*y; return true;
}

/* Part 2: two-parameter family with a genuine BT where fold & Hopf curves meet.
 *   x' = y
 *   y' = b1 + b2*y + x^2 + x*y     (parameters p=b1, q=b2)
 * Equilibria: y=0, b1 + x^2 = 0. J=[[0,1],[2x+y, b2+x]].
 * Fold: det J = -(2x+y) = 0 -> x=0 (at y=0) -> b1=0.  Hopf: tr J = b2+x = 0.
 * BT where both: x=0, b2=0, b1=0. So a BT sits at (b1,b2)=(0,0). */
static bool f_2p(const double*X,double b1,double b2,double*o,std::string*){
  double x=X[0],y=X[1]; o[0]=y; o[1]=b1 + b2*y + x*x + x*y; return true;
}

int main(){
  int fails=0;

  Model m; m.n=2; m.vector_field=f_dir;
  double a=0,b=0; std::string err;
  std::vector<double> x0={0.0,0.0};
  bool ok=bt_normal_form(m, x0, 0.0, &a, &b, &err);
  bool p1pass = ok && std::fabs(a-1.0)<0.2 && std::fabs(b-1.0)<0.3;
  printf("[P1] BT coeffs at origin: ok=%d a=%.4f b=%.4f (expect ~1,~1) %s\n",
         ok,a,b, p1pass?"PASS":"FAIL");
  if(!p1pass) fails++;

  /* Part 2: seed a fold point on the fold curve (x=0,y=0 at b1=0) for some b2,
   * then trace the fold curve and look for the BT at b2=0. Seed at b2=-0.6. */
  Model2 m2; m2.n=2; m2.vector_field=f_2p;
  TwoParamSettings s; s.p_min=-2; s.p_max=2; s.q_min=-2; s.q_max=2; s.h0=0.025; s.max_points=800;
  std::vector<double> xs={0.0,0.0};
  TwoParamCurve cur=two_param_curve(m2, TwoParamKind::Fold, xs, 0.0, -0.6, s);
  printf("[P2] fold curve: ok=%d points=%zu specials=%zu\n",cur.ok,cur.points.size(),cur.special_indices.size());
  bool foundBT=false;
  for(size_t idx:cur.special_indices){auto&pt=cur.points[idx];
    const char* kn = pt.special==SpecialPointKind::BogdanovTakens?"BT":
                     pt.special==SpecialPointKind::Cusp?"Cusp":
                     pt.special==SpecialPointKind::GeneralizedHopf?"GH":"other";
    printf("    %s at refined (p2=%.4f,q2=%.4f) nf=%d a=%.4f b=%.4f\n",
           kn,pt.p2,pt.q2,pt.has_codim2_nf,pt.bt_a,pt.bt_b);
    if(pt.special==SpecialPointKind::BogdanovTakens){
      foundBT=true;
      /* BT should be near (b1,b2)=(0,0); coeffs near a~1,b~1 */
      bool loc = std::fabs(pt.p2)<0.1 && std::fabs(pt.q2)<0.1;
      bool nf  = pt.has_codim2_nf && std::fabs(pt.bt_a-1.0)<0.4 && std::fabs(pt.bt_b-1.0)<0.5;
      printf("      -> location %s, normal-form %s\n", loc?"GOOD":"off", nf?"GOOD":"off");
      if(!loc) fails++;
    }
  }
  if(!foundBT){ printf("    (no BT detected on the fold curve)\n"); fails++; }

  if(fails){ printf("bt_codim2_smoke: %d FAIL\n",fails); return 1; }
  printf("bt_codim2_smoke: all checks pass\n");
  return 0;
}
