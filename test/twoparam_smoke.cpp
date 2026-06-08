#include "analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  // cusp normal form: x' = q + p x - x^3. Fold set: x'=0 and d/dx=0:
  //   p - 3x^2 = 0 -> x^2 = p/3 ; and q = x^3 - p x = x(x^2 - p) = x(p/3 - p)= -2/3 p x
  //   => 27 q^2 = 4 p^3  (cusp). We start near a fold and trace.
  Model2 m; m.n=1;
  m.vector_field=[](const double*x,double p,double q,double*f,std::string*)->bool{
    f[0]=q + p*x[0] - x[0]*x[0]*x[0]; return true; };
  // start: pick p=3 -> x^2=1 -> x=1, q = 1 - 3 = -2 (on the fold). check 27*4=4*27 ok
  TwoParamSettings s; s.p_min=0.2; s.p_max=8; s.q_min=-8; s.q_max=8; s.h0=0.05; s.max_points=600;
  TwoParamCurve c = two_param_curve(m, TwoParamKind::Fold, {1.0}, 3.0, -2.0, s);
  printf("fold curve: ok=%d points=%zu msg=%s\n", c.ok, c.points.size(), c.message.c_str());
  // verify points satisfy 27 q^2 ≈ 4 p^3
  int checked=0; double worst=0;
  for(auto&pt:c.points){ if(pt.p>0.3){ double lhs=27*pt.q*pt.q, rhs=4*pt.p*pt.p*pt.p;
      double rel=std::fabs(lhs-rhs)/(std::fabs(rhs)+1.0); worst=std::max(worst,rel); checked++; } }
  printf("   cusp relation 27q^2=4p^3: worst rel err=%.4f over %d pts\n", worst, checked);
  chk("fold curve traced", c.ok && c.points.size()>10);
  chk("points lie on the cusp 27q^2=4p^3", checked>0 && worst<0.05);
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
