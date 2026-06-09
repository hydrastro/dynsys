#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Cusp test. 1-D-in-x cusp normal form with params (p=alpha, q=beta) and a 3rd
 * parameter r that shifts beta:
 *   x' = p + (q + 0.5*r)*x - x^3
 * The fold set: x' = 0 and d/dx = (q+0.5r) - 3x^2 = 0.
 * The CUSP is where additionally d2/dx2 = -6x = 0 => x=0 => q+0.5r=0 and p=0.
 * So the cusp locus is { p=0, q=-0.5r, r free }, equilibrium x=0.
 * (To make it >=2-dim for the eigen machinery, add a trivial stable y.) */
static bool cusp3(const double*X,double p,double q,double r,double*o,std::string*){
  double x=X[0],y=X[1];
  o[0]=p + (q + 0.5*r)*x - x*x*x;
  o[1]=-y;   /* trivial decoupled stable direction */
  return true;
}
int main(){
  Model3 m; m.n=2; m.vector_field=cusp3;
  std::vector<double> x0={0.0,0.0};
  BTCurveSettings st; st.ds=0.1; st.max_points=120; st.r_min=-4; st.r_max=4; st.p_min=-4; st.p_max=4; st.q_min=-4; st.q_max=4;
  BTCurve C = cusp_curve(m, x0, 0.0, 0.0, 0.0, st);  // start at cusp (p=0,q=0,r=0)
  printf("cusp_curve: ok=%d points=%zu\n  msg: %s\n",C.ok,C.points.size(),C.message.c_str());
  if(!C.ok) return 1;
  double maxp=0,maxcond=0,rmin=1e9,rmax=-1e9,maxres=0;
  for(auto&P:C.points){ maxp=std::fmax(maxp,std::fabs(P.p)); maxcond=std::fmax(maxcond,std::fabs(P.q+0.5*P.r));
    rmin=std::fmin(rmin,P.r); rmax=std::fmax(rmax,P.r); maxres=std::fmax(maxres,P.residual); }
  printf("  max|p|=%.2e (fold/alpha) max|q+0.5r|=%.2e (cusp cond); r in [%.3f,%.3f]; max res=%.2e\n",maxp,maxcond,rmin,rmax,maxres);
  for(size_t i=0;i<C.points.size(); i+=std::max((size_t)1,C.points.size()/6)){auto&P=C.points[i];
    printf("    pt: p=%+.4f q=%+.4f r=%+.4f  (cusp c=%.3f)\n",P.p,P.q,P.r,P.a);}
  bool pass = C.ok && C.points.size()>=8 && maxp<1e-3 && maxcond<1e-3 && (rmax-rmin)>2.0 && maxres<1e-6;
  printf("%s\n", pass?"cusp_curve_smoke: all checks pass":"cusp_curve_smoke: FAIL");
  return pass?0:1;
}
