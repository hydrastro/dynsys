#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* x'=y, y'=b1 + b2*y + r*x^2 + x*y. params p=b1, q=b2, r=r.
 * BT curve solved analytically: {b1=0, b2=0, r free}, equilibrium (0,0).
 * (trace = b2+x = 0 and det = -2r*x = 0 => x=0 => b2=0, b1=0.) */
static bool bt3(const double*X,double b1,double b2,double r,double*o,std::string*){
  double x=X[0],y=X[1]; o[0]=y; o[1]=b1 + b2*y + r*x*x + x*y; return true;
}
int main(){
  Model3 m; m.n=2; m.vector_field=bt3;
  std::vector<double> x0={0.02,-0.01};   // near the BT point at r0=1
  BTCurveSettings st; st.ds=0.05; st.max_points=40; st.r_min=-3; st.r_max=3;
  BTCurve C = bt_curve(m, x0, 0.03, -0.02, 1.0, st);
  printf("bt_curve: ok=%d points=%zu\n  msg: %s\n",C.ok,C.points.size(),C.message.c_str());
  // check ALL points have b1~0, b2~0 (the known BT locus), r spread out, residual small
  double maxb1=0,maxb2=0,rmin=1e9,rmax=-1e9,maxres=0;
  for(auto&P:C.points){ maxb1=std::fmax(maxb1,std::fabs(P.p)); maxb2=std::fmax(maxb2,std::fabs(P.q));
    rmin=std::fmin(rmin,P.r); rmax=std::fmax(rmax,P.r); maxres=std::fmax(maxres,P.residual); }
  printf("  max|b1|=%.2e max|b2|=%.2e over curve; r in [%.3f,%.3f]; max residual=%.2e\n",maxb1,maxb2,rmin,rmax,maxres);
  // sample a few
  for(size_t i=0;i<C.points.size(); i+=std::max((size_t)1,C.points.size()/6)){auto&P=C.points[i];
    printf("    pt: b1=%+.4f b2=%+.4f r=%+.4f  (BT a=%.3f b=%.3f)\n",P.p,P.q,P.r,P.a,P.b);}
  bool pass = C.ok && C.points.size()>=8 && maxb1<1e-3 && maxb2<1e-3 && (rmax-rmin)>0.8 && maxres<1e-6;
  printf("%s\n", pass?"bt_curve_smoke: all checks pass":"bt_curve_smoke: FAIL");
  return pass?0:1;
}
