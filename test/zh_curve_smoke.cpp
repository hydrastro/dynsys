#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Zero-Hopf test in 3 params (p,q,r):
 *   x' = p*x - x^2            (fold: zero eigenvalue at x=0 when p=0)
 *   y' = -(1+0.2*r)*z + q*y   (Hopf pair at q=0, frequency 1+0.2r)
 *   z' =  (1+0.2*r)*y + q*z
 * Equilibrium (0,0,0). Jacobian eigenvalues: {p, q +- i(1+0.2r)}.
 * Zero-Hopf locus: p=0 AND q=0, r free. The third parameter r changes the
 * frequency but the locus is exactly {p=0,q=0}. */
static bool zh3(const double*X,double p,double q,double r,double*o,std::string*){
  double x=X[0],y=X[1],z=X[2]; double w=1.0+0.2*r;
  o[0]=(p + 0.3*r)*x - x*x;
  o[1]=-w*z + q*y;
  o[2]= w*y + q*z;
  return true;
}
int main(){
  Model3 m; m.n=3; m.vector_field=zh3;
  std::vector<double> x0={0.0,0.0,0.0};        // near the equilibrium
  BTCurveSettings st; st.ds=0.1; st.max_points=120; st.r_min=-4; st.r_max=4; st.p_min=-4; st.p_max=4; st.q_min=-4; st.q_max=4;
  BTCurve C = zh_curve(m, x0, -0.3, 0.0, 1.0, st);
  printf("zh_curve: ok=%d points=%zu\n  msg: %s\n",C.ok,C.points.size(),C.message.c_str());
  if(!C.ok){return 1;}
  double maxqcond=0,maxfold=0,maxq=0,rmin=1e9,rmax=-1e9,maxres=0;
  for(auto&P:C.points){ maxfold=std::fmax(maxfold,std::fabs(P.p+0.3*P.r)); maxq=std::fmax(maxq,std::fabs(P.q));
    rmin=std::fmin(rmin,P.r); rmax=std::fmax(rmax,P.r); maxres=std::fmax(maxres,P.residual); }
  (void)maxqcond;
  printf("  max|p+0.3r|=%.2e (fold cond) max|q|=%.2e (Hopf cond); r in [%.3f,%.3f]; max residual=%.2e\n",maxfold,maxq,rmin,rmax,maxres);
  for(size_t i=0;i<C.points.size(); i+=std::max((size_t)1,C.points.size()/6)){auto&P=C.points[i];
    printf("    pt: p=%+.4f q=%+.4f r=%+.4f  (ZH b=%.3f Re(c)=%.3f)\n",P.p,P.q,P.r,P.a,P.b);}
  bool pass = C.ok && C.points.size()>=8 && maxfold<1e-3 && maxq<1e-3 && (rmax-rmin)>2.0 && maxres<1e-6;
  printf("%s\n", pass?"zh_curve_smoke: all checks pass":"zh_curve_smoke: FAIL");
  return pass?0:1;
}
