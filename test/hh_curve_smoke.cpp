#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Hopf-Hopf test: two oscillators, dampings controlled by params, with a 3rd
 * parameter r shifting the second oscillator's damping.
 *   block 1 (x1,x2): damping p, frequency 1     -> pair p +- i
 *   block 2 (x3,x4): damping (q + 0.4 r), freq 2.3 -> pair (q+0.4r) +- 2.3i
 * Both pairs on the imaginary axis  <=>  p=0 AND q+0.4r=0.
 * Hopf-Hopf locus: { p=0, q=-0.4r, r free }, equilibrium origin. */
static bool hh4(const double*X,double p,double q,double r,double*o,std::string*){
  double x1=X[0],x2=X[1],x3=X[2],x4=X[3];
  double d2=q+0.4*r;
  o[0]= p*x1 - 1.0*x2;
  o[1]= 1.0*x1 + p*x2;
  o[2]= d2*x3 - 2.3*x4;
  o[3]= 2.3*x3 + d2*x4;
  return true;
}
int main(){
  Model3 m; m.n=4; m.vector_field=hh4;
  std::vector<double> x0={0,0,0,0};
  BTCurveSettings st; st.ds=0.1; st.max_points=120; st.r_min=-4; st.r_max=4; st.p_min=-4; st.p_max=4; st.q_min=-4; st.q_max=4;
  BTCurve C = hh_curve(m, x0, 0.0, 0.0, 0.0, st);  // start at HH (p=0,q=0,r=0)
  printf("hh_curve: ok=%d points=%zu\n  msg: %s\n",C.ok,C.points.size(),C.message.c_str());
  if(!C.ok) return 1;
  double maxp=0,maxcond=0,rmin=1e9,rmax=-1e9,maxres=0;
  for(auto&P:C.points){ maxp=std::fmax(maxp,std::fabs(P.p)); maxcond=std::fmax(maxcond,std::fabs(P.q+0.4*P.r));
    rmin=std::fmin(rmin,P.r); rmax=std::fmax(rmax,P.r); maxres=std::fmax(maxres,P.residual); }
  printf("  max|p|=%.2e (Hopf1) max|q+0.4r|=%.2e (Hopf2); r in [%.3f,%.3f]; max res=%.2e\n",maxp,maxcond,rmin,rmax,maxres);
  for(size_t i=0;i<C.points.size(); i+=std::max((size_t)1,C.points.size()/6)){auto&P=C.points[i];
    printf("    pt: p=%+.4f q=%+.4f r=%+.4f  (w1=%.3f w2=%.3f)\n",P.p,P.q,P.r,P.a,P.b);}
  bool pass = C.ok && C.points.size()>=8 && maxp<1e-3 && maxcond<1e-3 && (rmax-rmin)>2.0 && maxres<1e-6;
  printf("%s\n", pass?"hh_curve_smoke: all checks pass":"hh_curve_smoke: FAIL");
  return pass?0:1;
}
