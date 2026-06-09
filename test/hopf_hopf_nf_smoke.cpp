#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Two oscillators in polar-to-Cartesian form at the origin:
 *   mode 1 (x1,x2), freq W1, self-cubic A1 (r1^3 term)
 *   mode 2 (x3,x4), freq W2, self-cubic A2
 *   cross coupling: mode1 gets G12*r2^2 driving, mode2 gets G21*r1^2.
 * In Cartesian, r1^2 = x1^2+x2^2. The amplitude self term A1*r1^2 on r1 becomes
 * A1*(x1^2+x2^2) multiplying (x1,x2). So:
 *   x1' = -W1 x2 + (A1*(x1^2+x2^2) + G12*(x3^2+x4^2)) * x1
 *   x2' =  W1 x1 + (A1*(x1^2+x2^2) + G12*(x3^2+x4^2)) * x2
 *   x3' = -W2 x4 + (A2*(x3^2+x4^2) + G21*(x1^2+x2^2)) * x3
 *   x4' =  W2 x3 + (A2*(x3^2+x4^2) + G21*(x1^2+x2^2)) * x4
 * Then the amplitude normal form is EXACTLY r1'=r1(A1 r1^2 + G12 r2^2),
 * r2'=r2(A2 r2^2 + G21 r1^2), so p11=A1, p22=A2, p12=G12, p21=G21
 * (up to the overall convention factor the routine uses, which we check by
 * RATIO consistency + sign). */
static double W1=1.0, W2=1.7, A1=-0.3, A2=-0.5, G12=0.4, G21=0.25;
static bool hh(const double*X,double mu,double*o,std::string*){(void)mu;
  double x1=X[0],x2=X[1],x3=X[2],x4=X[3];
  double r1=x1*x1+x2*x2, r2=x3*x3+x4*x4;
  double c1=A1*r1+G12*r2, c2=A2*r2+G21*r1;
  o[0]=-W1*x2 + c1*x1;
  o[1]= W1*x1 + c1*x2;
  o[2]=-W2*x4 + c2*x3;
  o[3]= W2*x3 + c2*x4;
  return true;
}
int main(){
  Model m; m.n=4; m.vector_field=hh;
  std::vector<double> x={0,0,0,0};
  double p11,p12,p21,p22,om1,om2; std::string e;
  bool ok=hopf_hopf_normal_form(m,x,0.0,&p11,&p12,&p21,&p22,&om1,&om2,&e);
  printf("hopf_hopf_normal_form: ok=%d\n",ok); if(!ok){printf("err: %s\n",e.c_str());return 1;}
  printf("  omega1=%.4f omega2=%.4f (expect %.4f, %.4f)\n",om1,om2,W1,W2);
  printf("  p11=%.4f p22=%.4f  (prescribed A1=%.3f A2=%.3f)\n",p11,p22,A1,A2);
  printf("  p12=%.4f p21=%.4f  (prescribed G12=%.3f G21=%.3f)\n",p12,p21,G12,G21);
  // The routine's convention may scale all four by a common factor k. Check the
  // RATIOS match the prescribed ratios (scale-invariant) + signs.
  double k = p11/A1;  // inferred common scale
  printf("  inferred scale k=p11/A1=%.4f; k*A2=%.4f vs p22=%.4f; k*G12=%.4f vs p12=%.4f; k*G21=%.4f vs p21=%.4f\n",
         k, k*A2,p22, k*G12,p12, k*G21,p21);
  bool freqs = std::fabs(om1-W1)<1e-2 && std::fabs(om2-W2)<1e-2;
  bool ratios = std::fabs(k*A2-p22)<0.02*std::fabs(p22>1e-6?p22:1)+1e-3
             && std::fabs(k*G12-p12)<0.05*std::fabs(p12)+5e-3
             && std::fabs(k*G21-p21)<0.05*std::fabs(p21)+5e-3;
  bool signs = (p11*A1>0)&&(p22*A2>0)&&(p12*G12>0)&&(p21*G21>0);
  printf("  freqs=%d ratios=%d signs=%d\n",freqs,ratios,signs);
  bool pass = ok && freqs && signs && ratios;
  printf("%s\n", pass?"hopf_hopf_nf_smoke: all checks pass":"hopf_hopf_nf_smoke: FAIL");
  return pass?0:1;
}
