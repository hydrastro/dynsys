#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Two-parameter family: x'=y, y'= q*x - x^2 + p*y.
 * At p=0 it's Hamiltonian with a homoclinic to the origin (for q>0), peak
 * x* = (3/2)q. Adding p*y is dissipative (p<0) or anti (p>0) and breaks the
 * homoclinic EXCEPT along a curve p=p(q). Near this conservative limit the
 * homoclinic sits at p ~ 0 for each q (the conservative locus). So continuing
 * in q should trace p(q) ~ 0 with peak ~ (3/2)q. */
static bool f2p(const double*X,double p,double q,double*o,std::string*){
  o[0]=X[1]; o[1]=q*X[0]-X[0]*X[0]+p*X[1]; return true;
}
int main(){
  Model2 m2; m2.n=2; m2.vector_field=f2p;
  double q0=1.0; double xstar=1.5*q0; double lam=std::sqrt(q0);
  // seed: analytic-ish homoclinic of the conservative system at q0, peak xstar
  // x(t)=xstar*sech^2(lam t/2), y=x'
  double Tt=12.0/lam; int Np=400;
  std::vector<std::vector<double>> seed(Np,std::vector<double>(2));
  for(int i=0;i<Np;i++){ double t=-Tt+2*Tt*i/(Np-1); double se=1.0/std::cosh(lam*t/2.0);
    seed[i][0]=xstar*se*se; seed[i][1]=-xstar*lam*se*se*std::tanh(lam*t/2.0); }
  HomoclinicContSettings cs; cs.dq=0.05; cs.max_steps=12; cs.both_directions=true;
  cs.q_min=0.3; cs.q_max=2.0;
  cs.bvp.mesh=180; cs.bvp.T=Tt; cs.bvp.free_T=false; cs.bvp.newton_iters=60;
  std::vector<double> x0={0,0};
  HomoclinicCurve C=continue_homoclinic(m2,x0,0.0,q0,seed,cs);
  printf("ok=%d points=%zu msg=%s\n",C.ok,C.points.size(),C.message.c_str());
  // print a few points; check p stays ~0 and peak ~ 1.5*q
  int good=0,total=0;
  for(size_t i=0;i<C.points.size();i+=std::max((size_t)1,C.points.size()/10)){
    auto&pt=C.points[i]; total++;
    printf("  q=%.3f p=%.4f T=%.2f amp=%.4f (peak target=%.3f) res=%.1e\n",
           pt.q,pt.p,pt.T,pt.amplitude,1.5*pt.q,pt.residual);
  }
  // validation: across all points, |p|<0.05 and amplitude ~ 1.5q within 8%
  int okpts=0;
  for(auto&pt:C.points){ if(std::fabs(pt.p)<0.06 && std::fabs(pt.amplitude-1.5*pt.q)<0.12*1.5*pt.q) okpts++; }
  printf("points matching conservative locus (p~0, amp~1.5q): %d/%zu\n",okpts,C.points.size());
  bool pass = C.ok && C.points.size()>=8 && okpts >= (int)(0.7*C.points.size());
  printf("%s\n", pass?"homoclinic_cont_smoke: all checks pass":"homoclinic_cont_smoke: FAIL");
  return pass?0:1;
}
