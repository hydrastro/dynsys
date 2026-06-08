#include "analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
int main(){
  // same Bautin quintic; LPC curve in (mu, a): fold at mu = -a^2/4.
  Model2 m; m.n=2;
  m.vector_field=[](const double*X,double mu,double a,double*f,std::string*)->bool{
    double R=X[0]*X[0]+X[1]*X[1]; double g=mu+a*R-R*R;
    f[0]=-X[1]+X[0]*g; f[1]=X[0]+X[1]*g; return true; };
  double mu0=0.1, a0=1.0; double R=(a0+std::sqrt(a0*a0+4*mu0))/2; double rad=std::sqrt(R);
  int mesh=50; std::vector<std::vector<double>> guess;
  for(int i=0;i<mesh;i++){double th=2*M_PI*i/mesh; guess.push_back({rad*std::cos(th),rad*std::sin(th)});}
  TwoParamSettings s; s.p_min=-0.5; s.p_max=0.5; s.q_min=0.8; s.q_max=1.4; s.max_points=24;
  CycleSettings cs; cs.mesh=mesh; cs.dp=0.01; cs.max_steps=200;
  LPCCurve c = lpc_curve(m, guess, 2*M_PI, mu0, a0, s, cs);
  printf("LPC curve: ok=%d points=%zu msg=%s\n", c.ok, c.points.size(), c.message.c_str());
  int checked=0; double worst=0;
  for(auto&pt:c.points){ double expmu=-pt.q*pt.q/4; double err=std::fabs(pt.p-expmu);
    printf("   a=%.3f mu=%.4f (exp %.4f) err=%.4f\n", pt.q, pt.p, expmu, err);
    worst=std::max(worst,err); checked++; }
  printf("worst |mu - (-a^2/4)| = %.4f over %d pts\n", worst, checked);
  printf("=== %s ===\n", (c.ok && checked>=3 && worst<0.05)?"PASS":"FAIL");
  return (c.ok && checked>=3 && worst<0.05)?0:1;
}
