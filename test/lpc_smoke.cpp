#include "analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  // Subcritical-to-stable cycle via a quintic radial term (Bautin-like):
  // r' = mu r + a r^3 - r^5, theta'=1, with a>0 (so there's a fold of cycles).
  // In Cartesian: x'=-y + x*(mu + a*(x^2+y^2) - (x^2+y^2)^2), y'=x + y*(...).
  // Fold of cycles occurs where d mu/d(r^2) = 0 on mu = -a R^2 + R^4 (R=r^2):
  //   mu(R)= -a R + R^2 ; fold at dmu/dR=0 -> R = a/2 -> mu = -a^2/4.
  // With a=1: cycle fold at mu = -0.25. We trace the LPC in (mu, a) plane.
  double aProbe=1.0;
  Model2 m; m.n=2;
  m.vector_field=[](const double*X,double mu,double a,double*f,std::string*)->bool{
    double R=X[0]*X[0]+X[1]*X[1];
    double g=mu + a*R - R*R;
    f[0]=-X[1]+X[0]*g; f[1]=X[0]+X[1]*g; return true; };
  // seed: at mu=0.1, a=1 there's a stable big cycle. radius solves mu+aR-R^2=0 -> R=(a+sqrt(a^2+4mu))/2
  double mu0=0.1, a0=1.0; double R=(a0+std::sqrt(a0*a0+4*mu0))/2; double rad=std::sqrt(R);
  int mesh=60; std::vector<std::vector<double>> guess;
  for(int i=0;i<mesh;i++){double th=2*M_PI*i/mesh; guess.push_back({rad*std::cos(th),rad*std::sin(th)});}
  // First: 1-param continuation in mu at a=1, expect an LPC near mu=-0.25
  Model m1; m1.n=2;
  m1.vector_field=[aProbe](const double*X,double mu,double*f,std::string*)->bool{
    double R=X[0]*X[0]+X[1]*X[1]; double g=mu+aProbe*R-R*R;
    f[0]=-X[1]+X[0]*g; f[1]=X[0]+X[1]*g; return true; };
  CycleSettings cs; cs.mesh=mesh; cs.p_min=-0.4; cs.p_max=0.5; cs.dp=0.01; cs.max_steps=200;
  CycleBranch b=continue_limit_cycle(m1, guess, 2*M_PI, mu0, cs);
  printf("1-param cycle branch: ok=%d samples=%zu\n", b.ok, b.samples.size());
  int nfold=0; double mufold=999;
  for(auto&s:b.samples) if(s.is_fold){nfold++; if(std::fabs(s.p+0.25)<std::fabs(mufold+0.25)) mufold=s.p;}
  printf("   LPC samples=%d, nearest mu=%.4f (expect ~ -0.25)\n", nfold, mufold);
  chk("cycle branch traced", b.ok && b.samples.size()>10);
  chk("a fold-of-cycles detected near mu=-0.25", nfold>=1 && std::fabs(mufold+0.25)<0.08);
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
// (separate driver below would go here in a real suite)
