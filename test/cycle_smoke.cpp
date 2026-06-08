#include "analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  // supercritical Hopf normal form, param mu: cycle radius sqrt(mu), period 2pi
  Model m; m.n=2;
  m.vector_field=[](const double*X,double mu,double*f,std::string*)->bool{
    double r2=X[0]*X[0]+X[1]*X[1];
    f[0]=-X[1]+X[0]*(mu-r2); f[1]=X[0]+X[1]*(mu-r2); return true; };
  // initial guess at mu=1: a circle of radius 1, M points, period ~2pi
  const int M=60; std::vector<std::vector<double>> guess;
  for(int i=0;i<M;i++){ double th=2*M_PI*i/M; guess.push_back({std::cos(th), std::sin(th)}); }
  CycleSettings s; s.mesh=M; s.p_min=0.05; s.p_max=4.0; s.dp=0.05; s.max_steps=200;
  CycleBranch b = continue_limit_cycle(m, guess, 2*M_PI, 1.0, s);
  printf("cycle branch: ok=%d samples=%zu msg=%s\n", b.ok, b.samples.size(), b.message.c_str());
  // check the mu=1 sample: amplitude ~ 2*radius = 2 (peak-to-peak of x in [-1,1]), period ~ 2pi
  double worstA=0, worstT=0; int checked=0;
  for(auto&smp:b.samples){
    if(smp.p<0.1) continue;
    double expected_amp = 2.0*std::sqrt(smp.p);   // peak-to-peak of x = 2*sqrt(mu)
    double expected_T = 2*M_PI;
    worstA=std::max(worstA, std::fabs(smp.amplitude-expected_amp));
    worstT=std::max(worstT, std::fabs(smp.period-expected_T));
    checked++;
  }
  printf("   over %d samples: worst |amp-2sqrt(mu)|=%.4f, worst |T-2pi|=%.4f\n", checked, worstA, worstT);
  chk("cycle continuation traced", b.ok && b.samples.size()>20);
  chk("amplitude ~ 2 sqrt(mu) along the branch", checked>0 && worstA<0.05);
  chk("period ~ 2pi along the branch", checked>0 && worstT<0.05);
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
