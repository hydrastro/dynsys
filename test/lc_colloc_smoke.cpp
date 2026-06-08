#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  // Supercritical Hopf: x'=-y + x(mu - (x^2+y^2)), y'= x + y(mu - (x^2+y^2)).
  // For mu>0: stable circular cycle radius sqrt(mu), period ~ 2*pi.
  Model m; m.n=2;
  m.vector_field=[](const double*X,double mu,double*f,std::string*)->bool{
    double r2=X[0]*X[0]+X[1]*X[1];
    f[0]=-X[1]+X[0]*(mu-r2); f[1]=X[0]+X[1]*(mu-r2); return true; };
  // initial guess at mu=0.5: circle radius sqrt(0.5)=0.707, period 2pi
  double mu0=0.5, R=std::sqrt(mu0); int mesh=60;
  std::vector<std::vector<double>> guess;
  for(int i=0;i<mesh;i++){ double th=2*M_PI*i/mesh; guess.push_back({R*std::cos(th), R*std::sin(th)}); }
  CycleSettings s; s.mesh=mesh; s.p_min=0.05; s.p_max=1.5; s.dp=0.05; s.max_steps=60;
  CycleBranch b = continue_limit_cycle(m, guess, 2*M_PI, mu0, s);
  printf("cycle branch: ok=%d samples=%zu msg=%s\n", b.ok, b.samples.size(), b.message.c_str());
  if(b.ok && !b.samples.empty()){
    // check a sample near mu=0.5: amplitude (peak-to-peak) ~ 2*sqrt(mu)=1.414, period~6.283
    double worstA=1e9, worstT=1e9; 
    for(auto&smp:b.samples){
      if(std::fabs(smp.p-0.5)<0.06){
        double expA=2*std::sqrt(smp.p), expT=2*M_PI;
        printf("   mu=%.3f period=%.4f (exp %.4f) amp=%.4f (exp %.4f)\n", smp.p, smp.period, expT, smp.amplitude, expA);
        worstA=std::min(worstA, std::fabs(smp.amplitude-expA));
        worstT=std::min(worstT, std::fabs(smp.period-expT));
      }
    }
    chk("period ~ 2pi near mu=0.5", worstT<0.2);
    chk("amplitude ~ 2sqrt(mu) near mu=0.5", worstA<0.15);
    // amplitude should grow as sqrt(mu): compare two mus
    double aLo=0,aHi=0,pLo=0,pHi=0;
    for(auto&smp:b.samples){ if(std::fabs(smp.p-0.2)<0.06){aLo=smp.amplitude;pLo=smp.p;} if(std::fabs(smp.p-1.0)<0.06){aHi=smp.amplitude;pHi=smp.p;} }
    if(aLo>0&&aHi>0){ printf("   amp(mu=%.2f)=%.3f < amp(mu=%.2f)=%.3f ?\n",pLo,aLo,pHi,aHi);
      chk("amplitude grows with mu", aHi>aLo); }
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
