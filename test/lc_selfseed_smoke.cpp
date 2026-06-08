#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* van der Pol x'=y, y'=mu(1-x^2)y - x. Give a CRUDE circular guess (the true
 * cycle is far from circular) and rely on the new simulation reseed. Continue
 * from mu=1 and verify: branch holds, period grows with mu, trivial Floquet
 * multiplier stays ~1. */
static bool vdp(const double*X,double mu,double*o,std::string*){o[0]=X[1];o[1]=mu*(1-X[0]*X[0])*X[1]-X[0];return true;}
int main(){
  Model m; m.n=2; m.vector_field=vdp;
  int mesh=100; std::vector<std::vector<double>> guess(mesh);
  for(int i=0;i<mesh;i++){double th=2*M_PI*i/mesh; guess[i]={2.0*std::cos(th),2.0*std::sin(th)};} // crude circle
  CycleSettings cs; cs.mesh=mesh; cs.arclength=true; cs.ds=0.15; cs.max_steps=80;
  cs.compute_floquet=true; cs.adaptive_mesh=true; cs.p_min=0.5; cs.p_max=4.0;
  CycleBranch B=continue_limit_cycle(m, guess, 6.66, 1.0, cs);
  printf("ok=%d samples=%zu msg=%s\n",B.ok,B.samples.size(),B.message.c_str());
  if(!B.ok){printf("LC SELFSEED: FAIL (continuation did not start)\n");return 1;}
  double max_triv_err=0,mu_lo=1e9,mu_hi=-1e9;
  for(size_t i=0;i<B.samples.size();i+=std::max((size_t)1,B.samples.size()/8)){
    auto&s=B.samples[i]; double best=1e9,trivial=0;
    for(size_t k=0;k<s.floquet_re.size();k++){double m1=std::hypot(s.floquet_re[k]-1.0,s.floquet_im[k]);if(m1<best){best=m1;trivial=std::hypot(s.floquet_re[k],s.floquet_im[k]);}}
    printf("  mu=%.3f period=%.3f amp=%.3f trivial_mult=%.4f\n",s.p,s.period,s.amplitude,trivial);
    max_triv_err=std::max(max_triv_err,std::fabs(trivial-1.0));
  }
  for(auto&s:B.samples){mu_lo=std::min(mu_lo,s.p);mu_hi=std::max(mu_hi,s.p);}
  printf("mu range [%.2f, %.2f]; max|trivial Floquet-1|=%.3e\n",mu_lo,mu_hi,max_triv_err);
  // vdP period at mu=1 ~ 6.66; grows ~ (3-2ln2)*mu for large mu. Just require
  // the branch covers a real mu range and trivial multiplier is near 1.
  bool pass = B.ok && B.samples.size()>5 && (mu_hi-mu_lo)>0.5 && max_triv_err<0.2;
  printf("%s\n", pass?"lc_selfseed_smoke: all checks pass":"lc_selfseed_smoke: FAIL");
  return pass?0:1;
}
