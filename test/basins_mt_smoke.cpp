#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Duffing-like bistable map-ish advance with TWO attractors, pure lambdas
 * (thread-safe). Compare serial compute_basins vs parallel compute_basins_mt:
 * the cell_attractor arrays must be IDENTICAL. */
static bool adv(double x,double y,double*nx,double*ny){
  // damped pendulum-ish discrete: contract toward nearest of two fixed points
  // x' = x + 0.1*(-x + (x>0?1:-1)) ... simple bistable
  double tx = (x>=0)? 1.0 : -1.0;
  *nx = x + 0.15*(tx - x);
  *ny = y + 0.15*(0.0 - y);
  return true;
}
int main(){
  BasinOptions o; o.xmin=-2;o.xmax=2;o.ymin=-2;o.ymax=2;o.width=120;o.height=120;
  o.max_steps=300;o.settle_tol=1e-6;o.cluster_tol=0.05;o.diverge_r=1e6;o.max_attractors=8;
  BasinResult S=compute_basins(adv,o);
  auto mk=[&](int){ return AdvanceFn(adv); };
  BasinResult M=compute_basins_mt(mk,o);
  // compare
  int diff=0; size_t N=(size_t)o.width*o.height;
  for(size_t i=0;i<N;i++) if(S.cell_attractor[i]!=M.cell_attractor[i]) diff++;
  printf("serial attractors=%zu  mt attractors=%zu  cell diffs=%d / %zu\n",
         S.attractors.size(), M.attractors.size(), diff, N);
  printf("serial: conv=%ld div=%ld nonconv=%ld | mt: conv=%ld div=%ld nonconv=%ld\n",
         (long)S.n_converged,(long)S.n_diverged,(long)S.n_nonconvergent, (long)M.n_converged,(long)M.n_diverged,(long)M.n_nonconvergent);
  printf("%s\n", (diff==0 && S.attractors.size()==M.attractors.size())?"basins_mt_smoke: all checks pass":"basins_mt_smoke: FAIL");
  return (diff==0)?0:1;
}
