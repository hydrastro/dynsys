/* Locks the basin non-convergent fix (this iteration): a bistable system
 * gives 2 basins / 0 non-convergent; a chaotic Henon yields 0 false
 * attractors (non-convergent dominates). make test-basinchaos */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
int main(){
  int fails=0;
  // Bistable (2 fixed pts) — should converge, 2 basins, ~0 nonconvergent
  { const double h=0.01;
    auto adv=[&](double x,double y,double*nx,double*ny){
      auto f=[](double X,double Y,double*u,double*v){*u=X-X*X*X;*v=-Y;};
      double k1u,k1v,k2u,k2v,k3u,k3v,k4u,k4v;
      f(x,y,&k1u,&k1v);f(x+.5*h*k1u,y+.5*h*k1v,&k2u,&k2v);
      f(x+.5*h*k2u,y+.5*h*k2v,&k3u,&k3v);f(x+h*k3u,y+h*k3v,&k4u,&k4v);
      *nx=x+h/6*(k1u+2*k2u+2*k3u+k4u);*ny=y+h/6*(k1v+2*k2v+2*k3v+k4v);return true;};
    BasinOptions o;o.xmin=-2;o.xmax=2;o.ymin=-2;o.ymax=2;o.width=40;o.height=40;o.max_steps=4000;o.settle_tol=1e-5;o.cluster_tol=0.05;
    auto B=compute_basins(adv,o);
    printf("bistable: %zu basins, conv=%ld nonconv=%ld div=%ld\n",B.attractors.size(),B.n_converged,B.n_nonconvergent,B.n_diverged);
    if(B.attractors.size()!=2){printf("  FAIL basins\n");fails++;}
    if(B.n_nonconvergent > B.n_converged){printf("  FAIL too many nonconv\n");fails++;}
  }
  // Henon chaotic map — endpoints never settle -> should be mostly nonconvergent
  { auto adv=[&](double x,double y,double*nx,double*ny){*nx=1-1.4*x*x+y;*ny=0.3*x;return true;};
    BasinOptions o;o.xmin=-2;o.xmax=2;o.ymin=-2;o.ymax=2;o.width=40;o.height=40;o.max_steps=1500;o.settle_tol=1e-6;o.cluster_tol=0.05;
    auto B=compute_basins(adv,o);
    long total=B.n_converged+B.n_nonconvergent+B.n_diverged;
    printf("henon-chaos: conv=%ld nonconv=%ld div=%ld (nonconv should dominate the bounded part)\n",B.n_converged,B.n_nonconvergent,B.n_diverged);
    // the bounded cells should be mostly nonconvergent (chaotic), not falsely clustered
    if(B.n_converged > total/2){printf("  FAIL: chaotic endpoints falsely clustered as attractors\n");fails++;}
  }
  printf("=== %s ===\n",fails==0?"PASS":"FAIL");
  return fails;
}
