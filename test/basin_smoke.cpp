/* Locks compute_basins (PHASE B/C) against known basins: a bistable ODE
 * x'=x-x^3,y'=-y (two attractors at +/-1, boundary x=0) and Newton's
 * method for z^3-1 (three roots, fractal basins). make test-basin */
/* Validate compute_basins on systems with KNOWN basins. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;

int main(){
  int fails=0;

  /* (1) ODE x'=x-x^3, y'=-y : stable nodes at (+1,0) and (-1,0), saddle at
     origin. Basin boundary is x=0. RK4 small step as the advance fn. */
  {
    const double h=0.01;
    auto advance=[&](double x,double y,double*nx,double*ny)->bool{
      auto f=[](double X,double Y,double*u,double*v){ *u=X-X*X*X; *v=-Y; };
      double k1u,k1v,k2u,k2v,k3u,k3v,k4u,k4v;
      f(x,y,&k1u,&k1v);
      f(x+0.5*h*k1u,y+0.5*h*k1v,&k2u,&k2v);
      f(x+0.5*h*k2u,y+0.5*h*k2v,&k3u,&k3v);
      f(x+h*k3u,y+h*k3v,&k4u,&k4v);
      *nx=x+(h/6)*(k1u+2*k2u+2*k3u+k4u);
      *ny=y+(h/6)*(k1v+2*k2v+2*k3v+k4v);
      return true;
    };
    BasinOptions o; o.xmin=-2;o.xmax=2;o.ymin=-2;o.ymax=2;o.width=60;o.height=60;
    o.max_steps=4000; o.settle_tol=1e-5; o.cluster_tol=0.05;
    auto B=compute_basins(advance,o);
    printf("(1) bistable: found %zu attractors (expect 2)\n", B.attractors.size());
    if(B.attractors.size()!=2){printf("  FAIL attractor count\n");fails++;}
    // check: a point at x=+1.5 lands on the +1 attractor; x=-1.5 on the -1
    // find label at (x>0) vs (x<0): sample cell near (1.5,0) and (-1.5,0)
    auto labelAt=[&](double X,double Y){
      int i=(int)((X-o.xmin)/(o.xmax-o.xmin)*(o.width-1)+0.5);
      int j=(int)((Y-o.ymin)/(o.ymax-o.ymin)*(o.height-1)+0.5);
      return B.cell_attractor[(size_t)j*o.width+i];
    };
    int lp=labelAt(1.5,0), lm=labelAt(-1.5,0);
    printf("  label at (+1.5,0)=%d, (-1.5,0)=%d (should differ)\n", lp, lm);
    if(lp==lm||lp<0||lm<0){printf("  FAIL basin separation\n");fails++;}
    // check the attractor positions are near (+1,0) and (-1,0)
    bool foundp=false,foundm=false;
    for(auto&a:B.attractors){
      if(std::fabs(a.first-1)<0.1&&std::fabs(a.second)<0.1)foundp=true;
      if(std::fabs(a.first+1)<0.1&&std::fabs(a.second)<0.1)foundm=true;
    }
    if(!foundp||!foundm){printf("  FAIL attractor locations\n");fails++;}
  }

  /* (2) Newton's method for z^3-1: 3 roots (1, -1/2±i sqrt3/2), fractal
     basins. As a map: z -> z - (z^3-1)/(3z^2). Expect 3 attractors. */
  {
    auto advance=[&](double x,double y,double*nx,double*ny)->bool{
      // z = x+iy ; f=z^3-1 ; f'=3z^2 ; z' = z - f/f'
      double zr=x, zi=y;
      // z^2
      double z2r=zr*zr-zi*zi, z2i=2*zr*zi;
      // z^3
      double z3r=z2r*zr-z2i*zi, z3i=z2r*zi+z2i*zr;
      double fr=z3r-1, fi=z3i;
      double dpr=3*z2r, dpi=3*z2i; // f'
      double den=dpr*dpr+dpi*dpi;
      if(den<1e-300){*nx=x;*ny=y;return false;}
      // f/f'
      double qr=(fr*dpr+fi*dpi)/den, qi=(fi*dpr-fr*dpi)/den;
      *nx=zr-qr; *ny=zi-qi;
      return true;
    };
    BasinOptions o; o.xmin=-2;o.xmax=2;o.ymin=-2;o.ymax=2;o.width=80;o.height=80;
    o.max_steps=200; o.settle_tol=1e-6; o.cluster_tol=0.05;
    auto B=compute_basins(advance,o);
    printf("(2) Newton z^3-1: found %zu attractors (expect 3)\n", B.attractors.size());
    if(B.attractors.size()!=3){printf("  FAIL root count\n");fails++;}
    // roots: (1,0), (-0.5, 0.8660), (-0.5,-0.8660)
    int hits=0;
    for(auto&a:B.attractors){
      if(std::hypot(a.first-1,a.second)<0.05)hits++;
      if(std::hypot(a.first+0.5,a.second-0.86602540)<0.05)hits++;
      if(std::hypot(a.first+0.5,a.second+0.86602540)<0.05)hits++;
    }
    printf("  matched %d/3 known roots\n", hits);
    if(hits!=3){printf("  FAIL root locations\n");fails++;}
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
