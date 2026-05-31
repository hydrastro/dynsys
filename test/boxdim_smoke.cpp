/* Locks box_counting_dimension (this iteration) against known dimensions:
 * filled square ~2, line/circle ~1, Sierpinski ~1.585, Henon ~1.26.
 * make test-boxdim */
/* Validate box_counting_dimension against shapes of KNOWN dimension. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
using namespace dynsys::analysis;
int main(){
  int fails=0;
  std::mt19937 rng(12345);
  std::uniform_real_distribution<double> U(0,1);

  // (1) filled square -> D ~ 2
  { std::vector<double> xs,ys;
    for(int i=0;i<200000;i++){xs.push_back(U(rng));ys.push_back(U(rng));}
    auto R=box_counting_dimension(xs,ys,11);
    printf("(1) filled square: D=%.3f (expect ~2.0) R2=%.3f\n",R.dimension,R.r_squared);
    if(!(R.dimension>1.85 && R.dimension<2.05)){printf("  FAIL\n");fails++;}
  }
  // (2) line segment -> D ~ 1
  { std::vector<double> xs,ys;
    for(int i=0;i<100000;i++){double t=U(rng); xs.push_back(t); ys.push_back(0.37*t+0.1);}
    auto R=box_counting_dimension(xs,ys,12);
    printf("(2) line segment: D=%.3f (expect ~1.0) R2=%.3f\n",R.dimension,R.r_squared);
    if(!(R.dimension>0.9 && R.dimension<1.1)){printf("  FAIL\n");fails++;}
  }
  // (3) circle -> D ~ 1
  { std::vector<double> xs,ys;
    for(int i=0;i<100000;i++){double a=2*M_PI*U(rng); xs.push_back(std::cos(a)); ys.push_back(std::sin(a));}
    auto R=box_counting_dimension(xs,ys,12);
    printf("(3) circle: D=%.3f (expect ~1.0) R2=%.3f\n",R.dimension,R.r_squared);
    if(!(R.dimension>0.9 && R.dimension<1.15)){printf("  FAIL\n");fails++;}
  }
  // (4) Sierpinski triangle via chaos game -> D = log3/log2 ~ 1.585
  { std::vector<double> xs,ys;
    double vx[3]={0,1,0.5}, vy[3]={0,0,0.8660254};
    double x=0.1,y=0.1;
    for(int i=0;i<300000;i++){int k=(int)(U(rng)*3)%3; x=(x+vx[k])/2; y=(y+vy[k])/2; if(i>100){xs.push_back(x);ys.push_back(y);}}
    auto R=box_counting_dimension(xs,ys,12);
    printf("(4) Sierpinski: D=%.3f (expect ~1.585) R2=%.3f\n",R.dimension,R.r_squared);
    if(!(R.dimension>1.50 && R.dimension<1.66)){printf("  FAIL\n");fails++;}
  }
  // (5) Henon attractor -> D ~ 1.26
  { std::vector<double> xs,ys; double x=0.1,y=0;
    for(int i=0;i<300000;i++){double nx=1-1.4*x*x+y, ny=0.3*x; x=nx;y=ny; if(i>1000){xs.push_back(x);ys.push_back(y);}}
    auto R=box_counting_dimension(xs,ys,12);
    printf("(5) Henon: D=%.3f (expect ~1.26) R2=%.3f\n",R.dimension,R.r_squared);
    if(!(R.dimension>1.15 && R.dimension<1.40)){printf("  FAIL\n");fails++;}
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
