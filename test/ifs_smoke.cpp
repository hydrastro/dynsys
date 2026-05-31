/* Locks the chaos game / IFS (this iteration): Sierpinski attractor has
 * box-dim ~1.585, Barnsley fern ~1.8 and the right bounding box, and the
 * generator is deterministic per seed. Cross-checks box_counting_dimension.
 * make test-ifs */
/* Validate chaos_game: the attractors must have the right box-counting
 * dimension and stay bounded. Uses the box_counting_dimension we trust. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
int main(){
  int fails=0;

  // (1) Sierpinski triangle via 3 maps (scale 1/2 toward each vertex).
  //     Box-counting dimension = log3/log2 ~ 1.585.
  {
    std::vector<AffineMap> M(3);
    M[0]={0.5,0,0,0.5, 0.0, 0.0, 1.0/3};
    M[1]={0.5,0,0,0.5, 0.5, 0.0, 1.0/3};
    M[2]={0.5,0,0,0.5, 0.25, 0.5, 1.0/3};
    auto R=chaos_game(M, 300000, 7);
    printf("(1) Sierpinski: ok=%d, %zu pts, bbox[%.2f,%.2f]x[%.2f,%.2f]\n",
           R.ok, R.xs.size(), R.xmin,R.xmax,R.ymin,R.ymax);
    if(!R.ok){fails++;}
    std::vector<double> xs(R.xs.begin(),R.xs.end()), ys(R.ys.begin(),R.ys.end());
    auto D=box_counting_dimension(xs,ys,12);
    printf("    D_box=%.3f (expect ~1.585) R2=%.3f\n", D.dimension, D.r_squared);
    if(!(D.dimension>1.50 && D.dimension<1.66)){printf("    FAIL dimension\n");fails++;}
    // bounded in [0,1]x[0,~0.87]
    if(!(R.xmin>-0.01 && R.xmax<1.01 && R.ymin>-0.01)){printf("    FAIL bounds\n");fails++;}
  }

  // (2) Barnsley fern (4 maps). Should be bounded roughly x in [-2.2,2.7],
  //     y in [0,10], and box-dim ~1.8-1.9.
  {
    std::vector<AffineMap> M(4);
    M[0]={0,0,0,0.16, 0,0, 0.01};
    M[1]={0.85,0.04,-0.04,0.85, 0,1.6, 0.85};
    M[2]={0.20,-0.26,0.23,0.22, 0,1.6, 0.07};
    M[3]={-0.15,0.28,0.26,0.24, 0,0.44, 0.07};
    auto R=chaos_game(M, 400000, 3);
    printf("(2) Barnsley fern: ok=%d, %zu pts, bbox[%.2f,%.2f]x[%.2f,%.2f]\n",
           R.ok, R.xs.size(), R.xmin,R.xmax,R.ymin,R.ymax);
    if(!R.ok){fails++;}
    // sanity bounds for the classic fern
    if(!(R.xmin>-3 && R.xmax<3 && R.ymin>-0.1 && R.ymax<11)){printf("    FAIL fern bounds\n");fails++;}
    std::vector<double> xs(R.xs.begin(),R.xs.end()), ys(R.ys.begin(),R.ys.end());
    auto D=box_counting_dimension(xs,ys,12);
    printf("    D_box=%.3f (fern ~1.8) R2=%.3f\n", D.dimension, D.r_squared);
    if(!(D.dimension>1.6 && D.dimension<2.0)){printf("    FAIL fern dim\n");fails++;}
  }

  // (3) determinism: same seed -> identical first point
  {
    std::vector<AffineMap> M(2);
    M[0]={0.5,0,0,0.5,0,0,0.5}; M[1]={0.5,0,0,0.5,0.5,0.5,0.5};
    auto A=chaos_game(M,1000,42); auto B=chaos_game(M,1000,42);
    bool same = A.xs.size()==B.xs.size() && !A.xs.empty() && A.xs[0]==B.xs[0] && A.xs.back()==B.xs.back();
    printf("(3) determinism (same seed): %s\n", same?"identical":"DIFFER");
    if(!same)fails++;
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
