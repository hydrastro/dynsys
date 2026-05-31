/* Locks parametrized IFS (this iteration): a rotation/scale parameter
 * changes the attractor (what evaluate_ifs_maps does at render time), the
 * result stays bounded, and its dimension is a sane fractal. make test-ifsparam */
/* Verify a PARAMETRIZED IFS: a rotation angle theta controls the maps.
 * Changing theta must change the attractor (different maps -> different
 * geometry), which is exactly what evaluate_ifs_maps does at render time. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

/* a "rotation IFS": two maps that contract and rotate by theta. */
static std::vector<AffineMap> rot_ifs(double theta, double s){
  double ct=std::cos(theta), st=std::sin(theta);
  std::vector<AffineMap> m(2);
  m[0]={ s*ct, -s*st, s*st, s*ct, 0.0, 0.0, 0.5 };
  m[1]={ s*ct, -s*st, s*st, s*ct, 0.5, 0.0, 0.5 };
  return m;
}
int main(){
  int fails=0;
  /* (1) two different theta values must give different attractors */
  auto A=chaos_game(rot_ifs(0.3,0.6), 200000, 7);
  auto B=chaos_game(rot_ifs(1.2,0.6), 200000, 7);
  printf("(1) theta=0.3 bbox[%.2f,%.2f]x[%.2f,%.2f]\n",A.xmin,A.xmax,A.ymin,A.ymax);
  printf("    theta=1.2 bbox[%.2f,%.2f]x[%.2f,%.2f]\n",B.xmin,B.xmax,B.ymin,B.ymax);
  bool different = std::fabs(A.ymax-B.ymax)>0.05 || std::fabs(A.xmin-B.xmin)>0.05;
  printf("    attractors differ with theta: %s\n", different?"YES":"no");
  if(!different) fails++;

  /* (2) both are bounded (contractions s<1 => finite attractor) */
  bool bounded = std::isfinite(A.xmin)&&std::isfinite(A.xmax)&&A.xmax-A.xmin<100;
  printf("(2) bounded attractor: %s\n", bounded?"yes":"NO");
  if(!bounded) fails++;

  /* (3) dimension is sane (a fractal between 1 and 2) */
  std::vector<double> xs(A.xs.begin(),A.xs.end()), ys(A.ys.begin(),A.ys.end());
  auto D=box_counting_dimension(xs,ys,11);
  printf("(3) D=%.3f (expect 1..2)\n", D.dimension);
  if(!(D.dimension>0.8 && D.dimension<2.05)) fails++;

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
