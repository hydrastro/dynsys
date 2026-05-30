/* Locks scan_fixed_points_2d (PHASE6-UI): finds & classifies ALL
 * equilibria of a planar field. System x'=x-x^3, y'=-y has roots
 * (-1,0),(0,0),(1,0): two stable nodes and a saddle. make test-fp */
/* Verify scan_fixed_points_2d finds ALL equilibria and classifies them.
 * System (a classic with 3 fixed points):
 *   x' = x - x^3      (roots x = -1, 0, 1)
 *   y' = -y
 * Fixed points: (-1,0),(0,0),(1,0). Jacobian diag(1-3x^2, -1):
 *   at x=0:  eigs (+1,-1) -> saddle
 *   at x=+-1: eigs (-2,-1) -> stable node
 */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
int main(){
  PlanarField f;
  f.eval=[](double x,double y,double*u,double*v)->bool{ *u=x-x*x*x; *v=-y; return true; };
  auto pts=scan_fixed_points_2d(f,-3,3,-3,3,15);
  printf("found %zu fixed points\n", pts.size());
  int saddles=0, stable=0; bool gotm1=false,got0=false,gotp1=false;
  for(auto&p:pts){
    printf("  (% .4f,% .4f)  %s   eig: ", p.x,p.y,p.label.c_str());
    for(auto&e:p.eigenvalues) printf("%.3g%+.3gi ", e.real(),e.imag());
    printf(" dirs=%zu\n", p.directions.size());
    if(p.is_saddle) saddles++;
    if(p.label.find("stable")!=std::string::npos && p.label.find("un")==std::string::npos) stable++;
    if(std::fabs(p.x+1)<1e-3&&std::fabs(p.y)<1e-3) gotm1=true;
    if(std::fabs(p.x)<1e-3&&std::fabs(p.y)<1e-3) got0=true;
    if(std::fabs(p.x-1)<1e-3&&std::fabs(p.y)<1e-3) gotp1=true;
  }
  int fails=0;
  if(pts.size()!=3){printf("FAIL expected 3 points\n");fails++;}
  if(!(gotm1&&got0&&gotp1)){printf("FAIL missing a known root\n");fails++;}
  if(saddles!=1){printf("FAIL expected 1 saddle, got %d\n",saddles);fails++;}
  if(stable!=2){printf("FAIL expected 2 stable nodes, got %d\n",stable);fails++;}
  // saddle must have 2 real directions
  for(auto&p:pts) if(p.is_saddle && p.directions.size()!=2){printf("FAIL saddle missing eigendirections\n");fails++;}
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
