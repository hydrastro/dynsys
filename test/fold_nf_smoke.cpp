#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static int fails=0;
static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  // 1) scalar fold x' = p + x^2 ; at fold p=0,x=0 -> a = 1/2 * f'' = 1
  {
    Model m; m.n=1;
    m.vector_field=[](const double*x,double p,double*f,std::string*)->bool{ f[0]=p + x[0]*x[0]; return true; };
    std::vector<double> x={0.0}; double a=0,l0=0; std::string e;
    bool ok=fold_normal_form(m,x,0.0,&a,&l0,&e);
    printf("scalar fold: ok=%d a=%.4f lambda0=%.2e (expect a=1)\n", ok,a,l0);
    chk("scalar fold a == 1", ok && std::fabs(a-1.0)<1e-2);
    chk("lambda0 ~ 0 at fold", ok && l0<1e-3);
  }
  // 2) cusp-ish: x' = p + x^3 has a=0 at origin (no quadratic term) -> a ~ 0
  {
    Model m; m.n=1;
    m.vector_field=[](const double*x,double p,double*f,std::string*)->bool{ f[0]=p + x[0]*x[0]*x[0]; return true; };
    std::vector<double> x={0.0}; double a=0,l0=0; std::string e;
    bool ok=fold_normal_form(m,x,0.0,&a,&l0,&e);
    printf("cubic (cusp): ok=%d a=%.4f (expect ~0)\n", ok,a);
    chk("cubic a ~ 0 (cusp signal)", ok && std::fabs(a)<1e-2);
  }
  // 3) planar saddle-node: x'=p+x^2, y'=-y ; fold in x-direction, a=1
  {
    Model m; m.n=2;
    m.vector_field=[](const double*x,double p,double*f,std::string*)->bool{ f[0]=p+x[0]*x[0]; f[1]=-x[1]; return true; };
    std::vector<double> x={0.0,0.0}; double a=0,l0=0; std::string e;
    bool ok=fold_normal_form(m,x,0.0,&a,&l0,&e);
    printf("planar SN: ok=%d a=%.4f lambda0=%.2e (expect |a|=1)\n", ok,a,l0);
    chk("planar saddle-node |a| == 1", ok && std::fabs(std::fabs(a)-1.0)<1e-2);
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
