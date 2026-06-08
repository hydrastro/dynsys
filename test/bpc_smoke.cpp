#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static bool vdp(const double*X,double mu,double*o,std::string*){o[0]=X[1];o[1]=mu*(1-X[0]*X[0])*X[1]-X[0];return true;}
static bool tc(const double*X,double mu,double*o,std::string*){double x=X[0],y=X[1],z=X[2];double r2=x*x+y*y;o[0]=-y+x*(1-r2);o[1]=x+y*(1-r2);o[2]=(mu-x*x)*z;return true;}
int main(){
  // 1) van der Pol: NO branch point
  { Model m; m.n=2; m.vector_field=vdp;
    int mesh=100; std::vector<std::vector<double>> g(mesh);
    for(int i=0;i<mesh;i++){double th=2*M_PI*i/mesh; g[i]={2*std::cos(th),2*std::sin(th)};}
    CycleSettings cs; cs.mesh=mesh; cs.arclength=true; cs.ds=0.1; cs.max_steps=60; cs.compute_floquet=true; cs.adaptive_mesh=true; cs.p_min=0.5; cs.p_max=3.0;
    CycleBranch B=continue_limit_cycle(m,g,6.66,1.0,cs);
    int nbp=0; for(auto&s:B.samples) if(s.is_bp) nbp++;
    printf("van der Pol: BPC count=%d (expect 0)\n",nbp);
    if(!(B.ok && nbp==0)){printf("bpc_smoke: FAIL\n");return 1;}
  }
  // 2) transcritical of cycles: BPC at mu=0.5 (transverse multiplier through +1)
  { Model m; m.n=3; m.vector_field=tc;
    int mesh=80; std::vector<std::vector<double>> g(mesh);
    for(int i=0;i<mesh;i++){double th=2*M_PI*i/mesh; g[i]={std::cos(th),std::sin(th),0.0};}
    CycleSettings cs; cs.mesh=mesh; cs.arclength=true; cs.ds=0.03; cs.max_steps=200; cs.compute_floquet=true; cs.adaptive_mesh=false; cs.p_min=0.2; cs.p_max=0.8;
    CycleBranch B=continue_limit_cycle(m,g,2*M_PI,0.3,cs);
    int nbp=0; double bpc=-1; for(auto&s:B.samples) if(s.is_bp){nbp++; bpc=s.bp_p;}
    printf("transcritical-of-cycles: BPC count=%d at mu=%.4f (expect 1 near 0.5)\n",nbp,bpc);
    if(!(B.ok && nbp>=1 && std::fabs(bpc-0.5)<0.06)){printf("bpc_smoke: FAIL\n");return 1;}
  }
  printf("bpc_smoke: all checks pass\n");
  return 0;
}
