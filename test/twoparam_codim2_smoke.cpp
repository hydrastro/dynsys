#include "analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
static const char* kn(SpecialPointKind k){switch(k){case SpecialPointKind::Cusp:return"Cusp";case SpecialPointKind::BogdanovTakens:return"BT";case SpecialPointKind::GeneralizedHopf:return"GH";default:return"-";}}
int main(){
  // cusp normal form x'=q+px-x^3: fold curve is 27q^2=4p^3, with a CUSP at (0,0).
  Model2 m; m.n=1;
  m.vector_field=[](const double*x,double p,double q,double*f,std::string*)->bool{
    f[0]=q + p*x[0] - x[0]*x[0]*x[0]; return true; };
  TwoParamSettings s; s.p_min=-0.5; s.p_max=6; s.q_min=-6; s.q_max=6; s.h0=0.04; s.max_points=800;
  // start at p=3,x=1,q=-2 (on one fold branch) and trace toward the cusp
  TwoParamCurve c = two_param_curve(m, TwoParamKind::Fold, {1.0}, 3.0, -2.0, s);
  printf("fold curve: ok=%d points=%zu specials=%zu msg=%s\n", c.ok, c.points.size(), c.special_indices.size(), c.message.c_str());
  bool sawCusp=false;
  for(auto idx:c.special_indices){auto&sp=c.points[idx];
    printf("   codim-2 @ (p=%.4f,q=%.4f) kind=%s\n", sp.p, sp.q, kn(sp.special));
    if(sp.special==SpecialPointKind::Cusp && std::fabs(sp.p)<0.5 && std::fabs(sp.q)<0.5) sawCusp=true;
  }
  chk("fold curve traced", c.ok);
  chk("cusp detected near the origin", sawCusp);
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
