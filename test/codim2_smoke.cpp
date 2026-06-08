#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static int fails=0;
static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
static const char* kindname(SpecialPointKind k){
  switch(k){case SpecialPointKind::Fold:return"Fold";case SpecialPointKind::Hopf:return"Hopf";
  case SpecialPointKind::BogdanovTakens:return"BT";case SpecialPointKind::Cusp:return"Cusp";
  case SpecialPointKind::GeneralizedHopf:return"GH";case SpecialPointKind::BranchPoint:return"BP";
  case SpecialPointKind::EndOfBranch:return"End";default:return"-";}
}
int main(){
  // ---- pitchfork: x' = p x - x^3. Trivial branch x=0; BP at (0,0). ----
  {
    Model m; m.n=1;
    m.vector_field=[](const double*x,double p,double*f,std::string*)->bool{ f[0]=p*x[0]-x[0]*x[0]*x[0]; return true; };
    ContinuationSettings s; s.p_min=-2; s.p_max=2; s.h0=0.05; s.max_points=400; s.direction=1;
    s.detect_fold=true; s.detect_hopf=true;
    // start on trivial branch at p=-1, x=0
    Branch b=continue_equilibrium(m, {0.0}, -1.0, s);
    printf("pitchfork branch: ok=%d points=%zu specials=%zu\n", b.ok, b.points.size(), b.special_indices.size());
    int nbp=0, idx_bp=-1;
    for(size_t i=0;i<b.special_indices.size();i++){ auto& sp=b.points[b.special_indices[i]];
      printf("   special @p=%.4f kind=%s\n", sp.p, kindname(sp.special));
      if(sp.special==SpecialPointKind::BranchPoint){nbp++; idx_bp=(int)b.special_indices[i];}
    }
    chk("pitchfork: a branch point detected near p=0", nbp>=1);
    if(idx_bp>=0){
      // switch branch -> should land on x=±sqrt(p) branch
      Branch b2=switch_branch(m, b.points[idx_bp], s);
      printf("   switched: ok=%d points=%zu msg=%s\n", b2.ok, b2.points.size(), b2.message.c_str());
      chk("branch switch traces the crossing branch", b2.ok && b2.points.size()>5);
      if(b2.ok){
        // on the nontrivial branch at p=1, |x| should be ~1 (x^2=p)
        double worst=0; int checked=0;
        for(auto&pt:b2.points){ if(std::fabs(pt.p-1.0)<0.1){ checked++; worst=std::max(worst,std::fabs(pt.x[0]*pt.x[0]-pt.p)); } }
        if(checked) printf("   nontrivial branch x^2≈p check: worst |x^2-p|=%.3f (n=%d)\n", worst, checked);
        chk("crossing branch satisfies x^2≈p", checked==0 || worst<0.2);
      }
    }
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
