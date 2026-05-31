/* Locks continue_equilibrium (PHASE D): the fold x'=p-x^2 is traced
 * through its saddle-node at (0,0); the Hopf normal form is detected at
 * p=0 with the stability flip. Also guards the exact-Jacobian null-buffer
 * crash fixed this iteration. make test-continuation */
/* Verify continue_equilibrium on known bifurcations before building UI. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;

int main(){
  int fails=0;

  /* (1) Saddle-node/fold: x' = p - x^2. Equilibria x=±sqrt(p) for p>=0,
     none for p<0. A FOLD at (x,p)=(0,0). 1-D system. */
  {
    Model m; m.n=1;
    m.vector_field=[](const double*x,double p,double*f,std::string*){ f[0]=p - x[0]*x[0]; return true; };
    m.jacobian_x =[](const double*x,double,double*J,std::string*){ J[0]=-2.0*x[0]; return true; };
    m.dfdp      =[](const double*,double,double*d,std::string*){ d[0]=1.0; return true; };
    ContinuationSettings s; s.p_min=-1; s.p_max=4; s.max_points=400; s.direction=-1;
    // start on the upper branch at p=1 and trace toward the fold
    std::vector<double> x0={1.0};
    Branch b=continue_equilibrium(m,x0,1.0,s);
    printf("(1) fold x'=p-x^2: ok=%d, %zu points, %zu special\n", b.ok, b.points.size(), b.special_indices.size());
    // expect at least one fold detected near p=0,x=0
    bool found_fold=false;
    for(auto i:b.special_indices){ if(b.points[i].special==SpecialPointKind::Fold &&
        std::fabs(b.points[i].p)<0.05 && std::fabs(b.points[i].x[0])<0.1) found_fold=true; }
    printf("    fold near (0,0): %s\n", found_fold?"YES":"no");
    if(!found_fold) fails++;
    // verify points satisfy x^2≈p along the branch
    double maxerr=0; for(auto&pt:b.points){ double e=std::fabs(pt.x[0]*pt.x[0]-pt.p); if(e>maxerr)maxerr=e; }
    printf("    max |x^2 - p| on branch = %.2e\n", maxerr);
    if(maxerr>1e-6) fails++;
  }

  /* (2) Hopf: planar system x'=p x - y - x(x^2+y^2), y'=x + p y - y(x^2+y^2).
     Equilibrium at origin; eigenvalues p±i. Hopf at p=0. */
  {
    Model m; m.n=2;
    m.vector_field=[](const double*x,double p,double*f,std::string*){
      double r2=x[0]*x[0]+x[1]*x[1];
      f[0]=p*x[0]-x[1]-x[0]*r2; f[1]=x[0]+p*x[1]-x[1]*r2; return true; };
    m.jacobian_x=[](const double*x,double p,double*J,std::string*){
      double X=x[0],Y=x[1];
      J[0]=p-3*X*X-Y*Y; J[1]=-1-2*X*Y;
      J[2]=1-2*X*Y;     J[3]=p-X*X-3*Y*Y; return true; };
    m.dfdp=[](const double*x,double,double*d,std::string*){ d[0]=x[0]; d[1]=x[1]; return true; };
    ContinuationSettings s; s.p_min=-1; s.p_max=1; s.max_points=400; s.direction=+1;
    std::vector<double> x0={0.0,0.0};
    Branch b=continue_equilibrium(m,x0,-0.5,s);
    printf("(2) Hopf normal form: ok=%d, %zu points, %zu special\n", b.ok, b.points.size(), b.special_indices.size());
    bool found_hopf=false;
    for(auto i:b.special_indices){ if(b.points[i].special==SpecialPointKind::Hopf && std::fabs(b.points[i].p)<0.05) found_hopf=true; }
    printf("    Hopf near p=0: %s\n", found_hopf?"YES":"no");
    if(!found_hopf) fails++;
    // stability should flip across p=0 (origin stable for p<0, unstable for p>0)
    int stab_neg=-1, stab_pos=-1;
    for(auto&pt:b.points){ if(pt.p<-0.2) stab_neg=pt.stable?1:0; if(pt.p>0.2) stab_pos=pt.stable?1:0; }
    printf("    stable at p<0: %d, stable at p>0: %d (expect 1 then 0)\n", stab_neg, stab_pos);
    if(stab_neg!=1||stab_pos!=0) fails++;
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
