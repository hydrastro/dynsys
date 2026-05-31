/* Locks lyapunov_spectrum + Kaplan-Yorke (PHASE B) against known
 * systems: stable linear (exps -1,-2,-3, D=0), Henon map (sum=ln|detJ|,
 * D~1.26), Lorenz (0.906,0,-14.57, D~2.06). make test-lyap */
/* Verify lyapunov_spectrum + Kaplan-Yorke against known systems. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

static Model lorenz(){
  Model m; m.n=3;
  m.vector_field=[](const double*x,double,double*f,std::string*)->bool{
    const double s=10,r=28,b=8.0/3.0;
    f[0]=s*(x[1]-x[0]); f[1]=x[0]*(r-x[2])-x[1]; f[2]=x[0]*x[1]-b*x[2];
    return true;
  };
  m.jacobian_x=[](const double*x,double,double*J,std::string*)->bool{
    const double s=10,r=28,b=8.0/3.0;
    J[0]=-s;       J[1]=s;     J[2]=0;
    J[3]=r-x[2];   J[4]=-1;    J[5]=-x[0];
    J[6]=x[1];     J[7]=x[0];  J[8]=-b;
    return true;
  };
  return m;
}
static Model henon(){
  Model m; m.n=2;
  m.vector_field=[](const double*x,double,double*f,std::string*)->bool{
    const double a=1.4,b=0.3;
    f[0]=1-a*x[0]*x[0]+x[1]; f[1]=b*x[0]; return true;
  };
  m.jacobian_x=[](const double*x,double,double*J,std::string*)->bool{
    const double a=1.4,b=0.3;
    J[0]=-2*a*x[0]; J[1]=1; J[2]=b; J[3]=0; return true;
  };
  return m;
}
static Model stable_linear(){ // x'=-1 x, y'=-2 y, z'=-3 z -> exps -1,-2,-3
  Model m; m.n=3;
  m.vector_field=[](const double*x,double,double*f,std::string*)->bool{
    f[0]=-1*x[0]; f[1]=-2*x[1]; f[2]=-3*x[2]; return true; };
  m.jacobian_x=[](const double*,double,double*J,std::string*)->bool{
    for(int i=0;i<9;i++)J[i]=0; J[0]=-1;J[4]=-2;J[8]=-3; return true; };
  return m;
}

static bool close(double a,double b,double tol){ return std::fabs(a-b)<tol; }

int main(){
  int fails=0;
  printf("=== Lyapunov spectrum tests ===\n");

  // 1) stable linear: exact exponents -1,-2,-3, KY dimension 0
  {
    LyapunovOptions o; o.is_map=false; o.dt=0.01; o.transient=0; o.steps=20000; o.reorth_every=1;
    auto r=lyapunov_spectrum(stable_linear(), {1,1,1}, 0, o);
    printf("stable linear: [%.4f %.4f %.4f] KY=%.3f\n", r.exponents[0],r.exponents[1],r.exponents[2], r.kaplan_yorke);
    if(!close(r.exponents[0],-1,0.05)||!close(r.exponents[1],-2,0.05)||!close(r.exponents[2],-3,0.05)){printf("  FAIL exps\n");fails++;}
    if(r.kaplan_yorke!=0.0){printf("  FAIL KY (point attractor should be 0)\n");fails++;}
  }

  // 2) Henon map: lambda ~ (0.419, -1.62), sum = ln(0.3) = -1.204
  {
    LyapunovOptions o; o.is_map=true; o.transient=1000; o.steps=100000; o.reorth_every=1;
    auto r=lyapunov_spectrum(henon(), {0.1,0.0}, 0, o);
    printf("Henon: [%.4f %.4f] sum=%.4f (expect sum=ln0.3=%.4f) KY=%.4f\n",
           r.exponents[0],r.exponents[1], r.sum, std::log(0.3), r.kaplan_yorke);
    if(!close(r.exponents[0],0.419,0.02)){printf("  FAIL lambda1 (expect ~0.419)\n");fails++;}
    if(!close(r.sum,std::log(0.3),0.01)){printf("  FAIL sum != ln|det J|\n");fails++;}
    // KY for Henon ~ 1.26
    if(!close(r.kaplan_yorke,1.26,0.05)){printf("  FAIL KY (expect ~1.26)\n");fails++;}
  }

  // 3) Lorenz: lambda ~ (0.906, 0, -14.57), KY ~ 2.06
  {
    LyapunovOptions o; o.is_map=false; o.dt=0.005; o.transient=5000; o.steps=200000; o.reorth_every=1;
    auto r=lyapunov_spectrum(lorenz(), {1,1,1}, 0, o);
    printf("Lorenz: [%.4f %.4f %.4f] KY=%.4f (expect ~0.906,0,-14.57; KY~2.06)\n",
           r.exponents[0],r.exponents[1],r.exponents[2], r.kaplan_yorke);
    if(!close(r.exponents[0],0.906,0.05)){printf("  FAIL lambda1 (expect ~0.906)\n");fails++;}
    if(!close(r.exponents[1],0.0,0.03)){printf("  FAIL lambda2 (expect ~0)\n");fails++;}
    if(!close(r.exponents[2],-14.57,0.2)){printf("  FAIL lambda3 (expect ~-14.57)\n");fails++;}
    if(!close(r.kaplan_yorke,2.06,0.03)){printf("  FAIL KY (expect ~2.06)\n");fails++;}
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
