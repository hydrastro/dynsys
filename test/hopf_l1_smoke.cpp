/* Hopf first Lyapunov coefficient on known normal forms: supercritical (l1<0)
 * and subcritical (l1>0). The SIGN is the criticality verdict MatCont reports.
 * make test-hopfl1 */
/* Verify the first Lyapunov coefficient on known Hopf normal forms.
 * Supercritical:  x'=-y+x(mu-(x^2+y^2)),  y'= x+y(mu-(x^2+y^2))  -> l1<0
 * Subcritical:    x'=-y+x(mu+(x^2+y^2)),  y'= x+y(mu+(x^2+y^2))  -> l1>0
 * at mu=0 the equilibrium is the origin with eigenvalues +/- i. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
using namespace dynsys::analysis;
static int fails=0;
static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL"); if(!c)fails++;}

int main(){
  // supercritical
  Model sup; sup.n=2;
  sup.vector_field=[](const double*x,double,double*f,std::string*)->bool{
    double r2=x[0]*x[0]+x[1]*x[1];
    f[0]=-x[1]+x[0]*(0.0 - r2);
    f[1]= x[0]+x[1]*(0.0 - r2);
    return true; };
  double l1,w; std::string err;
  bool ok = hopf_first_lyapunov(sup, {0.0,0.0}, 0.0, &l1, &w, &err);
  printf("supercritical: ok=%d l1=%.4f omega=%.4f err=%s\n", ok, ok?l1:0, ok?w:0, err.c_str());
  chk("supercritical solved", ok);
  chk("supercritical l1 < 0", ok && l1 < -0.01);
  chk("omega ~ 1", ok && std::fabs(w-1.0)<0.05);

  // subcritical (sign of cubic flipped)
  Model sub; sub.n=2;
  sub.vector_field=[](const double*x,double,double*f,std::string*)->bool{
    double r2=x[0]*x[0]+x[1]*x[1];
    f[0]=-x[1]+x[0]*(0.0 + r2);
    f[1]= x[0]+x[1]*(0.0 + r2);
    return true; };
  ok = hopf_first_lyapunov(sub, {0.0,0.0}, 0.0, &l1, &w, &err);
  printf("subcritical: ok=%d l1=%.4f\n", ok, ok?l1:0);
  chk("subcritical l1 > 0", ok && l1 > 0.01);

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
