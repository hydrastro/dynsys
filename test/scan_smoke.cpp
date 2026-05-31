/* Locks the per-point largest-Lyapunov estimate used by the 2-parameter
 * scan (PHASE B): Henon is chaotic (lambda>0) at a=1.4,b=0.3 and
 * periodic (lambda<0) at a=1.0 and a=0.2. make test-scan */
/* Validate the per-point largest-Lyapunov estimate used by the 2-param
 * scan, on the Henon map at known (a,b) points. */
#include <cstdio>
#include <cmath>
#include <vector>
using V=std::vector<double>;
// Henon step
static void henon(const V&s,double a,double b,V&o){o.resize(2);o[0]=1-a*s[0]*s[0]+s[1];o[1]=b*s[0];}
// largest Lyapunov via shadow orbit (mirrors compute_scan_image)
static double lyap(double a,double b,int trans,int iters){
  double leps=1e-8; V s={0.1,0.0},sh;
  for(int k=0;k<trans;k++){V o;henon(s,a,b,o);s=o;}
  sh=s; sh[0]+=leps;
  double sum=0; long n=0;
  for(int k=0;k<iters;k++){
    V n1,n2; henon(s,a,b,n1); henon(sh,a,b,n2);
    double d2=0; for(int q=0;q<2;q++){double d=n2[q]-n1[q];d2+=d*d;}
    double dist=std::sqrt(d2);
    if(dist>1e-300&&std::isfinite(dist)){sum+=std::log(dist/leps);n++;
      double sc=leps/dist; sh.resize(2); for(int q=0;q<2;q++)sh[q]=n1[q]+(n2[q]-n1[q])*sc;}
    s=n1;
  }
  return n>0?sum/n:NAN;
}
int main(){
  int fails=0;
  // chaotic Henon a=1.4,b=0.3 -> lambda ~ +0.42
  double lc=lyap(1.4,0.3,1000,5000);
  printf("Henon a=1.4 b=0.3 (chaotic): lambda=%.4f (expect ~+0.42, must be >0)\n", lc);
  if(!(lc>0.2)){printf("  FAIL: should be positive\n");fails++;}
  // periodic window a=1.0,b=0.3 -> stable cycle, lambda < 0
  double lp=lyap(1.0,0.3,2000,5000);
  printf("Henon a=1.0 b=0.3 (periodic): lambda=%.4f (must be <0)\n", lp);
  if(!(lp<0)){printf("  FAIL: should be negative\n");fails++;}
  // another periodic point a=0.2 -> fixed point, strongly negative
  double lf=lyap(0.2,0.3,2000,5000);
  printf("Henon a=0.2 b=0.3 (fixed pt): lambda=%.4f (must be <0)\n", lf);
  if(!(lf<0)){printf("  FAIL: should be negative\n");fails++;}
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
