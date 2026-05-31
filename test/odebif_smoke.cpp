/* Locks the ODE bifurcation method (PHASE B): recording local maxima of
 * the observable reproduces the Rossler period-doubling route
 * (1 -> 2 -> 4 -> chaos peaks as c increases). make test-odebif */
/* Verify the ODE bifurcation "local maxima of observable" approach
 * produces a period-doubling structure for the Rossler system as c varies.
 * Rossler: x'=-y-z, y'=x+a y, z'=b+z(x-c). a=b=0.2.
 * Observable = x. As c increases, the number of distinct local maxima of x
 * should go 1 -> 2 -> 4 ... (period doubling). */
#include <cstdio>
#include <cmath>
#include <set>
#include <vector>
using V=std::vector<double>;
static void f(const V&s,double a,double b,double c,V&o){o.resize(3);o[0]=-s[1]-s[2];o[1]=s[0]+a*s[1];o[2]=b+s[2]*(s[0]-c);}
static void rk4(V&s,double a,double b,double c,double h){
  V k1,k2,k3,k4,t(3);
  f(s,a,b,c,k1);for(int i=0;i<3;i++)t[i]=s[i]+0.5*h*k1[i];
  f(t,a,b,c,k2);for(int i=0;i<3;i++)t[i]=s[i]+0.5*h*k2[i];
  f(t,a,b,c,k3);for(int i=0;i<3;i++)t[i]=s[i]+h*k3[i];
  f(t,a,b,c,k4);for(int i=0;i<3;i++)s[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}
static int count_maxima(double c){
  double a=0.2,b=0.2,h=0.02; V s={1,1,1};
  for(int i=0;i<200000;i++)rk4(s,a,b,c,h); // transient
  std::set<long> peaks;
  double w0=0,w1=0,w2=0; int have=0;
  for(int i=0;i<400000;i++){
    rk4(s,a,b,c,h);
    w0=w1;w1=w2;w2=s[0]; have=have<3?have+1:3;
    if(have==3 && w1>w0 && w1>=w2) peaks.insert((long)std::round(w1*100));
  }
  return (int)peaks.size();
}
int main(){
  int fails=0;
  // Known Rossler period-doubling route: c=2.5 (period-1), c=3.5 (period-2),
  // c=4.0 (period-4), c=5.0 (chaos). Counts are approximate (binned).
  struct T{double c; int lo; int hi; const char*lbl;};
  T t[]={
    {2.5, 1, 2, "period-1 (expect ~1 peak)"},
    {3.5, 2, 4, "period-2 (expect ~2)"},
    {4.0, 3, 8, "period-4 (expect ~4)"},
    {5.7, 10, 100000, "chaos (expect many)"},
  };
  for(auto&x:t){
    int n=count_maxima(x.c);
    bool ok = n>=x.lo && n<=x.hi;
    printf("Rossler c=%.1f: %d distinct maxima of x  [%s] %s\n", x.c, n, x.lbl, ok?"":"<-- FAIL");
    if(!ok)fails++;
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
