/* Locks the integrator tableaux (PHASE B) against exact solutions:
 * exponential decay (Euler 1st order; RK4/RK38 ~machine precision),
 * DOPRI45 embedded error ~ h^5 (32x per halving), RK4 energy
 * conservation on the harmonic oscillator. make test-solver */
#include <string>
/* Validate the integrator tableaux against known exact solutions.
 * Mirrors the math in step_ode_state / step_ode_adaptive exactly. */
#include <cstdio>
#include <cmath>
#include <vector>
#include <functional>
using V=std::vector<double>;
using F=std::function<void(const V&,double,V&)>;

static void euler(const V&y,double t,double h,const F&f,V&o){V k(y.size());f(y,t,k);o=y;for(size_t i=0;i<y.size();++i)o[i]+=h*k[i];}
static void rk4(const V&y,double t,double h,const F&f,V&o){
  size_t n=y.size();V k1(n),k2(n),k3(n),k4(n),tmp(n);
  f(y,t,k1);for(size_t i=0;i<n;++i)tmp[i]=y[i]+0.5*h*k1[i];
  f(tmp,t+0.5*h,k2);for(size_t i=0;i<n;++i)tmp[i]=y[i]+0.5*h*k2[i];
  f(tmp,t+0.5*h,k3);for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*k3[i];
  f(tmp,t+h,k4);o=y;for(size_t i=0;i<n;++i)o[i]+=h/6*(k1[i]+2*k2[i]+2*k3[i]+k4[i]);
}
static void rk38(const V&y,double t,double h,const F&f,V&o){
  size_t n=y.size();V k1(n),k2(n),k3(n),k4(n),tmp(n);
  f(y,t,k1);for(size_t i=0;i<n;++i)tmp[i]=y[i]+h/3*k1[i];
  f(tmp,t+h/3,k2);for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(-1.0/3*k1[i]+k2[i]);
  f(tmp,t+2*h/3,k3);for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(k1[i]-k2[i]+k3[i]);
  f(tmp,t+h,k4);o=y;for(size_t i=0;i<n;++i)o[i]+=h/8*(k1[i]+3*k2[i]+3*k3[i]+k4[i]);
}
// DOPRI45 single fixed step, returns 5th order solution + error estimate
static double dopri(const V&y,double t,double h,const F&f,V&o){
  size_t n=y.size();std::vector<V>k(7,V(n));V tmp(n),y5(n),y4(n);
  f(y,t,k[0]);
  for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(0.2*k[0][i]);f(tmp,t+0.2*h,k[1]);
  for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(3.0/40*k[0][i]+9.0/40*k[1][i]);f(tmp,t+0.3*h,k[2]);
  for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(44.0/45*k[0][i]-56.0/15*k[1][i]+32.0/9*k[2][i]);f(tmp,t+0.8*h,k[3]);
  for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(19372.0/6561*k[0][i]-25360.0/2187*k[1][i]+64448.0/6561*k[2][i]-212.0/729*k[3][i]);f(tmp,t+8.0/9*h,k[4]);
  for(size_t i=0;i<n;++i)tmp[i]=y[i]+h*(9017.0/3168*k[0][i]-355.0/33*k[1][i]+46732.0/5247*k[2][i]+49.0/176*k[3][i]-5103.0/18656*k[4][i]);f(tmp,t+h,k[5]);
  for(size_t i=0;i<n;++i)y5[i]=y[i]+h*(35.0/384*k[0][i]+500.0/1113*k[2][i]+125.0/192*k[3][i]-2187.0/6784*k[4][i]+11.0/84*k[5][i]);
  f(y5,t+h,k[6]);
  for(size_t i=0;i<n;++i)y4[i]=y[i]+h*(5179.0/57600*k[0][i]+7571.0/16695*k[2][i]+393.0/640*k[3][i]-92097.0/339200*k[4][i]+187.0/2100*k[5][i]+1.0/40*k[6][i]);
  o=y5;double e=0;for(size_t i=0;i<n;++i)e=std::max(e,std::fabs(y5[i]-y4[i]));return e;
}
static bool close(double a,double b,double t){return std::fabs(a-b)<t;}

int main(){
  int fails=0;
  // decay x'=-x, exact e^-t
  F decay=[](const V&y,double,V&k){k[0]=-y[0];};
  // integrate to t=2 with many small steps; check vs e^-2
  double target=std::exp(-2.0);
  for(const char*nm:{"euler","rk4","rk38"}){
    V y={1.0};double t=0,h=0.001;
    for(int i=0;i<2000;++i){V o;
      if(std::string(nm)=="euler")euler(y,t,h,decay,o);
      else if(std::string(nm)=="rk4")rk4(y,t,h,decay,o);
      else rk38(y,t,h,decay,o);
      y=o;t+=h;}
    double e=std::fabs(y[0]-target);
    printf("decay %-5s -> %.8f (exact %.8f) err=%.2e\n",nm,y[0],target,e);
    // RK4/RK38 should be very accurate; euler less so but converging
    if(std::string(nm)!="euler"&&e>1e-7){printf("  FAIL accuracy\n");fails++;}
    if(std::string(nm)=="euler"&&e>1e-2){printf("  FAIL euler sanity\n");fails++;}
  }
  // DOPRI embedded error must shrink ~h^5: halving h should cut error ~32x
  {
    V y={1.0},o1,o2; double e1=dopri(y,0,0.1,decay,o1);
    double e2=dopri(y,0,0.05,decay,o2);
    double ratio=e1/std::max(e2,1e-300);
    printf("DOPRI err(h=0.1)=%.3e err(h=0.05)=%.3e ratio=%.1f (expect ~32 for order-5 embedded)\n",e1,e2,ratio);
    if(ratio<16||ratio>64){printf("  FAIL order\n");fails++;}
  }
  // Harmonic oscillator energy: RK4 should nearly conserve over one period
  {
    F osc=[](const V&y,double,V&k){k[0]=y[1];k[1]=-y[0];};
    V y={1.0,0.0};double t=0,h=0.01;
    for(int i=0;i<629;++i){V o;rk4(y,t,h,osc,o);y=o;t+=h;} // ~2pi
    double E=y[0]*y[0]+y[1]*y[1];
    printf("oscillator energy after ~1 period (RK4): %.8f (expect ~1.0)\n",E);
    if(!close(E,1.0,1e-3)){printf("  FAIL energy\n");fails++;}
  }
  printf("=== %s ===\n",fails==0?"PASS":"FAIL");
  return fails;
}
