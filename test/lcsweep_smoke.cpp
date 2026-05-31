/* Locks the limit-cycle SWEEP (Phase D step 2): on the supercritical Hopf
 * normal form the measured amplitude tracks 2*sqrt(mu) and there's no cycle
 * for mu<=0; on van der Pol the period increases monotonically with mu.
 * This is the pipeline behind the period/amplitude-vs-parameter diagram.
 * make test-lcsweep */
/* Prove the limit-cycle SWEEP pipeline on systems with known behavior,
 * BEFORE building UI. We integrate each parameter value with plain RK4,
 * settle, then measure period+amplitude with the validated detector. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

/* tiny fixed-step RK4 for a 2D system, matching what the app does */
template<class F>
static void rk4_2d(F f, double&x, double&y, double p, double dt){
  double k1x,k1y,k2x,k2y,k3x,k3y,k4x,k4y;
  f(x,y,p,k1x,k1y);
  f(x+0.5*dt*k1x, y+0.5*dt*k1y, p, k2x,k2y);
  f(x+0.5*dt*k2x, y+0.5*dt*k2y, p, k3x,k3y);
  f(x+dt*k3x, y+dt*k3y, p, k4x,k4y);
  x += dt/6*(k1x+2*k2x+2*k3x+k4x);
  y += dt/6*(k1y+2*k2y+2*k3y+k4y);
}

int main(){
  int fails=0;
  const double dt=0.01;

  /* (1) Supercritical Hopf normal form:
     x' = mu x - omega y - x(x^2+y^2)
     y' = omega x + mu y - y(x^2+y^2)
     For mu<=0: spiral to origin (no cycle). For mu>0: stable limit cycle
     of radius sqrt(mu) -> peak-to-peak amplitude = 2*sqrt(mu); period ~ 2pi/omega. */
  {
    const double omega=1.0;
    printf("(1) Supercritical Hopf (omega=1): amplitude should ~ 2*sqrt(mu), period ~ 2pi=6.283\n");
    double mus[]={-0.2, 0.05, 0.2, 0.5, 1.0};
    for(double mu: mus){
      auto f=[&](double x,double y,double p,double&dx,double&dy){
        double r2=x*x+y*y; dx=p*x-omega*y-x*r2; dy=omega*x+p*y-y*r2; };
      double x=0.1,y=0.0;
      for(int i=0;i<40000;i++) rk4_2d(f,x,y,mu,dt);   // settle
      std::vector<double> sig; sig.reserve(40000);
      for(int i=0;i<40000;i++){ rk4_2d(f,x,y,mu,dt); sig.push_back(x); }
      auto R=limit_cycle_period_amplitude(sig,dt);
      double expect_amp = mu>0 ? 2*std::sqrt(mu) : 0.0;
      printf("   mu=%+.2f: ok=%d period=%.3f amp=%.4f (expect amp~%.4f)\n",
             mu, R.ok, R.period, R.amplitude, expect_amp);
      if(mu<=0){
        if(R.ok && R.amplitude>0.05){printf("     FAIL: should be no cycle\n");fails++;}
      } else {
        if(!R.ok){printf("     FAIL: cycle not detected\n");fails++;}
        else {
          if(std::fabs(R.amplitude-expect_amp)>0.15){printf("     FAIL amplitude\n");fails++;}
          if(std::fabs(R.period-2*M_PI)>0.3){printf("     FAIL period\n");fails++;}
        }
      }
    }
  }

  /* (2) van der Pol vs mu: period grows with mu (relaxation oscillation).
     mu=0.5 -> T~6.38, mu=1 -> ~6.66, mu=2 -> ~7.6, mu=3 -> ~8.86. monotone up. */
  {
    printf("(2) van der Pol: period should INCREASE with mu\n");
    double prevT=0; bool monotone=true;
    for(double mu: {0.5,1.0,2.0,3.0}){
      auto f=[&](double x,double y,double p,double&dx,double&dy){ dx=y; dy=p*(1-x*x)*y-x; };
      double x=2,y=0;
      for(int i=0;i<40000;i++) rk4_2d(f,x,y,mu,dt);
      std::vector<double> sig; for(int i=0;i<60000;i++){ rk4_2d(f,x,y,mu,dt); sig.push_back(x); }
      auto R=limit_cycle_period_amplitude(sig,dt);
      printf("   mu=%.1f: ok=%d period=%.3f amp=%.3f\n", mu, R.ok, R.period, R.amplitude);
      if(!R.ok){fails++;}
      if(prevT>0 && R.period < prevT-0.1){monotone=false;}
      prevT=R.period;
    }
    if(!monotone){printf("   FAIL: period not monotone increasing\n");fails++;}
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
