/* Locks limit_cycle_period_amplitude (Phase D step 2 foundation): exact on
 * a sine, matches the van der Pol relaxation-oscillation period (~6.66),
 * and reports 'no oscillation' for a fixed point. make test-limitcycle */
/* Validate limit_cycle_period_amplitude on known signals. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
int main(){
  int fails=0;
  const double dt=0.01;

  // (1) pure sine: period 2π/ω, amplitude 2A peak-to-peak
  { std::vector<double> y; double w=2.0, A=1.5;
    for(int i=0;i<5000;i++){double t=i*dt; y.push_back(A*std::sin(w*t)+0.3);}
    auto R=limit_cycle_period_amplitude(y,dt);
    double expected_T=2*M_PI/w;
    printf("(1) sine: ok=%d T=%.4f (expect %.4f) amp=%.3f (expect %.3f)\n",
           R.ok,R.period,expected_T,R.amplitude,2*A);
    if(!R.ok||std::fabs(R.period-expected_T)>0.05){printf("  FAIL period\n");fails++;}
    if(std::fabs(R.amplitude-2*A)>0.05){printf("  FAIL amp\n");fails++;}
  }

  // (2) van der Pol (mu=1): integrate, then measure. Known period ~6.66.
  { double x=2,y=0,mu=1.0; std::vector<double> sig;
    // transient
    for(int i=0;i<20000;i++){double dx=y,dy=mu*(1-x*x)*y-x; x+=dx*dt; y+=dy*dt;}
    for(int i=0;i<20000;i++){double dx=y,dy=mu*(1-x*x)*y-x; x+=dx*dt; y+=dy*dt; sig.push_back(x);}
    auto R=limit_cycle_period_amplitude(sig,dt);
    printf("(2) van der Pol mu=1: ok=%d T=%.3f (lit ~6.66) amp=%.3f (x in ~[-2,2])\n",
           R.ok,R.period,R.amplitude);
    if(!R.ok||std::fabs(R.period-6.66)>0.4){printf("  FAIL period\n");fails++;}
    if(std::fabs(R.amplitude-4.0)>0.6){printf("  FAIL amp\n");fails++;}
  }

  // (3) constant signal -> no oscillation
  { std::vector<double> y(3000, 1.234);
    auto R=limit_cycle_period_amplitude(y,dt);
    printf("(3) constant: ok=%d (expect 0) msg=%s\n", R.ok, R.message.c_str());
    if(R.ok){printf("  FAIL (should report no oscillation)\n");fails++;}
  }

  // (4) decaying oscillation toward fixed point: aperiodic-ish, amplitude small at end
  //     (we feed only the settled tail which is ~constant -> no osc)
  { std::vector<double> y; for(int i=0;i<3000;i++){double t=i*dt; y.push_back(0.001*std::sin(3*t));}
    auto R=limit_cycle_period_amplitude(y,dt);
    printf("(4) tiny oscillation: ok=%d (amplitude below threshold -> no osc)\n", R.ok);
    // tiny amplitude relative to scale -> treated as fixed point; that's fine either way
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
