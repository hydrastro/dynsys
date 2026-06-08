#include <cstdio>
#include <cmath>
// Verify the logistic<->Mandelbrot conjugacy c = r/2 - r^2/4 places the
// period-doubling cascade at the correct real-axis c-values.
// Known logistic bifurcation r-values and their Mandelbrot real-axis c:
//   r=3   (period 1->2)  -> c = -0.75   (boundary of main cardioid)
//   r=1+sqrt(6)=3.449    (2->4)         -> c = -1.25  (boundary of period-2 disk)
//   r=3.5699 (onset chaos)              -> c ~ -1.401 (Feigenbaum point, Myrberg)
//   r=4   (tip)                         -> c = -2.0
static double c_of_r(double r){ return r/2.0 - r*r/4.0; }
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}
int main(){
  printf("r=3.0    -> c=%.4f (expect -0.75)\n", c_of_r(3.0));
  chk("period-1->2 lands at c=-0.75", std::fabs(c_of_r(3.0)+0.75)<1e-6);
  double r2=1+std::sqrt(6.0);
  printf("r=%.4f -> c=%.4f (expect -1.25)\n", r2, c_of_r(r2));
  chk("period-2->4 lands at c=-1.25", std::fabs(c_of_r(r2)+1.25)<1e-6);
  printf("r=3.5699 -> c=%.4f (expect ~ -1.401, Feigenbaum)\n", c_of_r(3.5699));
  chk("chaos onset near c=-1.401", std::fabs(c_of_r(3.5699)+1.401)<0.01);
  printf("r=4.0    -> c=%.4f (expect -2.0)\n", c_of_r(4.0));
  chk("antenna tip at c=-2.0", std::fabs(c_of_r(4.0)+2.0)<1e-6);
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
