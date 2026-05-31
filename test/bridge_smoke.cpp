/* Verify the logistic <-> quadratic conjugacy that the 3D bridge relies on.
 * The logistic map x->r x(1-x) is conjugate to z->z^2+c with
 *   c = r/2 - r^2/4.
 * Known logistic bifurcation points and their Mandelbrot-needle images:
 *   r=3   (1->2 cycle)  => c=-0.75  (the period-1->2 boundary of M, c=-3/4)
 *   r=1+sqrt(6)=3.449.. (2->4)      => c=-1.25 (the period-2 disk boundary)
 *   r=4   (full chaos)              => c=-2    (tip of the needle)
 *   r=2   (fixed pt)                => c=0     (center of main cardioid)
 * We check c(r) lands on these. */
#include <cstdio>
#include <cmath>
static double c_of_r(double r){ return r/2.0 - r*r/4.0; }
static bool close(double a,double b,double t){ return std::fabs(a-b)<t; }
int main(){
  int fails=0;
  struct T{const char*name; double r; double c;};
  T t[]={
    {"r=2 (fixed pt) -> cardioid center c=0", 2.0, 0.0},
    {"r=3 (period 1->2) -> c=-0.75", 3.0, -0.75},
    {"r=1+sqrt(6) (period 2->4) -> c=-1.25", 1.0+std::sqrt(6.0), -1.25},
    {"r=4 (chaos, needle tip) -> c=-2", 4.0, -2.0},
  };
  for(auto&x:t){
    double c=c_of_r(x.r);
    bool ok=close(c,x.c,1e-9);
    printf("  %-44s c(r)=%.6f (expect %.6f) %s\n", x.name, c, x.c, ok?"":"<-- FAIL");
    if(!ok)fails++;
  }
  // Also verify the period-doubling onset matches between the two pictures:
  // logistic 2->4 at r=3.449490; quadratic 2->4 at c=-1.25. Already covered.
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
