/* Locks the escape-time logic used by the fractal view (PHASE C):
 * iterating z->z^2+c from 0 must classify known Mandelbrot points
 * (cardioid, period-2 bulb, Douady rabbit IN; c=1,0.5,-2.5,2i OUT).
 * make test-fractal */
/* Verify the escape-time logic computes the actual Mandelbrot set:
 * iterate z->z^2+c from z0=0, check membership for known points. */
#include <cmath>
#include <cstdio>
struct C{double re,im;};
// returns iterations to escape (maxit if bounded), and final r2
static int escape(C c,int maxit,double R2,double*r2out){
  double x=0,y=0; int it=0; double r2=0;
  for(;it<maxit;++it){
    double nx=x*x-y*y+c.re;
    double ny=2*x*y+c.im;
    x=nx; y=ny; r2=x*x+y*y;
    if(r2>R2){++it;break;}
  }
  *r2out=r2; return it;
}
int main(){
  const int M=500; const double R2=16.0;
  struct T{const char*name;C c;bool in;};
  T tests[]={
    {"c=0 (origin, in)",{0,0},true},
    {"c=-1 (period-2 bulb, in)",{-1,0},true},
    {"c=-0.5 (cardioid, in)",{-0.5,0},true},
    {"c=0.25 (cusp edge, in-ish)",{0.24,0},true},
    {"c=1 (out)",{1,0},false},
    {"c=0.5 (out)",{0.5,0},false},
    {"c=-2.5 (out)",{-2.5,0},false},
    {"c=2i (out)",{0,2},false},
    {"c=-0.123+0.745i Douady rabbit (in)",{-0.123,0.745},true},
  };
  int fails=0;
  for(auto&t:tests){
    double r2; int it=escape(t.c,M,R2,&r2);
    bool isin=(it>=M);
    printf("  %-40s it=%4d %s (expect %s) %s\n",t.name,it, isin?"IN ":"OUT",
           t.in?"IN":"OUT", isin==t.in?"":"<-- FAIL");
    if(isin!=t.in)fails++;
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
