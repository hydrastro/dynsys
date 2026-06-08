#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;

/* A self-contained validation of dynsys's analysis core against systems with
 * KNOWN analytic results — the kind of head-to-head one would run vs MatCont. */

static int pass=0, total=0;
static void check(const char* name, double got, double expect, double tol, const char* unit=""){
  total++; bool ok=std::fabs(got-expect)<tol; if(ok)pass++;
  printf("  [%s] %-42s got %.6g  expect %.6g %s  %s\n", ok?"PASS":"FAIL", name, got, expect, unit, ok?"":"  <-- MISMATCH");
}

int main(){
  printf("=== dynsys analysis validation against analytic results ===\n\n");

  /* 1. Lorenz: classic equilibria + eigenvalues at the origin (sigma=10,
   * rho=28, beta=8/3). Origin eigenvalues: beta and the roots of
   * lambda^2 + (sigma+1)lambda + sigma(1-rho) = 0. */
  printf("1. Lorenz system (sigma=10, rho=28, beta=8/3):\n");
  {
    double sig=10, rho=28, beta=8.0/3.0;
    std::vector<double> A={-sig,sig,0,  rho,-1,0,  0,0,-beta}; // Jacobian at origin
    std::vector<Complex> ev; eigenvalues(A,3,&ev);
    // expected: -beta=-2.667, and (-(sig+1)+-sqrt((sig+1)^2-4 sig(1-rho)))/2
    double disc=std::sqrt((sig+1)*(sig+1)-4*sig*(1-rho));
    double l1=(-(sig+1)+disc)/2, l2=(-(sig+1)-disc)/2;
    // find the largest real eigenvalue
    double maxre=-1e9,minre=1e9; for(auto&z:ev){maxre=std::max(maxre,z.real());minre=std::min(minre,z.real());}
    check("origin unstable eigenvalue", maxre, l1, 1e-3);
    check("origin most-negative eigenvalue", minre, l2, 1e-3);
  }

  /* 2. Hopf normal form first Lyapunov coefficient. x'=-y+x(mu-(x^2+y^2)),
   * y'=x+y(mu-(x^2+y^2)). At mu=0 this is a supercritical Hopf with
   * l1 = -1 (negative => supercritical). dynsys computes l1 sign/magnitude. */
  printf("\n2. Supercritical Hopf normal form (l1 should be < 0):\n");
  {
    Model m; m.n=2;
    m.vector_field=[](const double*X,double mu,double*o,std::string*){
      double x=X[0],y=X[1],r2=x*x+y*y; o[0]=-y+x*(mu-r2); o[1]=x+y*(mu-r2); return true; };
    std::vector<double> x0={0,0}; double l1=0,omega=0; std::string e;
    bool ok=hopf_first_lyapunov(m,x0,0.0,&l1,&omega,&e);
    printf("  hopf_first_lyapunov ok=%d l1=%.6g omega=%.6g (expect l1<0, omega~1)\n",ok,l1,omega);
    total++; if(ok && l1<0) pass++; else printf("    <-- expected supercritical (l1<0)\n");
    check("Hopf frequency omega", std::fabs(omega), 1.0, 1e-2);
  }

  /* 3. Van der Pol Lyapunov exponents (mu=1): one ~0 (along the cycle), the
   * spectrum sums to the average divergence. For a planar limit cycle the
   * largest exponent is ~0. */
  printf("\n3. Van der Pol (mu=1) Lyapunov spectrum (largest ~ 0 on the cycle):\n");
  {
    Model m; m.n=2;
    m.vector_field=[](const double*X,double mu,double*o,std::string*){
      o[0]=X[1]; o[1]=mu*(1-X[0]*X[0])*X[1]-X[0]; return true; };
    LyapunovOptions ls; ls.transient=20000; ls.steps=200000; ls.dt=0.01;
    std::vector<double> x0={2.0,0.0};
    LyapunovResult R=lyapunov_spectrum(m,x0,1.0,ls);
    printf("  exponents:"); for(double e:R.exponents) printf(" %.4f",e); printf("\n");
    if(!R.exponents.empty()) check("largest Lyapunov exponent", R.exponents[0], 0.0, 0.05);
  }

  /* 4. Homoclinic orbit (x'=y, y'=x-x^2): exact peak 1.5 (validated already,
   * repeated here as part of the benchmark suite). */
  printf("\n4. Homoclinic orbit x'=y,y'=x-x^2 (analytic peak 1.5):\n");
  {
    Model m; m.n=2; m.vector_field=[](const double*X,double,double*o,std::string*){o[0]=X[1];o[1]=X[0]-X[0]*X[0];return true;};
    double Tt=8; int Np=300; std::vector<std::vector<double>> seed(Np,std::vector<double>(2));
    for(int i=0;i<Np;i++){double t=-Tt+2*Tt*i/(Np-1);double s=1.0/std::cosh(t/2);seed[i][0]=1.3*std::exp(-0.35*t*t);seed[i][1]=1.3*(-0.7*t)*std::exp(-0.35*t*t);}
    HomoclinicSettings hs; hs.mesh=200; hs.T=Tt; hs.free_T=false; hs.newton_iters=120;
    std::vector<double> x0={0,0};
    HomoclinicResult R=solve_homoclinic(m,x0,0.0,seed,hs);
    double peak=0; for(auto&p:R.orbit)peak=std::max(peak,p[0]);
    check("homoclinic peak amplitude", peak, 1.5, 0.02);
  }

  printf("\nvalidation_smoke: %d/%d checks passed\n", pass, total);
  if(pass==total) printf("validation_smoke: all checks pass\n"); else printf("validation_smoke: FAIL\n");
  return (pass==total)?0:1;
}
