#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Construct a system with eigenvalues {0, +-i*omega} at the origin and KNOWN
 * quadratic coefficients. Linear part:
 *   x' = 0*x          (the zero direction; q0 = e1)
 *   y' = -omega*z
 *   z' = +omega*y     (the Hopf pair; eigenvalue +-i*omega)
 * Add quadratic terms in the x-equation and the (y,z) block so b and Re(c) are
 * predictable.
 *
 * With q0=(1,0,0): B(q0,q0) = second deriv in pure-x direction = the x^2-type
 * coefficients (times 2). The x-equation's x^2 coefficient is BX. The y,z eqs
 * have no x^2 term here, so B(q0,q0)=(2*BX, 0, 0).
 * Adjoint p0=(1,0,0) (since the x-row is decoupled at linear order). So
 *   b = <p0,B(q0,q0)> = 2*BX.
 * Set BX = 0.35 -> expect b = 0.70.
 *
 * For c = <p1,B(q0,q1)>: q1 is the complex eigenvector of the (y,z) block for
 * +i*omega: block [[0,-w],[w,0]], eigenvector for +i w is (1, -i)/... in (y,z).
 * So q1 = (0, 1, -i) (x-component 0). B(q0,q1) couples x with (y,z): we add a
 * term CX * x * y to the z-equation (mixed x-y). Then
 *   B(q0,q1)_z = d^2(z')/dx dy * q0_x * q1_y = CX * 1 * 1 = CX (real part).
 * p1 is adjoint of the (y,z) block; its z-component pairs with q1_z. We'll just
 * read Re(c) numerically and check it's CONSISTENT and STABLE, plus verify the
 * SIGN flips when we flip CX. (Exact p1 normalization makes the closed form
 * messy; the robust check is sign + reproducibility + b exact.) */
static double BX=0.35, CX=0.5, W=1.3;
static bool zh(const double*X,double mu,double*o,std::string*){(void)mu;
  double x=X[0],y=X[1],z=X[2];
  o[0]=0.0*x + BX*x*x;
  o[1]=-W*z + CX*x*y;
  o[2]= W*y;
  return true;
}
int main(){
  Model m; m.n=3; m.vector_field=zh;
  std::vector<double> x={0,0,0};
  double b=0,recc=0,om=0,s=0; std::string e;
  bool ok=zero_hopf_normal_form(m,x,0.0,&b,&recc,&om,&s,&e);
  printf("zero_hopf_normal_form: ok=%d\n",ok);
  if(!ok){printf("err: %s\n",e.c_str());return 1;}
  printf("  omega=%.4f (expect %.4f)\n",om,W);
  printf("  b=%.4f (expect 2*BX=%.4f)\n",b,2*BX);
  printf("  Re(c)=%.4f  s=b*Re(c)=%.4f\n",recc,s);
  // flip CX sign -> Re(c) should flip sign
  CX=-0.5; double b2,recc2,om2,s2;
  zero_hopf_normal_form(m,x,0.0,&b2,&recc2,&om2,&s2,&e);
  printf("  after CX->-CX: Re(c)=%.4f (should flip sign vs %.4f)\n",recc2,recc);
  bool pass = ok && std::fabs(om-W)<1e-3 && std::fabs(b-2*BX)<1e-2 && (recc*recc2<0);
  printf("%s\n", pass?"zero_hopf_nf_smoke: all checks pass":"zero_hopf_nf_smoke: FAIL");
  return pass?0:1;
}
