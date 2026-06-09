#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* BT normal form: x'=y, y'=b1 + b2*y + x^2 + x*y. p=b1, q=b2.
 * The Bogdanov-Takens point is EXACTLY at (b1,b2)=(0,0), equilibrium (0,0),
 * where J=[[0,1],[0,0]] has a double-zero eigenvalue. */
static bool bt(const double*X,double b1,double b2,double*o,std::string*){
  double x=X[0],y=X[1]; o[0]=y; o[1]=b1 + b2*y + x*x + x*y; return true;
}
int main(){
  Model2 m; m.n=2; m.vector_field=bt;
  // start from a guess NEAR but not at the BT point
  std::vector<double> x0={0.1,-0.05};
  Codim2Point R = locate_bogdanov_takens(m, x0, 0.08, -0.06);
  printf("locate_bogdanov_takens: ok=%d iters=%d residual=%.2e\n",R.ok,R.iters,R.residual);
  printf("  located: x=(%.5f,%.5f) p=b1=%.5f q=b2=%.5f (expect x=(0,0), b1=0, b2=0)\n",R.x[0],R.x[1],R.p,R.q);
  printf("  BT normal form: a=%.4f b=%.4f\n",R.a,R.b);
  printf("  msg: %s\n",R.message.c_str());
  bool pass = R.ok && std::fabs(R.p)<1e-4 && std::fabs(R.q)<1e-4 && std::hypot(R.x[0],R.x[1])<1e-4;
  printf("%s\n", pass?"bt_locate_smoke: all checks pass":"bt_locate_smoke: FAIL");
  return pass?0:1;
}
