#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static bool duff(const double*X,double p,double*o,std::string*){ (void)p; double x=X[0],y=X[1]; o[0]=y; o[1]=-x+x*x*x; return true; }
int main(){
  Model m; m.n=2; m.vector_field=duff;
  std::vector<double> x0={-1,0}, x1={1,0};
  int Np=300; std::vector<std::vector<double>> seed(Np);
  // seed slightly inside [-1,1] and stretch endpoints toward saddles along separatrix
  for(int i=0;i<Np;i++){ double x=-0.999+1.998*i/(Np-1); double y=(1-x*x)/std::sqrt(2.0); seed[i]={x,y}; }
  HomoclinicSettings st; st.mesh=160; st.newton_iters=80; st.free_T=true; st.T=4.0; // larger half-time => endpoints closer to saddles
  HeteroclinicResult R=solve_heteroclinic(m, x0, x1, 0.0, seed, st);
  printf("ok=%d residual=%.3e T=%.4f length=%.4f\n",R.ok,R.newton_residual,R.T,R.length);
  if(R.orbit.empty()){printf("HETERO: FAIL\n");return 1;}
  auto&a=R.orbit.front(); auto&b=R.orbit.back();
  double d0=std::hypot(a[0]-x0[0],a[1]-x0[1]), d1=std::hypot(b[0]-x1[0],b[1]-x1[1]);
  // unstable eigvec at (-1,0): J=[[0,1],[2,0]], lambda=+sqrt2, v=(1,sqrt2)/sqrt3
  double vx=1/std::sqrt(3.0), vy=std::sqrt(2.0)/std::sqrt(3.0);
  double sx=(a[0]-x0[0])/d0, sy=(a[1]-x0[1])/d0;
  double align0=std::fabs(sx*vx+sy*vy); // ~1 if start is along unstable eigvec
  // stable eigvec at (1,0): J=[[0,1],[2,0]], lambda=-sqrt2, v=(1,-sqrt2)/sqrt3
  double wx=1/std::sqrt(3.0), wy=-std::sqrt(2.0)/std::sqrt(3.0);
  double ex=(b[0]-x1[0])/d1, ey=(b[1]-x1[1])/d1;
  double align1=std::fabs(ex*wx+ey*wy);
  auto&mid=R.orbit[R.orbit.size()/2];
  printf("start(%.4f,%.4f) d0=%.4f align_unstable=%.4f\n",a[0],a[1],d0,align0);
  printf("end  (%.4f,%.4f) d1=%.4f align_stable=%.4f\n",b[0],b[1],d1,align1);
  printf("midpoint(%.4f,%.4f) expect (0,0.707)\n",mid[0],mid[1]);
  // correctness: converged + endpoints leave/arrive ALONG the right eigenvectors + correct midpoint
  bool pass = R.ok && align0>0.999 && align1>0.999 && std::fabs(mid[1]-0.7071)<0.02 && std::fabs(mid[0])<0.02;
  printf("%s\n", pass?"heteroclinic_smoke: all checks pass":"heteroclinic_smoke: FAIL");
  return pass?0:1;
}
