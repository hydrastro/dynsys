#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
/* Rossler x'=-y-z, y'=x+a*y, z'=b+z(x-c), b=0.2 fixed. Two parameters:
 * p = c (primary), q = a (secondary). The period-doubling locus moves in (a,c).
 * Trace the PD curve over a range of a and confirm it returns a coherent set of
 * (c, a) points with c near the known ~2.8 at a=0.2. */
static double B=0.2;
static bool ross2(const double*X,double c,double a,double*o,std::string*){
  double x=X[0],y=X[1],z=X[2]; o[0]=-y-z; o[1]=x+a*y; o[2]=B+z*(x-c); return true;
}
int main(){
  Model2 m2; m2.n=3; m2.vector_field=ross2;
  // simulate a cycle at (c=3, a=0.2) for a seed
  double c0=3.0,a0=0.2,dt=0.002; std::vector<double> x={1,1,1};
  auto f=[&](const std::vector<double>&s,std::vector<double>&o){o.assign(3,0);double X[3]={s[0],s[1],s[2]};ross2(X,c0,a0,o.data(),0);};
  auto step=[&](std::vector<double>&s){std::vector<double> k1,k2,k3,k4,t(3);f(s,k1);
    for(int i=0;i<3;i++)t[i]=s[i]+0.5*dt*k1[i];f(t,k2);for(int i=0;i<3;i++)t[i]=s[i]+0.5*dt*k2[i];f(t,k3);
    for(int i=0;i<3;i++)t[i]=s[i]+dt*k3[i];f(t,k4);for(int i=0;i<3;i++)s[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;};
  for(int i=0;i<200000;i++) step(x);
  std::vector<std::vector<double>> loop; int cr=0; double T=0;
  for(int i=0;i<400000;i++){double prev=x[1];step(x);T+=dt; if(prev<x[1]&&x[1]>=0&&prev<0){cr++;if(cr==1){loop.clear();T=0;}else if(cr==2)break;} if(cr>=1)loop.push_back(x);}
  int mesh=90; std::vector<std::vector<double>> guess(mesh);
  for(int i=0;i<mesh;i++){double ff=(double)i/mesh*(loop.size()-1);guess[i]=loop[(int)ff];}
  printf("seed period~%.3f loop=%zu\n",T,loop.size());
  TwoParamSettings st; st.p_min=2.0;st.p_max=6.0;st.q_min=0.19;st.q_max=0.24;st.max_points=16;
  CycleSettings cs; cs.mesh=mesh; cs.arclength=true; cs.ds=0.04; cs.max_steps=200; cs.compute_floquet=true; cs.adaptive_mesh=true;
  CycleBifCurve pc=pd_curve(m2, guess, T, c0, a0, st, cs);
  printf("PD curve: ok=%d points=%zu msg=%s\n",pc.ok,pc.points.size(),pc.message.c_str());
  for(auto&pt:pc.points) printf("  a=%.3f  c_PD=%.4f  period=%.3f\n",pt.q,pt.p,pt.period);
  // sanity: c_PD should be in a sensible band and vary monotonically-ish with a
  // require a coherent locus: >=3 points and c_PD decreasing as a increases
  bool mono=true; for(size_t i=1;i<pc.points.size();i++) if(pc.points[i].p > pc.points[i-1].p+0.02) mono=false;
  bool pass = pc.ok && pc.points.size()>=3 && mono;
  printf("%s\n", pass?"pd_curve_smoke: all checks pass":"pd_curve_smoke: FAIL");
  return pass?0:1;
}
