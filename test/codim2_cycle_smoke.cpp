#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static double B=0.2;
static bool ross2(const double*X,double c,double a,double*o,std::string*){double x=X[0],y=X[1],z=X[2];o[0]=-y-z;o[1]=x+a*y;o[2]=B+z*(x-c);return true;}

/* mirrors the library's codim-2 bracketing (cycle_bif_curve) for a direct logic
 * check on a synthetic PD curve with known degeneracies */
struct SP { double p,q,fold_test,ns_test,amplitude; };
static bool synthetic_scan_ok(){
  std::vector<SP> pts={{3.0,0.20,+0.40,-0.30,0.50},{2.9,0.25,+0.20,-0.25,0.45},{2.8,0.30,+0.05,-0.20,0.40},
    {2.7,0.35,-0.15,-0.15,0.35},{2.6,0.40,-0.30,-0.08,0.30},{2.5,0.45,-0.45,+0.05,0.25},
    {2.4,0.50,-0.55,+0.20,0.20},{2.3,0.55,-0.62,+0.30,0.08},{2.2,0.60,-0.66,+0.38,0.03}};
  auto iq=[&](const SP&a,const SP&b,double fa,double fb){double t=fa/(fa-fb);return a.q+t*(b.q-a.q);};
  int nff=0,npd=0,nd=0; double qff=-1,qpd=-1,qd=-1;
  for(size_t i=1;i<pts.size();i++){const SP&a=pts[i-1];const SP&b=pts[i];
    if(a.fold_test*b.fold_test<0){nff++;qff=iq(a,b,a.fold_test,b.fold_test);}
    if(a.ns_test*b.ns_test<0){npd++;qpd=iq(a,b,a.ns_test,b.ns_test);}
    if(b.amplitude<0.05&&a.amplitude>=0.05){nd++;qd=b.q;}}
  return nff==1&&std::fabs(qff-0.3125)<0.02&&npd==1&&std::fabs(qpd-0.4308)<0.02&&nd==1&&std::fabs(qd-0.60)<1e-9;
}
int main(){
  // (1) real Rossler PD curve: no false codim-2
  Model2 m2; m2.n=3; m2.vector_field=ross2;
  double c0=3.0,a0=0.2,dt=0.002; std::vector<double> x={1,1,1};
  auto f=[&](const std::vector<double>&s,std::vector<double>&o){o.assign(3,0);double X[3]={s[0],s[1],s[2]};ross2(X,c0,a0,o.data(),0);};
  auto step=[&](std::vector<double>&s){std::vector<double> k1,k2,k3,k4,t(3);f(s,k1);for(int i=0;i<3;i++)t[i]=s[i]+0.5*dt*k1[i];f(t,k2);for(int i=0;i<3;i++)t[i]=s[i]+0.5*dt*k2[i];f(t,k3);for(int i=0;i<3;i++)t[i]=s[i]+dt*k3[i];f(t,k4);for(int i=0;i<3;i++)s[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;};
  for(int i=0;i<200000;i++)step(x);
  std::vector<std::vector<double>> loop;int cr=0;double T=0;
  for(int i=0;i<400000;i++){double prev=x[1];step(x);T+=dt;if(prev<x[1]&&x[1]>=0&&prev<0){cr++;if(cr==1){loop.clear();T=0;}else if(cr==2)break;}if(cr>=1)loop.push_back(x);}
  int mesh=90; std::vector<std::vector<double>> g(mesh);for(int i=0;i<mesh;i++){double ff=(double)i/mesh*(loop.size()-1);g[i]=loop[(int)ff];}
  TwoParamSettings st; st.p_min=2.0;st.p_max=6.0;st.q_min=0.19;st.q_max=0.24;st.max_points=16;
  CycleSettings cs; cs.mesh=mesh; cs.arclength=true; cs.ds=0.04; cs.max_steps=200; cs.compute_floquet=true; cs.adaptive_mesh=true;
  CycleBifCurve pc=pd_curve(m2,g,T,c0,a0,st,cs);
  printf("real Rossler PD curve: points=%zu codim2=%zu (expect codim2=0)\n",pc.points.size(),pc.codim2.size());
  bool no_fp = pc.ok && pc.points.size()>=3 && pc.codim2.size()==0;
  // (2) synthetic detection
  bool detect = synthetic_scan_ok();
  printf("synthetic codim-2 detection (FoldFlip/PDNS/DegeneratePD): %s\n", detect?"all located":"FAILED");
  bool pass = no_fp && detect;
  printf("%s\n", pass?"codim2_cycle_smoke: all checks pass":"codim2_cycle_smoke: FAIL");
  return pass?0:1;
}
