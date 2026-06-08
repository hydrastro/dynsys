#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static double A=0.2,B=0.2,C;
static bool ross(const double*X,double c,double*o,std::string*){double x=X[0],y=X[1],z=X[2];o[0]=-y-z;o[1]=x+A*y;o[2]=B+z*(x-c);return true;}
int main(){
  Model m; m.n=3; m.vector_field=ross;
  C=3.0; double dt=0.002; std::vector<double> x={1,1,1};
  auto step=[&](std::vector<double>&s){double k1[3],k2[3],k3[3],k4[3],t[3];
    ross(s.data(),C,k1,0);for(int i=0;i<3;i++)t[i]=s[i]+0.5*dt*k1[i];
    ross(t,C,k2,0);for(int i=0;i<3;i++)t[i]=s[i]+0.5*dt*k2[i];
    ross(t,C,k3,0);for(int i=0;i<3;i++)t[i]=s[i]+dt*k3[i];
    ross(t,C,k4,0);for(int i=0;i<3;i++)s[i]+=dt*(k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;};
  for(int i=0;i<200000;i++) step(x); // settle onto cycle
  // record one period via z-section (z min crossing) 
  std::vector<std::vector<double>> loop; double pz=x[2]; int cr=0; double T=0; std::vector<double> xp=x;
  for(int i=0;i<400000;i++){ double prev=x[1]; step(x); T+=dt;
    if(prev<x[1] && x[1]>=0 && prev<0){ cr++; if(cr==1){loop.clear();T=0;} else if(cr==2) break; }
    if(cr>=1) loop.push_back(x);
  }
  printf("seed: period~%.3f, loop pts=%zu\n",T,loop.size());
  if(loop.size()<10){printf("seed failed\n");return 1;}
  int mesh=100; std::vector<std::vector<double>> guess(mesh);
  for(int i=0;i<mesh;i++){double f=(double)i/mesh*(loop.size()-1);guess[i]=loop[(int)f];}
  CycleSettings cs; cs.mesh=mesh; cs.arclength=true; cs.ds=0.03; cs.max_steps=300;
  cs.compute_floquet=true; cs.adaptive_mesh=true; cs.p_min=2.5; cs.p_max=6.0;
  CycleBranch b=continue_limit_cycle(m, guess, T, 3.0, cs);
  printf("ok=%d samples=%zu\n",b.ok,b.samples.size());
  int npd=0; double pdc=-1; for(auto&s:b.samples){ if(s.is_pd){npd++; pdc=s.pd_p;} }
  printf("  PD count=%d at c=%.4f (expect 1 near 2.83)\n", npd, pdc);
  bool pass = b.ok && npd==1 && std::fabs(pdc-2.83)<0.15;
  printf("%s\n", pass?"cycle_pdns_smoke: all checks pass":"cycle_pdns_smoke: FAIL");
  return pass?0:1;
}
