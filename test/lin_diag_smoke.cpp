#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static bool bt(const double*X,double p,double q,double*o,std::string*){double x=X[0],y=X[1];o[0]=y;o[1]=p+q*y+x*x+x*y;return true;}
int main(){
  Model2 m2; m2.n=2; m2.vector_field=bt;
  // q=-1, p around -0.15: W^u goes to the OTHER equilibrium (heteroclinic)
  std::vector<double> sad={std::sqrt(0.15),0};
  LinSettings st; st.eps=1e-5; st.dt=0.002; st.max_time=400; st.max_iter=30; st.tol=1e-8; st.p_lo=-0.5; st.p_hi=-0.02;
  LinResult R=lin_homoclinic(m2, sad, -0.15, -1.0, st);
  printf("ok=%d\nmessage: %s\n", R.ok, R.message.c_str());
  // Should now correctly say "heteroclinic, use solve_heteroclinic" or "did not return"
  bool informative = (R.message.find("heteroclinic")!=std::string::npos) || (R.message.find("did not return")!=std::string::npos);
  printf("%s\n", informative?"lin_diag_smoke: all checks pass":"lin_diag_smoke: FAIL");
  return informative?0:1;
}
