#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
using namespace dynsys::analysis;
static bool sys(const double*X,double p,double q,double*o,std::string*){(void)q;double x=X[0],y=X[1];o[0]=y;o[1]=x-x*x+p*y;return true;}
int main(){
  Model2 m2; m2.n=2; m2.vector_field=sys;
  std::vector<double> sad={0,0};
  HomoclinicSettings st; st.mesh=100; st.newton_iters=60; st.free_T=true; st.T=0;
  double p_out=-99;
  HomoclinicResult R=find_homoclinic(m2, sad, 0.0, -0.10, 0.10, st, &p_out);
  printf("find_homoclinic: ok=%d located p=%.5f (expect ~0.0) residual=%.3e T=%.3f amp=%.3f\n",R.ok,p_out,R.newton_residual,R.T,R.amplitude);
  printf("msg: %s\n",R.message.c_str());
  bool pass = std::fabs(p_out-0.0)<0.012;  // located near homoclinic at p=0
  printf("%s\n", pass?"find_homoclinic_smoke: all checks pass":"find_homoclinic_smoke: FAIL");
  return pass?0:1;
}
