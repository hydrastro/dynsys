/* Regression guard for pseudo-arclength limit-cycle continuation + Floquet.
 *
 * Two analytic checks:
 *  T1 (Floquet): supercritical Hopf normal form. For mu>0 there is a stable
 *      cycle of radius sqrt(mu), period ~2*pi. Its nontrivial Floquet multiplier
 *      must be inside the unit circle (stable); one multiplier is ~1 (trivial).
 *  T2 (fold of cycles): radial form  r' = r(mu + r^2 - r^4),  which has TWO
 *      positive cycles for mu in (-1/4, 0) meeting at a saddle-node of cycles at
 *      mu = -1/4. Arclength continuation seeded on the outer branch must TURN
 *      around the fold (reach p ~ -0.25 and set 'turned').
 */
#include "analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace dynsys::analysis;

static bool field_hopf(const double*x,double mu,double*f,std::string*){
  double r2=x[0]*x[0]+x[1]*x[1];
  f[0]=mu*x[0]-x[1]-x[0]*r2;
  f[1]=x[0]+mu*x[1]-x[1]*r2;
  return true;
}
static bool field_fold(const double*x,double mu,double*f,std::string*){
  double r2=x[0]*x[0]+x[1]*x[1];
  double g=mu + r2 - r2*r2;
  f[0]=x[0]*g - x[1];
  f[1]=x[1]*g + x[0];
  return true;
}

int main(){
  int fails=0;

  /* T1: Floquet on supercritical Hopf, mu=0.5 */
  {
    Model m; m.n=2; m.vector_field=field_hopf;
    double mu=0.5, R=std::sqrt(mu);
    int M=80; std::vector<std::vector<double>> g(M, std::vector<double>(2));
    for(int i=0;i<M;i++){double th=2*M_PI*i/M; g[i][0]=R*std::cos(th); g[i][1]=R*std::sin(th);}
    CycleSettings cs; cs.mesh=M; cs.arclength=true; cs.compute_floquet=true;
    cs.ds=0.05; cs.max_steps=20; cs.p_min=0.1; cs.p_max=0.9;
    CycleBranch b=continue_limit_cycle(m,g,2*M_PI,mu,cs);
    const CycleSample* s=nullptr; double bestd=1e9;
    for(auto&z:b.samples){double d=std::fabs(z.p-0.5); if(d<bestd){bestd=d;s=&z;}}
    bool ok = b.ok && s &&
              std::fabs(s->period-2*M_PI) < 0.2 &&
              std::fabs(s->amplitude-2*R) < 0.1 &&
              s->stable &&
              s->max_nontrivial_mult < 0.9 &&
              s->floquet_re.size()==2;
    printf("  T1 Floquet (Hopf cycle): ok=%d period=%.4f amp=%.4f maxmult=%.4f stable=%d %s\n",
           b.ok, s?s->period:0, s?s->amplitude:0, s?s->max_nontrivial_mult:-1, s?s->stable:0,
           ok?"PASS":"FAIL");
    if(!ok) fails++;
  }

  /* T2: fold of cycles, seed outer branch at mu=-0.1 */
  {
    Model m; m.n=2; m.vector_field=field_fold;
    double mu=-0.1;
    double r2=(1.0+std::sqrt(1.0+4.0*mu))/2.0, R=std::sqrt(r2);
    int M=80; std::vector<std::vector<double>> g(M, std::vector<double>(2));
    for(int i=0;i<M;i++){double th=2*M_PI*i/M; g[i][0]=R*std::cos(th); g[i][1]=R*std::sin(th);}
    CycleSettings cs; cs.mesh=M; cs.arclength=true; cs.compute_floquet=true;
    cs.ds=0.03; cs.max_steps=200; cs.p_min=-0.30; cs.p_max=0.30;
    CycleBranch b=continue_limit_cycle(m,g,2*M_PI,mu,cs);
    double pmin=1e9; for(auto&z:b.samples) pmin=std::min(pmin,z.p);
    /* the fold is at mu=-0.25; require we turned and reached near it */
    bool ok = b.ok && b.turned && pmin < -0.20 && pmin > -0.30;
    printf("  T2 fold-of-cycles: ok=%d turned=%d pmin=%.4f (fold@-0.25) samples=%zu %s\n",
           b.ok, b.turned, pmin, b.samples.size(), ok?"PASS":"FAIL");
    if(!ok) fails++;
  }


  /* T3: adaptive mesh keeps a known-stable cycle stable (it must not introduce
   * a spurious instability) and remeshing stays valid. Re-run the Hopf cycle
   * with adaptive_mesh on and require it is still classified stable. */
  {
    Model m; m.n=2; m.vector_field=field_hopf;
    double mu=0.5, R=std::sqrt(mu);
    int M=80; std::vector<std::vector<double>> g(M, std::vector<double>(2));
    for(int i=0;i<M;i++){double th=2*M_PI*i/M; g[i][0]=R*std::cos(th); g[i][1]=R*std::sin(th);}
    CycleSettings cs; cs.mesh=M; cs.arclength=true; cs.compute_floquet=true;
    cs.adaptive_mesh=true; cs.ds=0.05; cs.max_steps=20; cs.p_min=0.1; cs.p_max=0.9;
    CycleBranch b=continue_limit_cycle(m,g,2*M_PI,mu,cs);
    bool allstable=b.ok && !b.samples.empty();
    for(auto&z:b.samples) if(z.max_nontrivial_mult>1.05) allstable=false;
    printf("  T3 adaptive mesh (stable stays stable): ok=%d samples=%zu %s\n",
           b.ok,b.samples.size(), allstable?"PASS":"FAIL");
    if(!allstable) fails++;
  }

  if(fails){ printf("lpc_arclength_smoke: %d FAIL\n", fails); return 1; }
  printf("lpc_arclength_smoke: all checks pass\n");
  return 0;
}
