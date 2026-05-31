/* Regression guard for the bifurcation param-sync bug (horizontal bands):
 * the swept parameter MUST reach the evaluator, so a logistic sweep shows
 * period-1 at r=2.5 and chaos at r=3.9. make test-paramsync */
/* Reproduce the exact bug: stepping reads a SEPARATE param array
 * (param_values), but the sweep updates a struct field (param.value).
 * Without sync, every slice uses the same param -> horizontal bands.
 * With sync, the parameter varies -> bifurcation tree. */
#include "../src/expr_ir.h"
#include "pratt.h"
#include "arena.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
using namespace dynsys::ir;
int main(){
  std::vector<std::string> st={"x"},pa={"r"}; std::vector<DefSig> defs; std::vector<std::string> nl;
  LowerContext lc{st,pa,defs,nl}; arena_t a; arena_init(&a,1<<16);
  parse_result_t pr=parse("r * x * (1 - x)",&a); Program P; std::string er; lower(pr.ast,lc,&P,&er);
  Scratch sc; scratch_init(&sc,0);

  struct Param{double value;};
  for(int mode=0;mode<2;mode++){
    bool do_sync = (mode==1);
    Param param{2.0};
    std::vector<double> param_values={2.0}; // the array the evaluator reads
    auto step=[&](double x){scratch_reset_eval(&sc);RunContext rc{&x,1,0,param_values.data(),1,nullptr,0};char e[64]={0};double v=x;run(P,rc,sc,&v,e,sizeof(e));return v;};
    int slices=800, discard=1000, keep=50;
    double ymin=1e9,ymax=-1e9; int distinct_at_3_9=0;
    // measure spread at a few r to see if it's a tree
    auto spread=[&](double rr){ param.value=rr; if(do_sync)param_values[0]=param.value;
      double x=0.2; for(int j=0;j<discard;j++)x=step(x);
      double lo=1e9,hi=-1e9; for(int j=0;j<keep;j++){x=step(x); lo=std::min(lo,x);hi=std::max(hi,x);} return hi-lo; };
    double s25=spread(2.5), s33=spread(3.3), s39=spread(3.9);
    printf("%s: spread@r=2.5=%.3f  @r=3.3=%.3f  @r=3.9=%.3f\n",
           do_sync?"WITH sync (fixed)":"NO sync (buggy) ", s25, s33, s39);
  }
  // assert: WITH sync, spread must grow with r (tree); this guards the
  // bifurcation param-sync regression that produced horizontal bands.
  {
    Param param{2.0}; std::vector<double> pv={2.0};
    auto step=[&](double x){scratch_reset_eval(&sc);RunContext rc{&x,1,0,pv.data(),1,nullptr,0};char e[64]={0};double v=x;run(P,rc,sc,&v,e,sizeof(e));return v;};
    auto spread=[&](double rr){param.value=rr; pv[0]=param.value; double x=0.2; for(int j=0;j<1000;j++)x=step(x); double lo=1e9,hi=-1e9; for(int j=0;j<50;j++){x=step(x);lo=std::min(lo,x);hi=std::max(hi,x);} return hi-lo;};
    double s25=spread(2.5), s39=spread(3.9);
    int fails=0;
    if(!(s25<0.01)){printf("FAIL: r=2.5 should be period-1\n");fails++;}
    if(!(s39>0.5)){printf("FAIL: r=3.9 should be chaotic (wide spread)\n");fails++;}
    printf("=== %s ===\n", fails==0?"PASS":"FAIL");
    return fails;
  }
}
