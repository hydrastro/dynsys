/* Locks the bifurcation period-detection logic (this iteration): on the
 * logistic map it reports period 1,2,4 in the cascade, 3 in the period-3
 * window, and 0 (chaos) at r=3.9. Pure-double reference (no IR needed).
 * make test-period */
#include <cstdio>
#include <cmath>
#include <vector>
static int detect_period(const std::vector<double>&kept){
  const double tol=1e-4; const int maxP=16; const int m=(int)kept.size();
  for(int pgap=1;pgap<=maxP&&pgap<m;++pgap){ bool match=true;
    for(int q=0;q<pgap&&(m-1-q-pgap)>=0;++q) if(std::fabs(kept[m-1-q]-kept[m-1-q-pgap])>tol){match=false;break;}
    if(match)return pgap; }
  return 0;
}
int main(){
  struct T{double r;int expect;};
  T t[]={{2.9,1},{3.2,2},{3.5,4},{3.835,3},{3.9,0}};
  int fails=0;
  for(auto&x:t){
    double v=0.5;
    for(int i=0;i<2000;i++) v=x.r*v*(1-v);
    std::vector<double> kept;
    for(int i=0;i<200;i++){ v=x.r*v*(1-v); kept.push_back(v); }
    int per=detect_period(kept);
    printf("r=%.3f: detected period %d (expect %d) %s\n", x.r, per, x.expect, per==x.expect?"":"<-- FAIL");
    if(per!=x.expect)fails++;
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
