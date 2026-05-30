/* Locks the effective-dimension detection used by compile_system
 * (PHASE6-UI): a map/ODE carrying a trivial top dimension (z_next=z
 * or dz=0) detects with that dimension dropped, so 'embedded in 3D'
 * 2D systems are treated as 2D. make test-dim */
/* Verify the effective-dimension detection logic: a map carrying
 * z_next = z (identity 3rd dim) detects as 2D; a genuine 3D map stays
 * 3D. Mirrors the probe loop in compile_system. */
#include <cmath>
#include <cstdio>
#include <vector>
#include <functional>
using Vec=std::vector<double>;
// map: given state, return next
static int detect_dim_map(int dim, std::function<Vec(const Vec&)> mp){
  int detected=dim;
  auto trivial=[&](int k)->bool{
    const unsigned seeds[4]={1u,7u,13u,29u};
    for(unsigned s:seeds){
      Vec probe(dim); unsigned r=s;
      for(int i=0;i<dim;i++){r=r*1664525u+1013904223u; probe[i]=(double)(r%2000u)/1000.0-1.0;}
      Vec nx=mp(probe);
      if(std::fabs(nx[k]-probe[k])>1e-9) return false;
    }
    return true;
  };
  for(int k=dim;k-->1;){ if(trivial(k)) detected=k; else break; }
  return detected;
}
int main(){
  int fails=0;
  // Henon embedded in 3D: x'=1-1.4x^2+y, y'=0.3x, z'=z
  auto henon=[](const Vec&s)->Vec{ return {1-1.4*s[0]*s[0]+s[1], 0.3*s[0], s[2]}; };
  int d1=detect_dim_map(3,henon);
  printf("Henon-in-3D detected dim = %d (expect 2) %s\n", d1, d1==2?"OK":"FAIL"); fails+=d1!=2;
  // genuine 3D map: all coords mix
  auto g3=[](const Vec&s)->Vec{ return {s[0]+0.1*s[2], s[1]-0.2*s[0], 0.5*s[2]+s[1]}; };
  int d2=detect_dim_map(3,g3);
  printf("genuine 3D detected dim = %d (expect 3) %s\n", d2, d2==3?"OK":"FAIL"); fails+=d2!=3;
  // 2D map, no dummy: stays 2
  auto m2=[](const Vec&s)->Vec{ return {1-1.4*s[0]*s[0]+s[1], 0.3*s[0]}; };
  int d3=detect_dim_map(2,m2);
  printf("2D map detected dim = %d (expect 2) %s\n", d3, d3==2?"OK":"FAIL"); fails+=d3!=2;
  // Logistic embedded: x'=r x(1-x), y'=x, z'=z -> top dim z trivial, y depends on x (not identity) -> dim 2
  auto logi=[](const Vec&s)->Vec{ return {3.9*s[0]*(1-s[0]), s[0], s[2]}; };
  int d4=detect_dim_map(3,logi);
  printf("Logistic-in-3D detected dim = %d (expect 2) %s\n", d4, d4==2?"OK":"FAIL"); fails+=d4!=2;
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
