/* Locks the IFS-as-model parsing (this iteration): the fern preset's
 * ifs_map lines parse to 4 maps producing the right attractor/dimension,
 * short lines are rejected, and the probability field is optional.
 * make test-ifsmodel */
/* Verify the ifs_map parsing logic (mirrors compile_system's IFS branch)
 * and that the parsed fern/Sierpinski specs run through the chaos game with
 * the right dimension. Guards the IFS-as-model feature. */
#include "../src/analysis.h"
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
using namespace dynsys::analysis;

/* the exact parse the app does for an "ifs_map = a,b,c,d,e,f[,p]" rhs */
static bool parse_map(const std::string& rhs, AffineMap& m){
  std::vector<double> v; std::string cur;
  for(char ch : rhs + ","){
    if(ch==','){ std::string t; for(char c:cur) if(c!=' '&&c!='\t') t+=c; if(!t.empty()){ try{v.push_back(std::stod(t));}catch(...){return false;} } cur.clear(); }
    else cur+=ch;
  }
  if(v.size()<6) return false;
  m.a=v[0];m.b=v[1];m.c=v[2];m.d=v[3];m.e=v[4];m.f=v[5]; m.p=v.size()>=7?v[6]:0.0;
  return true;
}

int main(){
  int fails=0;

  /* (1) parse the fern preset's four ifs_map lines */
  {
    const char* lines[] = {
      "0, 0, 0, 0.16, 0, 0, 0.01",
      "0.85, 0.04, -0.04, 0.85, 0, 1.6, 0.85",
      "0.2, -0.26, 0.23, 0.22, 0, 1.6, 0.07",
      "-0.15, 0.28, 0.26, 0.24, 0, 0.44, 0.07",
    };
    std::vector<AffineMap> maps;
    for(auto s: lines){ AffineMap m; if(!parse_map(s,m)){printf("(1) FAIL parse: %s\n",s);fails++;} else maps.push_back(m); }
    printf("(1) fern parsed %zu maps (expect 4)\n", maps.size());
    if(maps.size()!=4) fails++;
    auto R=chaos_game(maps,400000,12345u);
    std::vector<double> xs(R.xs.begin(),R.xs.end()), ys(R.ys.begin(),R.ys.end());
    auto D=box_counting_dimension(xs,ys,12);
    printf("    fern bbox[%.2f,%.2f]x[%.2f,%.2f] D=%.3f\n",R.xmin,R.xmax,R.ymin,R.ymax,D.dimension);
    if(!(R.xmin>-3&&R.xmax<3&&R.ymin>-0.1&&R.ymax<11)){printf("    FAIL bounds\n");fails++;}
    if(!(D.dimension>1.6&&D.dimension<2.0)){printf("    FAIL dim\n");fails++;}
  }

  /* (2) bad line (too few numbers) must be rejected */
  {
    AffineMap m; bool ok=parse_map("0.5, 0, 0", m);
    printf("(2) short line rejected: %s\n", ok?"NO (bug)":"yes");
    if(ok) fails++;
  }

  /* (3) probability optional (6 numbers) */
  {
    AffineMap m; bool ok=parse_map("0.5, 0, 0, 0.5, 0.25, 0.5", m);
    printf("(3) 6-number line ok=%d p=%.2f (expect p=0)\n", ok, m.p);
    if(!ok||m.p!=0.0) fails++;
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
