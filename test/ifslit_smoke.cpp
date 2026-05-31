/* Locks IFS per-coefficient CONSTANT detection (this iteration): a
 * coefficient is editable iff it references no parameter. Negatives like
 * -0.04 (parsed as SUB(0,0.04)) count as constant; s*cos(theta) does not.
 * make test-ifslit */
/* Verify per-coefficient CONSTANT detection (editable iff no variable ref):
 * fern -> all constant (incl. negatives like -0.04); spiral -> e,f,p const,
 * the s*cos(theta) entries are parameter-driven. */
#include "expr_ir.h"
#include "pratt.h"
#include "ast.h"
#include "arena.h"
#include <cstdio>
#include <functional>
using namespace dynsys::ir;
static bool is_const_expr(node_t* n){
  if(!n) return true;
  if(n->kind==NODE_VAR) return false;
  if(n->kind==NODE_CONST) return true;
  if(n->kind==NODE_APP){ if(!is_const_expr(n->app.head)) return false;
    for(size_t i=0;i<n->app.argc;i++) if(!is_const_expr(n->app.args[i])) return false; return true; }
  return false;
}
static int classify(const char* e){ arena_t a; arena_init(&a,1<<16); parse_result_t pr=parse(e,&a);
  int c=pr.ast?(is_const_expr(pr.ast)?1:0):-1; arena_destroy(&a); return c; }
int main(){
  int fails=0;
  printf("=== fern (expect ALL constant, incl negatives) ===\n");
  const char* fern[]={"0","0.16","-0.04","0.85","1.6","-0.15","0.28","-0.26","0.23","0.44"};
  for(auto e: fern){ int c=classify(e); printf("  '%s' const=%d\n",e,c); if(c!=1)fails++; }
  printf("=== spiral (e,f,p const; a,b,c,d param) ===\n");
  struct{const char*e;int w;} sp[]={{"s*cos(theta)",0},{"0 - s*sin(theta)",0},{"s*sin(theta)",0},
    {"s*cos(theta)",0},{"0",1},{"0",1},{"0.5",1}};
  for(auto&s:sp){int c=classify(s.e); printf("  '%s' const=%d (want %d)\n",s.e,c,s.w); if(c!=s.w)fails++;}
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
