#include "expr_ir.h"
extern "C" {
#include "ast.h"
#include "arena.h"
}
#include <cstdio>
#include <cmath>
#include <vector>
#include <string>
using namespace dynsys::ir;
static int fails=0; static void chk(const char*l,bool c){printf("  %s : %s\n",l,c?"ok":"FAIL");if(!c)fails++;}

static bool eval_e(const std::vector<std::string>& params, const double* pvals, double* out){
  arena_t a; arena_init(&a, 1<<16);
  span_t sp = {0,0};
  node_t *e_node = ast_var(&a, "e", sp);
  std::vector<std::string> states;
  std::vector<DefSig> defs;
  std::vector<std::string> locals;
  LowerContext ctx{states, params, defs, locals};
  Program prog; std::string err;
  if(!lower(e_node, ctx, &prog, &err)){ printf("   lower failed: %s\n", err.c_str()); arena_destroy(&a); return false; }
  Scratch s; scratch_init(&s, 0); scratch_reset_eval(&s);
  RunContext rc{ nullptr, 0, 0.0, pvals, params.size(), nullptr, 0 };
  char eb[256]={0};
  bool ok = run(prog, rc, s, out, eb, sizeof(eb));
  if(!ok) printf("   run failed: %s\n", eb);
  arena_destroy(&a);
  return ok;
}

int main(){
  // Case 1: a param named "e" with value 0.25 -> e resolves to 0.25 (shadow)
  {
    std::vector<std::string> params = {"e"};
    double pv[1] = {0.25};
    double out = -1;
    bool ok = eval_e(params, pv, &out);
    printf("param e=0.25 -> e evaluates to %.6f\n", out);
    chk("param e shadows the constant (==0.25)", ok && std::fabs(out-0.25)<1e-12);
  }
  // Case 2: no params -> e is Euler's number 2.71828...
  {
    std::vector<std::string> params;
    double out = -1;
    bool ok = eval_e(params, nullptr, &out);
    printf("no params -> e evaluates to %.6f\n", out);
    chk("e falls back to Euler when not shadowed", ok && std::fabs(out-2.718281828)<1e-6);
  }
  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}
