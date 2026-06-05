/* Phase E: exercises the Sangaku CAS bridge (exact eigenvalues + Routh-Hurwitz
 * stability). GRACEFULLY SKIPS (passes) when LIZARD/SANGAKU_ROOT are unset or
 * the interpreter is absent, so CI without the CAS still goes green. When the
 * CAS is present it verifies saddle/focus/node/center classification exactly.
 * make test-cas */
/* Test the CAS bridge against the real Lizard+Sangaku build. */
#include "../src/cas_bridge.h"
#include <cstdio>
using namespace dynsys::cas;
static int fails=0;
static void check(const char* label, bool cond){ printf("  %s : %s\n", label, cond?"ok":"FAIL"); if(!cond)fails++; }

int main(){
  printf("CAS available: %s\n", is_available()?"yes":"no");
  if(!is_available()){ printf("(set LIZARD + SANGAKU_ROOT)\n"); return 0; }

  // saddle [[1,0],[0,-1]] -> eigenvalues 1, -1 -> 1 RHP, 1 LHP, unstable
  { auto R = eigen_report({{Rational(1),Rational(0)},{Rational(0),Rational(-1)}});
    printf("saddle: ok=%d charpoly_n=%zu eig=%zu rhp=%d lhp=%d stable=%d\n",
           R.ok, R.charpoly.size(), R.eigen_strings.size(), R.n_rhp, R.n_lhp, R.stable);
    check("saddle solved", R.ok);
    check("saddle has 1 RHP root", R.n_rhp==1);
    check("saddle unstable", !R.stable);
  }
  // stable focus [[-1,-1],[1,-1]] -> -1+/-i -> 0 RHP, stable
  { auto R = eigen_report({{Rational(-1),Rational(-1)},{Rational(1),Rational(-1)}});
    printf("stable focus: ok=%d rhp=%d lhp=%d stable=%d eig0=%s\n",
           R.ok, R.n_rhp, R.n_lhp, R.stable, R.eigen_strings.empty()?"":R.eigen_strings[0].c_str());
    check("focus solved", R.ok);
    check("focus stable (0 RHP)", R.stable && R.n_rhp==0);
  }
  // stable node 3x3 diag(-1,-2,-3) -> all LHP, stable
  { auto R = eigen_report({{Rational(-1),Rational(0),Rational(0)},
                           {Rational(0),Rational(-2),Rational(0)},
                           {Rational(0),Rational(0),Rational(-3)}});
    printf("stable node3: ok=%d rhp=%d stable=%d\n", R.ok, R.n_rhp, R.stable);
    check("node3 stable", R.ok && R.stable);
  }
  // unstable node diag(1,2) -> 2 RHP
  { auto R = eigen_report({{Rational(1),Rational(0)},{Rational(0),Rational(2)}});
    printf("unstable node: rhp=%d stable=%d\n", R.n_rhp, R.stable);
    check("unstable node has 2 RHP", R.n_rhp==2 && !R.stable);
  }
  // rational fraction entries: [[-1/2, 1],[-1, -1/2]] -> -1/2 +/- i, stable
  { auto R = eigen_report({{Rational(-1,2),Rational(1)},{Rational(-1),Rational(-1,2)}});
    printf("fractional: ok=%d charpoly: ", R.ok);
    for(auto&c:R.charpoly) printf("%s ", c.to_lisp().c_str());
    printf(" stable=%d\n", R.stable);
    check("fractional solved + stable", R.ok && R.stable);
  }

  printf("=== %s ===\n", fails==0?"PASS":"FAIL");
  return fails;
}

