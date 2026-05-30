/* Standalone smoke test for forward-mode AD over the IR.
 *
 * Builds tiny Programs by hand (no tpcas needed) and checks that
 * run_dual's value matches the plain evaluator and its derivative
 * matches the analytic derivative and a finite-difference estimate.
 *
 *   make test-ad
 */

#include "../src/expr_ir.h"
#include "../src/expr_ir_ad.h"

#include <cmath>
#include <cstdio>
#include <vector>

using namespace dynsys::ir;

static int g_fail = 0, g_checks = 0;
static void check(bool c, const char *what) {
  ++g_checks;
  if (!c) {
    ++g_fail;
    std::printf("  FAIL: %s\n", what);
  }
}
static bool close(double a, double b, double tol = 1e-6) {
  return std::fabs(a - b) <= tol * (1.0 + std::fabs(a) + std::fabs(b));
}

static Instr I(Op op, uint16_t a = 0, uint16_t b = 0) { return Instr{op, a, b}; }

/* Build a RunContext over a 2-state, 1-param system. */
static RunContext make_ctx(const double *state, double t, const double *params,
                           const Program *defs, size_t n_defs) {
  RunContext rc;
  rc.state = state;
  rc.n_state = 2;
  rc.t = t;
  rc.params = params;
  rc.n_params = 1;
  rc.defs = defs;
  rc.n_defs = n_defs;
  return rc;
}

int main() {
  std::printf("=== dynsys AD smoke test ===\n");

  double state[2] = {1.3, -0.7};
  double params[1] = {2.0};

  /* Program 1:  f = x0 * x1   (state0 * state1)
   * df/dx0 = x1, df/dx1 = x0 */
  {
    std::printf("AD: x*y\n");
    Program p;
    p.code = {I(Op::PushState, 0), I(Op::PushState, 1), I(Op::Mul)};
    DualScratch ds;
    dual_scratch_init(&ds, 0);
    RunContext rc = make_ctx(state, 0.0, params, nullptr, 0);
    double v = 0, d = 0;
    char err[128] = {0};
    DualSeed s0{DualSeed::Kind::State, 0};
    check(run_dual(p, rc, s0, ds, &v, &d, err, sizeof err), "run ok x0");
    check(close(v, state[0] * state[1]), "value x*y");
    check(close(d, state[1]), "df/dx0 = y");
    DualSeed s1{DualSeed::Kind::State, 1};
    run_dual(p, rc, s1, ds, &v, &d, err, sizeof err);
    check(close(d, state[0]), "df/dx1 = x");
  }

  /* Program 2:  f = sin(x0) + p0 * x1^2
   * use Pow for the square.
   * df/dx0 = cos(x0); df/dx1 = 2*p0*x1; df/dp0 = x1^2 */
  {
    std::printf("AD: sin(x) + p*y^2\n");
    Program p;
    p.constants = {2.0};
    p.code = {
        I(Op::PushState, 0), I(Op::CallBuiltin, (uint16_t)Builtin::Sin, 1),
        I(Op::PushParam, 0), I(Op::PushState, 1),
        I(Op::PushConst, 0), /* 2.0 */
        I(Op::CallBuiltin, (uint16_t)Builtin::Pow, 2), I(Op::Mul), I(Op::Add)};
    DualScratch ds;
    dual_scratch_init(&ds, 0);
    RunContext rc = make_ctx(state, 0.0, params, nullptr, 0);
    char err[128] = {0};
    double v = 0, d = 0;
    const double expect_v =
        std::sin(state[0]) + params[0] * state[1] * state[1];

    DualSeed s0{DualSeed::Kind::State, 0};
    check(run_dual(p, rc, s0, ds, &v, &d, err, sizeof err), "run ok");
    check(close(v, expect_v), "value");
    check(close(d, std::cos(state[0])), "df/dx0 = cos(x)");

    DualSeed s1{DualSeed::Kind::State, 1};
    run_dual(p, rc, s1, ds, &v, &d, err, sizeof err);
    check(close(d, 2.0 * params[0] * state[1]), "df/dx1 = 2 p y");

    DualSeed sp{DualSeed::Kind::Param, 0};
    run_dual(p, rc, sp, ds, &v, &d, err, sizeof err);
    check(close(d, state[1] * state[1]), "df/dp0 = y^2");
  }

  /* Program 3: AD vs central finite difference on a transcendental
   * mix:  f = exp(x0) / (1 + x1*x1) - log(p0 + 3)  */
  {
    std::printf("AD vs finite-difference\n");
    Program p;
    p.constants = {1.0, 3.0};
    p.code = {
        I(Op::PushState, 0), I(Op::CallBuiltin, (uint16_t)Builtin::Exp, 1),
        I(Op::PushConst, 0), /* 1 */
        I(Op::PushState, 1), I(Op::PushState, 1), I(Op::Mul), I(Op::Add),
        I(Op::Div),
        I(Op::PushParam, 0), I(Op::PushConst, 1), /* 3 */
        I(Op::Add), I(Op::CallBuiltin, (uint16_t)Builtin::Log, 1), I(Op::Sub)};
    DualScratch ds;
    dual_scratch_init(&ds, 0);
    char err[128] = {0};

    auto eval_plain = [&](const double *st, double pr) {
      Scratch sc;
      scratch_init(&sc, 0);
      double prm[1] = {pr};
      RunContext rc = make_ctx(st, 0.0, prm, nullptr, 0);
      double out = 0;
      run(p, rc, sc, &out, err, sizeof err);
      return out;
    };

    for (int which = 0; which < 2; ++which) {
      double v = 0, d = 0;
      RunContext rc = make_ctx(state, 0.0, params, nullptr, 0);
      DualSeed s{DualSeed::Kind::State, (size_t)which};
      run_dual(p, rc, s, ds, &v, &d, err, sizeof err);
      double h = 1e-6;
      double sp[2] = {state[0], state[1]};
      sp[which] += h;
      double fp = eval_plain(sp, params[0]);
      sp[which] -= 2 * h;
      double fm = eval_plain(sp, params[0]);
      double fd = (fp - fm) / (2 * h);
      char msg[64];
      std::snprintf(msg, sizeof msg, "AD vs FD on state %d", which);
      check(close(d, fd, 1e-4), msg);
    }
  }

  std::printf("=== %d/%d checks passed ===\n", g_checks - g_fail, g_checks);
  return g_fail == 0 ? 0 : 1;
}
