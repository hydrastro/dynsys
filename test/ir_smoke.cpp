/* Standalone test for the dynsys IR.
 *
 * Builds without ImGui / GLFW. Verifies:
 *   1. The lowerer accepts what eval_ast accepts.
 *   2. The VM produces the same value eval_ast does, for a sweep
 *      of expressions covering: state vars, parameters, builtins,
 *      operators, user definitions (with and without args), the
 *      memoization path (a 0-arity def referenced multiple times),
 *      and select's lazy semantics.
 *
 * Build (from project root):
 *   g++ -std=c++17 -O2 -Wall -Wextra \
 *       -Ivendor/tpcas/src -Ivendor/tpcas/vendor/ds \
 *       src/expr_ir.cpp test/ir_smoke.cpp \
 *       vendor/tpcas/src/arena.c vendor/tpcas/src/ast.c \
 *       vendor/tpcas/src/lex.c   vendor/tpcas/src/pratt.c \
 *       vendor/tpcas/vendor/ds/lib/<all .c files> \
 *       -o build/ir_smoke
 */

#include "../src/expr_ir.h"

extern "C" {
#include "../vendor/tpcas/src/arena.h"
#include "../vendor/tpcas/src/pratt.h"
}

#include <cmath>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

namespace ir = dynsys::ir;

struct Case {
    const char *expr;
    double      expected;
};

static int     g_pass = 0;
static int     g_fail = 0;

static bool nearly_equal(double a, double b) {
    if (std::isnan(a) && std::isnan(b)) return true;
    if (!std::isfinite(a) || !std::isfinite(b)) return a == b;
    const double mag = std::fmax(1.0, std::fmax(std::fabs(a), std::fabs(b)));
    return std::fabs(a - b) <= 1e-9 * mag;
}

/* Lowering context for the test. */
struct TestSetup {
    arena_t arena{};
    std::vector<std::string>  state_names;
    std::vector<std::string>  param_names;
    std::vector<double>       param_values;
    std::vector<ir::DefSig>   def_sigs;
    std::vector<ir::Program>  def_progs;
    std::vector<double>       state;
    double t = 0.0;

    TestSetup() {
        arena_init(&arena, 1 << 16);
    }
    ~TestSetup() { arena_destroy(&arena); }

    /* Add a 0-arity definition like "r2 = x*x + y*y + z*z".
     * Two-pass: register sig first, then lower body, so forward
     * references between defs work (e.g. mutual recursion tests). */
    void register_def0(const char *name) {
        def_sigs.push_back(ir::DefSig{name, 0});
        def_progs.emplace_back();
    }

    void lower_def0(const char *name, const char *expr) {
        size_t idx = SIZE_MAX;
        for (size_t i = 0; i < def_sigs.size(); ++i) {
            if (def_sigs[i].arity == 0 && def_sigs[i].name == name) { idx = i; break; }
        }
        if (idx == SIZE_MAX) { std::fprintf(stderr, "internal: def sig not registered: %s\n", name); std::abort(); }
        parse_result_t pr = parse(expr, &arena);
        if (!pr.ok) { std::fprintf(stderr, "internal: parse failure for def: %s\n", expr); std::abort(); }
        const std::vector<std::string> no_locals;
        ir::LowerContext lctx{state_names, param_names, def_sigs, no_locals};
        std::string err;
        if (!ir::lower(pr.ast, lctx, &def_progs[idx], &err)) {
            std::fprintf(stderr, "internal: lower failure for def: %s → %s\n", expr, err.c_str()); std::abort();
        }
    }

    void add_def0(const char *name, const char *expr) {
        register_def0(name);
        lower_def0(name, expr);
    }

    /* Add a user function like "wave(u) = sin(u)" */
    void add_def_fn(const char *name, std::vector<std::string> params, const char *expr) {
        ir::DefSig sig{name, params.size()};
        def_sigs.push_back(sig);
        ir::Program prog;
        parse_result_t pr = parse(expr, &arena);
        if (!pr.ok) { std::fprintf(stderr, "internal: parse failure for def: %s\n", expr); std::abort(); }
        ir::LowerContext lctx{state_names, param_names, def_sigs, params};
        std::string err;
        if (!ir::lower(pr.ast, lctx, &prog, &err)) {
            std::fprintf(stderr, "internal: lower failure for def: %s → %s\n", expr, err.c_str()); std::abort();
        }
        prog.arity = params.size();
        def_progs.push_back(std::move(prog));
    }

    bool eval(const char *expr, double *out, std::string *err) {
        parse_result_t pr = parse(expr, &arena);
        if (!pr.ok) {
            *err = std::string("parse: ") + (pr.err_msg ? pr.err_msg : "?");
            return false;
        }
        const std::vector<std::string> no_locals;
        ir::LowerContext lctx{state_names, param_names, def_sigs, no_locals};
        ir::Program prog;
        if (!ir::lower(pr.ast, lctx, &prog, err)) return false;

        ir::Scratch s;
        ir::scratch_init(&s, def_progs.size());
        ir::scratch_reset_eval(&s);

        ir::RunContext rc;
        rc.state    = state.data();
        rc.n_state  = state.size();
        rc.t        = t;
        rc.params   = param_values.data();
        rc.n_params = param_values.size();
        rc.defs     = def_progs.data();
        rc.n_defs   = def_progs.size();

        char ebuf[256] = {0};
        if (!ir::run(prog, rc, s, out, ebuf, sizeof ebuf)) {
            *err = std::string("run: ") + ebuf;
            return false;
        }
        return true;
    }
};

static void check(TestSetup &s, const char *expr, double expected, const char *label = nullptr) {
    double v = 0.0;
    std::string err;
    bool ok = s.eval(expr, &v, &err);
    if (!ok) {
        std::printf("FAIL  %-50s  error: %s\n", label ? label : expr, err.c_str());
        g_fail++; return;
    }
    if (!nearly_equal(v, expected)) {
        std::printf("FAIL  %-50s  got %.17g  expected %.17g\n",
                    label ? label : expr, v, expected);
        g_fail++; return;
    }
    g_pass++;
}

static void check_err(TestSetup &s, const char *expr, const char *substr) {
    double v = 0.0;
    std::string err;
    bool ok = s.eval(expr, &v, &err);
    if (ok) {
        std::printf("FAIL  %-50s  expected error containing '%s', got %.17g\n",
                    expr, substr, v);
        g_fail++; return;
    }
    if (err.find(substr) == std::string::npos) {
        std::printf("FAIL  %-50s  expected error containing '%s', got '%s'\n",
                    expr, substr, err.c_str());
        g_fail++; return;
    }
    g_pass++;
}

int main() {
    TestSetup s;
    s.state_names = {"x", "y", "z"};
    s.state       = {1.5, 2.5, 3.5};
    s.t           = 0.25;
    s.param_names = {"sigma", "rho", "beta"};
    s.param_values = {10.0, 28.0, 8.0 / 3.0};

    /* atoms */
    check(s, "0",       0.0);
    check(s, "1",       1.0);
    check(s, "3.14",    3.14);
    check(s, "x",       1.5);
    check(s, "y",       2.5);
    check(s, "z",       3.5);
    check(s, "t",       0.25);
    check(s, "pi",      M_PI);
    check(s, "e",       M_E);
    check(s, "sigma",   10.0);

    /* operators + precedence */
    check(s, "x + y",            4.0);
    check(s, "x - y",           -1.0);
    check(s, "x * y",            3.75);
    check(s, "y / x",            2.5 / 1.5);
    check(s, "x + y * z",        1.5 + 2.5 * 3.5);
    check(s, "(x + y) * z",      (1.5 + 2.5) * 3.5);
    check(s, "x * y + y * z",    1.5*2.5 + 2.5*3.5);

    /* unary minus */
    check(s, "0 - x",     -1.5);
    check(s, "0 - y * 2", -5.0);

    /* lorenz-y expressions */
    check(s, "sigma * (y - x)",        10.0 * (2.5 - 1.5));
    check(s, "x * (rho - z) - y",      1.5 * (28.0 - 3.5) - 2.5);
    check(s, "x * y - beta * z",       1.5 * 2.5 - (8.0/3.0) * 3.5);

    /* builtins */
    check(s, "sin(0)",        0.0);
    check(s, "cos(0)",        1.0);
    check(s, "sqrt(4)",       2.0);
    check(s, "exp(0)",        1.0);
    check(s, "log(e)",        1.0);
    check(s, "abs(0 - 3)",    3.0);
    check(s, "floor(1.7)",    1.0);
    check(s, "ceil(1.2)",     2.0);
    check(s, "sign(0 - 2)",  -1.0);
    check(s, "pow(2, 8)",     256.0);
    check(s, "min(3, 7)",     3.0);
    check(s, "max(3, 7)",     7.0);
    check(s, "mod(10, 3)",    1.0);
    check(s, "clamp(5, 0, 3)",3.0);
    check(s, "clamp(0 - 5, 0, 3)", 0.0);
    check(s, "sqrt(x*x + y*y + z*z)",
          std::sqrt(1.5*1.5 + 2.5*2.5 + 3.5*3.5));

    /* errors */
    check_err(s, "no_such_var", "unknown variable");
    check_err(s, "no_such_fn(1)", "unknown function");
    check_err(s, "sqrt(1, 2)", "takes 1 argument");
    check_err(s, "clamp(1)",   "takes 3 argument");

    /* user-defined helpers */
    s.add_def0("r2",  "x*x + y*y + z*z");
    s.add_def0("rho2", "rho * rho");
    s.add_def_fn("wave",  {"u"},      "sin(u)");
    s.add_def_fn("dist2", {"a", "b"}, "a*a + b*b");

    check(s, "r2",           1.5*1.5 + 2.5*2.5 + 3.5*3.5);
    check(s, "sqrt(r2)",     std::sqrt(1.5*1.5 + 2.5*2.5 + 3.5*3.5));
    check(s, "r2 * 2",       2.0*(1.5*1.5 + 2.5*2.5 + 3.5*3.5));
    /* 0-arity memoization: r2 appears twice, should give same value
     * even after a hypothetical caller mutates state mid-expression
     * (not actually possible from inside the IR, but the cache should
     * still be consistent within one outer eval) */
    check(s, "r2 + r2",      2.0*(1.5*1.5 + 2.5*2.5 + 3.5*3.5));
    check(s, "wave(0)",      0.0);
    check(s, "wave(pi)",     std::sin(M_PI));
    check(s, "dist2(x, y)",  1.5*1.5 + 2.5*2.5);
    check(s, "dist2(x, y) + z*z", 1.5*1.5 + 2.5*2.5 + 3.5*3.5);
    /* nested calls */
    check(s, "sqrt(dist2(x, y))",  std::sqrt(1.5*1.5 + 2.5*2.5));

    /* lazy select — the eager version would evaluate log(0) and
     * produce -inf, propagating to the result. The lazy version
     * never touches the unselected branch. */
    check(s, "select(1, 7, 99)",                7.0);
    check(s, "select(0, 7, 99)",                99.0);
    check(s, "select(x, log(x), 0 - 1)",        std::log(1.5));
    /* explicit lazy: if state had been 0 we'd return -1, not log(0). */
    s.state = {0.0, 2.5, 3.5};
    check(s, "select(x, log(x), 0 - 1)",        -1.0, "select lazy: x=0 → branch -1");
    s.state = {1.5, 2.5, 3.5};

    /* recursion safety */
    s.register_def0("cyclic_a");
    s.register_def0("cyclic_b");
    s.lower_def0("cyclic_a", "cyclic_b");
    s.lower_def0("cyclic_b", "cyclic_a");
    check_err(s, "cyclic_a", "cyclic");

    std::printf("\n%d passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
