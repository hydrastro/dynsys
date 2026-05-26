#pragma once

/* ============================================================
 * dynsys expression IR.
 *
 * compile_system parses each .dyn expression with tpcas into a
 * node_t AST. The hot integration loop used to walk that AST and
 * resolve every variable / function name by strcmp on every
 * evaluation. That ran ~tens of thousands of strcmps per frame for
 * a 3D system at 60fps × 96 substeps × 4 RK stages.
 *
 * This module lowers each AST to a flat opcode array once at
 * compile time. The hot eval becomes a tight switch with no string
 * comparisons, no name lookups, and no per-call vector::assign.
 *
 * The module is intentionally free of AppState/imgui dependencies
 * so it can be unit-tested on its own.
 * ============================================================ */

#include <cstdint>
#include <cstddef>
#include <string>
#include <vector>

extern "C" {
#include "../vendor/tpcas/src/ast.h"
}

namespace dynsys::ir {

enum class Op : uint8_t {
    PushConst,     /* a = constant index                       */
    PushState,     /* a = state-variable index                 */
    PushParam,     /* a = parameter index                      */
    PushLocal,     /* a = local (user-def argument) slot       */
    PushT, PushPi, PushE,
    Neg,
    Add, Sub, Mul, Div,
    CallBuiltin,   /* a = Builtin id, b = arity                */
    CallDef,       /* a = def index, b = arity                 */
    /* lazy conditional: pop top; if !=0 fall through to TRUE
     * branch; if ==0 jump forward by `a` bytes-of-Instr. The
     * TRUE branch ends with Jump(b) skipping past the FALSE
     * branch. Used by select(). */
    BrIfZero,      /* a = relative offset (skip distance)      */
    Jump,          /* a = relative offset                      */
};

enum class Builtin : uint8_t {
    Sin, Cos, Tan, Asin, Acos, Atan,
    Exp, Log, Log10, Sqrt, Abs, Floor, Ceil, Sign,
    Pow, Min, Max, Mod,
    Clamp,
    /* Select is lowered to BrIfZero/Jump, not CallBuiltin. */
    Unknown,
};

struct Instr {
    Op       op;
    uint16_t a;
    uint16_t b;
};

struct Program {
    std::vector<Instr>  code;
    std::vector<double> constants;
    /* For user-def programs: arity of the def. 0 for top-level
     * expressions (equations, observables, plot3d axes, ...). */
    size_t              arity = 0;
};

struct DefSig {
    std::string name;
    size_t      arity;
};

/* Resolution context for lowering. The lowerer doesn't know
 * AppState — it only needs the names that resolve, in order. */
struct LowerContext {
    const std::vector<std::string>            &state_names;
    const std::vector<std::string>            &param_names;
    const std::vector<DefSig>                 &defs;
    /* Non-empty only when lowering a user-def body; gives the
     * argument names for that def in declaration order. */
    const std::vector<std::string>            &locals;
};

/* Returns true on success. On failure `err` is populated and
 * `out` may be partially built; callers should discard it. */
bool lower(const node_t *ast,
           const LowerContext &ctx,
           Program *out,
           std::string *err);

/* Per-thread / per-eval reusable working memory. Owned by the
 * caller so allocations don't churn across the hot path. */
struct Scratch {
    std::vector<double>  stack;     /* value stack                       */
    std::vector<double>  locals;    /* nested call-frame args            */
    std::vector<uint8_t> active_def;/* cycle guard, size = n_defs        */
    std::vector<uint8_t> cached_def;/* memo flag for 0-arity defs        */
    std::vector<double>  cache_def; /* memoized 0-arity def results      */
    int                  depth = 0;
};

/* Configure scratch storage for a given def count. Must be called
 * (or re-called) whenever the system is re-compiled. */
void scratch_init(Scratch *s, size_t n_defs);

/* Clears the per-eval memo so a fresh "outer" run starts with no
 * cached 0-arity-def values. The caller invokes this before each
 * top-level run() so semantics match the old eval_expr_at, where
 * the cache was reset on every call. */
void scratch_reset_eval(Scratch *s);

struct RunContext {
    const double  *state;     /* state vector,    size = n_state         */
    size_t         n_state;
    double         t;
    const double  *params;    /* parameter values, size = n_params       */
    size_t         n_params;
    const Program *defs;      /* user-def programs, indexed by def index */
    size_t         n_defs;
};

bool run(const Program  &program,
         const RunContext &ctx,
         Scratch        &scratch,
         double         *out,
         char           *err_buf,
         size_t          err_cap);

/* Helpers, exposed for the test driver. */
const char *builtin_name(Builtin b);
Builtin     builtin_from_name(const char *name);  /* returns Unknown if absent */

}  /* namespace dynsys::ir */
