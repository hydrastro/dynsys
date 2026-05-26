#include "expr_ir.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstring>

namespace dynsys::ir {

namespace {

void set_err(std::string *err, const char *fmt, ...) {
    if (!err) return;
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    *err = buf;
}

void set_err_buf(char *err_buf, size_t err_cap, const char *fmt, ...) {
    if (!err_buf || err_cap == 0) return;
    va_list ap;
    va_start(ap, fmt);
    std::vsnprintf(err_buf, err_cap, fmt, ap);
    va_end(ap);
}

struct BuiltinSpec {
    const char *name;
    Builtin     id;
    int         arity;   /* -1 = variable, used by clamp(3) */
};

const BuiltinSpec kBuiltins[] = {
    {"sin",   Builtin::Sin,   1},
    {"cos",   Builtin::Cos,   1},
    {"tan",   Builtin::Tan,   1},
    {"asin",  Builtin::Asin,  1},
    {"acos",  Builtin::Acos,  1},
    {"atan",  Builtin::Atan,  1},
    {"exp",   Builtin::Exp,   1},
    {"log",   Builtin::Log,   1},
    {"log10", Builtin::Log10, 1},
    {"sqrt",  Builtin::Sqrt,  1},
    {"abs",   Builtin::Abs,   1},
    {"floor", Builtin::Floor, 1},
    {"ceil",  Builtin::Ceil,  1},
    {"sign",  Builtin::Sign,  1},
    {"pow",   Builtin::Pow,   2},
    {"min",   Builtin::Min,   2},
    {"max",   Builtin::Max,   2},
    {"mod",   Builtin::Mod,   2},
    {"clamp", Builtin::Clamp, 3},
};

const BuiltinSpec *find_builtin(const char *name) {
    for (const auto &b : kBuiltins) {
        if (std::strcmp(b.name, name) == 0) return &b;
    }
    return nullptr;
}

bool parse_double(const char *s, double *out) {
    if (!s) return false;
    char *end = nullptr;
    const double v = std::strtod(s, &end);
    if (end == s || (end && *end != '\0')) return false;
    *out = v;
    return true;
}

void emit(Program *p, Op op, uint16_t a = 0, uint16_t b = 0) {
    p->code.push_back({op, a, b});
}

uint16_t intern_const(Program *p, double v) {
    for (size_t i = 0; i < p->constants.size(); ++i) {
        if (p->constants[i] == v) return static_cast<uint16_t>(i);
    }
    p->constants.push_back(v);
    return static_cast<uint16_t>(p->constants.size() - 1);
}

int local_index(const LowerContext &ctx, const char *name) {
    for (size_t i = 0; i < ctx.locals.size(); ++i) {
        if (ctx.locals[i] == name) return static_cast<int>(i);
    }
    return -1;
}

int state_index(const LowerContext &ctx, const char *name) {
    for (size_t i = 0; i < ctx.state_names.size(); ++i) {
        if (ctx.state_names[i] == name) return static_cast<int>(i);
    }
    return -1;
}

int param_index(const LowerContext &ctx, const char *name) {
    for (size_t i = 0; i < ctx.param_names.size(); ++i) {
        if (ctx.param_names[i] == name) return static_cast<int>(i);
    }
    return -1;
}

int def_index_for(const LowerContext &ctx, const char *name, size_t arity) {
    for (size_t i = 0; i < ctx.defs.size(); ++i) {
        if (ctx.defs[i].arity == arity && ctx.defs[i].name == name) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

bool lower_node(const node_t *ast, const LowerContext &ctx,
                Program *out, std::string *err);

/* select(c, a, b) — lazy. Emit:
 *   [lower c]
 *   BrIfZero  off_to_FALSE
 *   [lower a]
 *   Jump      off_to_END
 *   [FALSE]: [lower b]
 *   [END]
 *
 * Offsets are filled in after the targets are known. */
bool lower_select(const node_t *ast, const LowerContext &ctx,
                  Program *out, std::string *err) {
    if (ast->app.argc != 3) {
        set_err(err, "select takes exactly 3 arguments (cond, a, b)");
        return false;
    }
    if (!lower_node(ast->app.args[0], ctx, out, err)) return false;

    const size_t br_idx = out->code.size();
    emit(out, Op::BrIfZero, 0);

    if (!lower_node(ast->app.args[1], ctx, out, err)) return false;

    const size_t jmp_idx = out->code.size();
    emit(out, Op::Jump, 0);

    const size_t false_target = out->code.size();
    if (!lower_node(ast->app.args[2], ctx, out, err)) return false;
    const size_t end_target = out->code.size();

    /* Branch offsets are in units of "instructions from end of branch
     * instruction" — i.e. number of instrs to skip ahead. */
    out->code[br_idx].a  = static_cast<uint16_t>(false_target - (br_idx + 1));
    out->code[jmp_idx].a = static_cast<uint16_t>(end_target   - (jmp_idx + 1));
    return true;
}

bool lower_application(const node_t *ast, const LowerContext &ctx,
                       Program *out, std::string *err) {
    /* Operator application: head is a NODE_CONST that carries an op_info_t. */
    if (ast->app.head->kind == NODE_CONST && ast->app.head->cnst.op != nullptr) {
        const op_info_t *op = ast->app.head->cnst.op;

        if (op == &OP_SUB && ast->app.argc == 1) {
            if (!lower_node(ast->app.args[0], ctx, out, err)) return false;
            emit(out, Op::Neg);
            return true;
        }
        if (ast->app.argc != 2) {
            set_err(err, "operator %s expects 2 arguments (got %zu)",
                    op->name, ast->app.argc);
            return false;
        }
        if (!lower_node(ast->app.args[0], ctx, out, err)) return false;
        if (!lower_node(ast->app.args[1], ctx, out, err)) return false;
        if (op == &OP_ADD) { emit(out, Op::Add); return true; }
        if (op == &OP_SUB) { emit(out, Op::Sub); return true; }
        if (op == &OP_MUL) { emit(out, Op::Mul); return true; }
        if (op == &OP_DIV) { emit(out, Op::Div); return true; }
        set_err(err, "operator %s is not allowed inside numeric expressions", op->name);
        return false;
    }

    /* Function call: head is a NODE_VAR with a name. */
    if (ast->app.head->kind == NODE_VAR) {
        const char *name = ast->app.head->var.name;
        const size_t argc = ast->app.argc;

        if (std::strcmp(name, "select") == 0) {
            return lower_select(ast, ctx, out, err);
        }

        const BuiltinSpec *bi = find_builtin(name);
        if (bi != nullptr) {
            if (bi->arity != static_cast<int>(argc)) {
                set_err(err, "builtin '%s' takes %d argument(s), got %zu",
                        bi->name, bi->arity, argc);
                return false;
            }
            for (size_t i = 0; i < argc; ++i) {
                if (!lower_node(ast->app.args[i], ctx, out, err)) return false;
            }
            emit(out, Op::CallBuiltin,
                 static_cast<uint16_t>(bi->id),
                 static_cast<uint16_t>(argc));
            return true;
        }

        const int di = def_index_for(ctx, name, argc);
        if (di < 0) {
            set_err(err, "unknown function: %s/%zu", name, argc);
            return false;
        }
        for (size_t i = 0; i < argc; ++i) {
            if (!lower_node(ast->app.args[i], ctx, out, err)) return false;
        }
        emit(out, Op::CallDef,
             static_cast<uint16_t>(di),
             static_cast<uint16_t>(argc));
        return true;
    }

    set_err(err, "unsupported application form");
    return false;
}

bool lower_node(const node_t *ast, const LowerContext &ctx,
                Program *out, std::string *err) {
    if (ast == nullptr) {
        set_err(err, "null AST node");
        return false;
    }

    switch (ast->kind) {
    case NODE_CONST: {
        double v = 0.0;
        if (!parse_double(ast->cnst.name, &v)) {
            set_err(err, "non-numeric constant: %s", ast->cnst.name);
            return false;
        }
        emit(out, Op::PushConst, intern_const(out, v));
        return true;
    }

    case NODE_VAR: {
        const char *name = ast->var.name;

        if (int i = local_index(ctx, name); i >= 0) {
            emit(out, Op::PushLocal, static_cast<uint16_t>(i));
            return true;
        }
        if (int i = state_index(ctx, name); i >= 0) {
            emit(out, Op::PushState, static_cast<uint16_t>(i));
            return true;
        }
        if (std::strcmp(name, "t")  == 0) { emit(out, Op::PushT);  return true; }
        if (std::strcmp(name, "pi") == 0) { emit(out, Op::PushPi); return true; }
        if (std::strcmp(name, "e")  == 0) { emit(out, Op::PushE);  return true; }
        if (int i = param_index(ctx, name); i >= 0) {
            emit(out, Op::PushParam, static_cast<uint16_t>(i));
            return true;
        }
        if (int i = def_index_for(ctx, name, 0); i >= 0) {
            emit(out, Op::CallDef, static_cast<uint16_t>(i), 0);
            return true;
        }
        set_err(err, "unknown variable: %s", name);
        return false;
    }

    case NODE_APP:
        return lower_application(ast, ctx, out, err);

    case NODE_LAMBDA:
    case NODE_FORALL:
    case NODE_EXISTS:
        set_err(err, "binders are not allowed inside numeric expressions");
        return false;
    }

    set_err(err, "unsupported AST node kind");
    return false;
}

}  /* namespace */

bool lower(const node_t *ast, const LowerContext &ctx,
           Program *out, std::string *err) {
    out->code.clear();
    out->constants.clear();
    return lower_node(ast, ctx, out, err);
}

void scratch_init(Scratch *s, size_t n_defs) {
    s->stack.clear();
    s->stack.reserve(64);
    s->locals.clear();
    s->locals.reserve(16);
    s->active_def.assign(n_defs, 0);
    s->cached_def.assign(n_defs, 0);
    s->cache_def.assign(n_defs, 0.0);
    s->depth = 0;
}

void scratch_reset_eval(Scratch *s) {
    std::fill(s->cached_def.begin(), s->cached_def.end(), 0);
    /* active_def is reset by every call returning; if a previous
     * run failed mid-stream we want to clear too. */
    std::fill(s->active_def.begin(), s->active_def.end(), 0);
    s->stack.clear();
    s->locals.clear();
    s->depth = 0;
}

namespace {

/* The hot path. `frame_base` indexes the first arg slot of the
 * current call frame within scratch.locals; the top-level run()
 * starts with frame_base = 0. */
bool exec(const Program &program, const RunContext &ctx, Scratch &scratch,
         size_t frame_base, char *err_buf, size_t err_cap) {
    const Instr *code = program.code.data();
    const size_t n    = program.code.size();
    const double *constants = program.constants.data();
    auto       &stack = scratch.stack;

    for (size_t pc = 0; pc < n; ++pc) {
        const Instr ins = code[pc];
        switch (ins.op) {
        case Op::PushConst: stack.push_back(constants[ins.a]); break;
        case Op::PushState:
            if (ins.a >= ctx.n_state) {
                set_err_buf(err_buf, err_cap, "state index out of range");
                return false;
            }
            stack.push_back(ctx.state[ins.a]);
            break;
        case Op::PushParam:
            if (ins.a >= ctx.n_params) {
                set_err_buf(err_buf, err_cap, "param index out of range");
                return false;
            }
            stack.push_back(ctx.params[ins.a]);
            break;
        case Op::PushLocal: {
            const size_t idx = frame_base + ins.a;
            if (idx >= scratch.locals.size()) {
                set_err_buf(err_buf, err_cap, "local index out of range");
                return false;
            }
            stack.push_back(scratch.locals[idx]);
            break;
        }
        case Op::PushT:  stack.push_back(ctx.t);  break;
        case Op::PushPi: stack.push_back(M_PI);   break;
        case Op::PushE:  stack.push_back(M_E);    break;

        case Op::Neg:
            stack.back() = -stack.back();
            break;
        case Op::Add: {
            const double b = stack.back(); stack.pop_back();
            stack.back() += b;
            break;
        }
        case Op::Sub: {
            const double b = stack.back(); stack.pop_back();
            stack.back() -= b;
            break;
        }
        case Op::Mul: {
            const double b = stack.back(); stack.pop_back();
            stack.back() *= b;
            break;
        }
        case Op::Div: {
            const double b = stack.back(); stack.pop_back();
            stack.back() /= b;
            break;
        }

        case Op::CallBuiltin: {
            const Builtin id = static_cast<Builtin>(ins.a);
            const size_t  argc = ins.b;
            if (stack.size() < argc) {
                set_err_buf(err_buf, err_cap, "stack underflow in builtin call");
                return false;
            }
            double *args = stack.data() + (stack.size() - argc);
            double r = 0.0;
            switch (id) {
            case Builtin::Sin:   r = std::sin(args[0]); break;
            case Builtin::Cos:   r = std::cos(args[0]); break;
            case Builtin::Tan:   r = std::tan(args[0]); break;
            case Builtin::Asin:  r = std::asin(args[0]); break;
            case Builtin::Acos:  r = std::acos(args[0]); break;
            case Builtin::Atan:  r = std::atan(args[0]); break;
            case Builtin::Exp:   r = std::exp(args[0]); break;
            case Builtin::Log:   r = std::log(args[0]); break;
            case Builtin::Log10: r = std::log10(args[0]); break;
            case Builtin::Sqrt:  r = std::sqrt(args[0]); break;
            case Builtin::Abs:   r = std::fabs(args[0]); break;
            case Builtin::Floor: r = std::floor(args[0]); break;
            case Builtin::Ceil:  r = std::ceil(args[0]); break;
            case Builtin::Sign:  r = (args[0] > 0.0) - (args[0] < 0.0); break;
            case Builtin::Pow:   r = std::pow(args[0], args[1]); break;
            case Builtin::Min:   r = std::fmin(args[0], args[1]); break;
            case Builtin::Max:   r = std::fmax(args[0], args[1]); break;
            case Builtin::Mod:   r = std::fmod(args[0], args[1]); break;
            case Builtin::Clamp: r = std::fmax(args[1], std::fmin(args[2], args[0])); break;
            case Builtin::Unknown:
                set_err_buf(err_buf, err_cap, "unknown builtin id %u",
                            static_cast<unsigned>(ins.a));
                return false;
            }
            stack.resize(stack.size() - argc);
            stack.push_back(r);
            break;
        }

        case Op::CallDef: {
            const uint16_t def_idx = ins.a;
            const uint16_t argc    = ins.b;
            if (def_idx >= ctx.n_defs) {
                set_err_buf(err_buf, err_cap, "def index out of range");
                return false;
            }
            if (stack.size() < argc) {
                set_err_buf(err_buf, err_cap, "stack underflow in def call");
                return false;
            }

            /* Memoized 0-arity defs: hit cache before entering recursion. */
            if (argc == 0 && scratch.cached_def[def_idx]) {
                stack.push_back(scratch.cache_def[def_idx]);
                break;
            }
            if (scratch.active_def[def_idx]) {
                set_err_buf(err_buf, err_cap,
                            "cyclic definition involving def#%u", def_idx);
                return false;
            }
            if (scratch.depth > 64) {
                set_err_buf(err_buf, err_cap,
                            "call depth exceeded (def#%u)", def_idx);
                return false;
            }
            /* Move args from value stack onto locals stack in order:
             * stack ...,a,b,c   ->   locals ...,a,b,c; stack ... */
            const size_t callee_frame_base = scratch.locals.size();
            for (uint16_t i = 0; i < argc; ++i) {
                scratch.locals.push_back(stack[stack.size() - argc + i]);
            }
            stack.resize(stack.size() - argc);

            scratch.active_def[def_idx] = 1;
            scratch.depth += 1;
            const bool ok = exec(ctx.defs[def_idx], ctx, scratch,
                                 callee_frame_base, err_buf, err_cap);
            scratch.depth -= 1;
            scratch.active_def[def_idx] = 0;
            scratch.locals.resize(scratch.locals.size() - argc);
            if (!ok) return false;

            if (argc == 0) {
                scratch.cached_def[def_idx] = 1;
                scratch.cache_def[def_idx]  = stack.back();
            }
            break;
        }

        case Op::BrIfZero: {
            const double v = stack.back(); stack.pop_back();
            if (v == 0.0) pc += ins.a;
            break;
        }
        case Op::Jump:
            pc += ins.a;
            break;
        }
    }
    return true;
}

}  /* namespace */

bool run(const Program  &program,
         const RunContext &ctx,
         Scratch          &scratch,
         double           *out,
         char             *err_buf,
         size_t            err_cap) {
    const size_t sp_before = scratch.stack.size();
    if (!exec(program, ctx, scratch, /*frame_base=*/0, err_buf, err_cap)) return false;
    if (scratch.stack.size() != sp_before + 1) {
        set_err_buf(err_buf, err_cap,
                    "internal: stack imbalance (was %zu now %zu)",
                    sp_before, scratch.stack.size());
        /* Roll back to a clean state for next call. */
        scratch.stack.resize(sp_before);
        return false;
    }
    *out = scratch.stack.back();
    scratch.stack.pop_back();
    return true;
}

const char *builtin_name(Builtin b) {
    for (const auto &spec : kBuiltins) {
        if (spec.id == b) return spec.name;
    }
    return "?";
}

Builtin builtin_from_name(const char *name) {
    const BuiltinSpec *s = find_builtin(name);
    return s ? s->id : Builtin::Unknown;
}

}  /* namespace dynsys::ir */
