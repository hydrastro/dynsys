/* ============================================================
 * Forward-mode AD executor over the dynsys IR. See expr_ir_ad.h.
 *
 * Mirrors expr_ir.cpp's exec() exactly in its value component and
 * propagates derivatives by the standard rules. Kept structurally
 * parallel to exec() so the two are easy to diff and keep in sync.
 * ============================================================ */

#include "expr_ir_ad.h"

#include <cmath>
#include <cstdarg>
#include <cstdio>

namespace dynsys::ir {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kE = 2.71828182845904523536;

void set_err(char *buf, std::size_t cap, const char *fmt, ...) {
  if (!buf || cap == 0) return;
  va_list ap;
  va_start(ap, fmt);
  std::vsnprintf(buf, cap, fmt, ap);
  va_end(ap);
}

bool exec_dual(const Program &program, const RunContext &ctx,
               const DualSeed &seed, DualScratch &scratch,
               std::size_t frame_base, char *err, std::size_t cap) {
  const Instr *code = program.code.data();
  const std::size_t n = program.code.size();
  const double *constants = program.constants.data();
  auto &stack = scratch.stack;

  auto push = [&](double v, double d) { stack.push_back(Dual{v, d}); };

  for (std::size_t pc = 0; pc < n; ++pc) {
    const Instr ins = code[pc];
    switch (ins.op) {
      case Op::PushConst:
        push(constants[ins.a], 0.0);
        break;
      case Op::PushState: {
        if (ins.a >= ctx.n_state) {
          set_err(err, cap, "state index out of range");
          return false;
        }
        const double d = (seed.kind == DualSeed::Kind::State &&
                          seed.index == ins.a)
                             ? 1.0
                             : 0.0;
        push(ctx.state[ins.a], d);
        break;
      }
      case Op::PushParam: {
        if (ins.a >= ctx.n_params) {
          set_err(err, cap, "param index out of range");
          return false;
        }
        const double d = (seed.kind == DualSeed::Kind::Param &&
                          seed.index == ins.a)
                             ? 1.0
                             : 0.0;
        push(ctx.params[ins.a], d);
        break;
      }
      case Op::PushLocal: {
        const std::size_t idx = frame_base + ins.a;
        if (idx >= scratch.locals.size()) {
          set_err(err, cap, "local index out of range");
          return false;
        }
        stack.push_back(scratch.locals[idx]);
        break;
      }
      case Op::PushT:
        push(ctx.t, seed.kind == DualSeed::Kind::Time ? 1.0 : 0.0);
        break;
      case Op::PushPi:
        push(kPi, 0.0);
        break;
      case Op::PushE:
        push(kE, 0.0);
        break;

      case Op::Neg:
        stack.back().v = -stack.back().v;
        stack.back().d = -stack.back().d;
        break;
      case Op::Add: {
        const Dual b = stack.back();
        stack.pop_back();
        stack.back().v += b.v;
        stack.back().d += b.d;
        break;
      }
      case Op::Sub: {
        const Dual b = stack.back();
        stack.pop_back();
        stack.back().v -= b.v;
        stack.back().d -= b.d;
        break;
      }
      case Op::Mul: {
        const Dual b = stack.back();
        stack.pop_back();
        Dual &a = stack.back();
        /* (uv)' = u'v + uv' */
        a.d = a.d * b.v + a.v * b.d;
        a.v = a.v * b.v;
        break;
      }
      case Op::Div: {
        const Dual b = stack.back();
        stack.pop_back();
        Dual &a = stack.back();
        /* (u/v)' = (u'v - uv') / v^2 */
        const double inv = 1.0 / b.v;
        a.d = (a.d * b.v - a.v * b.d) * inv * inv;
        a.v = a.v * inv;
        break;
      }

      case Op::CallBuiltin: {
        const Builtin id = static_cast<Builtin>(ins.a);
        const std::size_t argc = ins.b;
        if (stack.size() < argc) {
          set_err(err, cap, "stack underflow in builtin call");
          return false;
        }
        Dual *a = stack.data() + (stack.size() - argc);
        double rv = 0.0, rd = 0.0;
        switch (id) {
          case Builtin::Sin:
            rv = std::sin(a[0].v);
            rd = std::cos(a[0].v) * a[0].d;
            break;
          case Builtin::Cos:
            rv = std::cos(a[0].v);
            rd = -std::sin(a[0].v) * a[0].d;
            break;
          case Builtin::Tan: {
            rv = std::tan(a[0].v);
            const double sec = 1.0 / std::cos(a[0].v);
            rd = sec * sec * a[0].d;
            break;
          }
          case Builtin::Asin:
            rv = std::asin(a[0].v);
            rd = a[0].d / std::sqrt(1.0 - a[0].v * a[0].v);
            break;
          case Builtin::Acos:
            rv = std::acos(a[0].v);
            rd = -a[0].d / std::sqrt(1.0 - a[0].v * a[0].v);
            break;
          case Builtin::Atan:
            rv = std::atan(a[0].v);
            rd = a[0].d / (1.0 + a[0].v * a[0].v);
            break;
          case Builtin::Exp:
            rv = std::exp(a[0].v);
            rd = rv * a[0].d;
            break;
          case Builtin::Log:
            rv = std::log(a[0].v);
            rd = a[0].d / a[0].v;
            break;
          case Builtin::Log10:
            rv = std::log10(a[0].v);
            rd = a[0].d / (a[0].v * std::log(10.0));
            break;
          case Builtin::Sqrt:
            rv = std::sqrt(a[0].v);
            rd = rv != 0.0 ? 0.5 * a[0].d / rv : 0.0;
            break;
          case Builtin::Abs:
            rv = std::fabs(a[0].v);
            /* subgradient 0 at the kink */
            rd = (a[0].v > 0.0) ? a[0].d : (a[0].v < 0.0 ? -a[0].d : 0.0);
            break;
          case Builtin::Floor:
          case Builtin::Ceil:
          case Builtin::Sign:
            /* piecewise-constant: derivative 0 a.e. */
            rv = id == Builtin::Floor ? std::floor(a[0].v)
                 : id == Builtin::Ceil
                     ? std::ceil(a[0].v)
                     : static_cast<double>((a[0].v > 0.0) - (a[0].v < 0.0));
            rd = 0.0;
            break;
          case Builtin::Pow: {
            /* d/dx u^v = v u^{v-1} u' + u^v ln(u) v' */
            rv = std::pow(a[0].v, a[1].v);
            const double term1 =
                a[1].v * std::pow(a[0].v, a[1].v - 1.0) * a[0].d;
            const double term2 =
                (a[0].v > 0.0) ? rv * std::log(a[0].v) * a[1].d : 0.0;
            rd = term1 + term2;
            break;
          }
          case Builtin::Min:
            /* derivative follows the selected argument */
            if (a[0].v <= a[1].v) {
              rv = a[0].v;
              rd = a[0].d;
            } else {
              rv = a[1].v;
              rd = a[1].d;
            }
            break;
          case Builtin::Max:
            if (a[0].v >= a[1].v) {
              rv = a[0].v;
              rd = a[0].d;
            } else {
              rv = a[1].v;
              rd = a[1].d;
            }
            break;
          case Builtin::Mod:
            /* fmod(u,v): away from wrap points, d/dx = u' (v held). */
            rv = std::fmod(a[0].v, a[1].v);
            rd = a[0].d;
            break;
          case Builtin::Clamp: {
            /* clamp(x, lo, hi) = max(lo, min(hi, x)); follow the
             * active branch. args: a[0]=x, a[1]=lo, a[2]=hi. */
            if (a[0].v < a[1].v) {
              rv = a[1].v;
              rd = a[1].d;
            } else if (a[0].v > a[2].v) {
              rv = a[2].v;
              rd = a[2].d;
            } else {
              rv = a[0].v;
              rd = a[0].d;
            }
            break;
          }
          case Builtin::Unknown:
            set_err(err, cap, "unknown builtin id %u",
                    static_cast<unsigned>(ins.a));
            return false;
        }
        stack.resize(stack.size() - argc);
        push(rv, rd);
        break;
      }

      case Op::CallDef: {
        const std::uint16_t def_idx = ins.a;
        const std::uint16_t argc = ins.b;
        if (def_idx >= ctx.n_defs) {
          set_err(err, cap, "def index out of range");
          return false;
        }
        if (stack.size() < argc) {
          set_err(err, cap, "stack underflow in def call");
          return false;
        }
        if (argc == 0 && scratch.cached_def[def_idx]) {
          stack.push_back(scratch.cache_def[def_idx]);
          break;
        }
        if (scratch.active_def[def_idx]) {
          set_err(err, cap, "cyclic definition involving def#%u", def_idx);
          return false;
        }
        if (scratch.depth > 64) {
          set_err(err, cap, "call depth exceeded (def#%u)", def_idx);
          return false;
        }
        const std::size_t callee_frame_base = scratch.locals.size();
        for (std::uint16_t i = 0; i < argc; ++i)
          scratch.locals.push_back(stack[stack.size() - argc + i]);
        stack.resize(stack.size() - argc);

        scratch.active_def[def_idx] = 1;
        scratch.depth += 1;
        const bool ok = exec_dual(ctx.defs[def_idx], ctx, seed, scratch,
                                  callee_frame_base, err, cap);
        scratch.depth -= 1;
        scratch.active_def[def_idx] = 0;
        scratch.locals.resize(scratch.locals.size() - argc);
        if (!ok) return false;

        if (argc == 0) {
          scratch.cached_def[def_idx] = 1;
          scratch.cache_def[def_idx] = stack.back();
        }
        break;
      }

      case Op::BrIfZero: {
        const double v = stack.back().v;
        stack.pop_back();
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

}  // namespace

void dual_scratch_init(DualScratch *s, std::size_t n_defs) {
  s->stack.clear();
  s->stack.reserve(64);
  s->locals.clear();
  s->locals.reserve(16);
  s->active_def.assign(n_defs, 0);
  s->cached_def.assign(n_defs, 0);
  s->cache_def.assign(n_defs, Dual{});
  s->depth = 0;
}

void dual_scratch_reset_eval(DualScratch *s) {
  std::fill(s->cached_def.begin(), s->cached_def.end(),
            static_cast<std::uint8_t>(0));
  s->stack.clear();
  s->locals.clear();
  s->depth = 0;
}

bool run_dual(const Program &program, const RunContext &ctx,
              const DualSeed &seed, DualScratch &scratch, double *out_value,
              double *out_deriv, char *err_buf, std::size_t err_cap) {
  dual_scratch_reset_eval(&scratch);
  const std::size_t sp_before = scratch.stack.size();
  if (!exec_dual(program, ctx, seed, scratch, 0, err_buf, err_cap))
    return false;
  if (scratch.stack.size() != sp_before + 1) {
    set_err(err_buf, err_cap, "internal: dual stack imbalance");
    scratch.stack.resize(sp_before);
    return false;
  }
  const Dual top = scratch.stack.back();
  scratch.stack.pop_back();
  if (out_value) *out_value = top.v;
  if (out_deriv) *out_deriv = top.d;
  return true;
}

}  // namespace dynsys::ir
