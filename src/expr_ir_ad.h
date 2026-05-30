#pragma once

/* ============================================================
 * dynsys forward-mode automatic differentiation over the IR
 * (roadmap phase 3).
 *
 * Evaluates a lowered Program carrying dual numbers (value,
 * derivative) so that one pass yields both f and the directional
 * derivative df in the direction of a chosen seed. Running it once
 * per state variable (seed = e_i on that state) gives an exact
 * Jacobian column; seeding the parameter gives df/dp exactly.
 *
 * This is a separate executor from expr_ir's run(): the hot
 * integration path stays the plain double evaluator, untouched.
 * The AD path is only used by the analysis layer, where exactness
 * matters more than per-step speed. It mirrors exec() opcode for
 * opcode, including CallDef recursion and the select() branch ops,
 * so its value component is identical to run()'s.
 *
 * Derivative semantics for non-smooth builtins (abs, floor, ceil,
 * sign, min, max, mod, clamp) follow the usual AD conventions:
 * subgradient 0 at kinks, derivative of the selected branch for
 * min/max/clamp, derivative of the dividend for mod. These match
 * what a finite-difference Jacobian converges to away from the
 * measure-zero kink set, and are documented at each site.
 * ============================================================ */

#include <cstddef>

#include "expr_ir.h"

namespace dynsys::ir {

/* A value/derivative pair flowing through the dual-number VM. */
struct Dual {
  double v = 0.0;  /* primal value      */
  double d = 0.0;  /* derivative w.r.t. the active seed */
};

/* Seed selecting which input carries derivative 1.0.
 *   kind == State: d/d(state[index])
 *   kind == Param: d/d(param[index])
 *   kind == Time : d/dt
 * All other inputs are seeded with derivative 0. */
struct DualSeed {
  enum class Kind { State, Param, Time } kind = Kind::State;
  std::size_t index = 0;
};

/* Per-eval scratch for the dual executor, analogous to ir::Scratch
 * but holding Dual values. Owned by the caller so allocations don't
 * churn across repeated Jacobian-column evaluations. */
struct DualScratch {
  std::vector<Dual> stack;
  std::vector<Dual> locals;
  std::vector<std::uint8_t> active_def;
  std::vector<std::uint8_t> cached_def;
  std::vector<Dual> cache_def;
  int depth = 0;
};

void dual_scratch_init(DualScratch *s, std::size_t n_defs);
void dual_scratch_reset_eval(DualScratch *s);

/* Evaluate `program` in dual-number arithmetic. Returns the
 * primal in out_value and the directional derivative (w.r.t. the
 * seed) in out_deriv. The RunContext is the same one used by the
 * plain evaluator; the seed picks the differentiation direction. */
bool run_dual(const Program &program, const RunContext &ctx,
              const DualSeed &seed, DualScratch &scratch, double *out_value,
              double *out_deriv, char *err_buf, std::size_t err_cap);

}  // namespace dynsys::ir
