# Comparison against the previous dynsys v4 package

This tree is a cleaned merge of the uploaded dynsys package and the earlier `dynsys_laboratory_v4` package.

## Uploaded package strengths

- Contains the full v4 pplane-style interaction layer.
- Adds `src/expr_ir.{h,cpp}`, a lowered expression bytecode/VM for hot-path evaluation.
- Adds `--headless` differential/benchmark mode.
- Adds `test/ir_smoke.cpp`, which validates operator precedence, builtins, user functions, lazy `select`, error paths, and recursion detection.
- Preallocates RK scratch states, removing vector allocation churn in the RK4 hot path.

## Uploaded package problems

- Still contained legacy root-level `dynsys.c` using the old Lizard/FreeType path.
- Still contained old artifacts: `FSEX302.ttf`, `result.png`, `imgui.ini`.
- README roadmap still listed some v4 features as future work.
- Makefile did not expose the IR smoke test through `make test`.

## Earlier v4 package strengths

- Cleaner package shape: no legacy FreeType/Lizard root implementation.
- Smaller, easier-to-audit tree.
- v4 pplane workbench was coherent and scoped correctly.

## Earlier v4 package limitations

- Used the AST evaluator in hot loops.
- No standalone IR test target.
- No headless benchmark/differential testing path.

## Best-of-both result

This package keeps the uploaded optimized implementation and removes the stale garbage. It adds Makefile test/benchmark convenience targets and updates docs so the tree reflects the current architecture.
