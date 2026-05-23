# vendored data structures subset

This directory contains the pieces of the custom data-structures library used by TPCAS:

- `context` + `allocators` for the arena-backed allocation layer
- `str` for dynamic error-message construction in the parser-combinator runner
- supporting `common`, `status`, `error`, and `diagnostic` modules

Build artifacts, tests, Unicode data files, and unrelated containers from the original library archive were intentionally omitted from this standalone TPCAS package.
