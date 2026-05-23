#ifndef AST_H
#define AST_H

#include <stdbool.h>
#include <stddef.h>
#include "arena.h"

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Source positions — every token and AST node carries one,
 * so errors can be reported with caret diagnostics.
 * ============================================================ */
typedef struct {
    size_t start;   /* byte offset in source */
    size_t end;
    size_t line;    /* 1-based */
    size_t col;     /* 1-based */
} span_t;

/* ============================================================
 * Types (for HOL).  In simple typed lambda calculus you need
 *   base types (Bool, Nat, ...) and function types  T -> U.
 * Type variables are kept around in case we want polymorphism
 * later. Types are stored as a separate small tree; for now
 * they're parsed and printed but not enforced (no type checker).
 * ============================================================ */
typedef struct type type_t;

typedef enum {
    TYPE_BASE,
    TYPE_ARROW,
    TYPE_VAR
} type_kind_t;

struct type {
    type_kind_t kind;
    union {
        struct { const char *name; } base;
        struct { type_t *dom, *cod; } arrow;
        struct { const char *name; } var;
    };
};

type_t *type_base (arena_t *a, const char *name);
type_t *type_arrow(arena_t *a, type_t *dom, type_t *cod);
type_t *type_var  (arena_t *a, const char *name);

/* ============================================================
 * Operator metadata.  Single source of truth for syntax,
 * precedence and associativity.  Adding XOR is a one-liner.
 * ============================================================ */
typedef enum { ASSOC_LEFT, ASSOC_RIGHT, ASSOC_NONE } assoc_t;
typedef enum { FIXITY_PREFIX, FIXITY_INFIX, FIXITY_POSTFIX } fixity_t;

typedef struct {
    const char *name;        /* canonical name in AST: "AND", "NOT", ... */
    const char *syntax;      /* surface syntax: "&&" */
    int         arity;
    int         precedence;  /* higher = binds tighter */
    assoc_t     assoc;
    fixity_t    fixity;
} op_info_t;

/* Pre-defined operators (extern so transforms can compare by &OP_AND etc.) */
extern const op_info_t OP_NOT;
extern const op_info_t OP_AND;
extern const op_info_t OP_OR;
extern const op_info_t OP_IMP;
extern const op_info_t OP_IFF;
extern const op_info_t OP_EQ;     /* equality, useful for FOL/HOL */
extern const op_info_t OP_ADD;    /* arithmetic, useful for embedding TPCAS as an expression parser */
extern const op_info_t OP_SUB;
extern const op_info_t OP_MUL;
extern const op_info_t OP_DIV;

const op_info_t *op_lookup(const char *syntax);
const op_info_t *const *op_all(size_t *count);

/* ============================================================
 * AST.  Unified for PL, FOL and HOL:
 *   - NODE_CONST    : named constant (true, false, 0, MyConst)
 *                     If it's a built-in operator, .cnst.op != NULL.
 *   - NODE_VAR      : variable (free or bound — we don't pre-resolve)
 *   - NODE_APP      : application f(a, b, ...). Covers operator
 *                     applications (head is a const with op != NULL),
 *                     function calls (FOL), predicate calls (FOL),
 *                     and higher-order application (HOL).
 *   - NODE_LAMBDA   : \x[:T]. body          (HOL)
 *   - NODE_FORALL   : forall x[:T]. body    (FOL/HOL)
 *   - NODE_EXISTS   : exists x[:T]. body    (FOL/HOL)
 *
 * Binders hold ONE variable.  Multi-variable binders sugar to
 * nested binders at parse time.
 * ============================================================ */
typedef struct node node_t;

typedef enum {
    NODE_CONST,
    NODE_VAR,
    NODE_APP,
    NODE_LAMBDA,
    NODE_FORALL,
    NODE_EXISTS
} node_kind_t;

struct node {
    node_kind_t kind;
    span_t      span;
    type_t     *type;          /* optional; for HOL annotations */
    union {
        struct { const char *name; const op_info_t *op; } cnst;
        struct { const char *name; } var;
        struct { node_t *head; size_t argc; node_t **args; } app;
        struct { const char *bvar; type_t *btype; node_t *body; } bind;
    };
};

/* ---- constructors (arena-allocated) ---- */
node_t *ast_const  (arena_t *a, const char *name, const op_info_t *op, span_t span);
node_t *ast_var    (arena_t *a, const char *name, span_t span);
node_t *ast_app    (arena_t *a, node_t *head, size_t argc, node_t **args, span_t span);
node_t *ast_app1   (arena_t *a, node_t *head, node_t *a1, span_t span);
node_t *ast_app2   (arena_t *a, node_t *head, node_t *a1, node_t *a2, span_t span);
node_t *ast_lambda (arena_t *a, const char *bvar, type_t *btype, node_t *body, span_t span);
node_t *ast_forall (arena_t *a, const char *bvar, type_t *btype, node_t *body, span_t span);
node_t *ast_exists (arena_t *a, const char *bvar, type_t *btype, node_t *body, span_t span);
node_t *ast_binder (arena_t *a, node_kind_t kind, const char *bvar, type_t *btype, node_t *body, span_t span);

/* ---- query helpers ---- */
bool ast_is_op_app(const node_t *n, const op_info_t *op, int argc);
bool ast_equal     (const node_t *a, const node_t *b);  /* alpha-eq */
node_t *ast_copy   (arena_t *a, const node_t *n);

/* ---- generic rewriters.  Return non-NULL to replace, NULL to keep. ---- */
typedef node_t *(*ast_rewrite_fn)(arena_t *a, node_t *n, void *ctx);

node_t *ast_rewrite_bottomup(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx);
node_t *ast_rewrite_topdown (arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx);
node_t *ast_rewrite_fixpoint(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx, int max_iters);

/* free-variable check */
bool ast_free_in(const node_t *n, const char *name);


#ifdef __cplusplus
}
#endif
#endif
