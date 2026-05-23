#include "ast.h"
#include <string.h>

/* ============================================================
 * Operator table — higher precedence binds tighter.
 * Numbers chosen so there's room to insert new operators between
 * existing ones without a global renumbering.
 * ============================================================ */

const op_info_t OP_NOT = {"NOT", "!",   1, 90, ASSOC_NONE,  FIXITY_PREFIX};
const op_info_t OP_MUL = {"MUL", "*",   2, 85, ASSOC_LEFT,  FIXITY_INFIX };
const op_info_t OP_DIV = {"DIV", "/",   2, 85, ASSOC_LEFT,  FIXITY_INFIX };
const op_info_t OP_ADD = {"ADD", "+",   2, 80, ASSOC_LEFT,  FIXITY_INFIX };
const op_info_t OP_SUB = {"SUB", "-",   2, 80, ASSOC_LEFT,  FIXITY_INFIX };
const op_info_t OP_AND = {"AND", "&&",  2, 70, ASSOC_LEFT,  FIXITY_INFIX };
const op_info_t OP_OR  = {"OR",  "||",  2, 60, ASSOC_LEFT,  FIXITY_INFIX };
const op_info_t OP_EQ  = {"EQ",  "=",   2, 50, ASSOC_NONE,  FIXITY_INFIX };
const op_info_t OP_IMP = {"IMP", "=>",  2, 40, ASSOC_RIGHT, FIXITY_INFIX };
const op_info_t OP_IFF = {"IFF", "<=>", 2, 30, ASSOC_NONE,  FIXITY_INFIX };

static const op_info_t *const OPS[] = {
    &OP_NOT, &OP_MUL, &OP_DIV, &OP_ADD, &OP_SUB,
    &OP_AND, &OP_OR, &OP_EQ, &OP_IMP, &OP_IFF
};

const op_info_t *op_lookup(const char *syntax) {
    size_t n = sizeof(OPS) / sizeof(OPS[0]);
    for (size_t i = 0; i < n; i++)
        if (strcmp(OPS[i]->syntax, syntax) == 0) return OPS[i];
    return NULL;
}

const op_info_t *const *op_all(size_t *count) {
    *count = sizeof(OPS) / sizeof(OPS[0]);
    return OPS;
}

/* ============================================================ */
/* Type constructors                                            */
/* ============================================================ */

type_t *type_base(arena_t *a, const char *name) {
    type_t *t = arena_alloc(a, sizeof(type_t));
    t->kind = TYPE_BASE;
    t->base.name = name;
    return t;
}
type_t *type_arrow(arena_t *a, type_t *dom, type_t *cod) {
    type_t *t = arena_alloc(a, sizeof(type_t));
    t->kind = TYPE_ARROW;
    t->arrow.dom = dom;
    t->arrow.cod = cod;
    return t;
}
type_t *type_var(arena_t *a, const char *name) {
    type_t *t = arena_alloc(a, sizeof(type_t));
    t->kind = TYPE_VAR;
    t->var.name = name;
    return t;
}

/* ============================================================ */
/* AST constructors                                             */
/* ============================================================ */

static node_t *node_alloc(arena_t *a, node_kind_t kind, span_t span) {
    node_t *n = arena_alloc(a, sizeof(node_t));
    n->kind = kind;
    n->span = span;
    n->type = NULL;
    return n;
}

node_t *ast_const(arena_t *a, const char *name, const op_info_t *op, span_t span) {
    node_t *n = node_alloc(a, NODE_CONST, span);
    n->cnst.name = name;
    n->cnst.op   = op;
    return n;
}

node_t *ast_var(arena_t *a, const char *name, span_t span) {
    node_t *n = node_alloc(a, NODE_VAR, span);
    n->var.name = name;
    return n;
}

node_t *ast_app(arena_t *a, node_t *head, size_t argc, node_t **args, span_t span) {
    node_t *n = node_alloc(a, NODE_APP, span);
    n->app.head = head;
    n->app.argc = argc;
    n->app.args = arena_alloc(a, sizeof(node_t *) * (argc ? argc : 1));
    for (size_t i = 0; i < argc; i++) n->app.args[i] = args[i];
    return n;
}

node_t *ast_app1(arena_t *a, node_t *head, node_t *a1, span_t span) {
    node_t *args[1] = {a1};
    return ast_app(a, head, 1, args, span);
}

node_t *ast_app2(arena_t *a, node_t *head, node_t *a1, node_t *a2, span_t span) {
    node_t *args[2] = {a1, a2};
    return ast_app(a, head, 2, args, span);
}

node_t *ast_binder(arena_t *a, node_kind_t kind, const char *bvar, type_t *btype,
                   node_t *body, span_t span) {
    node_t *n = node_alloc(a, kind, span);
    n->bind.bvar  = bvar;
    n->bind.btype = btype;
    n->bind.body  = body;
    return n;
}

node_t *ast_lambda(arena_t *a, const char *bvar, type_t *btype, node_t *body, span_t span) {
    return ast_binder(a, NODE_LAMBDA, bvar, btype, body, span);
}
node_t *ast_forall(arena_t *a, const char *bvar, type_t *btype, node_t *body, span_t span) {
    return ast_binder(a, NODE_FORALL, bvar, btype, body, span);
}
node_t *ast_exists(arena_t *a, const char *bvar, type_t *btype, node_t *body, span_t span) {
    return ast_binder(a, NODE_EXISTS, bvar, btype, body, span);
}

/* ============================================================ */
/* Query helpers                                                */
/* ============================================================ */

bool ast_is_op_app(const node_t *n, const op_info_t *op, int argc) {
    return n && n->kind == NODE_APP
        && n->app.head && n->app.head->kind == NODE_CONST
        && n->app.head->cnst.op == op
        && (argc < 0 || n->app.argc == (size_t)argc);
}

/* ============================================================ */
/* Alpha equivalence: two binders are equal if their bodies are  */
/* equal under a position-aware rename.                          */
/* ============================================================ */

#define EQ_MAX_DEPTH 128

static int find_bind_idx(const char **binds, int depth, const char *name) {
    for (int i = depth - 1; i >= 0; i--)
        if (binds[i] && strcmp(binds[i], name) == 0) return i;
    return -1;
}

static bool eq_helper(const node_t *a, const node_t *b,
                      const char **lb, const char **rb, int depth)
{
    if (!a || !b) return a == b;
    if (a->kind != b->kind) return false;
    switch (a->kind) {
        case NODE_CONST:
            return strcmp(a->cnst.name, b->cnst.name) == 0;
        case NODE_VAR: {
            int li = find_bind_idx(lb, depth, a->var.name);
            int ri = find_bind_idx(rb, depth, b->var.name);
            if (li != ri) return false;
            if (li == -1) return strcmp(a->var.name, b->var.name) == 0;
            return true;
        }
        case NODE_APP:
            if (a->app.argc != b->app.argc) return false;
            if (!eq_helper(a->app.head, b->app.head, lb, rb, depth)) return false;
            for (size_t i = 0; i < a->app.argc; i++)
                if (!eq_helper(a->app.args[i], b->app.args[i], lb, rb, depth)) return false;
            return true;
        case NODE_LAMBDA:
        case NODE_FORALL:
        case NODE_EXISTS:
            if (depth >= EQ_MAX_DEPTH) return false;
            lb[depth] = a->bind.bvar;
            rb[depth] = b->bind.bvar;
            bool ok = eq_helper(a->bind.body, b->bind.body, lb, rb, depth + 1);
            lb[depth] = NULL;
            rb[depth] = NULL;
            return ok;
    }
    return false;
}

bool ast_equal(const node_t *a, const node_t *b) {
    const char *lb[EQ_MAX_DEPTH] = {0}, *rb[EQ_MAX_DEPTH] = {0};
    return eq_helper(a, b, lb, rb, 0);
}

/* ============================================================ */
/* Deep copy                                                    */
/* ============================================================ */

node_t *ast_copy(arena_t *a, const node_t *n) {
    if (!n) return NULL;
    switch (n->kind) {
        case NODE_CONST: return ast_const(a, n->cnst.name, n->cnst.op, n->span);
        case NODE_VAR:   return ast_var  (a, n->var.name,  n->span);
        case NODE_APP: {
            node_t **args = arena_alloc(a, sizeof(node_t *) * (n->app.argc ? n->app.argc : 1));
            for (size_t i = 0; i < n->app.argc; i++)
                args[i] = ast_copy(a, n->app.args[i]);
            return ast_app(a, ast_copy(a, n->app.head), n->app.argc, args, n->span);
        }
        case NODE_LAMBDA:
        case NODE_FORALL:
        case NODE_EXISTS:
            return ast_binder(a, n->kind, n->bind.bvar, n->bind.btype,
                              ast_copy(a, n->bind.body), n->span);
    }
    return NULL;
}

/* ============================================================ */
/* free-variable test                                           */
/* ============================================================ */

bool ast_free_in(const node_t *n, const char *name) {
    if (!n) return false;
    switch (n->kind) {
        case NODE_CONST: return false;
        case NODE_VAR:   return strcmp(n->var.name, name) == 0;
        case NODE_APP:
            if (ast_free_in(n->app.head, name)) return true;
            for (size_t i = 0; i < n->app.argc; i++)
                if (ast_free_in(n->app.args[i], name)) return true;
            return false;
        case NODE_LAMBDA:
        case NODE_FORALL:
        case NODE_EXISTS:
            if (strcmp(n->bind.bvar, name) == 0) return false; /* shadowed */
            return ast_free_in(n->bind.body, name);
    }
    return false;
}

/* ============================================================ */
/* Rewriters                                                    */
/*                                                              */
/* fn returns NULL  ⇒  "no match here, keep the node"           */
/* fn returns node  ⇒  "replace with this"                      */
/*                                                              */
/* recurse_children rebuilds parents only when a child changed; */
/* this gives the fixpoint driver a cheap pointer-equality      */
/* convergence test.                                            */
/* ============================================================ */

static node_t *recurse_children(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx, bool topdown);

node_t *ast_rewrite_bottomup(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx) {
    if (!n) return NULL;
    n = recurse_children(a, n, fn, ctx, false);
    node_t *r = fn(a, n, ctx);
    return r ? r : n;
}

node_t *ast_rewrite_topdown(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx) {
    if (!n) return NULL;
    node_t *r = fn(a, n, ctx);
    if (r) n = r;
    return recurse_children(a, n, fn, ctx, true);
}

static node_t *recurse_children(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx, bool topdown) {
    node_t *(*rec)(arena_t *, node_t *, ast_rewrite_fn, void *) =
        topdown ? ast_rewrite_topdown : ast_rewrite_bottomup;
    switch (n->kind) {
        case NODE_CONST:
        case NODE_VAR:
            return n;
        case NODE_APP: {
            node_t *new_head = rec(a, n->app.head, fn, ctx);
            bool changed = new_head != n->app.head;
            node_t **new_args = NULL;
            if (n->app.argc) {
                new_args = arena_alloc(a, sizeof(node_t *) * n->app.argc);
                for (size_t i = 0; i < n->app.argc; i++) {
                    new_args[i] = rec(a, n->app.args[i], fn, ctx);
                    if (new_args[i] != n->app.args[i]) changed = true;
                }
            }
            if (!changed) return n;
            return ast_app(a, new_head, n->app.argc, new_args, n->span);
        }
        case NODE_LAMBDA:
        case NODE_FORALL:
        case NODE_EXISTS: {
            node_t *new_body = rec(a, n->bind.body, fn, ctx);
            if (new_body == n->bind.body) return n;
            return ast_binder(a, n->kind, n->bind.bvar, n->bind.btype, new_body, n->span);
        }
    }
    return n;
}

node_t *ast_rewrite_fixpoint(arena_t *a, node_t *n, ast_rewrite_fn fn, void *ctx, int max_iters) {
    for (int i = 0; i < max_iters; i++) {
        node_t *next = ast_rewrite_bottomup(a, n, fn, ctx);
        if (next == n) return n;
        n = next;
    }
    return n;
}
