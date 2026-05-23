#include "print.h"
#include <stdio.h>

/* ============================================================
 * Precedence-aware printing.
 *   ctx_prec = the precedence the surrounding context expects.
 *   If this node's outer precedence < ctx_prec → parenthesise.
 *
 * For an infix op at precedence P with associativity A:
 *   left child:  ctx = (A == LEFT)  ? P    : P + 1
 *   right child: ctx = (A == RIGHT) ? P    : P + 1
 * This is the standard trick: equal precedence is allowed on
 * the "tight" side, requires parens on the other.
 * ============================================================ */

#define PREC_BINDER 10
#define PREC_ATOM   1000

void print_type(FILE *out, const type_t *t) {
    if (!t) { fputs("?", out); return; }
    switch (t->kind) {
        case TYPE_BASE: fputs(t->base.name, out); return;
        case TYPE_VAR:  fputs(t->var.name,  out); return;
        case TYPE_ARROW:
            /* arrow is right-assoc; parenthesise non-atomic domain */
            if (t->arrow.dom->kind == TYPE_ARROW) {
                fputc('(', out);
                print_type(out, t->arrow.dom);
                fputc(')', out);
            } else {
                print_type(out, t->arrow.dom);
            }
            fputs(" -> ", out);
            print_type(out, t->arrow.cod);
            return;
    }
}

static void rec(FILE *out, const node_t *n, int ctx_prec);

static void print_app_args(FILE *out, const node_t *n) {
    fputc('(', out);
    for (size_t i = 0; i < n->app.argc; i++) {
        if (i) fputs(", ", out);
        rec(out, n->app.args[i], 0);
    }
    fputc(')', out);
}

static void rec(FILE *out, const node_t *n, int ctx_prec) {
    if (!n) { fputs("<null>", out); return; }
    switch (n->kind) {
        case NODE_CONST:
            fputs(n->cnst.name, out);
            return;

        case NODE_VAR:
            fputs(n->var.name, out);
            return;

        case NODE_APP:
            if (n->app.head->kind == NODE_CONST && n->app.head->cnst.op) {
                const op_info_t *op = n->app.head->cnst.op;
                if (op->fixity == FIXITY_PREFIX && n->app.argc == 1) {
                    bool wrap = op->precedence < ctx_prec;
                    if (wrap) fputc('(', out);
                    fputs(op->syntax, out);
                    rec(out, n->app.args[0], op->precedence);
                    if (wrap) fputc(')', out);
                    return;
                }
                if (op->fixity == FIXITY_INFIX && n->app.argc == 2) {
                    bool wrap = op->precedence < ctx_prec;
                    if (wrap) fputc('(', out);
                    int lprec = (op->assoc == ASSOC_LEFT)  ? op->precedence : op->precedence + 1;
                    int rprec = (op->assoc == ASSOC_RIGHT) ? op->precedence : op->precedence + 1;
                    rec(out, n->app.args[0], lprec);
                    fputc(' ', out); fputs(op->syntax, out); fputc(' ', out);
                    rec(out, n->app.args[1], rprec);
                    if (wrap) fputc(')', out);
                    return;
                }
            }
            /* generic application f(a, b, ...). Wrap head if non-atomic */
            if (n->app.head->kind == NODE_VAR || n->app.head->kind == NODE_CONST) {
                rec(out, n->app.head, PREC_ATOM);
            } else {
                fputc('(', out);
                rec(out, n->app.head, 0);
                fputc(')', out);
            }
            print_app_args(out, n);
            return;

        case NODE_LAMBDA:
        case NODE_FORALL:
        case NODE_EXISTS: {
            const char *kw = (n->kind == NODE_LAMBDA) ? "\\" :
                             (n->kind == NODE_FORALL) ? "forall" : "exists";
            bool wrap = PREC_BINDER < ctx_prec;
            if (wrap) fputc('(', out);
            fputs(kw, out);
            if (n->kind != NODE_LAMBDA) fputc(' ', out);
            fputs(n->bind.bvar, out);
            if (n->bind.btype) { fputs(":", out); print_type(out, n->bind.btype); }
            fputs(". ", out);
            rec(out, n->bind.body, 0);
            if (wrap) fputc(')', out);
            return;
        }
    }
}

void print_expr(FILE *out, const node_t *n) {
    rec(out, n, 0);
}

/* ----------- indented tree view ----------- */
static void indent(FILE *out, int depth) {
    for (int i = 0; i < depth; i++) fputs("  ", out);
}

static void tree_rec(FILE *out, const node_t *n, int depth) {
    indent(out, depth);
    if (!n) { fputs("<null>\n", out); return; }
    switch (n->kind) {
        case NODE_CONST:
            fprintf(out, "CONST %s%s\n", n->cnst.name, n->cnst.op ? "  [op]" : "");
            return;
        case NODE_VAR:
            fprintf(out, "VAR %s\n", n->var.name);
            return;
        case NODE_APP:
            fprintf(out, "APP (argc=%zu)\n", n->app.argc);
            tree_rec(out, n->app.head, depth + 1);
            for (size_t i = 0; i < n->app.argc; i++)
                tree_rec(out, n->app.args[i], depth + 1);
            return;
        case NODE_LAMBDA:
        case NODE_FORALL:
        case NODE_EXISTS: {
            const char *kw = (n->kind == NODE_LAMBDA) ? "LAMBDA" :
                             (n->kind == NODE_FORALL) ? "FORALL" : "EXISTS";
            fprintf(out, "%s %s", kw, n->bind.bvar);
            if (n->bind.btype) { fputs(" : ", out); print_type(out, n->bind.btype); }
            fputc('\n', out);
            tree_rec(out, n->bind.body, depth + 1);
            return;
        }
    }
}

void print_expr_indent(FILE *out, const node_t *n) {
    tree_rec(out, n, 0);
}
