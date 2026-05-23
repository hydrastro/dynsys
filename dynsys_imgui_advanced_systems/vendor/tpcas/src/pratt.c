/* ============================================================
 * Pratt parser (top-down operator precedence).
 *
 * Each token has a "null denotation" (nud) — how it's parsed
 * as a prefix — and a "left denotation" (led) — how it's parsed
 * after a left-hand-side has been built.
 *
 * The main loop:
 *
 *   parse_expr(rbp):
 *       t = consume();  left = nud(t)
 *       while infix_lbp(cur) > rbp:
 *           t = consume();  left = led(t, left)
 *       return left
 *
 * Associativity-as-direction (the trick from the original
 * tpcas) is encoded by choosing the recursive call's rbp:
 *   - LEFT-assoc:  recurse with rbp = op.precedence
 *                  (same-precedence sibling exits → folds left)
 *   - RIGHT-assoc: recurse with rbp = op.precedence - 1
 *                  (same-precedence sibling re-enters → folds right)
 * Same idea as the old multi-pass sweep, one pass.
 *
 * Function application f(args...) is handled as a postfix
 * "(" with very high binding power, so juxtaposition of an
 * identifier and "(" parses naturally. This is what made
 * extending to FOL/HOL impossible before.
 * ============================================================ */

#include "pratt.h"
#include "lex.h"
#include <string.h>
#include <stdio.h>

#define BP_APPLY 100   /* binding power of "(" in led position */

typedef struct {
    lexer_t    *lex;
    token_t     cur;
    arena_t    *arena;
    bool        had_error;
    span_t      err_span;
    const char *err_msg;
} parser_t;

static void p_advance(parser_t *p) { p->cur = lex_next(p->lex); }

static void p_error(parser_t *p, span_t s, const char *msg) {
    if (!p->had_error) { p->had_error = true; p->err_span = s; p->err_msg = msg; }
}

static span_t span_join(span_t a, span_t b) {
    span_t r = a;
    r.end = b.end > a.end ? b.end : a.end;
    return r;
}

/* forward decls */
static node_t  *parse_expr(parser_t *p, int rbp);
static node_t  *parse_nud (parser_t *p, token_t t);
static node_t  *parse_led (parser_t *p, token_t t, node_t *left);
static node_t  *parse_binder_intro(parser_t *p, node_kind_t kind, span_t start_span);
static type_t  *parse_type(parser_t *p);
static type_t  *parse_type_atom(parser_t *p);

static int infix_lbp(token_t t) {
    if (t.kind == TOK_LPAREN) return BP_APPLY;
    if (t.kind == TOK_OP && t.op->fixity == FIXITY_INFIX) return t.op->precedence;
    return 0;
}

/* ============================================================ */
/* nud — prefix / atom dispatch                                 */
/* ============================================================ */
static node_t *parse_nud(parser_t *p, token_t t) {
    switch (t.kind) {
        case TOK_OP:
            /* Arithmetic negation is syntactic sugar for 0 - expr.  The AST
             * still has only the ordinary binary subtraction operator. */
            if (t.op == &OP_SUB) {
                node_t *arg = parse_expr(p, OP_NOT.precedence);
                if (!arg) return NULL;
                node_t *zero = ast_const(p->arena, "0", NULL, t.span);
                node_t *head = ast_const(p->arena, OP_SUB.name, &OP_SUB, t.span);
                return ast_app2(p->arena, head, zero, arg, span_join(t.span, arg->span));
            }
            if (t.op->fixity != FIXITY_PREFIX) {
                p_error(p, t.span, "operator used in prefix position is not prefix");
                return NULL;
            }
            {
                node_t *arg = parse_expr(p, t.op->precedence);
                if (!arg) return NULL;
                node_t *head = ast_const(p->arena, t.op->name, t.op, t.span);
                return ast_app1(p->arena, head, arg, span_join(t.span, arg->span));
            }

        case TOK_LPAREN: {
            node_t *e = parse_expr(p, 0);
            if (!e) return NULL;
            if (p->cur.kind != TOK_RPAREN) {
                p_error(p, p->cur.span, "expected ')'");
                return NULL;
            }
            p_advance(p);
            return e;
        }

        case TOK_BACKSLASH:
            return parse_binder_intro(p, NODE_LAMBDA, t.span);

        case TOK_IDENT:
            if (kw_is_true(t.lexeme) || kw_is_false(t.lexeme))
                return ast_const(p->arena, t.lexeme, NULL, t.span);
            if (kw_is_forall(t.lexeme)) return parse_binder_intro(p, NODE_FORALL, t.span);
            if (kw_is_exists(t.lexeme)) return parse_binder_intro(p, NODE_EXISTS, t.span);
            if (kw_is_lambda(t.lexeme)) return parse_binder_intro(p, NODE_LAMBDA, t.span);
            return ast_var(p->arena, t.lexeme, t.span);

        case TOK_NUMBER:
            return ast_const(p->arena, t.lexeme, NULL, t.span);

        case TOK_EOF:
            p_error(p, t.span, "unexpected end of input");
            return NULL;
        case TOK_ERROR:
            p_error(p, t.span, t.err_msg ? t.err_msg : "lexer error");
            return NULL;
        default:
            p_error(p, t.span, "unexpected token");
            return NULL;
    }
}

/* ============================================================ */
/* led — infix / postfix dispatch                                */
/* ============================================================ */
static node_t *parse_led(parser_t *p, token_t t, node_t *left) {
    if (t.kind == TOK_LPAREN) {
        /* function application: left(arg, arg, ...) */
        node_t **args = NULL;
        size_t cap = 0, count = 0;
        if (p->cur.kind != TOK_RPAREN) {
            for (;;) {
                node_t *a = parse_expr(p, 0);
                if (!a) return NULL;
                if (count >= cap) {
                    size_t ncap = cap ? cap * 2 : 4;
                    node_t **na = arena_alloc(p->arena, sizeof(node_t *) * ncap);
                    for (size_t i = 0; i < count; i++) na[i] = args[i];
                    args = na;
                    cap  = ncap;
                }
                args[count++] = a;
                if (p->cur.kind == TOK_COMMA) { p_advance(p); continue; }
                break;
            }
        }
        if (p->cur.kind != TOK_RPAREN) {
            p_error(p, p->cur.span, "expected ',' or ')' in call");
            return NULL;
        }
        span_t end_span = p->cur.span;
        p_advance(p);
        return ast_app(p->arena, left, count, args, span_join(left->span, end_span));
    }

    if (t.kind == TOK_OP && t.op->fixity == FIXITY_INFIX) {
        int next_rbp = (t.op->assoc == ASSOC_RIGHT)
                         ? t.op->precedence - 1
                         : t.op->precedence;
        node_t *right = parse_expr(p, next_rbp);
        if (!right) return NULL;
        node_t *head = ast_const(p->arena, t.op->name, t.op, t.span);
        return ast_app2(p->arena, head, left, right,
                        span_join(left->span, right->span));
    }

    p_error(p, t.span, "unexpected token in infix position");
    return NULL;
}

/* ============================================================ */
/* main expression loop                                          */
/* ============================================================ */
static node_t *parse_expr(parser_t *p, int rbp) {
    token_t t = p->cur;
    p_advance(p);
    node_t *left = parse_nud(p, t);
    if (!left) return NULL;
    while (!p->had_error && infix_lbp(p->cur) > rbp) {
        token_t op = p->cur;
        p_advance(p);
        left = parse_led(p, op, left);
        if (!left) return NULL;
    }
    return left;
}

/* ============================================================ */
/* binders: forall x[:T] [y[:U] ...] . body                      */
/* Multiple vars sugar to nested binders (right-folded).        */
/* ============================================================ */
static node_t *parse_binder_intro(parser_t *p, node_kind_t kind, span_t start_span) {
    /* gather: at least one identifier, each optionally typed */
    enum { MAX_VARS = 32 };
    const char *names[MAX_VARS];
    type_t     *types[MAX_VARS];
    size_t      n = 0;

    while (p->cur.kind == TOK_IDENT && !kw_is_any(p->cur.lexeme)) {
        if (n >= MAX_VARS) {
            p_error(p, p->cur.span, "too many bound variables");
            return NULL;
        }
        names[n] = p->cur.lexeme;
        types[n] = NULL;
        p_advance(p);
        if (p->cur.kind == TOK_COLON) {
            p_advance(p);
            types[n] = parse_type(p);
            if (!types[n]) return NULL;
        }
        n++;
    }
    if (n == 0) {
        p_error(p, p->cur.span, "expected variable name after binder");
        return NULL;
    }
    if (p->cur.kind != TOK_DOT) {
        p_error(p, p->cur.span, "expected '.' after binder vars");
        return NULL;
    }
    p_advance(p);

    node_t *body = parse_expr(p, 0);
    if (!body) return NULL;

    /* right-fold:  ∀ x y. P  ≡  ∀ x. ∀ y. P */
    node_t *acc = body;
    for (size_t i = n; i > 0; i--) {
        acc = ast_binder(p->arena, kind, names[i - 1], types[i - 1],
                         acc, span_join(start_span, acc->span));
    }
    return acc;
}

/* ============================================================ */
/* types:  atom -> type    (right-assoc)                         */
/*         atom = IDENT | ( type )                               */
/* ============================================================ */
static type_t *parse_type_atom(parser_t *p) {
    if (p->cur.kind == TOK_LPAREN) {
        p_advance(p);
        type_t *t = parse_type(p);
        if (!t) return NULL;
        if (p->cur.kind != TOK_RPAREN) {
            p_error(p, p->cur.span, "expected ')' in type");
            return NULL;
        }
        p_advance(p);
        return t;
    }
    if (p->cur.kind == TOK_IDENT) {
        type_t *t = type_base(p->arena, p->cur.lexeme);
        p_advance(p);
        return t;
    }
    p_error(p, p->cur.span, "expected type");
    return NULL;
}

static type_t *parse_type(parser_t *p) {
    type_t *lhs = parse_type_atom(p);
    if (!lhs) return NULL;
    if (p->cur.kind == TOK_ARROW_T) {
        p_advance(p);
        type_t *rhs = parse_type(p);
        if (!rhs) return NULL;
        return type_arrow(p->arena, lhs, rhs);
    }
    return lhs;
}

/* ============================================================ */
/* entry point                                                   */
/* ============================================================ */
parse_result_t parse(const char *src, arena_t *arena) {
    lexer_t lex;
    lex_init(&lex, src, arena);
    parser_t p;
    memset(&p, 0, sizeof p);
    p.lex   = &lex;
    p.arena = arena;
    p.cur   = lex_next(&lex);

    parse_result_t r;
    memset(&r, 0, sizeof r);

    node_t *ast = parse_expr(&p, 0);
    if (p.had_error) {
        r.ok = false;
        r.err_span = p.err_span;
        r.err_msg  = p.err_msg;
        return r;
    }
    if (p.cur.kind != TOK_EOF) {
        r.ok = false;
        r.err_span = p.cur.span;
        r.err_msg  = "unexpected token after expression";
        return r;
    }
    r.ok  = true;
    r.ast = ast;
    return r;
}
