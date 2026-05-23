/* ============================================================
 * tpcas grammar — built with the combinator library.
 *
 * Read this side-by-side with pratt.c. They produce the same
 * AST. The shape of the code is almost the inverse:
 *
 *   Pratt:        one loop driven by a precedence argument.
 *   Combinators:  one parser per precedence level, stacked.
 *
 * Direction-of-parsing reads off the wrapping:
 *   chainl1  for left-associative operators
 *   chainr1  for right-associative operators
 *   chain_none for non-associative ones
 *
 * Atom layer is one `pc_alts`; function application is a
 * postfix `atom "(" args ")"` parsed via pc_pair + pc_opt_or.
 * ============================================================ */

#include "combo.h"
#include "pc.h"
#include "ast.h"
#include <string.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================
 * Character predicates
 * ============================================================ */
static bool is_ws(int c)          { return c==' '||c=='\t'||c=='\r'||c=='\n'; }
static bool is_alpha_start(int c) { return (c>='a'&&c<='z')||(c>='A'&&c<='Z')||c=='_'; }
static bool is_alnum_more (int c) { return is_alpha_start(c) || (c>='0'&&c<='9'); }
static bool is_digit_c    (int c) { return c>='0' && c<='9'; }
static bool not_newline   (int c) { return c != '\n'; }

/* ============================================================
 * Reserved words
 * ============================================================ */
static const char *RESERVED[] = {
    "forall", "exists", "lambda", "fun", "true", "false", NULL
};
static bool is_reserved(const char *s) {
    for (size_t i = 0; RESERVED[i]; i++)
        if (strcmp(RESERVED[i], s) == 0) return true;
    return false;
}

/* ============================================================
 * Custom primitives (use the parser_fn protocol directly)
 *
 *   p_ident      : alpha (alnum)*, NOT in RESERVED, return char*
 *   p_keyword(s) : matches exactly `s` followed by a non-alnum,
 *                  so "forall" doesn't accidentally match in "forallow"
 *   p_natural    : digit+ or digit+ '.' digit+, return char*
 *
 * Implemented as small parser_fns to keep the library general.
 * They use the public pc_note_err / pc_advance helpers so error
 * accumulation works uniformly.
 * ============================================================ */

static presult_t fn_ident(const parser_t *self, pstate_t *st) {
    (void)self;
    if (st->pos >= st->len || !is_alpha_start((unsigned char)st->input[st->pos])) {
        pc_note_err(st, "identifier");
        return (presult_t){false, false, NULL};
    }
    size_t spos = st->pos, sline = st->line, scol = st->col;
    size_t start = st->pos;
    while (st->pos < st->len && is_alnum_more((unsigned char)st->input[st->pos]))
        pc_advance(st);
    char *word = arena_strndup(st->arena, st->input + start, st->pos - start);
    if (is_reserved(word)) {
        st->pos = spos; st->line = sline; st->col = scol;
        pc_note_err(st, "identifier");
        return (presult_t){false, false, NULL};
    }
    return (presult_t){true, true, word};
}
static parser_t *p_ident(arena_t *a) {
    parser_t *p = arena_alloc(a, sizeof *p);
    p->fn = fn_ident; p->state = NULL; p->desc = "identifier";
    return p;
}

typedef struct { const char *lit; size_t n; } kw_state_t;
static presult_t fn_kw(const parser_t *self, pstate_t *st) {
    kw_state_t *s = self->state;
    if (st->pos + s->n > st->len)                                   { pc_note_err(st, self->desc); return (presult_t){false,false,NULL}; }
    if (memcmp(st->input + st->pos, s->lit, s->n) != 0)             { pc_note_err(st, self->desc); return (presult_t){false,false,NULL}; }
    /* the char following must NOT be ident-continuing */
    if (st->pos + s->n < st->len
        && is_alnum_more((unsigned char)st->input[st->pos + s->n])) { pc_note_err(st, self->desc); return (presult_t){false,false,NULL}; }
    for (size_t i = 0; i < s->n; i++) pc_advance(st);
    return (presult_t){true, true, (void*)s->lit};
}
static parser_t *p_keyword(arena_t *a, const char *kw) {
    parser_t *p = arena_alloc(a, sizeof *p);
    kw_state_t *s = arena_alloc(a, sizeof *s);
    s->lit = kw; s->n = strlen(kw);
    p->fn = fn_kw; p->state = s; p->desc = kw;
    return p;
}

static presult_t fn_natural(const parser_t *self, pstate_t *st) {
    (void)self;
    if (st->pos >= st->len || !is_digit_c((unsigned char)st->input[st->pos])) {
        pc_note_err(st, "number");
        return (presult_t){false, false, NULL};
    }
    size_t start = st->pos;
    while (st->pos < st->len && is_digit_c((unsigned char)st->input[st->pos]))
        pc_advance(st);
    if (st->pos + 1 < st->len && st->input[st->pos] == '.'
        && is_digit_c((unsigned char)st->input[st->pos + 1])) {
        pc_advance(st);
        while (st->pos < st->len && is_digit_c((unsigned char)st->input[st->pos]))
            pc_advance(st);
    }
    return (presult_t){true, true, arena_strndup(st->arena, st->input + start, st->pos - start)};
}
static parser_t *p_natural(arena_t *a) {
    parser_t *p = arena_alloc(a, sizeof *p);
    p->fn = fn_natural; p->state = NULL; p->desc = "number";
    return p;
}

/* ============================================================
 * Punctuation with negative lookahead.
 *
 * The combinator parser doesn't tokenize, so operator parsers
 * run pc_string at the current position. Operators that are
 * prefixes of others ("=" prefix of "=>") need to refuse to
 * match when the forbidden continuation is present, otherwise
 * the inner-precedence operator eats the start of the outer
 * one. forbid_next is a string of single chars; if any of them
 * follows the literal in input, the match is rejected.
 * ============================================================ */
typedef struct { const char *lit; size_t n; const char *forbid_next; } punct_state_t;
static presult_t fn_punct(const parser_t *self, pstate_t *st) {
    punct_state_t *s = self->state;
    if (st->pos + s->n > st->len)                                  { pc_note_err(st, self->desc); return (presult_t){false,false,NULL}; }
    if (memcmp(st->input + st->pos, s->lit, s->n) != 0)            { pc_note_err(st, self->desc); return (presult_t){false,false,NULL}; }
    if (s->forbid_next && st->pos + s->n < st->len) {
        char nx = st->input[st->pos + s->n];
        if (strchr(s->forbid_next, nx))                            { pc_note_err(st, self->desc); return (presult_t){false,false,NULL}; }
    }
    for (size_t i = 0; i < s->n; i++) pc_advance(st);
    return (presult_t){true, true, (void*)s->lit};
}
static parser_t *p_punct(arena_t *a, const char *lit, const char *forbid_next) {
    parser_t *p = arena_alloc(a, sizeof *p);
    punct_state_t *s = arena_alloc(a, sizeof *s);
    s->lit = lit; s->n = strlen(lit); s->forbid_next = forbid_next;
    p->fn = fn_punct; p->state = s; p->desc = lit;
    return p;
}

/* ============================================================
 * Whitespace / lexeme helpers
 * ============================================================ */
static parser_t *make_ws(arena_t *a) {
    parser_t *ws_chr  = pc_satisfy(a, is_ws, "whitespace");
    parser_t *hash    = pc_char(a, '#');
    parser_t *comment_body = pc_many(a, pc_satisfy(a, not_newline, "comment"));
    parser_t *comment = pc_seq_r(a, hash, comment_body);
    return pc_many(a, pc_alt(a, ws_chr, comment));
}
static parser_t *lex(arena_t *a, parser_t *p, parser_t *ws) {
    return pc_seq_l(a, p, ws);
}

/* ============================================================
 * Semantic actions — build AST nodes
 *
 * Each map/combine function takes void* and returns void*.
 * Casts at the boundaries. Context structs carry the data we
 * need (the operator metadata, the binder kind, etc.).
 * ============================================================ */

/* ---- operator combiner: builds (op left right) ---- */
typedef struct { const op_info_t *op; } op_ctx_t;
static op_ctx_t *mk_op_ctx(arena_t *a, const op_info_t *op) {
    op_ctx_t *c = arena_alloc(a, sizeof *c); c->op = op; return c;
}
static void *combine_binop(void *l, void *o, void *r, arena_t *a, void *ctx) {
    (void)o;
    op_ctx_t *c = ctx;
    node_t *L = l, *R = r;
    span_t s = L->span; s.end = R->span.end;
    node_t *head = ast_const(a, c->op->name, c->op, s);
    return ast_app2(a, head, L, R, s);
}

static void *map_op_value(void *v, arena_t *a, void *ctx) {
    (void)v; (void)a;
    return ctx;
}

static void *combine_binop_from_value(void *l, void *o, void *r, arena_t *a, void *ctx) {
    (void)ctx;
    const op_info_t *op = o;
    node_t *L = l, *R = r;
    span_t s = L->span; s.end = R->span.end;
    node_t *head = ast_const(a, op->name, op, s);
    return ast_app2(a, head, L, R, s);
}

/* ---- prefix NOT: (\!  expr) → !expr ---- */
static void *map_wrap_not(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    pc_pair_t *p = v;  /* (bang_span, inner) where bang_span is pc_spanned_t */
    pc_spanned_t *bs = p->first;
    node_t *inner = p->second;
    span_t s = bs->span; s.end = inner->span.end;
    node_t *head = ast_const(a, OP_NOT.name, &OP_NOT, bs->span);
    return ast_app1(a, head, inner, s);
}

/* ---- prefix arithmetic negation: -expr → 0 - expr ---- */
static void *map_wrap_neg(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    pc_pair_t *p = v;
    pc_spanned_t *ms = p->first;
    node_t *inner = p->second;
    span_t s = ms->span; s.end = inner->span.end;
    node_t *zero = ast_const(a, "0", NULL, ms->span);
    node_t *head = ast_const(a, OP_SUB.name, &OP_SUB, ms->span);
    return ast_app2(a, head, zero, inner, s);
}

/* ---- var node from a spanned ident ---- */
static void *map_var(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    pc_spanned_t *sp = v;
    return ast_var(a, (const char *)sp->value, sp->span);
}
/* ---- const node (true / false / number) from a spanned value ---- */
static void *map_const_from_spanned(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    pc_spanned_t *sp = v;
    return ast_const(a, (const char *)sp->value, NULL, sp->span);
}

/* ---- application: atom [args]?  →  if args present, build app ---- */
static void *map_apply(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    pc_pair_t *p = v;
    node_t *head = p->first;
    pc_list_t *args = p->second;   /* may be NULL (no parens) */
    if (!args) return head;
    span_t s = head->span;
    if (args->count > 0) {
        node_t *last = args->items[args->count - 1];
        s.end = last->span.end;
    }
    node_t **arr = arena_alloc(a, sizeof(node_t *) * (args->count ? args->count : 1));
    for (size_t i = 0; i < args->count; i++) arr[i] = args->items[i];
    return ast_app(a, head, args->count, arr, s);
}

/* ---- types: parsed as right-assoc -> ---- */
typedef struct { void *type_parser; } type_box_t;  /* unused but keeps shape */

static void *combine_arrow_type(void *l, void *o, void *r, arena_t *a, void *ctx) {
    (void)o; (void)ctx;
    return type_arrow(a, (type_t *)l, (type_t *)r);
}
static void *map_type_atom_ident(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    return type_base(a, (const char *)v);
}

/* ---- binders ---- */
/* var decl is (name, type-or-NULL) */
typedef struct { const char *name; type_t *type; span_t span; } vardecl_t;

static void *map_vardecl(void *v, arena_t *a, void *ctx) {
    (void)ctx;
    pc_pair_t *p = v;
    pc_spanned_t *idsp = p->first;
    type_t *t = p->second;       /* NULL if no type annotation */
    vardecl_t *vd = arena_alloc(a, sizeof *vd);
    vd->name = (const char *)idsp->value;
    vd->type = t;
    vd->span = idsp->span;
    return vd;
}

typedef struct { node_kind_t kind; } binder_kind_ctx_t;
static binder_kind_ctx_t *mk_bk(arena_t *a, node_kind_t k) {
    binder_kind_ctx_t *c = arena_alloc(a, sizeof *c); c->kind = k; return c;
}

/* binder data: kind  +  list of vardecls  +  body. We assemble the
 * nested binder chain right-folded. */
static void *map_binder(void *v, arena_t *a, void *ctx) {
    binder_kind_ctx_t *bk = ctx;
    pc_pair_t *p = v;
    pc_pair_t *intro = p->first;             /* (kw_span, varlist) */
    node_t   *body  = p->second;
    pc_spanned_t *kws = intro->first;
    pc_list_t *vars = intro->second;
    node_t *acc = body;
    for (size_t i = vars->count; i > 0; i--) {
        vardecl_t *vd = vars->items[i - 1];
        span_t s = (i == 1) ? kws->span : vd->span;
        s.end = acc->span.end;
        acc = ast_binder(a, bk->kind, vd->name, vd->type, acc, s);
    }
    return acc;
}

/* ============================================================
 * Grammar
 *
 * top      ::= ws iff_expr eof
 * iff_expr ::= imp_expr ("<=>" imp_expr)?       non-assoc
 * imp_expr ::= eq_expr  ("=>"  imp_expr)?       right-assoc
 * eq_expr  ::= or_expr  ("="   or_expr)?        non-assoc
 * or_expr  ::= and_expr ("||"  and_expr)*       left-assoc
 * and_expr ::= add_expr ("&&"  add_expr)*       left-assoc
 * add_expr ::= mul_expr (("+"|"-") mul_expr)*  left-assoc
 * mul_expr ::= prefix   (("*"|"/") prefix)*    left-assoc
 * prefix   ::= "!" prefix | "-" prefix | app_expr
 * app_expr ::= atom ("(" expr ("," expr)* ")")? postfix
 * atom     ::= "(" expr ")" | binder | "true" | "false"
 *            | ident | number
 * binder   ::= ("forall"|"exists"|"lambda"|"fun"|"\\") vdecl+ "." expr
 * vdecl    ::= ident (":" type)?
 * type     ::= type_atom ("->" type)?
 * type_atom::= ident | "(" type ")"
 * ============================================================ */

combo_result_t combo_parse(const char *src, arena_t *a) {
    parser_t *ws = make_ws(a);

    /* tokens / lexemes */
    parser_t *t_ident   = lex(a, p_ident(a),    ws);  /* → char* */
    parser_t *t_num     = lex(a, p_natural(a),  ws);  /* → char* */
    parser_t *t_true    = lex(a, p_keyword(a, "true"),   ws);
    parser_t *t_false   = lex(a, p_keyword(a, "false"),  ws);
    parser_t *t_forall  = lex(a, p_keyword(a, "forall"), ws);
    parser_t *t_exists  = lex(a, p_keyword(a, "exists"), ws);
    parser_t *t_lambda  = lex(a, p_keyword(a, "lambda"), ws);
    parser_t *t_fun     = lex(a, p_keyword(a, "fun"),    ws);
    parser_t *t_bslash  = lex(a, pc_char(a, '\\'),       ws);

    parser_t *s_lp      = lex(a, pc_string(a, "("),  ws);
    parser_t *s_rp      = lex(a, pc_string(a, ")"),  ws);
    parser_t *s_comma   = lex(a, pc_string(a, ","),  ws);
    parser_t *s_dot     = lex(a, pc_string(a, "."),  ws);
    parser_t *s_colon   = lex(a, pc_string(a, ":"),  ws);
    parser_t *s_arrow   = lex(a, pc_string(a, "->"), ws);

    /* operator lexemes.
     * "=" must not match when followed by ">"  ("=>" lives at a higher level).
     * The others have no prefix conflicts. */
    parser_t *s_not     = lex(a, p_punct(a, "!",   NULL), ws);
    parser_t *s_mul     = lex(a, p_punct(a, "*",   NULL), ws);
    parser_t *s_div     = lex(a, p_punct(a, "/",   NULL), ws);
    parser_t *s_plus    = lex(a, p_punct(a, "+",   NULL), ws);
    parser_t *s_minus   = lex(a, p_punct(a, "-",   ">" ), ws);
    parser_t *s_and     = lex(a, p_punct(a, "&&",  NULL), ws);
    parser_t *s_or      = lex(a, p_punct(a, "||",  NULL), ws);
    parser_t *s_iff     = lex(a, p_punct(a, "<=>", NULL), ws);
    parser_t *s_imp     = lex(a, p_punct(a, "=>",  NULL), ws);
    parser_t *s_eq      = lex(a, p_punct(a, "=",   ">" ), ws);

    /* expression refs for mutual recursion */
    parser_t *expr_ref       = pc_ref(a);
    parser_t *prefix_expr_ref = pc_ref(a);

    /* === type subgrammar === */
    parser_t *type_ref = pc_ref(a);
    parser_t *type_atom = pc_alt(a,
        pc_map(a, t_ident, map_type_atom_ident, NULL),
        pc_between(a, s_lp, type_ref, s_rp));
    parser_t *type_full = pc_chainr1(a, type_atom, s_arrow, combine_arrow_type, NULL);
    pc_set(type_ref, type_full);

    /* === binders === */
    /* vdecl = ident (":" type)? */
    parser_t *vdecl = pc_map(a,
        pc_pair(a,
                pc_with_span(a, t_ident),
                pc_opt_or(a, pc_seq_r(a, s_colon, type_ref), NULL)),
        map_vardecl, NULL);
    parser_t *vdecls = pc_many1(a, vdecl);

    /* binder intro: (kw_with_span,  vdecls)  followed by  "."  body */
    /* Reused for each kind. Helper to build one binder parser: */
    #define BINDER_OF(KW, KIND_CTX)                                              \
        pc_map(a,                                                                \
            pc_pair(a,                                                           \
                pc_pair(a, pc_with_span(a, (KW)), vdecls),                       \
                pc_seq_r(a, s_dot, expr_ref)),                                   \
            map_binder, (KIND_CTX))

    binder_kind_ctx_t *ctx_forall = mk_bk(a, NODE_FORALL);
    binder_kind_ctx_t *ctx_exists = mk_bk(a, NODE_EXISTS);
    binder_kind_ctx_t *ctx_lambda = mk_bk(a, NODE_LAMBDA);

    parser_t *binder_alts[] = {
        BINDER_OF(t_forall, ctx_forall),
        BINDER_OF(t_exists, ctx_exists),
        BINDER_OF(t_lambda, ctx_lambda),
        BINDER_OF(t_fun,    ctx_lambda),
        BINDER_OF(t_bslash, ctx_lambda),
    };
    parser_t *binder = pc_alts(a, binder_alts, 5);

    /* === atom === */
    parser_t *paren_expr   = pc_between(a, s_lp, expr_ref, s_rp);
    parser_t *atom_true    = pc_map(a, pc_with_span(a, t_true),  map_const_from_spanned, NULL);
    parser_t *atom_false   = pc_map(a, pc_with_span(a, t_false), map_const_from_spanned, NULL);
    parser_t *atom_num     = pc_map(a, pc_with_span(a, t_num),   map_const_from_spanned, NULL);
    parser_t *atom_var     = pc_map(a, pc_with_span(a, t_ident), map_var, NULL);

    parser_t *atom_alts[] = {
        atom_true, atom_false, binder, atom_var, atom_num, paren_expr
    };
    parser_t *atom = pc_alts(a, atom_alts, 6);

    /* === application: atom (LPAREN expr ("," expr)* RPAREN)? === */
    parser_t *arglist = pc_sepby1(a, expr_ref, s_comma);
    parser_t *call_tail = pc_between(a, s_lp, arglist, s_rp);  /* → pc_list_t* */
    parser_t *app_expr = pc_map(a,
        pc_pair(a, atom, pc_opt_or(a, call_tail, NULL)),
        map_apply, NULL);

    /* === prefix operators === */
    parser_t *not_then_inner = pc_pair(a, pc_with_span(a, s_not), prefix_expr_ref);
    parser_t *neg_then_inner = pc_pair(a, pc_with_span(a, s_minus), prefix_expr_ref);
    parser_t *prefix_alts[] = {
        pc_map(a, not_then_inner, map_wrap_not, NULL),
        pc_map(a, neg_then_inner, map_wrap_neg, NULL),
        app_expr,
    };
    parser_t *prefix_expr = pc_alts(a, prefix_alts, 3);
    pc_set(prefix_expr_ref, prefix_expr);

    /* === precedence stack ===
     *   "direction" of each level is encoded by the chain combinator */
    parser_t *mul_op = pc_alt(a,
        pc_map(a, s_mul, map_op_value, (void *)&OP_MUL),
        pc_map(a, s_div, map_op_value, (void *)&OP_DIV));
    parser_t *add_op = pc_alt(a,
        pc_map(a, s_plus,  map_op_value, (void *)&OP_ADD),
        pc_map(a, s_minus, map_op_value, (void *)&OP_SUB));

    parser_t *mul_expr = pc_chainl1   (a, prefix_expr, mul_op, combine_binop_from_value, NULL);
    parser_t *add_expr = pc_chainl1   (a, mul_expr,   add_op, combine_binop_from_value, NULL);
    parser_t *and_expr = pc_chainl1   (a, add_expr,   s_and,  combine_binop, mk_op_ctx(a, &OP_AND));
    parser_t *or_expr  = pc_chainl1   (a, and_expr,   s_or,   combine_binop, mk_op_ctx(a, &OP_OR));
    parser_t *eq_expr  = pc_chain_none(a, or_expr,    s_eq,   combine_binop, mk_op_ctx(a, &OP_EQ));
    parser_t *imp_expr = pc_chainr1   (a, eq_expr,    s_imp,  combine_binop, mk_op_ctx(a, &OP_IMP));
    parser_t *iff_expr = pc_chain_none(a, imp_expr,   s_iff,  combine_binop, mk_op_ctx(a, &OP_IFF));

    pc_set(expr_ref, iff_expr);

    /* === top-level: leading ws, expr, EOF === */
    parser_t *top = pc_seq_l(a,
                              pc_seq_r(a, ws, expr_ref),
                              pc_eof(a));

    pc_result_t r = pc_run(top, src, a);
    combo_result_t out;
    memset(&out, 0, sizeof out);
    if (!r.ok) {
        out.err_pos     = r.err_pos;
        out.err_line    = r.err_line;
        out.err_col     = r.err_col;
        out.err_msg     = r.err_expected;
        return out;
    }
    out.ok  = true;
    out.ast = (node_t *)r.value;
    return out;

    #undef BINDER_OF
}
