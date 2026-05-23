#include "lex.h"
#include <ctype.h>
#include <string.h>
#include <stdlib.h>

void lex_init(lexer_t *l, const char *src, arena_t *arena) {
    l->src   = src;
    l->len   = strlen(src);
    l->pos   = 0;
    l->line  = 1;
    l->col   = 1;
    l->arena = arena;
}

static void advance(lexer_t *l) {
    if (l->pos >= l->len) return;
    if (l->src[l->pos] == '\n') { l->line++; l->col = 1; }
    else l->col++;
    l->pos++;
}

static char peek_ch(const lexer_t *l, size_t k) {
    if (l->pos + k >= l->len) return '\0';
    return l->src[l->pos + k];
}

static void skip_ws(lexer_t *l) {
    while (l->pos < l->len) {
        char c = l->src[l->pos];
        if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
            advance(l);
        } else if (c == '#') {
            while (l->pos < l->len && l->src[l->pos] != '\n') advance(l);
        } else break;
    }
}

static span_t make_span(const lexer_t *l, size_t sp, size_t sl, size_t sc) {
    span_t s = {sp, l->pos, sl, sc};
    return s;
}

static token_t mk_tok(tok_kind_t k, span_t s, const char *lex) {
    token_t t;
    memset(&t, 0, sizeof t);
    t.kind = k;
    t.span = s;
    t.lexeme = lex;
    return t;
}

static token_t mk_err(span_t s, const char *msg) {
    token_t t;
    memset(&t, 0, sizeof t);
    t.kind = TOK_ERROR;
    t.span = s;
    t.err_msg = msg;
    return t;
}

bool kw_is_forall(const char *s) { return strcmp(s, "forall") == 0; }
bool kw_is_exists(const char *s) { return strcmp(s, "exists") == 0; }
bool kw_is_lambda(const char *s) { return strcmp(s, "lambda") == 0 || strcmp(s, "fun") == 0; }
bool kw_is_true  (const char *s) { return strcmp(s, "true")   == 0; }
bool kw_is_false (const char *s) { return strcmp(s, "false")  == 0; }
bool kw_is_any(const char *s) {
    return kw_is_forall(s) || kw_is_exists(s) || kw_is_lambda(s)
        || kw_is_true(s)   || kw_is_false(s);
}

/* Match the longest operator at p. Longest-match disambiguates <=> from =>. */
static const op_info_t *match_op(const char *p, size_t remaining, size_t *out_len) {
    size_t count;
    const op_info_t *const *ops = op_all(&count);
    const op_info_t *best = NULL;
    size_t best_len = 0;
    for (size_t i = 0; i < count; i++) {
        size_t l = strlen(ops[i]->syntax);
        if (l > remaining || l <= best_len) continue;
        if (strncmp(p, ops[i]->syntax, l) == 0) {
            best = ops[i];
            best_len = l;
        }
    }
    if (best) { *out_len = best_len; return best; }
    return NULL;
}

token_t lex_next(lexer_t *l) {
    skip_ws(l);
    size_t start = l->pos, sline = l->line, scol = l->col;

    if (l->pos >= l->len)
        return mk_tok(TOK_EOF, make_span(l, start, sline, scol), "");

    char c = l->src[l->pos];

    /* single-char punctuation */
    if (c == '(') { advance(l); return mk_tok(TOK_LPAREN,    make_span(l,start,sline,scol), "("); }
    if (c == ')') { advance(l); return mk_tok(TOK_RPAREN,    make_span(l,start,sline,scol), ")"); }
    if (c == ',') { advance(l); return mk_tok(TOK_COMMA,     make_span(l,start,sline,scol), ","); }
    if (c == '.') { advance(l); return mk_tok(TOK_DOT,       make_span(l,start,sline,scol), "."); }
    if (c == ':') { advance(l); return mk_tok(TOK_COLON,     make_span(l,start,sline,scol), ":"); }
    if (c == '\\'){ advance(l); return mk_tok(TOK_BACKSLASH, make_span(l,start,sline,scol), "\\"); }

    /* type arrow */
    if (c == '-' && peek_ch(l, 1) == '>') {
        advance(l); advance(l);
        return mk_tok(TOK_ARROW_T, make_span(l, start, sline, scol), "->");
    }

    /* operator */
    size_t oplen;
    const op_info_t *op = match_op(l->src + l->pos, l->len - l->pos, &oplen);
    if (op) {
        for (size_t i = 0; i < oplen; i++) advance(l);
        token_t t = mk_tok(TOK_OP, make_span(l, start, sline, scol), op->syntax);
        t.op = op;
        return t;
    }

    /* number: integer or decimal literal. Keep this intentionally small;
     * exponent notation can be layered in later without changing the AST. */
    if (isdigit((unsigned char)c)) {
        while (l->pos < l->len && isdigit((unsigned char)l->src[l->pos])) advance(l);
        if (l->pos < l->len && l->src[l->pos] == '.'
            && isdigit((unsigned char)peek_ch(l, 1))) {
            advance(l);
            while (l->pos < l->len && isdigit((unsigned char)l->src[l->pos])) advance(l);
        }
        size_t n = l->pos - start;
        char *lex = arena_strndup(l->arena, l->src + start, n);
        return mk_tok(TOK_NUMBER, make_span(l, start, sline, scol), lex);
    }

    /* identifier */
    if (isalpha((unsigned char)c) || c == '_') {
        while (l->pos < l->len &&
               (isalnum((unsigned char)l->src[l->pos]) || l->src[l->pos] == '_'))
            advance(l);
        size_t n = l->pos - start;
        char *lex = arena_strndup(l->arena, l->src + start, n);
        return mk_tok(TOK_IDENT, make_span(l, start, sline, scol), lex);
    }

    /* unknown */
    advance(l);
    return mk_err(make_span(l, start, sline, scol), "unexpected character");
}

token_t lex_peek(lexer_t *l) {
    size_t sp = l->pos, sl = l->line, sc = l->col;
    token_t t = lex_next(l);
    l->pos = sp; l->line = sl; l->col = sc;
    return t;
}
