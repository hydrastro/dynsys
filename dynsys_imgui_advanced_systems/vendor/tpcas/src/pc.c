#include "pc.h"
#include "lib/str.h"
#include <string.h>
#include <ctype.h>

/* ============================================================
 * Helpers
 * ============================================================ */

void pc_advance(pstate_t *st) {
    if (st->pos >= st->len) return;
    if (st->input[st->pos] == '\n') { st->line++; st->col = 1; }
    else st->col++;
    st->pos++;
}
static void advance_one(pstate_t *st) { pc_advance(st); }

/* Record "expected X" at the current position.
 *
 * We keep ALL the things expected at the FURTHEST position reached.
 * When something is expected at a strictly greater position, the
 * old list is cleared. When at the same position, append (dedup
 * by pointer identity — the descs are static or arena strings).
 * pc_run formats the array into one human-readable message. */
void pc_note_err(pstate_t *st, const char *what) {
    if (what == NULL) return;
    if (st->pos > st->err_pos) {
        st->err_pos             = st->pos;
        st->err_line            = st->line;
        st->err_col             = st->col;
        st->err_expected_count  = 0;
    } else if (st->pos < st->err_pos) {
        return;
    }
    /* at err_pos: append with dedup */
    for (size_t i = 0; i < st->err_expected_count; i++)
        if (st->err_expected_arr[i] == what) return;
    if (st->err_expected_count >= st->err_expected_cap) {
        size_t ncap = st->err_expected_cap ? st->err_expected_cap * 2 : 8;
        const char **na = arena_alloc(st->arena, sizeof(const char *) * ncap);
        for (size_t i = 0; i < st->err_expected_count; i++)
            na[i] = st->err_expected_arr[i];
        st->err_expected_arr = na;
        st->err_expected_cap = ncap;
    }
    st->err_expected_arr[st->err_expected_count++] = what;
}

pc_snapshot_t pc_snapshot(const pstate_t *st) {
    pc_snapshot_t s;
    s.pos       = st->pos;     s.line      = st->line;    s.col       = st->col;
    s.err_pos   = st->err_pos; s.err_line  = st->err_line;s.err_col   = st->err_col;
    s.err_arr   = st->err_expected_arr;
    s.err_count = st->err_expected_count;
    s.err_cap   = st->err_expected_cap;
    return s;
}
void pc_restore_err(pstate_t *st, pc_snapshot_t snap) {
    st->err_pos             = snap.err_pos;
    st->err_line            = snap.err_line;
    st->err_col             = snap.err_col;
    st->err_expected_arr    = snap.err_arr;
    st->err_expected_count  = snap.err_count;
    st->err_expected_cap    = snap.err_cap;
}

/* internal alias kept for the rest of this file */
static void note_err(pstate_t *st, const char *what) { pc_note_err(st, what); }

static parser_t *mkpar(arena_t *a, parser_fn fn, void *state, const char *desc) {
    parser_t *p = arena_alloc(a, sizeof(parser_t));
    p->fn = fn;
    p->state = state;
    p->desc = desc;
    return p;
}

#define FAIL(c)      ((presult_t){.ok=false, .consumed=(c), .value=NULL})
#define OK(c, v)     ((presult_t){.ok=true,  .consumed=(c), .value=(v)})

/* ============================================================
 * Growable list
 * ============================================================ */
pc_list_t *pc_list_new(arena_t *a) {
    pc_list_t *l = arena_alloc(a, sizeof *l);
    return l;
}
void pc_list_push(pc_list_t *l, void *v, arena_t *a) {
    if (l->count >= l->cap) {
        size_t ncap = l->cap ? l->cap * 2 : 4;
        void **na = arena_alloc(a, sizeof(void *) * ncap);
        for (size_t i = 0; i < l->count; i++) na[i] = l->items[i];
        l->items = na;
        l->cap   = ncap;
    }
    l->items[l->count++] = v;
}

/* ============================================================
 * Primitives
 * ============================================================ */

typedef struct { bool (*pred)(int); const char *desc; } st_sat_t;
static presult_t fn_satisfy(const parser_t *self, pstate_t *st) {
    st_sat_t *s = self->state;
    if (st->pos >= st->len || !s->pred((unsigned char)st->input[st->pos])) {
        note_err(st, s->desc);
        return FAIL(false);
    }
    advance_one(st);
    return OK(true, NULL);
}
parser_t *pc_satisfy(arena_t *a, bool (*pred)(int), const char *desc) {
    st_sat_t *s = arena_alloc(a, sizeof *s);
    s->pred = pred; s->desc = desc;
    return mkpar(a, fn_satisfy, s, desc);
}

typedef struct { char c; } st_char_t;
static presult_t fn_char(const parser_t *self, pstate_t *st) {
    st_char_t *s = self->state;
    if (st->pos >= st->len || st->input[st->pos] != s->c) {
        note_err(st, self->desc);
        return FAIL(false);
    }
    advance_one(st);
    return OK(true, NULL);
}
parser_t *pc_char(arena_t *a, char c) {
    st_char_t *s = arena_alloc(a, sizeof *s);
    s->c = c;
    char *d = arena_alloc(a, 4);
    d[0] = '\''; d[1] = c; d[2] = '\''; d[3] = '\0';
    return mkpar(a, fn_char, s, d);
}

typedef struct { const char *s; size_t n; } st_string_t;
static presult_t fn_string(const parser_t *self, pstate_t *st) {
    st_string_t *s = self->state;
    if (st->pos + s->n > st->len || memcmp(st->input + st->pos, s->s, s->n) != 0) {
        note_err(st, self->desc);
        return FAIL(false);  /* atomic — never partial */
    }
    for (size_t i = 0; i < s->n; i++) advance_one(st);
    return OK(s->n > 0, (void *)s->s);
}
parser_t *pc_string(arena_t *a, const char *str) {
    st_string_t *s = arena_alloc(a, sizeof *s);
    s->s = str; s->n = strlen(str);
    size_t dlen = s->n + 3;
    char *d = arena_alloc(a, dlen);
    d[0] = '"'; memcpy(d + 1, str, s->n); d[1 + s->n] = '"'; d[2 + s->n] = '\0';
    return mkpar(a, fn_string, s, d);
}

typedef struct { bool (*pred)(int); const char *desc; } st_tw_t;
static presult_t fn_take_while1(const parser_t *self, pstate_t *st) {
    st_tw_t *s = self->state;
    size_t start = st->pos;
    while (st->pos < st->len && s->pred((unsigned char)st->input[st->pos]))
        advance_one(st);
    if (st->pos == start) {
        note_err(st, s->desc);
        return FAIL(false);
    }
    char *out = arena_strndup(st->arena, st->input + start, st->pos - start);
    return OK(true, out);
}
parser_t *pc_take_while1(arena_t *a, bool (*pred)(int), const char *desc) {
    st_tw_t *s = arena_alloc(a, sizeof *s);
    s->pred = pred; s->desc = desc;
    return mkpar(a, fn_take_while1, s, desc);
}

static presult_t fn_eof(const parser_t *self, pstate_t *st) {
    (void)self;
    if (st->pos >= st->len) return OK(false, NULL);
    note_err(st, "end of input");
    return FAIL(false);
}
parser_t *pc_eof(arena_t *a) { return mkpar(a, fn_eof, NULL, "end of input"); }

typedef struct { void *v; } st_pure_t;
static presult_t fn_pure(const parser_t *self, pstate_t *st) {
    (void)st;
    return OK(false, ((st_pure_t *)self->state)->v);
}
parser_t *pc_pure(arena_t *a, void *value) {
    st_pure_t *s = arena_alloc(a, sizeof *s);
    s->v = value;
    return mkpar(a, fn_pure, s, "pure");
}

static presult_t fn_fail(const parser_t *self, pstate_t *st) {
    note_err(st, self->desc);
    return FAIL(false);
}
parser_t *pc_fail(arena_t *a, const char *desc) {
    return mkpar(a, fn_fail, NULL, desc);
}

/* ============================================================
 * Combinators
 * ============================================================ */

typedef struct { parser_t *p; const char *desc; } st_label_t;
static presult_t fn_label(const parser_t *self, pstate_t *st) {
    st_label_t *s = self->state;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok && !r.consumed) note_err(st, s->desc);
    return r;
}
parser_t *pc_label(arena_t *a, parser_t *p, const char *desc) {
    st_label_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->desc = desc;
    return mkpar(a, fn_label, s, desc);
}

typedef struct { parser_t *p; } st_un_t;
static presult_t fn_try(const parser_t *self, pstate_t *st) {
    st_un_t *s = self->state;
    size_t pos = st->pos, line = st->line, col = st->col;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) {
        st->pos = pos; st->line = line; st->col = col;
        r.consumed = false;
    }
    return r;
}
parser_t *pc_try(arena_t *a, parser_t *p) {
    st_un_t *s = arena_alloc(a, sizeof *s);
    s->p = p;
    return mkpar(a, fn_try, s, p->desc);
}

typedef struct { parser_t *p1; parser_t *p2; } st_bin_t;
static presult_t fn_alt(const parser_t *self, pstate_t *st) {
    st_bin_t *s = self->state;
    presult_t r = s->p1->fn(s->p1, st);
    if (r.ok)        return r;
    if (r.consumed)  return r;   /* committed — don't try p2 */
    return s->p2->fn(s->p2, st);
}
parser_t *pc_alt(arena_t *a, parser_t *p1, parser_t *p2) {
    st_bin_t *s = arena_alloc(a, sizeof *s);
    s->p1 = p1; s->p2 = p2;
    return mkpar(a, fn_alt, s, p1->desc);
}
parser_t *pc_alts(arena_t *a, parser_t **arr, size_t n) {
    if (n == 0) return pc_fail(a, "no alternatives");
    parser_t *acc = arr[n - 1];
    for (size_t i = n - 1; i > 0; i--) acc = pc_alt(a, arr[i - 1], acc);
    return acc;
}

typedef enum { SEQ_LEFT, SEQ_RIGHT, SEQ_PAIR } seq_mode_t;
typedef struct { parser_t *p1; parser_t *p2; seq_mode_t mode; } st_seq_t;
static presult_t fn_seq(const parser_t *self, pstate_t *st) {
    st_seq_t *s = self->state;
    presult_t r1 = s->p1->fn(s->p1, st);
    if (!r1.ok) return r1;
    presult_t r2 = s->p2->fn(s->p2, st);
    if (!r2.ok)
        return (presult_t){.ok=false, .consumed = r1.consumed || r2.consumed, .value=NULL};
    void *v = NULL;
    switch (s->mode) {
        case SEQ_LEFT:  v = r1.value; break;
        case SEQ_RIGHT: v = r2.value; break;
        case SEQ_PAIR: {
            pc_pair_t *p = arena_alloc(st->arena, sizeof *p);
            p->first = r1.value; p->second = r2.value; v = p;
            break;
        }
    }
    return OK(r1.consumed || r2.consumed, v);
}
static parser_t *mk_seq(arena_t *a, parser_t *p1, parser_t *p2, seq_mode_t mode) {
    st_seq_t *s = arena_alloc(a, sizeof *s);
    s->p1 = p1; s->p2 = p2; s->mode = mode;
    return mkpar(a, fn_seq, s, p1->desc);
}
parser_t *pc_seq_l(arena_t *a, parser_t *p1, parser_t *p2) { return mk_seq(a, p1, p2, SEQ_LEFT); }
parser_t *pc_seq_r(arena_t *a, parser_t *p1, parser_t *p2) { return mk_seq(a, p1, p2, SEQ_RIGHT); }
parser_t *pc_pair (arena_t *a, parser_t *p1, parser_t *p2) { return mk_seq(a, p1, p2, SEQ_PAIR); }

typedef struct { parser_t *p; bool at_least_one; } st_many_t;
static presult_t fn_many(const parser_t *self, pstate_t *st) {
    st_many_t *s = self->state;
    pc_list_t *list = pc_list_new(st->arena);
    bool any = false;
    for (;;) {
        size_t before = st->pos;
        pc_snapshot_t snap = pc_snapshot(st);
        presult_t r = s->p->fn(s->p, st);
        if (!r.ok) {
            if (r.consumed) return FAIL(true);
            /* iteration's body errors are not the user's concern at this
             * point — many is happy with 0+. Restore. */
            pc_restore_err(st, snap);
            break;
        }
        if (st->pos == before) {
            pc_list_push(list, r.value, st->arena);
            break;
        }
        pc_list_push(list, r.value, st->arena);
        any = true;
    }
    if (s->at_least_one && list->count == 0) {
        note_err(st, self->desc);
        return FAIL(false);
    }
    return OK(any, list);
}
parser_t *pc_many(arena_t *a, parser_t *p) {
    st_many_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->at_least_one = false;
    return mkpar(a, fn_many, s, p->desc);
}
parser_t *pc_many1(arena_t *a, parser_t *p) {
    st_many_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->at_least_one = true;
    return mkpar(a, fn_many, s, p->desc);
}

typedef struct { parser_t *p; parser_t *sep; } st_sepby_t;
static presult_t fn_sepby1(const parser_t *self, pstate_t *st) {
    st_sepby_t *s = self->state;
    pc_list_t *list = pc_list_new(st->arena);
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) return r;
    pc_list_push(list, r.value, st->arena);
    bool consumed = r.consumed;
    for (;;) {
        presult_t rs = s->sep->fn(s->sep, st);
        if (!rs.ok) {
            if (rs.consumed) return FAIL(true);
            break;
        }
        presult_t rv = s->p->fn(s->p, st);
        if (!rv.ok) return FAIL(true); /* sep matched, body must follow */
        pc_list_push(list, rv.value, st->arena);
        consumed = true;
    }
    return OK(consumed, list);
}
parser_t *pc_sepby1(arena_t *a, parser_t *p, parser_t *sep) {
    st_sepby_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->sep = sep;
    return mkpar(a, fn_sepby1, s, p->desc);
}

parser_t *pc_between(arena_t *a, parser_t *open, parser_t *body, parser_t *close) {
    return pc_seq_l(a, pc_seq_r(a, open, body), close);
}

typedef struct { parser_t *p; void *default_v; } st_opt_t;
static presult_t fn_opt(const parser_t *self, pstate_t *st) {
    st_opt_t *s = self->state;
    pc_snapshot_t snap = pc_snapshot(st);
    presult_t r = s->p->fn(s->p, st);
    if (r.ok)        return r;
    if (r.consumed)  return r;
    /* p was optional and didn't fire — its expectations aren't real
     * errors at this position. Restore. */
    pc_restore_err(st, snap);
    return OK(false, s->default_v);
}
parser_t *pc_opt_or(arena_t *a, parser_t *p, void *default_value) {
    st_opt_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->default_v = default_value;
    return mkpar(a, fn_opt, s, p->desc);
}

static presult_t fn_with_span(const parser_t *self, pstate_t *st) {
    st_un_t *s = self->state;
    size_t spos = st->pos, sline = st->line, scol = st->col;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) return r;
    pc_spanned_t *sp = arena_alloc(st->arena, sizeof *sp);
    sp->value = r.value;
    sp->span.start = spos;
    sp->span.end   = st->pos;
    sp->span.line  = sline;
    sp->span.col   = scol;
    return OK(r.consumed, sp);
}
parser_t *pc_with_span(arena_t *a, parser_t *p) {
    st_un_t *s = arena_alloc(a, sizeof *s);
    s->p = p;
    return mkpar(a, fn_with_span, s, p->desc);
}

typedef struct { parser_t *p; pc_map_fn fn; void *ctx; } st_map_t;
static presult_t fn_map(const parser_t *self, pstate_t *st) {
    st_map_t *s = self->state;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) return r;
    return OK(r.consumed, s->fn(r.value, st->arena, s->ctx));
}
parser_t *pc_map(arena_t *a, parser_t *p, pc_map_fn fn, void *ctx) {
    st_map_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->fn = fn; s->ctx = ctx;
    return mkpar(a, fn_map, s, p->desc);
}

/* ============================================================
 * Expression chains — the heart of the precedence story.
 *
 *   chainl1  parses  p (op p)*   folding LEFT.
 *   chainr1  parses  p (op p)*   folding RIGHT.
 *   chain_none parses p (op p)?  exactly one application.
 *
 * Left/right choice is exactly your "direction" idea.
 * ============================================================ */
typedef struct { parser_t *p; parser_t *op; pc_combine_fn fn; void *ctx; } st_chain_t;

static presult_t fn_chainl1(const parser_t *self, pstate_t *st) {
    st_chain_t *s = self->state;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) return r;
    void *acc = r.value;
    bool any = r.consumed;
    for (;;) {
        presult_t ro = s->op->fn(s->op, st);
        if (!ro.ok) {
            if (ro.consumed) return FAIL(true);
            break;
        }
        presult_t rr = s->p->fn(s->p, st);
        if (!rr.ok) return FAIL(true);
        acc = s->fn(acc, ro.value, rr.value, st->arena, s->ctx);
        any = true;
    }
    return OK(any, acc);
}
parser_t *pc_chainl1(arena_t *a, parser_t *p, parser_t *op, pc_combine_fn fn, void *ctx) {
    st_chain_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->op = op; s->fn = fn; s->ctx = ctx;
    return mkpar(a, fn_chainl1, s, p->desc);
}

static presult_t fn_chainr1(const parser_t *self, pstate_t *st) {
    st_chain_t *s = self->state;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) return r;
    /* collect operands and operators into parallel lists, then fold right */
    pc_list_t *operands = pc_list_new(st->arena);
    pc_list_t *operators = pc_list_new(st->arena);
    pc_list_push(operands, r.value, st->arena);
    bool any = r.consumed;
    for (;;) {
        presult_t ro = s->op->fn(s->op, st);
        if (!ro.ok) {
            if (ro.consumed) return FAIL(true);
            break;
        }
        presult_t rr = s->p->fn(s->p, st);
        if (!rr.ok) return FAIL(true);
        pc_list_push(operators, ro.value, st->arena);
        pc_list_push(operands,  rr.value, st->arena);
        any = true;
    }
    void *acc = operands->items[operands->count - 1];
    for (size_t i = operands->count - 1; i > 0; i--)
        acc = s->fn(operands->items[i - 1], operators->items[i - 1],
                    acc, st->arena, s->ctx);
    return OK(any, acc);
}
parser_t *pc_chainr1(arena_t *a, parser_t *p, parser_t *op, pc_combine_fn fn, void *ctx) {
    st_chain_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->op = op; s->fn = fn; s->ctx = ctx;
    return mkpar(a, fn_chainr1, s, p->desc);
}

static presult_t fn_chain_none(const parser_t *self, pstate_t *st) {
    st_chain_t *s = self->state;
    presult_t r = s->p->fn(s->p, st);
    if (!r.ok) return r;
    presult_t ro = s->op->fn(s->op, st);
    if (!ro.ok) {
        if (ro.consumed) return FAIL(true);
        return r;
    }
    presult_t rr = s->p->fn(s->p, st);
    if (!rr.ok) return FAIL(true);
    void *combined = s->fn(r.value, ro.value, rr.value, st->arena, s->ctx);
    return OK(true, combined);
}
parser_t *pc_chain_none(arena_t *a, parser_t *p, parser_t *op, pc_combine_fn fn, void *ctx) {
    st_chain_t *s = arena_alloc(a, sizeof *s);
    s->p = p; s->op = op; s->fn = fn; s->ctx = ctx;
    return mkpar(a, fn_chain_none, s, p->desc);
}

/* ============================================================
 * Forward references
 * ============================================================ */
typedef struct { parser_t *target; } st_ref_t;
static presult_t fn_ref(const parser_t *self, pstate_t *st) {
    st_ref_t *s = self->state;
    if (!s->target) { note_err(st, "unresolved reference"); return FAIL(false); }
    return s->target->fn(s->target, st);
}
parser_t *pc_ref(arena_t *a) {
    st_ref_t *s = arena_alloc(a, sizeof *s);
    s->target = NULL;
    return mkpar(a, fn_ref, s, "<ref>");
}
void pc_set(parser_t *ref, parser_t *target) {
    ((st_ref_t *)ref->state)->target = target;
}

/* ============================================================
 * Run
 * ============================================================ */
static const char *format_expected(arena_t *a,
                                   const char **arr, size_t n) {
    if (n == 0) return "syntax error";
    if (n == 1) return arr[0];

    ds_str_t *msg = str_create();
    if (!msg) return "syntax error";

    for (size_t i = 0; i < n; i++) {
        if (i > 0) {
            if (n == 2) {
                str_pushc(msg, ' ');
            } else {
                str_append_cstr(msg, ", ");
            }
            if (i == n - 1) str_append_cstr(msg, "or ");
        }
        str_append_cstr(msg, arr[i]);
    }

    const char *out = arena_strndup(a, FUNC_str_cstr(msg), FUNC_str_len(msg));
    str_destroy(msg);
    return out ? out : "syntax error";
}

pc_result_t pc_run(parser_t *p, const char *input, arena_t *arena) {
    pstate_t st;
    memset(&st, 0, sizeof st);
    st.input = input;
    st.len   = strlen(input);
    st.line  = 1;
    st.col   = 1;
    st.arena = arena;

    presult_t r = p->fn(p, &st);

    pc_result_t out;
    memset(&out, 0, sizeof out);
    if (!r.ok) {
        out.err_pos      = st.err_pos;
        out.err_line     = st.err_line;
        out.err_col      = st.err_col;
        out.err_expected = format_expected(arena,
                                           st.err_expected_arr,
                                           st.err_expected_count);
    } else {
        out.ok    = true;
        out.value = r.value;
    }
    return out;
}
