#ifndef PC_H
#define PC_H

#include <stdbool.h>
#include <stddef.h>
#include "arena.h"
#include "ast.h"   /* for span_t */

#ifdef __cplusplus
extern "C" {
#endif

/* ============================================================
 * Parser combinators — focused, Parsec-style.
 *
 * A parser is a function-pointer + closure-state pair, both
 * arena-allocated. It runs against a mutable parser_state_t
 * and returns a parse_result_t.
 *
 * Failures live in the state, not the result. The state tracks
 * the FURTHEST position any parser reached and what was expected
 * there. On overall failure we report that — that's almost
 * always the most useful location for the user.
 *
 * Backtracking: pc_alt(p1, p2) tries p2 only if p1 fails WITHOUT
 * consuming input. pc_try(p) wraps p so its failures look
 * "unconsumed" — that's how you opt into backtracking.
 * ============================================================ */

typedef struct parser parser_t;

typedef struct {
    const char *input;
    size_t      len;
    size_t      pos;
    size_t      line, col;
    arena_t    *arena;
    /* furthest-error tracking — list of "expected X" descs at err_pos */
    size_t      err_pos;
    size_t      err_line, err_col;
    const char **err_expected_arr;   /* arena-owned */
    size_t      err_expected_count;
    size_t      err_expected_cap;
    /* compatibility: pc_run formats the array into one string here */
    const char *err_expected;
} pstate_t;

typedef struct {
    bool  ok;
    bool  consumed;
    void *value;
} presult_t;

typedef presult_t (*parser_fn)(const parser_t *self, pstate_t *st);

struct parser {
    parser_fn   fn;
    void       *state;
    const char *desc;
};

/* Public helpers — call these from custom primitives so error
 * accumulation works uniformly. */
void pc_note_err(pstate_t *st, const char *what);
void pc_advance (pstate_t *st);

/* For pc_many / pc_opt_or style suppression of body errors. */
typedef struct {
    size_t pos, line, col;
    size_t err_pos, err_line, err_col;
    const char **err_arr;
    size_t err_count, err_cap;
} pc_snapshot_t;
pc_snapshot_t pc_snapshot   (const pstate_t *st);
void          pc_restore_err(pstate_t *st, pc_snapshot_t snap);

/* helpers used in semantic actions */
typedef struct { void *first; void *second; } pc_pair_t;
typedef struct { void *value; span_t span; } pc_spanned_t;
typedef struct { void **items; size_t count; size_t cap; } pc_list_t;

pc_list_t *pc_list_new (arena_t *a);
void       pc_list_push(pc_list_t *l, void *v, arena_t *a);

/* ============================================================
 * Primitives
 * ============================================================ */
parser_t *pc_satisfy   (arena_t *a, bool (*pred)(int), const char *desc);
parser_t *pc_char      (arena_t *a, char c);
parser_t *pc_string    (arena_t *a, const char *s);    /* atomic — never partial */
parser_t *pc_take_while1(arena_t *a, bool (*pred)(int), const char *desc); /* → char* */
parser_t *pc_eof       (arena_t *a);
parser_t *pc_pure      (arena_t *a, void *value);
parser_t *pc_fail      (arena_t *a, const char *desc);

/* ============================================================
 * Combinators
 * ============================================================ */
parser_t *pc_label     (arena_t *a, parser_t *p, const char *desc);
parser_t *pc_try       (arena_t *a, parser_t *p);
parser_t *pc_alt       (arena_t *a, parser_t *p1, parser_t *p2);
parser_t *pc_alts      (arena_t *a, parser_t **arr, size_t n);

parser_t *pc_seq_l     (arena_t *a, parser_t *p1, parser_t *p2);  /* → p1's value */
parser_t *pc_seq_r     (arena_t *a, parser_t *p1, parser_t *p2);  /* → p2's value */
parser_t *pc_pair      (arena_t *a, parser_t *p1, parser_t *p2);  /* → pc_pair_t* */

parser_t *pc_many      (arena_t *a, parser_t *p);                 /* → pc_list_t* */
parser_t *pc_many1     (arena_t *a, parser_t *p);                 /* → pc_list_t* */
parser_t *pc_sepby1    (arena_t *a, parser_t *p, parser_t *sep);  /* → pc_list_t* */
parser_t *pc_between   (arena_t *a, parser_t *open, parser_t *body, parser_t *close);
parser_t *pc_opt_or    (arena_t *a, parser_t *p, void *default_value);
parser_t *pc_with_span (arena_t *a, parser_t *p);                 /* → pc_spanned_t* */

typedef void *(*pc_map_fn)(void *v, arena_t *a, void *ctx);
parser_t *pc_map       (arena_t *a, parser_t *p, pc_map_fn fn, void *ctx);

/* Expression-level combinators — encode "direction of parsing" directly. */
typedef void *(*pc_combine_fn)(void *left, void *op_v, void *right, arena_t *a, void *ctx);
parser_t *pc_chainl1   (arena_t *a, parser_t *p, parser_t *op, pc_combine_fn fn, void *ctx);
parser_t *pc_chainr1   (arena_t *a, parser_t *p, parser_t *op, pc_combine_fn fn, void *ctx);
parser_t *pc_chain_none(arena_t *a, parser_t *p, parser_t *op, pc_combine_fn fn, void *ctx);

/* Forward references for recursive grammars. */
parser_t *pc_ref (arena_t *a);
void      pc_set (parser_t *ref, parser_t *target);

/* ============================================================
 * Run
 * ============================================================ */
typedef struct {
    bool        ok;
    void       *value;
    size_t      err_pos, err_line, err_col;
    const char *err_expected;
} pc_result_t;

pc_result_t pc_run(parser_t *p, const char *input, arena_t *arena);


#ifdef __cplusplus
}
#endif
#endif
