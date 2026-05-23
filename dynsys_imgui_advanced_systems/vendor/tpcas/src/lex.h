#ifndef LEX_H
#define LEX_H

#include <stdbool.h>
#include <stddef.h>
#include "ast.h"
#include "arena.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    TOK_EOF,
    TOK_LPAREN,     /* (        */
    TOK_RPAREN,     /* )        */
    TOK_COMMA,      /* ,        */
    TOK_DOT,        /* .        */
    TOK_COLON,      /* :        */
    TOK_ARROW_T,    /* ->  (type arrow) */
    TOK_BACKSLASH,  /* \  (alternate lambda intro) */
    TOK_OP,         /* operator, see .op */
    TOK_IDENT,      /* identifier or keyword */
    TOK_NUMBER,
    TOK_ERROR
} tok_kind_t;

typedef struct {
    tok_kind_t   kind;
    span_t       span;
    const char  *lexeme;   /* arena-owned, NUL-terminated */
    const op_info_t *op;   /* only when kind == TOK_OP */
    const char  *err_msg;  /* only when kind == TOK_ERROR */
} token_t;

typedef struct {
    const char *src;
    size_t      len;
    size_t      pos;
    size_t      line, col;
    arena_t    *arena;
} lexer_t;

void    lex_init(lexer_t *l, const char *src, arena_t *arena);
token_t lex_next(lexer_t *l);
token_t lex_peek(lexer_t *l);

bool kw_is_forall(const char *s);
bool kw_is_exists(const char *s);
bool kw_is_lambda(const char *s);
bool kw_is_true  (const char *s);
bool kw_is_false (const char *s);
bool kw_is_any   (const char *s);


#ifdef __cplusplus
}
#endif
#endif
