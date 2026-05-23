#ifndef PARSE_H
#define PARSE_H

#include "ast.h"
#include "arena.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    bool        ok;
    node_t     *ast;
    span_t      err_span;
    const char *err_msg;
} parse_result_t;

/* Parse the whole input as one expression. */
parse_result_t parse(const char *src, arena_t *arena);


#ifdef __cplusplus
}
#endif
#endif
