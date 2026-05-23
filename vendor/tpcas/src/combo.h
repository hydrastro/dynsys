#ifndef COMBO_H
#define COMBO_H

#include "ast.h"
#include "arena.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    bool        ok;
    node_t     *ast;
    size_t      err_pos, err_line, err_col;
    const char *err_msg;
} combo_result_t;

/* Same input/output shape as the Pratt parser, so they can be
 * driven with the same test harness and their ASTs compared via
 * ast_equal. */
combo_result_t combo_parse(const char *src, arena_t *arena);


#ifdef __cplusplus
}
#endif
#endif
