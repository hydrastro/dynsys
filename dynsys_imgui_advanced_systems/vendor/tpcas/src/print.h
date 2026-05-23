#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include "ast.h"

#ifdef __cplusplus
extern "C" {
#endif

/* Pretty-print an expression. Round-trips through the parser. */
void print_expr(FILE *out, const node_t *n);
void print_expr_indent(FILE *out, const node_t *n);  /* multi-line tree view */
void print_type(FILE *out, const type_t *t);


#ifdef __cplusplus
}
#endif
#endif
