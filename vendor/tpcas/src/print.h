#ifndef PRINT_H
#define PRINT_H

#include <stdio.h>
#include "ast.h"

/* Pretty-print an expression. Round-trips through the parser. */
void print_expr(FILE *out, const node_t *n);
void print_expr_indent(FILE *out, const node_t *n);  /* multi-line tree view */
void print_type(FILE *out, const type_t *t);

#endif
