/* ============================================================
 * tpcas test harness — runs the same inputs through both
 * parsers and verifies the ASTs are equal (alpha-equivalent
 * via ast_equal).
 *
 * Pass:  both parsers produce ASTs that ast_equal returns true.
 * Fail:  parse error, or mismatched ASTs.
 *
 * Usage:
 *   ./build/tpcas           run the built-in suite
 *   ./build/tpcas "EXPR"    parse one expression with both
 *                          parsers and dump both ASTs
 * ============================================================ */

#define _POSIX_C_SOURCE 200809L
#include "arena.h"
#include "ast.h"
#include "print.h"
#include "pratt.h"
#include "combo.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const char *TESTS[] = {
    /* propositional core */
    "A",
    "!A",
    "!!A",
    "A && B",
    "A || B",
    "A && B && C",
    "A || B || C",
    "A && B || C",
    "A || B && C",
    "!A && B",
    "!(A && B)",
    "A => B",
    "A => B => C",
    "(A => B) => C",
    "A <=> B",
    "A = B",
    "A && B => C || D",
    /* first-order */
    "P(x)",
    "P(x, y)",
    "f(g(x))",
    "f(g(x), h(y, z))",
    "forall x. P(x)",
    "exists x. P(x)",
    "forall x y. P(x, y)",
    "forall x y z. P(x, y, z)",
    "forall x. P(x) && Q(x)",
    "forall x. exists y. R(x, y)",
    /* higher-order with types */
    "forall x:Nat. P(x)",
    "forall P:Nat -> Bool. forall x:Nat. P(x) || !P(x)",
    "lambda x. f(x)",
    "lambda x:Nat. x",
    "(lambda x. x)(a)",
    "lambda f. lambda x. f(f(x))",
    "\\x. f(x)",
    "fun x y. P(x, y)",
    /* literals */
    "true",
    "false",
    "true && false",
    "forall x. x = x",
    /* arithmetic extension, useful for embedding TPCAS as an expression parser */
    "x + y",
    "x + y * z",
    "(x + y) * z",
    "x * (28 - z) - y",
    "x * y - (8 / 3) * z",
    "-x + 1.5",
};

static int run_one(const char *src, bool verbose) {
    arena_t arena;
    arena_init(&arena, 65536);

    parse_result_t  pr = parse(src, &arena);
    combo_result_t  cr = combo_parse(src, &arena);

    int rc = 0;

    if (!pr.ok) {
        fprintf(stderr, "  PRATT FAILED:    %s\n", pr.err_msg ? pr.err_msg : "?");
        rc = 1;
    }
    if (!cr.ok) {
        fprintf(stderr, "  COMBO FAILED at %zu:%zu: expected %s\n",
                cr.err_line, cr.err_col, cr.err_msg ? cr.err_msg : "?");
        rc = 1;
    }
    if (pr.ok && cr.ok) {
        if (!ast_equal(pr.ast, cr.ast)) {
            fprintf(stderr, "  AST MISMATCH:\n    pratt: ");
            print_expr(stderr, pr.ast); fprintf(stderr, "\n    combo: ");
            print_expr(stderr, cr.ast); fprintf(stderr, "\n");
            rc = 1;
        }
    }

    if (verbose) {
        if (pr.ok) { fprintf(stderr, "  pratt: "); print_expr(stderr, pr.ast); fprintf(stderr, "\n"); }
        if (cr.ok) { fprintf(stderr, "  combo: "); print_expr(stderr, cr.ast); fprintf(stderr, "\n"); }
    }

    arena_destroy(&arena);
    return rc;
}

int main(int argc, char **argv) {
    if (argc > 1) {
        printf("=== %s ===\n", argv[1]);
        int rc = run_one(argv[1], true);
        return rc;
    }
    size_t n = sizeof(TESTS) / sizeof(TESTS[0]);
    size_t pass = 0, fail = 0;
    for (size_t i = 0; i < n; i++) {
        printf("[%2zu/%2zu] %-50s ", i + 1, n, TESTS[i]);
        fflush(stdout);
        int rc = run_one(TESTS[i], false);
        if (rc == 0) { printf("OK\n");   pass++; }
        else         { printf("FAIL\n"); fail++; }
    }
    printf("\n%zu/%zu passed, %zu failed\n", pass, n, fail);
    return fail == 0 ? 0 : 1;
}
