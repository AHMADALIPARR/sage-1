/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * C port of the pushdown grammar in prolog/sage1_grammar.pl (legal/2,
 * advance/3) and a seeded stochastic generator for its language.
 * bench/check_grammar.pl checks both against the Prolog original.
 *
 * Tokens: 0 <s>, 1 (, 2 ), 3 [, 4 ], 5 a, 6 </s>.
 */
#ifndef UOR_SAGE1_BENCH_GRAMMAR_H
#define UOR_SAGE1_BENCH_GRAMMAR_H
#include <stddef.h>
#include <stdint.h>

#define G_TOKENS 7
#define G_MAX_DEPTH_LIMIT 16

/* Parser state p(Stack, Last, MaxDepth); stack[0] is the bottom. */
typedef struct {
    uint32_t stack[G_MAX_DEPTH_LIMIT];
    size_t depth, max_depth;
    uint32_t last;
} g_parser;

void g_parser_init(g_parser *p, size_t max_depth);
/* legal/2 as a 0/1 mask of length G_TOKENS. */
void g_legal_mask(const g_parser *p, uint8_t mask[G_TOKENS]);
/* advance/3. */
void g_advance(g_parser *p, uint32_t tok);

/* Generator state: the parser plus, per nesting level, whether the last
 * item was an atom and how many items the level holds. */
typedef struct {
    g_parser parser;
    uint8_t prev_atom[G_MAX_DEPTH_LIMIT + 1];
    uint32_t count[G_MAX_DEPTH_LIMIT + 1];
    int done;
} g_gen;

void g_gen_init(g_gen *g, size_t max_depth);
/* Exact next-token distribution of the generator (sums to 1). */
void g_gen_dist(const g_gen *g, double dist[G_TOKENS]);
void g_gen_advance(g_gen *g, uint32_t tok);

/* splitmix64 */
uint64_t g_rng_next(uint64_t *s);
/* Sample one sentence "<s> ... </s>" into toks (capacity cap).
 * Returns its length, or 0 if it would not fit. */
size_t g_sample_sentence(uint64_t *rng, size_t max_depth, uint32_t *toks, size_t cap);

#endif
