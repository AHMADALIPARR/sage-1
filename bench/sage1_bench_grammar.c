/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Grammar port and corpus generator for the SAGE-1 perplexity benchmark.
 *
 * Generator (a probabilistic version of sentence//0 in sage1_grammar.pl,
 * bounded by MaxDepth like legal/2):
 *   A level at depth d either stops or adds one more item.
 *     stop probability   0.35 at d = 0 (never before the first item),
 *                        0.30 + 0.15·d inside a group;
 *                        1 when no item is allowed (atom after atom at max depth).
 *     stop token         </s> at d = 0, else the closer of the open bracket.
 *   An item is a, ( ... ) or [ ... ] with weights that depend on the
 *   enclosing bracket (none / ( / [):
 *                        a 0.40 0.60 0.30,  ( 0.35 0.25 0.30,  [ 0.25 0.15 0.40
 *   restricted to what is allowed (no a after a; groups only below MaxDepth)
 *   and renormalized.
 * The next token therefore depends on the stack, not only on the previous
 * token, which is the information a bigram model plus the mask cannot see.
 */
#include "sage1_bench_grammar.h"
#include <string.h>

enum { T_BOS = 0, T_LP = 1, T_RP = 2, T_LB = 3, T_RB = 4, T_A = 5, T_EOS = 6 };

static uint32_t closer_of(uint32_t open) { return open == T_LP ? T_RP : T_RB; }
static int is_opener(uint32_t t) { return t == T_LP || t == T_LB; }

void g_parser_init(g_parser *p, size_t max_depth) {
    memset(p, 0, sizeof *p);
    p->max_depth = max_depth > G_MAX_DEPTH_LIMIT ? G_MAX_DEPTH_LIMIT : max_depth;
    p->last = T_BOS;
}

void g_legal_mask(const g_parser *p, uint8_t mask[G_TOKENS]) {
    memset(mask, 0, G_TOKENS);
    if (p->depth < p->max_depth) { mask[T_LP] = 1; mask[T_LB] = 1; }
    if (p->depth > 0) mask[closer_of(p->stack[p->depth - 1])] = 1;
    if (p->last != T_A) mask[T_A] = 1;
    if (p->depth == 0 && p->last != T_BOS) mask[T_EOS] = 1;
}

void g_advance(g_parser *p, uint32_t tok) {
    if (is_opener(tok) && p->depth < G_MAX_DEPTH_LIMIT) {
        p->stack[p->depth++] = tok;
    } else if (p->depth > 0 && tok == closer_of(p->stack[p->depth - 1])) {
        p->depth--;
    }
    p->last = tok;
}

void g_gen_init(g_gen *g, size_t max_depth) {
    memset(g, 0, sizeof *g);
    g_parser_init(&g->parser, max_depth);
}

void g_gen_dist(const g_gen *g, double dist[G_TOKENS]) {
    memset(dist, 0, G_TOKENS * sizeof *dist);
    if (g->done) return;
    size_t d = g->parser.depth;
    int top = d == 0 ? 0 : (g->parser.stack[d - 1] == T_LP ? 1 : 2);
    static const double w_atom[3] = {0.40, 0.60, 0.30};
    static const double w_lp[3]   = {0.35, 0.25, 0.30};
    static const double w_lb[3]   = {0.25, 0.15, 0.40};
    double wa = g->prev_atom[d] ? 0.0 : w_atom[top];
    double wp = d < g->parser.max_depth ? w_lp[top] : 0.0;
    double wb = d < g->parser.max_depth ? w_lb[top] : 0.0;
    double wsum = wa + wp + wb;
    double stop;
    if (d == 0 && g->count[0] == 0) stop = 0.0;
    else if (wsum == 0.0) stop = 1.0;
    else stop = d == 0 ? 0.35 : 0.30 + 0.15 * (double)d;
    if (stop > 1.0) stop = 1.0;
    uint32_t stop_tok = d == 0 ? T_EOS : closer_of(g->parser.stack[d - 1]);
    dist[stop_tok] += stop;
    if (wsum > 0.0) {
        dist[T_A]  += (1.0 - stop) * wa / wsum;
        dist[T_LP] += (1.0 - stop) * wp / wsum;
        dist[T_LB] += (1.0 - stop) * wb / wsum;
    }
}

void g_gen_advance(g_gen *g, uint32_t tok) {
    size_t d = g->parser.depth;
    if (tok == T_EOS) {
        g->done = 1;
    } else if (tok == T_A) {
        g->prev_atom[d] = 1;
        g->count[d]++;
    } else if (is_opener(tok)) {
        g->prev_atom[d] = 0;
        g->count[d]++;
        g->prev_atom[d + 1] = 0;
        g->count[d + 1] = 0;
    }
    /* closers only pop the stack: the group already counted in its parent */
    g_advance(&g->parser, tok);
}

uint64_t g_rng_next(uint64_t *s) {
    uint64_t z = (*s += 0x9E3779B97F4A7C15u);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
    return z ^ (z >> 31);
}

size_t g_sample_sentence(uint64_t *rng, size_t max_depth, uint32_t *toks, size_t cap) {
    g_gen g;
    g_gen_init(&g, max_depth);
    size_t len = 0;
    if (cap == 0) return 0;
    toks[len++] = T_BOS;
    while (!g.done) {
        double dist[G_TOKENS], u = (double)(g_rng_next(rng) >> 11) * 0x1.0p-53, acc = 0.0;
        g_gen_dist(&g, dist);
        uint32_t pick = G_TOKENS;
        for (uint32_t t = 0; t < G_TOKENS; t++)
            if (dist[t] > 0.0) { acc += dist[t]; pick = t; if (u < acc) break; }
        if (len >= cap) return 0;
        toks[len++] = pick;
        g_gen_advance(&g, pick);
    }
    return len;
}
