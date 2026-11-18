/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Held-out perplexity benchmark: masked Markov softmax vs the SAGE-1 SLERP
 * step vs variant A (log-ratio step, τ sweep) vs variant B (low-rank
 * residual, rank sweep), all under the same grammar mask.
 *
 * Corpus: sentences sampled from the stochastic grammar in
 * sage1_bench_grammar.c (seeded splitmix64), split into train / held-out.
 * Every model sees only bigram statistics of the training split:
 *   counts  C[k][j] + 1 for all 7 × 7 pairs (add-one), as the CSR matrix P.
 * At each held-out position the state is the one-hot e_k of the previous
 * token (sage1_condition, as in the Prolog decoder), the model produces
 * amplitudes x, and p_j = x_j² over the legal set, renormalized — exactly
 * what sage1_select scores. The softmax baseline is
 *   p_j ∝ exp( log P[k][j] / T )  over the legal set.
 *
 * Variant B's target rows T come from training data with the mask: for each
 * previous token k, the Luce choice model w_k that maximizes the masked
 * likelihood  Σ_events log( w_k[j] / Σ_{l legal} w_k[l] )  plus one
 * add-one pseudo event per (k, j) over the full vocabulary, fitted by the
 * MM iteration  w_j ← (C[k][j] + 1) / Σ_e [j ∈ L_e] / W_e(L_e)  (2000
 * iterations, from P[k]). T[k] = sqrt(w_k / Σ w_k). U, V by
 * sage1_lowrank_fit(S = sqrt(P), T, rank).
 *
 * Per method: held-out perplexity exp(−mean log p(next)), top-1 accuracy,
 * mean margin (top-1 minus top-2 probability, positions with ≥ 2 legal
 * tokens), and mask violations: positions where sage1_select picks an
 * illegal token or the held-out token is illegal under the mask (must be 0).
 *
 *   sage1_bench [--sentences N] [--seed S] [--max-depth D] [--out FILE]
 *   sage1_bench --dump-masks D        all parser states, Prolog terms
 *   sage1_bench --dump-corpus N S D   N sampled sentences, Prolog terms
 */
#include "sage1_bench_grammar.h"
#include "sage1_variants.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define V G_TOKENS
#define CAP 4096

typedef struct { uint32_t prev, next; uint8_t mask[V]; double oracle; } event;

typedef struct { event *e; size_t n, cap; } events;

static void push_event(events *ev, event x) {
    if (ev->n == ev->cap) {
        ev->cap = ev->cap ? 2 * ev->cap : 1024;
        ev->e = realloc(ev->e, ev->cap * sizeof *ev->e);
        if (!ev->e) { perror("realloc"); exit(2); }
    }
    ev->e[ev->n++] = x;
}

static void collect(uint64_t *rng, size_t sentences, size_t max_depth, events *ev,
                    size_t *tokens) {
    uint32_t toks[CAP];
    for (size_t s = 0; s < sentences; s++) {
        size_t len = 0;
        while (!(len = g_sample_sentence(rng, max_depth, toks, CAP))) {}
        g_gen gen;
        g_gen_init(&gen, max_depth);
        for (size_t i = 1; i < len; i++) {
            event x;
            x.prev = toks[i - 1];
            x.next = toks[i];
            g_legal_mask(&gen.parser, x.mask);
            double dist[V];
            g_gen_dist(&gen, dist);
            x.oracle = dist[x.next];
            push_event(ev, x);
            g_gen_advance(&gen, toks[i]);
        }
        *tokens += len - 1;
    }
}

typedef struct {
    double nll, margin, acc;
    size_t n, n_margin, violations;
} score;

/* Score one position from amplitudes x (unit, nonnegative). */
static void score_amp(score *sc, const double *x, const event *e) {
    size_t sel;
    double prob;
    sage1_status st = sage1_select(x, e->mask, V, &sel, &prob);
    if (st != SAGE1_OK || !e->mask[sel]) sc->violations++;
    if (!e->mask[e->next]) sc->violations++;
    double tot = 0.0, p[V];
    size_t legal = 0;
    for (size_t j = 0; j < V; j++) {
        p[j] = e->mask[j] ? x[j] * x[j] : 0.0;
        tot += p[j];
        legal += e->mask[j];
    }
    double pn = p[e->next] / tot;
    sc->nll += -log(pn);
    double t1 = -1.0, t2 = -1.0;
    for (size_t j = 0; j < V; j++) if (e->mask[j]) {
        double v = p[j] / tot;
        if (v > t1) { t2 = t1; t1 = v; } else if (v > t2) t2 = v;
    }
    if (st == SAGE1_OK && fabs(prob - t1) > 1e-9) sc->violations++;
    if (legal >= 2) { sc->margin += t1 - t2; sc->n_margin++; }
    if (st == SAGE1_OK && sel == e->next) sc->acc += 1.0;
    sc->n++;
}

static double ppl(const score *s) { return exp(s->nll / (double)s->n); }

static FILE *out;
static void row(const char *method, const char *param, const score *s) {
    fprintf(out, "| %s | %s | %.4f | %.4f | %.2f%% | %zu |\n", method, param, ppl(s),
            s->margin / (double)s->n_margin, 100.0 * s->acc / (double)s->n, s->violations);
}

static int dump_masks(size_t max_depth) {
    /* every stack over {(, [} of depth ≤ max_depth, every last token */
    printf("max_depth(%zu).\n", max_depth);
    for (size_t d = 0; d <= max_depth; d++)
        for (uint32_t bits = 0; bits < (1u << d); bits++)
            for (uint32_t last = 0; last < V; last++) {
                g_parser p;
                g_parser_init(&p, max_depth);
                for (size_t i = 0; i < d; i++) p.stack[i] = (bits >> i) & 1 ? 3 : 1;
                p.depth = d;
                p.last = last;
                uint8_t m[V];
                g_legal_mask(&p, m);
                printf("m([");
                for (size_t i = 0; i < d; i++)   /* Prolog stack: top first */
                    printf("%s%u", i ? "," : "", p.stack[d - 1 - i]);
                printf("],%u,[", last);
                for (size_t j = 0; j < V; j++) printf("%s%u", j ? "," : "", m[j]);
                printf("]).\n");
            }
    return 0;
}

static int dump_corpus(size_t n, uint64_t seed, size_t max_depth) {
    uint32_t toks[CAP];
    printf("max_depth(%zu).\n", max_depth);
    for (size_t s = 0; s < n; s++) {
        size_t len = 0;
        while (!(len = g_sample_sentence(&seed, max_depth, toks, CAP))) {}
        printf("s([");
        for (size_t i = 0; i < len; i++) printf("%s%u", i ? "," : "", toks[i]);
        printf("]).\n");
    }
    return 0;
}

int main(int argc, char **argv) {
    size_t sentences = 20000, max_depth = 3;
    uint64_t seed = 20261118u;
    const char *out_path = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--dump-masks") && i + 1 < argc)
            return dump_masks(strtoul(argv[i + 1], NULL, 10));
        if (!strcmp(argv[i], "--dump-corpus") && i + 3 < argc)
            return dump_corpus(strtoul(argv[i + 1], NULL, 10),
                               strtoull(argv[i + 2], NULL, 10),
                               strtoul(argv[i + 3], NULL, 10));
        if (!strcmp(argv[i], "--sentences") && i + 1 < argc) sentences = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--max-depth") && i + 1 < argc) max_depth = strtoul(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--out") && i + 1 < argc) out_path = argv[++i];
        else { fprintf(stderr, "usage: see header of sage1_bench.c\n"); return 2; }
    }
    if (max_depth == 0 || max_depth > G_MAX_DEPTH_LIMIT || sentences < 5) return 2;
    out = stdout;
    if (out_path && !(out = fopen(out_path, "w"))) { perror(out_path); return 2; }

    size_t n_train = sentences * 4 / 5, n_test = sentences - n_train;
    size_t train_tok = 0, test_tok = 0;
    uint64_t rng = seed;
    events train = {0}, test = {0};
    collect(&rng, n_train, max_depth, &train, &train_tok);
    collect(&rng, n_test, max_depth, &test, &test_tok);

    /* add-one bigram counts → dense CSR */
    uint64_t C[V][V];
    for (size_t k = 0; k < V; k++) for (size_t j = 0; j < V; j++) C[k][j] = 1;
    for (size_t i = 0; i < train.n; i++) C[train.e[i].prev][train.e[i].next]++;
    size_t row_ptr[V + 1];
    uint32_t col_idx[V * V];
    uint64_t counts[V * V];
    double P[V][V], S[V * V], T[V * V];
    for (size_t k = 0; k < V; k++) {
        row_ptr[k] = k * V;
        uint64_t tot = 0;
        for (size_t j = 0; j < V; j++) tot += C[k][j];
        for (size_t j = 0; j < V; j++) {
            col_idx[k * V + j] = (uint32_t)j;
            counts[k * V + j] = C[k][j];
            P[k][j] = (double)C[k][j] / (double)tot;
            S[k * V + j] = sqrt(P[k][j]);
        }
    }
    row_ptr[V] = V * V;
    sage1_csr m = {V, V * V, row_ptr, col_idx, counts};
    if (sage1_validate(&m) != SAGE1_OK) { fprintf(stderr, "bad CSR\n"); return 1; }

    /* masked Luce MLE per previous token (MM iteration) → T rows.
     * Events are grouped by (previous token, mask): group[k][bits] counts. */
    static double group[V][1u << V];
    for (size_t i = 0; i < train.n; i++) {
        unsigned bits = 0;
        for (size_t j = 0; j < V; j++) bits |= (unsigned)train.e[i].mask[j] << j;
        group[train.e[i].prev][bits] += 1.0;
    }
    for (size_t k = 0; k < V; k++) {
        double w[V];
        for (size_t j = 0; j < V; j++) w[j] = P[k][j];
        for (int it = 0; it < 2000; it++) {
            double den[V], wall = 0.0;
            for (size_t j = 0; j < V; j++) wall += w[j];
            for (size_t j = 0; j < V; j++) den[j] = (double)V / wall;
            for (unsigned bits = 1; bits < (1u << V); bits++) {
                if (group[k][bits] == 0.0) continue;
                double W = 0.0;
                for (size_t j = 0; j < V; j++) if (bits >> j & 1) W += w[j];
                for (size_t j = 0; j < V; j++) if (bits >> j & 1) den[j] += group[k][bits] / W;
            }
            double sum = 0.0;
            for (size_t j = 0; j < V; j++) { w[j] = (double)C[k][j] / den[j]; sum += w[j]; }
            for (size_t j = 0; j < V; j++) w[j] /= sum;
        }
        for (size_t j = 0; j < V; j++) T[k * V + j] = sqrt(w[j]);
    }

    fprintf(out, "Corpus: %zu sentences from the SAGE-1 grammar (max depth %zu, seed %llu); "
            "train %zu sentences / %zu tokens, held-out %zu sentences / %zu tokens.\n\n",
            sentences, max_depth, (unsigned long long)seed, n_train, train_tok, n_test, test_tok);
    fprintf(out, "| Method | Params | Held-out PPL | Margin | Top-1 | Mask violations |\n");
    fprintf(out, "|---|---|---:|---:|---:|---:|\n");

    double x[V], amp[V], q[V], y[V], z[V], w[V], U[V * V], Vm[V * V];
    size_t total_viol = 0;

    /* generator's own next-token distribution: the floor for any model */
    {
        score s = {0};
        for (size_t i = 0; i < test.n; i++) s.nll += -log(test.e[i].oracle);
        s.n = test.n;
        fprintf(out, "| Generator (true distribution) | — | %.4f | — | — | — |\n", ppl(&s));
    }

    static const double temps[] = {0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 2.5};
    static const double lambdas[] = {0.25, 0.3, 0.5, 0.9, 1.0};
    char param[64];

    for (size_t t = 0; t < sizeof temps / sizeof *temps; t++) {
        score s = {0};
        for (size_t i = 0; i < test.n; i++) {
            const event *e = &test.e[i];
            double mx = -INFINITY, tot = 0.0;
            for (size_t j = 0; j < V; j++) if (e->mask[j]) {
                double l = log(P[e->prev][j]) / temps[t];
                if (l > mx) mx = l;
            }
            for (size_t j = 0; j < V; j++) {
                amp[j] = e->mask[j] ? exp(log(P[e->prev][j]) / temps[t] - mx) : 0.0;
                tot += amp[j];
            }
            for (size_t j = 0; j < V; j++) amp[j] = sqrt(amp[j] / tot);
            score_amp(&s, amp, e);
        }
        snprintf(param, sizeof param, "T = %.2f", temps[t]);
        row("Masked Markov softmax", param, &s);
        total_viol += s.violations;
    }
    for (size_t l = 0; l < sizeof lambdas / sizeof *lambdas; l++) {
        score s = {0};
        for (size_t i = 0; i < test.n; i++) {
            const event *e = &test.e[i];
            sage1_condition(e->prev, V, x);
            if (sage1_step(&m, x, lambdas[l], amp, q, y) != SAGE1_OK) return 1;
            score_amp(&s, amp, e);
        }
        snprintf(param, sizeof param, "λ = %.2f", lambdas[l]);
        row("SAGE-1 SLERP step", param, &s);
        total_viol += s.violations;
    }
    for (size_t t = 0; t < sizeof temps / sizeof *temps; t++) {
        score s = {0};
        for (size_t i = 0; i < test.n; i++) {
            const event *e = &test.e[i];
            sage1_condition(e->prev, V, x);
            if (sage1_logratio_step(&m, x, temps[t], amp, q) != SAGE1_OK) return 1;
            score_amp(&s, amp, e);
        }
        snprintf(param, sizeof param, "τ = %.2f", temps[t]);
        row("A: log-ratio step", param, &s);
        total_viol += s.violations;
    }
    for (size_t r = 0; r <= V; r++) {
        if (sage1_lowrank_fit(S, T, V, r, U, Vm) != SAGE1_OK) {
            fprintf(stderr, "low-rank fit failed at rank %zu\n", r);
            return 1;
        }
        score s = {0};
        for (size_t i = 0; i < test.n; i++) {
            const event *e = &test.e[i];
            sage1_condition(e->prev, V, x);
            if (sage1_lowrank_step(&m, x, 1.0, U, Vm, r, amp, y, z, w) != SAGE1_OK) return 1;
            score_amp(&s, amp, e);
        }
        snprintf(param, sizeof param, "rank %zu, λ = 1", r);
        row("B: low-rank residual", param, &s);
        total_viol += s.violations;
    }
    if (out != stdout) fclose(out);
    free(train.e);
    free(test.e);
    if (total_viol) { fprintf(stderr, "mask violations: %zu\n", total_viol); return 1; }
    return 0;
}
