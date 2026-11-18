/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Extended checks for the SAGE-1 C11 kernel: input validation, aliasing
 * rules, Markov equivalence at lambda = 1, unit-sphere invariants on
 * pseudo-random CSR matrices, and decoder tie-breaking.
 */
#include "sage1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { \
    fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
    failures++; } } while (0)

static uint64_t lcg_state = 0x5A6E1u;
static uint64_t lcg(void) {
    lcg_state = lcg_state * 6364136223846793005ULL + 1442695040888963407ULL;
    return lcg_state >> 33;
}

typedef struct {
    size_t n, e;
    size_t *row;
    uint32_t *col;
    uint64_t *cnt;
} owned_csr;

/* Random CSR: every row gets 1..max_deg distinct successors. */
static owned_csr random_csr(size_t n, size_t max_deg) {
    owned_csr m;
    m.n = n;
    m.row = calloc(n + 1, sizeof *m.row);
    m.col = calloc(n * max_deg, sizeof *m.col);
    m.cnt = calloc(n * max_deg, sizeof *m.cnt);
    size_t e = 0;
    for (size_t i = 0; i < n; i++) {
        m.row[i] = e;
        size_t deg = 1 + (size_t)(lcg() % max_deg);
        size_t start = (size_t)(lcg() % n);
        for (size_t d = 0; d < deg && d < n; d++) {
            m.col[e] = (uint32_t)((start + d * 7) % n);
            /* skip duplicates produced by the stride wrapping */
            int dup = 0;
            for (size_t k = m.row[i]; k < e; k++) if (m.col[k] == m.col[e]) dup = 1;
            if (dup) continue;
            m.cnt[e] = 1 + lcg() % 9;
            e++;
        }
    }
    m.row[n] = e;
    m.e = e;
    return m;
}

static sage1_csr view(const owned_csr *m) {
    sage1_csr v = {m->n, m->e, m->row, m->col, m->cnt};
    return v;
}

static void free_csr(owned_csr *m) { free(m->row); free(m->col); free(m->cnt); }

static double norm2(const double *v, size_t n) {
    double s = 0; for (size_t i = 0; i < n; i++) s += v[i] * v[i]; return s;
}

static void test_validation(void) {
    const size_t rows[] = {0, 2, 3, 4};
    const uint32_t cols[] = {1, 2, 2, 0};
    const uint64_t counts[] = {3, 1, 1, 1};
    sage1_csr m = {3, 4, rows, cols, counts};
    CHECK(sage1_validate(&m) == SAGE1_OK);
    CHECK(sage1_validate(NULL) == SAGE1_BAD_INPUT);

    const uint32_t bad_cols[] = {1, 3, 2, 0};
    sage1_csr bc = {3, 4, rows, bad_cols, counts};
    CHECK(sage1_validate(&bc) == SAGE1_BAD_MATRIX);

    const uint64_t zero_counts[] = {3, 0, 1, 1};
    sage1_csr zc = {3, 4, rows, cols, zero_counts};
    CHECK(sage1_validate(&zc) == SAGE1_BAD_MATRIX);

    const uint64_t big[] = {UINT64_MAX, 1, 1, 1};
    sage1_csr ov = {3, 4, rows, cols, big};
    CHECK(sage1_validate(&ov) == SAGE1_BAD_MATRIX);

    const size_t bad_end[] = {0, 2, 3, 3};
    sage1_csr be = {3, 4, bad_end, cols, counts};
    CHECK(sage1_validate(&be) == SAGE1_BAD_INPUT);

    const size_t bad_start[] = {1, 2, 3, 4};
    sage1_csr bs = {3, 4, bad_start, cols, counts};
    CHECK(sage1_validate(&bs) == SAGE1_BAD_INPUT);
}

static void test_inputs_and_aliasing(void) {
    const size_t rows[] = {0, 2, 3, 4};
    const uint32_t cols[] = {1, 2, 2, 0};
    const uint64_t counts[] = {3, 1, 1, 1};
    sage1_csr m = {3, 4, rows, cols, counts};
    double p[] = {1, 0, 0}, q[3], y[3], out[3];
    double x[] = {1, 0, 0};

    CHECK(sage1_propagate(&m, p, p) == SAGE1_BAD_INPUT);          /* aliasing */
    double half[] = {0.5, 0, 0};
    CHECK(sage1_propagate(&m, half, q) == SAGE1_BAD_INPUT);       /* mass != 1 */
    double neg[] = {1.5, -0.5, 0};
    CHECK(sage1_propagate(&m, neg, q) == SAGE1_BAD_INPUT);        /* negative */
    double nan_p[] = {NAN, 0, 0};
    CHECK(sage1_propagate(&m, nan_p, q) == SAGE1_BAD_INPUT);

    CHECK(sage1_step(&m, x, 1.0, x, q, y) == SAGE1_BAD_INPUT);    /* out == x */
    CHECK(sage1_step(&m, x, 1.0, out, q, q) == SAGE1_BAD_INPUT);  /* q == y */
    CHECK(sage1_step(&m, x, 1.5, out, q, y) == SAGE1_BAD_INPUT);  /* lambda > 1 */
    CHECK(sage1_step(&m, x, -0.1, out, q, y) == SAGE1_BAD_INPUT);
    CHECK(sage1_step(&m, x, NAN, out, q, y) == SAGE1_BAD_INPUT);
    double not_unit[] = {0.5, 0.5, 0};
    CHECK(sage1_step(&m, not_unit, 0.5, out, q, y) == SAGE1_BAD_INPUT);

    double a[] = {1, 0, 0}, b[] = {0, 1, 0};
    CHECK(sage1_geodesic(a, b, 3, 0.5, a) == SAGE1_BAD_INPUT);
    CHECK(sage1_geodesic(a, b, 3, 0.5, out) == SAGE1_OK);
    CHECK(fabs(out[0] - sqrt(0.5)) < 1e-14 && fabs(out[1] - sqrt(0.5)) < 1e-14);
    /* quarter of the way along a 90 degree arc */
    CHECK(sage1_geodesic(a, b, 3, 0.25, out) == SAGE1_OK);
    CHECK(fabs(out[0] - cos(acos(-1.0) / 8)) < 1e-14 && fabs(out[1] - sin(acos(-1.0) / 8)) < 1e-14);

    CHECK(sage1_condition(3, 3, out) == SAGE1_BAD_INPUT);
    CHECK(sage1_condition(0, 0, out) == SAGE1_BAD_INPUT);
}

static void test_select(void) {
    const double s = sqrt(0.25);
    double x[] = {s, s, s, s};            /* four-way tie */
    const uint8_t all[] = {1, 1, 1, 1};
    const uint8_t back[] = {0, 0, 1, 1};
    size_t sel = 99; double prob = 0;
    CHECK(sage1_select(x, all, 4, &sel, &prob) == SAGE1_OK);
    CHECK(sel == 0 && fabs(prob - 0.25) < 1e-15);
    CHECK(sage1_select(x, back, 4, &sel, &prob) == SAGE1_OK);
    CHECK(sel == 2 && fabs(prob - 0.5) < 1e-15);

    double y[] = {0, sqrt(0.9), sqrt(0.1), 0};
    const uint8_t mask_zero_mass[] = {1, 0, 0, 1};
    CHECK(sage1_select(y, mask_zero_mass, 4, &sel, &prob) == SAGE1_DEAD_END);
    const uint8_t mask_minor[] = {1, 0, 1, 1};
    CHECK(sage1_select(y, mask_minor, 4, &sel, &prob) == SAGE1_OK);
    CHECK(sel == 2 && fabs(prob - 1.0) < 1e-15);
}

/* lambda = 1 reproduces the Markov chain p_{t+1} = P^T p_t after squaring. */
static void test_markov_equivalence(void) {
    for (int trial = 0; trial < 40; trial++) {
        size_t n = 2 + (size_t)(lcg() % 60);
        owned_csr om = random_csr(n, 5);
        sage1_csr m = view(&om);
        CHECK(sage1_validate(&m) == SAGE1_OK);
        double *x = calloc(n, sizeof *x), *out = calloc(n, sizeof *out);
        double *q = calloc(n, sizeof *q), *y = calloc(n, sizeof *y);
        double *p = calloc(n, sizeof *p), *pn = calloc(n, sizeof *pn);
        x[lcg() % n] = 1.0;
        for (size_t i = 0; i < n; i++) p[i] = x[i] * x[i];
        for (int t = 0; t < 25; t++) {
            CHECK(sage1_step(&m, x, 1.0, out, q, y) == SAGE1_OK);
            CHECK(sage1_propagate(&m, p, pn) == SAGE1_OK);
            double err = 0;
            for (size_t i = 0; i < n; i++) err = fmax(err, fabs(out[i] * out[i] - pn[i]));
            CHECK(err < 1e-12);
            CHECK(fabs(norm2(out, n) - 1.0) < 1e-12);
            memcpy(x, out, n * sizeof *x);
            memcpy(p, pn, n * sizeof *p);
            /* re-normalize the Markov reference so propagate's mass check holds */
            double mass = 0; for (size_t i = 0; i < n; i++) mass += p[i];
            for (size_t i = 0; i < n; i++) p[i] /= mass;
        }
        free(x); free(out); free(q); free(y); free(p); free(pn);
        free_csr(&om);
    }
}

/* Intermediate lambda stays on the sphere and in the nonnegative orthant,
 * and moves monotonically along the arc from x towards the target. */
static void test_geodesic_invariants(void) {
    for (int trial = 0; trial < 40; trial++) {
        size_t n = 2 + (size_t)(lcg() % 40);
        owned_csr om = random_csr(n, 4);
        sage1_csr m = view(&om);
        double *x = calloc(n, sizeof *x), *out = calloc(n, sizeof *out);
        double *q = calloc(n, sizeof *q), *y = calloc(n, sizeof *y);
        double *tgt = calloc(n, sizeof *tgt);
        x[lcg() % n] = 1.0;
        CHECK(sage1_step(&m, x, 1.0, tgt, q, y) == SAGE1_OK);
        double prev = -1.0;
        for (int k = 0; k <= 10; k++) {
            double lambda = k / 10.0;
            CHECK(sage1_step(&m, x, lambda, out, q, y) == SAGE1_OK);
            CHECK(fabs(norm2(out, n) - 1.0) < 1e-12);
            double d = 0;
            for (size_t i = 0; i < n; i++) { CHECK(out[i] >= 0.0); d += out[i] * tgt[i]; }
            CHECK(d >= prev - 1e-12);
            prev = d;
        }
        free(x); free(out); free(q); free(y); free(tgt);
        free_csr(&om);
    }
}

int main(void) {
    test_validation();
    test_inputs_and_aliasing();
    test_select();
    test_markov_equivalence();
    test_geodesic_invariants();
    if (failures) { fprintf(stderr, "SAGE-1 extended: %d failure(s)\n", failures); return 1; }
    puts("SAGE-1 extended: all checks passed");
    return 0;
}
