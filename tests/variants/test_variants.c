/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Checks for variant A (sage1_logratio_step) and variant B
 * (sage1_lowrank_step / sage1_lowrank_fit): unit-sphere invariants,
 * reduction to sage1_step, the 1/τ log-ratio law, exact full-rank fit,
 * and the Jacobi eigen-solver.
 */
#include "sage1_variants.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { \
    fprintf(stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, #cond); \
    failures++; } } while (0)

static uint64_t rng = 0x5A6E1B0Bu;
static uint64_t next_u64(void) {
    uint64_t z = (rng += 0x9E3779B97F4A7C15u);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9u;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBu;
    return z ^ (z >> 31);
}
static double unif(void) { return (double)(next_u64() >> 11) * 0x1.0p-53; }

#define N 9
static size_t row_ptr[N + 1];
static uint32_t col_idx[N * N];
static uint64_t counts[N * N];

static sage1_csr random_csr(void) {
    size_t e = 0;
    row_ptr[0] = 0;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < N; j++)
            if (j == (i + 1) % N || next_u64() % 3 == 0) {
                col_idx[e] = (uint32_t)j;
                counts[e] = 1 + next_u64() % 20;
                e++;
            }
        row_ptr[i + 1] = e;
    }
    sage1_csr m = {N, e, row_ptr, col_idx, counts};
    return m;
}

static void random_unit(double *x) {
    double s = 0.0;
    for (size_t i = 0; i < N; i++) { x[i] = unif() + 1e-3; s += x[i] * x[i]; }
    for (size_t i = 0; i < N; i++) x[i] /= sqrt(s);
}

static int is_unit_nonneg(const double *v, size_t n) {
    double s = 0.0;
    for (size_t i = 0; i < n; i++) { if (!(v[i] >= 0.0)) return 0; s += v[i] * v[i]; }
    return fabs(s - 1.0) < 1e-12;
}

static void test_logratio(void) {
    double x[N], a[N], b[N], q[N], y[N], tau_set[] = {0.25, 0.5, 1.0, 2.0, 4.0};
    for (int trial = 0; trial < 200; trial++) {
        sage1_csr m = random_csr();
        random_unit(x);
        CHECK(sage1_logratio_step(&m, x, 1.0, a, q) == SAGE1_OK);
        CHECK(sage1_step(&m, x, 1.0, b, q, y) == SAGE1_OK);
        for (size_t i = 0; i < N; i++) CHECK(fabs(a[i] - b[i]) < 1e-12);
        for (size_t t = 0; t < sizeof tau_set / sizeof *tau_set; t++) {
            double tau = tau_set[t];
            CHECK(sage1_logratio_step(&m, x, tau, a, q) == SAGE1_OK);
            CHECK(is_unit_nonneg(a, N));
            /* log(a_i²/a_j²) = (1/τ) log(b_i²/b_j²) on the support */
            for (size_t i = 0; i < N; i++)
                for (size_t j = 0; j < N; j++)
                    if (b[i] > 1e-6 && b[j] > 1e-6) {
                        double lhs = log(a[i] * a[i] / (a[j] * a[j]));
                        double rhs = log(b[i] * b[i] / (b[j] * b[j])) / tau;
                        CHECK(fabs(lhs - rhs) < 1e-9);
                    }
            for (size_t i = 0; i < N; i++) CHECK((a[i] == 0.0) == (b[i] == 0.0));
        }
    }
    /* one-hot state: support of the output is the CSR row */
    sage1_csr m = random_csr();
    double e0[N] = {1.0};
    CHECK(sage1_logratio_step(&m, e0, 0.5, a, q) == SAGE1_OK);
    for (size_t j = 0; j < N; j++) {
        int in_row = 0;
        for (size_t k = row_ptr[0]; k < row_ptr[1]; k++) in_row |= col_idx[k] == j;
        CHECK((a[j] > 0.0) == in_row);
    }
    /* input validation */
    CHECK(sage1_logratio_step(&m, x, 0.0, a, q) == SAGE1_BAD_INPUT);
    CHECK(sage1_logratio_step(&m, x, -1.0, a, q) == SAGE1_BAD_INPUT);
    CHECK(sage1_logratio_step(&m, x, NAN, a, q) == SAGE1_BAD_INPUT);
    CHECK(sage1_logratio_step(&m, x, INFINITY, a, q) == SAGE1_BAD_INPUT);
    CHECK(sage1_logratio_step(&m, x, 1.0, x, q) == SAGE1_BAD_INPUT);
    CHECK(sage1_logratio_step(&m, x, 1.0, a, a) == SAGE1_BAD_INPUT);
    double bad[N] = {0.5};
    CHECK(sage1_logratio_step(&m, bad, 1.0, a, q) == SAGE1_BAD_INPUT);
    /* extreme temperature: no overflow, still a unit vector */
    random_unit(x);
    CHECK(sage1_logratio_step(&m, x, 1e-4, a, q) == SAGE1_OK);
    CHECK(is_unit_nonneg(a, N));
}

static void dense_rows(const sage1_csr *m, double *S) {
    memset(S, 0, N * N * sizeof *S);
    for (size_t i = 0; i < N; i++) {
        double tot = 0.0;
        for (size_t k = m->row_ptr[i]; k < m->row_ptr[i + 1]; k++) tot += (double)m->counts[k];
        for (size_t k = m->row_ptr[i]; k < m->row_ptr[i + 1]; k++)
            S[i * N + m->col_idx[k]] = sqrt((double)m->counts[k] / tot);
    }
}

static void test_lowrank(void) {
    double x[N], a[N], b[N], q[N], y[N], z[N], w[N];
    double S[N * N], T[N * N], U[N * N], V[N * N];
    for (int trial = 0; trial < 100; trial++) {
        sage1_csr m = random_csr();
        random_unit(x);
        /* rank 0 is sage1_step */
        CHECK(sage1_lowrank_step(&m, x, 0.6, NULL, NULL, 0, a, y, z, NULL) == SAGE1_OK);
        CHECK(sage1_step(&m, x, 0.6, b, q, y) == SAGE1_OK);
        for (size_t i = 0; i < N; i++) CHECK(fabs(a[i] - b[i]) < 1e-12);
        /* dense amplitude rows and random unit nonnegative targets */
        dense_rows(&m, S);
        for (size_t k = 0; k < N; k++) random_unit(T + k * N);
        sage1_status st = sage1_lowrank_fit(S, T, N, N, U, V);
        if (st == SAGE1_NUMERIC) continue;   /* singular S: legitimately refused */
        CHECK(st == SAGE1_OK);
        for (size_t k = 0; k < N; k++) {
            double ek[N] = {0};
            ek[k] = 1.0;
            CHECK(sage1_lowrank_step(&m, ek, 1.0, U, V, N, a, y, z, w) == SAGE1_OK);
            CHECK(is_unit_nonneg(a, N));
            for (size_t j = 0; j < N; j++) CHECK(fabs(a[j] - T[k * N + j]) < 1e-8);
        }
        /* every rank keeps the sphere invariant; error is monotone in rank */
        double prev = INFINITY;
        for (size_t r = 1; r <= N; r++) {
            CHECK(sage1_lowrank_fit(S, T, N, r, U, V) == SAGE1_OK);
            double err = 0.0;
            for (size_t k = 0; k < N; k++) {
                /* raw residual fit error ||S_k + U V^T S_k − T_k||² */
                for (size_t i = 0; i < N; i++) {
                    double s = S[k * N + i];
                    for (size_t rr = 0; rr < r; rr++) {
                        double v = 0.0;
                        for (size_t j = 0; j < N; j++) v += V[j * r + rr] * S[k * N + j];
                        s += U[i * r + rr] * v;
                    }
                    err += (s - T[k * N + i]) * (s - T[k * N + i]);
                }
                CHECK(sage1_lowrank_step(&m, x, 0.5, U, V, r, a, y, z, w) == SAGE1_OK);
                CHECK(is_unit_nonneg(a, N));
            }
            CHECK(err <= prev + 1e-10);
            prev = err;
        }
        CHECK(prev < 1e-14);
    }
    sage1_csr m = random_csr();
    random_unit(x);
    CHECK(sage1_lowrank_step(&m, x, 0.5, NULL, V, 1, a, y, z, w) == SAGE1_BAD_INPUT);
    CHECK(sage1_lowrank_step(&m, x, 0.5, U, V, N + 1, a, y, z, w) == SAGE1_BAD_INPUT);
    CHECK(sage1_lowrank_step(&m, x, 1.5, NULL, NULL, 0, a, y, z, w) == SAGE1_BAD_INPUT);
    CHECK(sage1_lowrank_step(&m, x, 0.5, NULL, NULL, 0, a, a, z, w) == SAGE1_BAD_INPUT);
    double sing[4] = {1, 1, 1, 1}, tt[4] = {1, 0, 0, 1}, u2[4], v2[4];
    CHECK(sage1_lowrank_fit(sing, tt, 2, 1, u2, v2) == SAGE1_NUMERIC);
    CHECK(sage1_lowrank_fit(sing, tt, 2, 3, u2, v2) == SAGE1_BAD_INPUT);
}

static void test_eig(void) {
    double A[N * N], B[N * N], ev[N], E[N * N];
    for (int trial = 0; trial < 50; trial++) {
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j <= i; j++) A[i * N + j] = A[j * N + i] = unif() - 0.5;
        memcpy(B, A, sizeof A);
        CHECK(sage1_sym_eig(B, N, ev, E) == SAGE1_OK);
        for (size_t i = 1; i < N; i++) CHECK(ev[i - 1] >= ev[i]);
        for (size_t i = 0; i < N; i++)
            for (size_t j = 0; j < N; j++) {
                double s = 0.0, o = 0.0;
                for (size_t k = 0; k < N; k++) {
                    s += E[i * N + k] * ev[k] * E[j * N + k];
                    o += E[k * N + i] * E[k * N + j];
                }
                CHECK(fabs(s - A[i * N + j]) < 1e-12);
                CHECK(fabs(o - (i == j ? 1.0 : 0.0)) < 1e-12);
            }
    }
}

int main(void) {
    test_logratio();
    test_lowrank();
    test_eig();
    if (failures) { fprintf(stderr, "%d check(s) failed\n", failures); return 1; }
    puts("variants: all checks passed");
    return 0;
}
