/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Cross-checks the Futhark backend (futhark/sage1.fut, compiled with
 * `futhark <backend> --library`) against the C11 kernel on the same inputs:
 * propagate, step, multi-step evolve, geodesic, select, condition, and error
 * statuses, over pseudo-random CSR matrices up to SAGE1_FUT_MAX_N states.
 *
 * Build: cc -std=c11 -DFUT_HEADER='"sage1_c.h"' test_futhark_vs_c.c sage1_c.c sage1.c -lm
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "sage1.h"
#include FUT_HEADER

#ifndef FUT_BACKEND
#define FUT_BACKEND "unknown"
#endif

static const double TOL = 1e-12;
static int failures = 0;
static double max_err_prop = 0, max_err_step = 0, max_err_evolve = 0, max_err_geo = 0, max_err_sel = 0;
static long n_cases = 0;

#define CHECK(cond, ...) do { if (!(cond)) { failures++; \
    fprintf(stderr, "%s:%d: ", __FILE__, __LINE__); fprintf(stderr, __VA_ARGS__); \
    fputc('\n', stderr); } } while (0)
#define FUT(call) do { int fut_rc_ = (call); if (fut_rc_ != 0) { \
    char *e_ = futhark_context_get_error(ctx); \
    fprintf(stderr, "futhark error %d at %s:%d: %s\n", fut_rc_, __FILE__, __LINE__, e_ ? e_ : ""); \
    free(e_); exit(2); } } while (0)

static struct futhark_context *ctx;

static uint64_t rng = 0xC0FFEEu;
static uint64_t next_u64(void) {
    rng ^= rng << 13; rng ^= rng >> 7; rng ^= rng << 17; return rng;
}
static double next_unit(void) { return (double)(next_u64() >> 11) / 9007199254740992.0; }

typedef struct {
    size_t n, e;
    size_t *row; uint32_t *col; uint64_t *cnt;
    int64_t *row64, *col64;
} csr_t;

static csr_t make_csr(size_t n, size_t max_deg) {
    csr_t m;
    m.n = n;
    m.row = malloc((n + 1) * sizeof *m.row);
    m.col = malloc(n * max_deg * sizeof *m.col);
    m.cnt = malloc(n * max_deg * sizeof *m.cnt);
    unsigned char *seen = calloc(n, 1);
    size_t e = 0;
    for (size_t i = 0; i < n; i++) {
        m.row[i] = e;
        size_t deg = 1 + next_u64() % max_deg;
        for (size_t d = 0; d < deg; d++) {
            uint32_t c = (uint32_t)(next_u64() % n);
            if (seen[c]) continue;
            seen[c] = 1;
            m.col[e] = c;
            m.cnt[e] = 1 + next_u64() % 1000;
            e++;
        }
        for (size_t k = m.row[i]; k < e; k++) seen[m.col[k]] = 0;
    }
    m.row[n] = e;
    m.e = e;
    free(seen);
    m.row64 = malloc((n + 1) * sizeof *m.row64);
    m.col64 = malloc(e * sizeof *m.col64);
    for (size_t i = 0; i <= n; i++) m.row64[i] = (int64_t)m.row[i];
    for (size_t k = 0; k < e; k++) m.col64[k] = m.col[k];
    return m;
}

static void free_csr(csr_t *m) {
    free(m->row); free(m->col); free(m->cnt); free(m->row64); free(m->col64);
}

static sage1_csr view(const csr_t *m) {
    sage1_csr v = {m->n, m->e, m->row, m->col, m->cnt};
    return v;
}

typedef struct {
    struct futhark_i64_1d *row, *col;
    struct futhark_u64_1d *cnt;
} fut_csr;

static fut_csr upload(const csr_t *m) {
    fut_csr f;
    f.row = futhark_new_i64_1d(ctx, m->row64, (int64_t)m->n + 1);
    f.col = futhark_new_i64_1d(ctx, m->col64, (int64_t)m->e);
    f.cnt = futhark_new_u64_1d(ctx, m->cnt, (int64_t)m->e);
    if (!f.row || !f.col || !f.cnt) { fprintf(stderr, "upload failed\n"); exit(2); }
    return f;
}

static void free_fut(fut_csr *f) {
    futhark_free_i64_1d(ctx, f->row); futhark_free_i64_1d(ctx, f->col); futhark_free_u64_1d(ctx, f->cnt);
}

static double max_abs_diff(const double *a, const double *b, size_t n) {
    double d = 0;
    for (size_t i = 0; i < n; i++) d = fmax(d, fabs(a[i] - b[i]));
    return d;
}

typedef struct futhark_opaque_tup2_i32_arr1d_f64 tup2;
typedef struct futhark_opaque_tup3_i32_i64_f64 tup3;

/* Unpack an (i32 status, []f64) result into a host buffer. */
static int32_t take(tup2 *t, double *out) {
    int32_t st = 0;
    struct futhark_f64_1d *a = NULL;
    FUT(futhark_project_opaque_tup2_i32_arr1d_f64_0(ctx, &st, t));
    FUT(futhark_project_opaque_tup2_i32_arr1d_f64_1(ctx, &a, t));
    FUT(futhark_values_f64_1d(ctx, a, out));
    FUT(futhark_context_sync(ctx));
    futhark_free_f64_1d(ctx, a);
    futhark_free_opaque_tup2_i32_arr1d_f64(ctx, t);
    return st;
}

static int32_t take_select(int (*entry)(struct futhark_context *, tup3 **, const struct futhark_f64_1d *,
                                         const struct futhark_u8_1d *),
                           const struct futhark_f64_1d *x, const struct futhark_u8_1d *mask,
                           int64_t *sel, double *prob) {
    tup3 *t = NULL;
    int32_t st = 0;
    FUT(entry(ctx, &t, x, mask));
    FUT(futhark_project_opaque_tup3_i32_i64_f64_0(ctx, &st, t));
    FUT(futhark_project_opaque_tup3_i32_i64_f64_1(ctx, sel, t));
    FUT(futhark_project_opaque_tup3_i32_i64_f64_2(ctx, prob, t));
    FUT(futhark_context_sync(ctx));
    futhark_free_opaque_tup3_i32_i64_f64(ctx, t);
    return st;
}

static void random_unit(double *x, size_t n, int sparse) {
    double s = 0;
    for (size_t i = 0; i < n; i++) {
        x[i] = (sparse && (next_u64() % 4)) ? 0.0 : next_unit();
        s += x[i] * x[i];
    }
    if (s == 0) { x[0] = 1; s = 1; }
    s = sqrt(s);
    for (size_t i = 0; i < n; i++) x[i] /= s;
}

static void check_matrix(size_t n, size_t max_deg) {
    csr_t m = make_csr(n, max_deg);
    sage1_csr cm = view(&m);
    CHECK(sage1_validate(&cm) == SAGE1_OK, "fixture invalid");
    fut_csr f = upload(&m);
    double *x = malloc(n * sizeof *x), *p = malloc(n * sizeof *p);
    double *qc = malloc(n * sizeof *qc), *qf = malloc(n * sizeof *qf);
    double *outc = malloc(n * sizeof *outc), *outf = malloc(n * sizeof *outf);
    double *s1 = malloc(n * sizeof *s1), *s2 = malloc(n * sizeof *s2);
    double *y = malloc(n * sizeof *y);
    uint8_t *mask = malloc(n);
    int32_t st;

    /* propagate */
    random_unit(x, n, 1);
    for (size_t i = 0; i < n; i++) p[i] = x[i] * x[i];
    double mass = 0; for (size_t i = 0; i < n; i++) mass += p[i];
    for (size_t i = 0; i < n; i++) p[i] /= mass;
    CHECK(sage1_propagate(&cm, p, qc) == SAGE1_OK, "C propagate");
    struct futhark_f64_1d *pin = futhark_new_f64_1d(ctx, p, (int64_t)n);
    { tup2 *t_; FUT(futhark_entry_propagate(ctx, &t_, f.row, f.col, f.cnt, pin)); st = take(t_, qf); }
    futhark_free_f64_1d(ctx, pin);
    CHECK(st == SAGE1_OK, "futhark propagate status %d", st);
    double d = max_abs_diff(qc, qf, n);
    max_err_prop = fmax(max_err_prop, d);
    CHECK(d < TOL, "propagate n=%zu diff %.3g", n, d);

    /* step at several lambdas */
    const double lambdas[] = {0.0, 0.25, 0.5, 0.9, 1.0};
    for (size_t li = 0; li < sizeof lambdas / sizeof *lambdas; li++) {
        random_unit(x, n, li % 2);
        CHECK(sage1_step(&cm, x, lambdas[li], outc, s1, s2) == SAGE1_OK, "C step");
        struct futhark_f64_1d *xin = futhark_new_f64_1d(ctx, x, (int64_t)n);
        { tup2 *t_; FUT(futhark_entry_step(ctx, &t_, f.row, f.col, f.cnt, xin, lambdas[li])); st = take(t_, outf); }
        futhark_free_f64_1d(ctx, xin);
        CHECK(st == SAGE1_OK, "futhark step status %d", st);
        d = max_abs_diff(outc, outf, n);
        max_err_step = fmax(max_err_step, d);
        CHECK(d < TOL, "step n=%zu lambda=%g diff %.3g", n, lambdas[li], d);

        /* select on the identical (C) state with a random mask */
        for (size_t i = 0; i < n; i++) mask[i] = (uint8_t)(next_u64() % 3 == 0);
        size_t selc = 0; double probc = 0;
        sage1_status sc = sage1_select(outc, mask, n, &selc, &probc);
        struct futhark_f64_1d *xs = futhark_new_f64_1d(ctx, outc, (int64_t)n);
        struct futhark_u8_1d *mk = futhark_new_u8_1d(ctx, mask, (int64_t)n);
        int64_t self = 0; double probf = 0;
        st = take_select(futhark_entry_select, xs, mk, &self, &probf);
        futhark_free_f64_1d(ctx, xs); futhark_free_u8_1d(ctx, mk);
        CHECK(st == (int32_t)sc, "select status C %d futhark %d", sc, st);
        if (sc == SAGE1_OK) {
            CHECK((size_t)self == selc, "select index C %zu futhark %lld", selc, (long long)self);
            max_err_sel = fmax(max_err_sel, fabs(probc - probf));
            CHECK(fabs(probc - probf) < TOL, "select prob diff %.3g", fabs(probc - probf));
        }
        n_cases++;
    }

    /* 20-step evolution: C loop vs fused Futhark loop */
    random_unit(x, n, 1);
    memcpy(y, x, n * sizeof *y);
    for (int t = 0; t < 20; t++) {
        CHECK(sage1_step(&cm, y, 0.6, outc, s1, s2) == SAGE1_OK, "C evolve");
        memcpy(y, outc, n * sizeof *y);
    }
    struct futhark_f64_1d *x0 = futhark_new_f64_1d(ctx, x, (int64_t)n);
    { tup2 *t_; FUT(futhark_entry_evolve(ctx, &t_, f.row, f.col, f.cnt, x0, 0.6, 20)); st = take(t_, outf); }
    futhark_free_f64_1d(ctx, x0);
    CHECK(st == SAGE1_OK, "futhark evolve status %d", st);
    d = max_abs_diff(y, outf, n);
    max_err_evolve = fmax(max_err_evolve, d);
    CHECK(d < 1e-11, "evolve n=%zu diff %.3g", n, d);

    /* geodesic between two random unit vectors */
    random_unit(x, n, 0); random_unit(y, n, 1);
    CHECK(sage1_geodesic(x, y, n, 0.3, outc) == SAGE1_OK, "C geodesic");
    struct futhark_f64_1d *gx = futhark_new_f64_1d(ctx, x, (int64_t)n);
    struct futhark_f64_1d *gy = futhark_new_f64_1d(ctx, y, (int64_t)n);
    { tup2 *t_; FUT(futhark_entry_geodesic(ctx, &t_, gx, gy, 0.3)); st = take(t_, outf); }
    futhark_free_f64_1d(ctx, gx); futhark_free_f64_1d(ctx, gy);
    d = max_abs_diff(outc, outf, n);
    max_err_geo = fmax(max_err_geo, d);
    CHECK(st == SAGE1_OK && d < TOL, "geodesic n=%zu diff %.3g", n, d);

    free(x); free(p); free(qc); free(qf); free(outc); free(outf);
    free(s1); free(s2); free(y); free(mask);
    free_fut(&f);
    free_csr(&m);
}

static void check_statuses(void) {
    /* bad lambda, non-unit state, dead end, condition */
    const size_t rows[] = {0, 2, 3, 4};
    const uint32_t cols[] = {1, 2, 2, 0};
    const uint64_t counts[] = {3, 1, 1, 1};
    csr_t m = {3, 4, (size_t *)rows, (uint32_t *)cols, (uint64_t *)counts, NULL, NULL};
    int64_t r64[] = {0, 2, 3, 4}, c64[] = {1, 2, 2, 0};
    m.row64 = r64; m.col64 = c64;
    fut_csr f = upload(&m);
    int32_t st;
    double x[] = {1, 0, 0}, half[] = {0.5, 0.5, 0}, out[3];
    struct futhark_f64_1d *xi = futhark_new_f64_1d(ctx, x, 3);
    { tup2 *t_; FUT(futhark_entry_step(ctx, &t_, f.row, f.col, f.cnt, xi, 1.5)); st = take(t_, out); }
    CHECK(st == SAGE1_BAD_INPUT, "bad lambda status %d", st);
    { tup2 *t_; FUT(futhark_entry_step(ctx, &t_, f.row, f.col, f.cnt, xi, 1.0)); st = take(t_, out); }
    CHECK(st == SAGE1_OK && fabs(out[1] * out[1] - 0.75) < 1e-12 && fabs(out[2] * out[2] - 0.25) < 1e-12,
          "original fixture step");
    struct futhark_f64_1d *hi = futhark_new_f64_1d(ctx, half, 3);
    { tup2 *t_; FUT(futhark_entry_step(ctx, &t_, f.row, f.col, f.cnt, hi, 0.5)); st = take(t_, out); }
    CHECK(st == SAGE1_BAD_INPUT, "non-unit status %d", st);
    uint8_t none[] = {0, 1, 1};
    struct futhark_u8_1d *mk = futhark_new_u8_1d(ctx, none, 3);
    int64_t sel; double prob;
    st = take_select(futhark_entry_select, xi, mk, &sel, &prob);
    CHECK(st == SAGE1_DEAD_END, "dead end status %d", st);
    { tup2 *t_; FUT(futhark_entry_condition(ctx, &t_, 1, 3)); st = take(t_, out); }
    CHECK(st == SAGE1_OK && out[0] == 0 && out[1] == 1 && out[2] == 0, "condition");
    { tup2 *t_; FUT(futhark_entry_condition(ctx, &t_, 3, 3)); st = take(t_, out); }
    CHECK(st == SAGE1_BAD_INPUT, "condition range status %d", st);
    futhark_free_f64_1d(ctx, xi); futhark_free_f64_1d(ctx, hi); futhark_free_u8_1d(ctx, mk);
    free_fut(&f);
}

int main(void) {
    struct futhark_context_config *cfg = futhark_context_config_new();
    ctx = futhark_context_new(cfg);
    char *err = futhark_context_get_error(ctx);
    if (err) { fprintf(stderr, "context: %s\n", err); return 2; }

    size_t max_n = 20000;
    const char *env = getenv("SAGE1_FUT_MAX_N");
    if (env) max_n = (size_t)strtoull(env, NULL, 10);

    check_statuses();
    const size_t sizes[] = {2, 3, 5, 8, 17, 64, 257, 1000, 4096, 20000};
    for (size_t s = 0; s < sizeof sizes / sizeof *sizes; s++) {
        if (sizes[s] > max_n) continue;
        for (int rep = 0; rep < 4; rep++) check_matrix(sizes[s], 1 + (size_t)rep * 4);
    }

    printf("SAGE-1 Futhark (%s) vs C: %ld step/select cases, max |diff| "
           "propagate %.3g, step %.3g, evolve(20) %.3g, geodesic %.3g, select prob %.3g\n",
           FUT_BACKEND, n_cases, max_err_prop, max_err_step, max_err_evolve, max_err_geo, max_err_sel);
    futhark_context_free(ctx);
    futhark_context_config_free(cfg);
    if (failures) { fprintf(stderr, "%d failure(s)\n", failures); return 1; }
    puts("SAGE-1 Futhark cross-check: all checks passed");
    return 0;
}
