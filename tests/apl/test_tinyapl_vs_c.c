/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Checks the TinyAPL array engine (apl/sage1_engine_checked.tinyapl) against
 * the C11 kernel on its 5-vertex self-test patch. Usage:
 *   tinyapl apl/sage1_engine_checked.tinyapl > out.txt
 *   test_tinyapl_vs_c out.txt
 * The TinyAPL output is x1, x2, ||x1||, ||x2||, DetSelect x2.
 */
#include "sage1.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 5

/* Same patch as the TinyAPL self-test: rows are count-normalized, so the
 * weights 0.5/0.5, 0.3/0.4/0.3 become counts 1/1, 3/4/3. */
static const size_t rows[] = {0, 2, 5, 7, 9, 11};
static const uint32_t cols[] = {0, 1, 0, 1, 2, 1, 3, 2, 4, 3, 4};
static const uint64_t counts[] = {1, 1, 3, 4, 3, 1, 1, 1, 1, 1, 1};

/* Read all numbers from TinyAPL's printed output (¯ is minus, ⏨ is exponent). */
static size_t read_numbers(const char *path, double *out, size_t cap) {
    FILE *f = fopen(path, "rb");
    if (!f) { perror(path); exit(2); }
    static char buf[1 << 16], norm[1 << 16];
    size_t len = fread(buf, 1, sizeof buf - 1, f);
    fclose(f);
    buf[len] = 0;
    size_t j = 0;
    for (size_t i = 0; i < len && j < sizeof norm - 2;) {
        unsigned char c = (unsigned char)buf[i];
        if (c == 0xC2 && (unsigned char)buf[i + 1] == 0xAF) { norm[j++] = '-'; i += 2; }
        else if (c == 0xE2 && (unsigned char)buf[i + 1] == 0x8F && (unsigned char)buf[i + 2] == 0xA8) {
            norm[j++] = 'e'; i += 3;
        } else if ((c >= '0' && c <= '9') || c == '.') { norm[j++] = (char)c; i++; }
        else { norm[j++] = ' '; i++; }
    }
    norm[j] = 0;
    size_t n = 0;
    char *p = norm;
    while (*p && n < cap) {
        char *end;
        double v = strtod(p, &end);
        if (end == p) { p++; continue; }
        out[n++] = v;
        p = end;
    }
    return n;
}

static void unit(double *x) {
    double s = 0;
    for (int i = 0; i < N; i++) s += x[i] * x[i];
    s = sqrt(s);
    for (int i = 0; i < N; i++) x[i] /= s;
}

int main(int argc, char **argv) {
    if (argc != 2) { fprintf(stderr, "usage: %s tinyapl-output.txt\n", argv[0]); return 2; }
    double got[64];
    size_t k = read_numbers(argv[1], got, 64);
    if (k != 2 * N + 3) { fprintf(stderr, "expected %d numbers from TinyAPL, got %zu\n", 2 * N + 3, k); return 1; }

    const sage1_csr m = {N, 11, rows, cols, counts};
    double x0[N], x1[N], xb[N], x2[N], q[N], y[N];
    for (int i = 0; i < N; i++) x0[i] = 1.0;
    unit(x0);
    if (sage1_step(&m, x0, 0.35, x1, q, y) != SAGE1_OK) return 1;
    const int peak[N] = {1, 0, 1, 0, 1};               /* peaks 1 3 5 (1-based) */
    for (int i = 0; i < N; i++) xb[i] = x1[i] * (peak[i] ? 1.25 : 1.0);
    unit(xb);
    if (sage1_step(&m, xb, 0.35, x2, q, y) != SAGE1_OK) return 1;
    const uint8_t all[N] = {1, 1, 1, 1, 1};
    size_t sel; double prob;
    if (sage1_select(x2, all, N, &sel, &prob) != SAGE1_OK) return 1;

    double err1 = 0, err2 = 0, n1 = 0, n2 = 0;
    for (int i = 0; i < N; i++) {
        err1 = fmax(err1, fabs(got[i] - x1[i]));
        err2 = fmax(err2, fabs(got[N + i] - x2[i]));
        n1 += x1[i] * x1[i]; n2 += x2[i] * x2[i];
    }
    double errn = fmax(fabs(got[2 * N] - sqrt(n1)), fabs(got[2 * N + 1] - sqrt(n2)));
    printf("TinyAPL vs C: max |x1 diff| %.3g, max |x2 diff| %.3g, norm diff %.3g, "
           "DetSelect TinyAPL %g C %zu\n", err1, err2, errn, got[2 * N + 2], sel);
    int ok = err1 < 1e-12 && err2 < 1e-12 && errn < 1e-12 && (size_t)got[2 * N + 2] == sel;
    puts(ok ? "SAGE-1 TinyAPL engine matches the C kernel" : "MISMATCH");
    return ok ? 0 : 1;
}
