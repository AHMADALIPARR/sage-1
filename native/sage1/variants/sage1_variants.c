/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * SAGE-1 variants A (temperature-scaled log-ratio step) and B (low-rank
 * residual projection). Definitions in sage1_variants.h.
 */
#include "sage1_variants.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>

static int finite_unit(const double *v, size_t n) {
    if (!v || n == 0) return 0;
    double s = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(v[i]) || v[i] < 0.0) return 0;
        s += v[i] * v[i];
    }
    return isfinite(s) && fabs(s - 1.0) <= 1e-8;
}

sage1_status sage1_logratio_step(const sage1_csr *m, const double *x, double tau,
                                 double *out, double *q) {
    if (!m || !x || !out || !q || out == x || q == x || q == out ||
        !isfinite(tau) || !(tau > 0.0) || !finite_unit(x, m->states))
        return SAGE1_BAD_INPUT;
    size_t n = m->states;
    for (size_t i = 0; i < n; i++) q[i] = x[i] * x[i];
    sage1_status st = sage1_propagate(m, q, out);   /* out holds P^T(x⊙x) */
    if (st != SAGE1_OK) return st;
    double lmax = -INFINITY;
    for (size_t i = 0; i < n; i++)
        if (out[i] > 0.0) { double l = log(out[i]); if (l > lmax) lmax = l; }
    if (!isfinite(lmax)) return SAGE1_NUMERIC;
    double norm2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double v = out[i] > 0.0 ? exp((log(out[i]) - lmax) / (2.0 * tau)) : 0.0;
        if (!isfinite(v)) return SAGE1_NUMERIC;
        out[i] = v;
        norm2 += v * v;
    }
    if (!(norm2 > 0.0) || !isfinite(norm2)) return SAGE1_NUMERIC;
    double inv = 1.0 / sqrt(norm2);
    for (size_t i = 0; i < n; i++) out[i] *= inv;
    return SAGE1_OK;
}

sage1_status sage1_lowrank_step(const sage1_csr *m, const double *x, double lambda,
                                const double *U, const double *V, size_t rank,
                                double *out, double *y, double *z, double *w) {
    if (!m || !x || !out || !y || !z || out == x || out == y || out == z ||
        y == x || z == x || y == z || !finite_unit(x, m->states) ||
        (rank > 0 && (!U || !V || !w)))
        return SAGE1_BAD_INPUT;
    size_t n = m->states;
    if (rank > n) return SAGE1_BAD_INPUT;
    for (size_t i = 0; i < n; i++) z[i] = x[i] * x[i];
    sage1_status st = sage1_propagate(m, z, y);
    if (st != SAGE1_OK) return st;
    for (size_t i = 0; i < n; i++) y[i] = sqrt(y[i]);
    for (size_t r = 0; r < rank; r++) {
        double s = 0.0;
        for (size_t i = 0; i < n; i++) s += V[i * rank + r] * y[i];
        w[r] = s;
    }
    double norm2 = 0.0;
    for (size_t i = 0; i < n; i++) {
        double s = y[i];
        for (size_t r = 0; r < rank; r++) s += U[i * rank + r] * w[r];
        s = fabs(s);
        if (!isfinite(s)) return SAGE1_NUMERIC;
        z[i] = s;
        norm2 += s * s;
    }
    if (!(norm2 > 0.0) || !isfinite(norm2)) return SAGE1_NUMERIC;
    double inv = 1.0 / sqrt(norm2);
    for (size_t i = 0; i < n; i++) z[i] *= inv;
    return sage1_geodesic(x, z, n, lambda, out);
}

sage1_status sage1_sym_eig(double *A, size_t n, double *evals, double *evecs) {
    if (!A || !evals || !evecs || n == 0) return SAGE1_BAD_INPUT;
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) evecs[i * n + j] = (i == j) ? 1.0 : 0.0;
    for (int sweep = 0; sweep < 100; sweep++) {
        double off = 0.0, diag = 0.0;
        for (size_t i = 0; i < n; i++) {
            diag += A[i * n + i] * A[i * n + i];
            for (size_t j = i + 1; j < n; j++) off += A[i * n + j] * A[i * n + j];
        }
        if (off <= 1e-30 * (diag > 0.0 ? diag : 1.0)) break;
        for (size_t p = 0; p < n; p++)
            for (size_t r = p + 1; r < n; r++) {
                double apr = A[p * n + r];
                if (apr == 0.0) continue;
                double theta = (A[r * n + r] - A[p * n + p]) / (2.0 * apr);
                double t = (theta >= 0.0 ? 1.0 : -1.0) /
                           (fabs(theta) + sqrt(theta * theta + 1.0));
                double c = 1.0 / sqrt(t * t + 1.0), s = t * c;
                for (size_t k = 0; k < n; k++) {
                    double akp = A[k * n + p], akr = A[k * n + r];
                    A[k * n + p] = c * akp - s * akr;
                    A[k * n + r] = s * akp + c * akr;
                }
                for (size_t k = 0; k < n; k++) {
                    double apk = A[p * n + k], ark = A[r * n + k];
                    A[p * n + k] = c * apk - s * ark;
                    A[r * n + k] = s * apk + c * ark;
                }
                for (size_t k = 0; k < n; k++) {
                    double vkp = evecs[k * n + p], vkr = evecs[k * n + r];
                    evecs[k * n + p] = c * vkp - s * vkr;
                    evecs[k * n + r] = s * vkp + c * vkr;
                }
            }
    }
    for (size_t i = 0; i < n; i++) {
        evals[i] = A[i * n + i];
        if (!isfinite(evals[i])) return SAGE1_NUMERIC;
    }
    /* selection sort, descending; ties keep the lower index first */
    for (size_t i = 0; i < n; i++) {
        size_t best = i;
        for (size_t j = i + 1; j < n; j++) if (evals[j] > evals[best]) best = j;
        if (best != i) {
            double t = evals[i]; evals[i] = evals[best]; evals[best] = t;
            for (size_t k = 0; k < n; k++) {
                t = evecs[k * n + i];
                evecs[k * n + i] = evecs[k * n + best];
                evecs[k * n + best] = t;
            }
        }
    }
    for (size_t j = 0; j < n; j++) {
        size_t arg = 0;
        for (size_t k = 1; k < n; k++)
            if (fabs(evecs[k * n + j]) > fabs(evecs[arg * n + j]) + 1e-12) arg = k;
        if (evecs[arg * n + j] < 0.0)
            for (size_t k = 0; k < n; k++) evecs[k * n + j] = -evecs[k * n + j];
    }
    return SAGE1_OK;
}

sage1_status sage1_lowrank_fit(const double *S, const double *T, size_t n,
                               size_t rank, double *U, double *V) {
    if (!S || !T || n == 0 || rank > n || (rank > 0 && (!U || !V)))
        return SAGE1_BAD_INPUT;
    if (rank == 0) return SAGE1_OK;
    sage1_status st = SAGE1_OK;
    double *A = malloc(n * n * sizeof *A), *X = malloc(n * n * sizeof *X);
    double *G = malloc(n * n * sizeof *G), *E = malloc(n * n * sizeof *E);
    double *ev = malloc(n * sizeof *ev);
    if (!A || !X || !G || !E || !ev) { st = SAGE1_NUMERIC; goto done; }
    /* Solve S X = R, R = T − S. */
    memcpy(A, S, n * n * sizeof *A);
    for (size_t i = 0; i < n * n; i++) {
        if (!isfinite(S[i]) || !isfinite(T[i])) { st = SAGE1_BAD_INPUT; goto done; }
        X[i] = T[i] - S[i];
    }
    for (size_t c = 0; c < n; c++) {
        size_t piv = c;
        for (size_t r = c + 1; r < n; r++)
            if (fabs(A[r * n + c]) > fabs(A[piv * n + c])) piv = r;
        if (!(fabs(A[piv * n + c]) > 1e-12)) { st = SAGE1_NUMERIC; goto done; }
        if (piv != c)
            for (size_t k = 0; k < n; k++) {
                double t = A[c * n + k]; A[c * n + k] = A[piv * n + k]; A[piv * n + k] = t;
                t = X[c * n + k]; X[c * n + k] = X[piv * n + k]; X[piv * n + k] = t;
            }
        for (size_t r = 0; r < n; r++) {
            if (r == c) continue;
            double f = A[r * n + c] / A[c * n + c];
            if (f == 0.0) continue;
            for (size_t k = 0; k < n; k++) {
                A[r * n + k] -= f * A[c * n + k];
                X[r * n + k] -= f * X[c * n + k];
            }
        }
    }
    for (size_t r = 0; r < n; r++) {
        double d = A[r * n + r];
        for (size_t k = 0; k < n; k++) X[r * n + k] /= d;
    }
    /* Left singular vectors of R^T = eigenvectors of R^T R. */
    for (size_t i = 0; i < n; i++)
        for (size_t j = 0; j < n; j++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++)
                s += (T[k * n + i] - S[k * n + i]) * (T[k * n + j] - S[k * n + j]);
            G[i * n + j] = s;
        }
    st = sage1_sym_eig(G, n, ev, E);
    if (st != SAGE1_OK) goto done;
    for (size_t i = 0; i < n; i++)
        for (size_t r = 0; r < rank; r++) U[i * rank + r] = E[i * n + r];
    for (size_t i = 0; i < n; i++)
        for (size_t r = 0; r < rank; r++) {
            double s = 0.0;
            for (size_t k = 0; k < n; k++) s += X[i * n + k] * E[k * n + r];
            V[i * rank + r] = s;
        }
done:
    free(A); free(X); free(G); free(E); free(ev);
    return st;
}
