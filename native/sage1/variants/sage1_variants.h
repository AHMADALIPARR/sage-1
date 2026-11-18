/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * Two SAGE-1 step variants built on the unchanged kernel in ../sage1.c.
 *
 * Variant A, temperature-scaled log-ratio step (replaces the SLERP angle):
 *     q   = P^T (x ⊙ x)                       (sage1_propagate)
 *     out = q^(1/(2τ)) / || q^(1/(2τ)) ||      componentwise, 0 stays 0
 *   so log(out_i² / out_j²) = (1/τ) · log(q_i / q_j) for every pair i, j:
 *   the target's log-ratios are scaled by 1/τ and the state moves all the
 *   way there. τ = 1 is sage1_step at λ = 1; τ < 1 sharpens, τ > 1 flattens.
 *   Output is a nonnegative unit vector (computed with a max shift in the
 *   log domain so no power overflows).
 *
 * Variant B, low-rank residual projection after CSR propagation:
 *     y   = sqrt( P^T (x ⊙ x) )
 *     z   = y + U (V^T y)                      U, V: states × rank
 *     ŷ   = |z| / || z ||                      (|·| keeps the nonnegative orthant;
 *                                               probabilities z_i² are unchanged)
 *     out = SLERP(x, ŷ, λ)                     (sage1_geodesic)
 *   The grammar mask is applied afterwards by sage1_select, as for sage1_step.
 */
#ifndef UOR_SAGE1_VARIANTS_H
#define UOR_SAGE1_VARIANTS_H
#include "sage1.h"
#ifdef __cplusplus
extern "C" {
#endif

/* Variant A. tau finite and > 0; q scratch length m->states; out distinct. */
sage1_status sage1_logratio_step(const sage1_csr *m, const double *x, double tau,
                                 double *out, double *q);

/* Variant B. U and V are row-major states × rank; rank 0 is sage1_step.
 * Scratch y and z length m->states, w length rank (may be NULL if rank 0). */
sage1_status sage1_lowrank_step(const sage1_csr *m, const double *x, double lambda,
                                const double *U, const double *V, size_t rank,
                                double *out, double *y, double *z, double *w);

/* Fit U, V for variant B from amplitude rows.
 *   S: n × n, row k = sqrt(P[k,:]) (the propagated amplitude from state e_k)
 *   T: n × n, row k = target amplitude for state e_k
 * Residual R = T − S. Solve S X = R (Gaussian elimination, partial pivoting),
 * so M = X^T maps every S row to its residual: M S[k]^T = R[k]^T. The rank-r
 * reduced-rank regression of R^T on S^T keeps the top-r left singular
 * vectors U_r of R^T (cyclic Jacobi on R^T R, deterministic sweep order,
 * sign fixed so the largest-|·| entry is positive):
 *     U = U_r,  V = X U_r      so  U V^T = U_r U_r^T M.
 * rank = n reproduces every T row exactly from its S row.
 * U and V outputs: row-major n × rank. Returns SAGE1_NUMERIC if S is singular. */
sage1_status sage1_lowrank_fit(const double *S, const double *T, size_t n,
                               size_t rank, double *U, double *V);

/* Symmetric eigen-decomposition by cyclic Jacobi (n × n, row-major).
 * A is overwritten; eigenvalues descending in evals, eigenvectors in the
 * columns of evecs, each with its largest-|·| entry made positive. */
sage1_status sage1_sym_eig(double *A, size_t n, double *evals, double *evecs);

#ifdef __cplusplus
}
#endif
#endif
