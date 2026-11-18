#ifndef UOR_SAGE1_H
#define UOR_SAGE1_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* CSR row i contains transitions i -> col_idx[k]. */
typedef struct {
    size_t states, edges;
    const size_t *row_ptr;      /* length states + 1 */
    const uint32_t *col_idx;    /* length edges */
    const uint64_t *counts;     /* length edges */
} sage1_csr;
typedef enum {
    SAGE1_OK=0, SAGE1_BAD_INPUT=1, SAGE1_BAD_MATRIX=2,
    SAGE1_DEAD_END=3, SAGE1_NUMERIC=4
} sage1_status;
/* Check indices, row layout, nonzero row totals, and finite positive dimensions. */
sage1_status sage1_validate(const sage1_csr *m);
/* p -> q with row-normalized count matrix; all caller buffers length m->states. */
sage1_status sage1_propagate(const sage1_csr *m, const double *p, double *q);
/* x and y must be unit vectors in the nonnegative orthant; output normalized. */
sage1_status sage1_geodesic(const double *x, const double *y, size_t n,
                             double lambda, double *out);
/* Target sqrt(P^T (x^2)), then geodesic motion; scratch q and y length n. */
sage1_status sage1_step(const sage1_csr *m, const double *x, double lambda,
                         double *out, double *q, double *y);
/* Decode legal transitions with mask of length states; ties -> smallest index. */
sage1_status sage1_select(const double *x, const uint8_t *legal, size_t n,
                           size_t *selected, double *probability);
/* Collapse onto selected branch; use with explicit parser-stack updates. */
sage1_status sage1_condition(size_t selected, size_t n, double *x);
#ifdef __cplusplus
}
#endif
#endif
