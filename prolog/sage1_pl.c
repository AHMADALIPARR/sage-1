/* SPDX-License-Identifier: AGPL-3.0-only OR LicenseRef-SnapKitty-Opaque-1.0
 * Copyright (C) 2026 Ahmad Ali Parr / SnapKitty Research Lab
 *
 * SWI-Prolog foreign library exposing the SAGE-1 C11 kernel.
 *
 *   sage1_validate(+Csr, -Status)
 *   sage1_propagate(+Csr, +P, -Q)
 *   sage1_step(+Csr, +X, +Lambda, -Y)
 *   sage1_geodesic(+X, +Y, +Lambda, -Out)
 *   sage1_select(+X, +Mask, -Result)      Result = selected(Index, Prob) ; dead_end
 *   sage1_condition(+Selected, +N, -X)
 *
 * Csr = csr(States, RowPtr, ColIdx, Counts) with Prolog lists. Vectors are
 * lists of numbers, masks are lists of 0/1. Kernel errors other than a dead
 * end raise error(sage1_error(Status), Context).
 */
#include <SWI-Prolog.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "sage1.h"

static const char *status_name(sage1_status s) {
    switch (s) {
    case SAGE1_OK: return "ok";
    case SAGE1_BAD_INPUT: return "bad_input";
    case SAGE1_BAD_MATRIX: return "bad_matrix";
    case SAGE1_DEAD_END: return "dead_end";
    case SAGE1_NUMERIC: return "numeric";
    }
    return "unknown";
}

static int raise_status(sage1_status s, const char *pred) {
    term_t ex = PL_new_term_ref();
    if (!PL_unify_term(ex,
            PL_FUNCTOR_CHARS, "error", 2,
              PL_FUNCTOR_CHARS, "sage1_error", 1, PL_CHARS, status_name(s),
              PL_FUNCTOR_CHARS, "context", 2, PL_CHARS, pred, PL_VARIABLE))
        return FALSE;
    return PL_raise_exception(ex);
}

static int list_length(term_t list, size_t *n) {
    size_t len = 0;
    if (PL_skip_list(list, 0, &len) != PL_LIST) return PL_type_error("list", list);
    *n = len;
    return TRUE;
}

static int get_doubles(term_t list, size_t expect, double **out) {
    size_t n = 0;
    if (!list_length(list, &n)) return FALSE;
    if (n != expect) return PL_domain_error("vector_of_states_length", list);
    double *v = malloc((n ? n : 1) * sizeof *v);
    if (!v) return PL_resource_error("memory");
    term_t head = PL_new_term_ref(), tail = PL_copy_term_ref(list);
    for (size_t i = 0; PL_get_list(tail, head, tail); i++) {
        if (!PL_get_float(head, &v[i])) { free(v); return PL_type_error("number", head); }
    }
    *out = v;
    return TRUE;
}

static int unify_doubles(term_t t, const double *v, size_t n) {
    term_t l = PL_new_term_ref();
    PL_put_nil(l);
    for (size_t i = n; i-- > 0;) {
        term_t h = PL_new_term_ref();
        if (!PL_put_float(h, v[i]) || !PL_cons_list(l, h, l)) return FALSE;
    }
    return PL_unify(t, l);
}

typedef struct {
    sage1_csr m;
    size_t *row;
    uint32_t *col;
    uint64_t *cnt;
} pl_csr;

static void free_csr(pl_csr *c) { free(c->row); free(c->col); free(c->cnt); }

static int get_u64_list(term_t list, size_t *n, uint64_t **out) {
    if (!list_length(list, n)) return FALSE;
    uint64_t *v = malloc((*n ? *n : 1) * sizeof *v);
    if (!v) return PL_resource_error("memory");
    term_t head = PL_new_term_ref(), tail = PL_copy_term_ref(list);
    for (size_t i = 0; PL_get_list(tail, head, tail); i++) {
        int64_t x;
        if (!PL_get_int64_ex(head, &x)) { free(v); return FALSE; }
        if (x < 0) { free(v); return PL_domain_error("not_less_than_zero", head); }
        v[i] = (uint64_t)x;
    }
    *out = v;
    return TRUE;
}

static int get_csr(term_t t, pl_csr *c) {
    memset(c, 0, sizeof *c);
    if (!PL_is_functor(t, PL_new_functor(PL_new_atom("csr"), 4)))
        return PL_type_error("csr", t);
    term_t a = PL_new_term_ref();
    int64_t states;
    size_t nrow, ncol, ncnt;
    uint64_t *rows = NULL, *cols = NULL;
    _PL_get_arg(1, t, a);
    if (!PL_get_int64_ex(a, &states)) return FALSE;
    if (states <= 0) return PL_domain_error("positive_integer", a);
    _PL_get_arg(2, t, a);
    if (!get_u64_list(a, &nrow, &rows)) return FALSE;
    _PL_get_arg(3, t, a);
    if (!get_u64_list(a, &ncol, &cols)) { free(rows); return FALSE; }
    _PL_get_arg(4, t, a);
    if (!get_u64_list(a, &ncnt, &c->cnt)) { free(rows); free(cols); return FALSE; }
    int ok = nrow == (size_t)states + 1 && ncol == ncnt;
    c->row = malloc(nrow * sizeof *c->row);
    c->col = malloc((ncol ? ncol : 1) * sizeof *c->col);
    if (!c->row || !c->col) ok = 0;
    for (size_t i = 0; ok && i < nrow; i++) {
        if (rows[i] > ncol) ok = 0; else c->row[i] = (size_t)rows[i];
    }
    for (size_t i = 0; ok && i < ncol; i++) {
        if (cols[i] > UINT32_MAX) ok = 0; else c->col[i] = (uint32_t)cols[i];
    }
    free(rows); free(cols);
    if (!ok) { free_csr(c); return raise_status(SAGE1_BAD_MATRIX, "sage1_csr/4"); }
    c->m.states = (size_t)states;
    c->m.edges = ncol;
    c->m.row_ptr = c->row;
    c->m.col_idx = c->col;
    c->m.counts = c->cnt;
    return TRUE;
}

static foreign_t pl_validate(term_t csr, term_t status) {
    pl_csr c;
    if (!get_csr(csr, &c)) return FALSE;
    sage1_status s = sage1_validate(&c.m);
    free_csr(&c);
    return PL_unify_atom_chars(status, status_name(s));
}

static foreign_t pl_propagate(term_t csr, term_t p, term_t q) {
    pl_csr c;
    double *pv = NULL, *qv = NULL;
    if (!get_csr(csr, &c)) return FALSE;
    int rc = FALSE;
    if (get_doubles(p, c.m.states, &pv) && (qv = malloc(c.m.states * sizeof *qv))) {
        sage1_status s = sage1_propagate(&c.m, pv, qv);
        rc = s == SAGE1_OK ? unify_doubles(q, qv, c.m.states) : raise_status(s, "sage1_propagate/3");
    }
    free(pv); free(qv); free_csr(&c);
    return rc;
}

static foreign_t pl_step(term_t csr, term_t x, term_t lambda, term_t y) {
    pl_csr c;
    double *xv = NULL, *out = NULL, *q = NULL, *yy = NULL, lam;
    if (!PL_get_float_ex(lambda, &lam)) return FALSE;
    if (!get_csr(csr, &c)) return FALSE;
    int rc = FALSE;
    size_t n = c.m.states;
    if (get_doubles(x, n, &xv) && (out = malloc(n * sizeof *out)) &&
        (q = malloc(n * sizeof *q)) && (yy = malloc(n * sizeof *yy))) {
        sage1_status s = sage1_step(&c.m, xv, lam, out, q, yy);
        rc = s == SAGE1_OK ? unify_doubles(y, out, n) : raise_status(s, "sage1_step/4");
    }
    free(xv); free(out); free(q); free(yy); free_csr(&c);
    return rc;
}

static foreign_t pl_geodesic(term_t x, term_t y, term_t lambda, term_t out) {
    size_t n = 0;
    double *xv = NULL, *yv = NULL, *ov = NULL, lam;
    if (!PL_get_float_ex(lambda, &lam) || !list_length(x, &n)) return FALSE;
    int rc = FALSE;
    if (get_doubles(x, n, &xv) && get_doubles(y, n, &yv) &&
        (ov = malloc((n ? n : 1) * sizeof *ov))) {
        sage1_status s = sage1_geodesic(xv, yv, n, lam, ov);
        rc = s == SAGE1_OK ? unify_doubles(out, ov, n) : raise_status(s, "sage1_geodesic/4");
    }
    free(xv); free(yv); free(ov);
    return rc;
}

static foreign_t pl_select(term_t x, term_t mask, term_t result) {
    size_t n = 0, nm = 0;
    double *xv = NULL;
    if (!list_length(x, &n) || !list_length(mask, &nm)) return FALSE;
    if (n != nm) return PL_domain_error("mask_length", mask);
    uint8_t *legal = malloc(n ? n : 1);
    if (!legal) return PL_resource_error("memory");
    term_t head = PL_new_term_ref(), tail = PL_copy_term_ref(mask);
    for (size_t i = 0; PL_get_list(tail, head, tail); i++) {
        int b;
        if (!PL_get_integer_ex(head, &b)) { free(legal); return FALSE; }
        legal[i] = b != 0;
    }
    int rc = FALSE;
    if (get_doubles(x, n, &xv)) {
        size_t sel = 0;
        double prob = 0;
        sage1_status s = sage1_select(xv, legal, n, &sel, &prob);
        if (s == SAGE1_OK)
            rc = PL_unify_term(result, PL_FUNCTOR_CHARS, "selected", 2,
                               PL_INT64, (int64_t)sel, PL_FLOAT, prob);
        else if (s == SAGE1_DEAD_END)
            rc = PL_unify_atom_chars(result, "dead_end");
        else
            rc = raise_status(s, "sage1_select/3");
    }
    free(xv); free(legal);
    return rc;
}

static foreign_t pl_condition(term_t selected, term_t n_t, term_t x) {
    int64_t sel, n;
    if (!PL_get_int64_ex(selected, &sel) || !PL_get_int64_ex(n_t, &n)) return FALSE;
    if (sel < 0 || n <= 0) return raise_status(SAGE1_BAD_INPUT, "sage1_condition/3");
    double *v = malloc((size_t)n * sizeof *v);
    if (!v) return PL_resource_error("memory");
    sage1_status s = sage1_condition((size_t)sel, (size_t)n, v);
    int rc = s == SAGE1_OK ? unify_doubles(x, v, (size_t)n) : raise_status(s, "sage1_condition/3");
    free(v);
    return rc;
}

install_t install_sage1_pl(void) {
    PL_register_foreign_in_module("sage1", "sage1_validate", 2, pl_validate, 0);
    PL_register_foreign_in_module("sage1", "sage1_propagate", 3, pl_propagate, 0);
    PL_register_foreign_in_module("sage1", "sage1_step", 4, pl_step, 0);
    PL_register_foreign_in_module("sage1", "sage1_geodesic", 4, pl_geodesic, 0);
    PL_register_foreign_in_module("sage1", "sage1_select", 3, pl_select, 0);
    PL_register_foreign_in_module("sage1", "sage1_condition", 3, pl_condition, 0);
}
