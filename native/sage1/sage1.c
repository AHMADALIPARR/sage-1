#include "sage1.h"
#include <math.h>
#include <float.h>
#include <string.h>

static int finite_unit(const double *v, size_t n) {
    if (!v || n == 0) return 0;
    double s=0.0;
    for (size_t i=0;i<n;i++) {
        if (!isfinite(v[i]) || v[i]<0.0) return 0;
        s+=v[i]*v[i];
    }
    return isfinite(s) && fabs(s-1.0)<=1e-8;
}

sage1_status sage1_validate(const sage1_csr *m) {
    if (!m || !m->row_ptr || !m->col_idx || !m->counts || !m->states || !m->edges ||
        m->states>UINT32_MAX || m->row_ptr[0]!=0 || m->row_ptr[m->states]!=m->edges)
        return SAGE1_BAD_INPUT;
    for (size_t i=0;i<m->states;i++) {
        size_t a=m->row_ptr[i], b=m->row_ptr[i+1];
        if (a>=b || b>m->edges) return SAGE1_BAD_MATRIX;
        uint64_t total=0;
        for (size_t k=a;k<b;k++) {
            if (m->col_idx[k]>=m->states || m->counts[k]==0 ||
                UINT64_MAX-total<m->counts[k]) return SAGE1_BAD_MATRIX;
            total+=m->counts[k];
        }
    }
    return SAGE1_OK;
}

sage1_status sage1_propagate(const sage1_csr *m, const double *p, double *q) {
    sage1_status status=sage1_validate(m);
    if (status!=SAGE1_OK) return status;
    if (!p || !q || p==q) return SAGE1_BAD_INPUT;
    double mass=0.0;
    for (size_t i=0;i<m->states;i++) {
        if (!isfinite(p[i]) || p[i]<0.0) return SAGE1_BAD_INPUT;
        mass+=p[i];
    }
    if (!isfinite(mass) || fabs(mass-1.0)>1e-8) return SAGE1_BAD_INPUT;
    memset(q,0,m->states*sizeof(double));
    for (size_t i=0;i<m->states;i++) {
        uint64_t total=0;
        for (size_t k=m->row_ptr[i];k<m->row_ptr[i+1];k++) total+=m->counts[k];
        for (size_t k=m->row_ptr[i];k<m->row_ptr[i+1];k++)
            q[m->col_idx[k]]+=p[i]*((double)m->counts[k]/(double)total);
    }
    double sum=0.0;
    for (size_t j=0;j<m->states;j++) sum+=q[j];
    if (!(sum>0.0) || !isfinite(sum)) return SAGE1_NUMERIC;
    for (size_t j=0;j<m->states;j++) q[j]/=sum;
    return SAGE1_OK;
}

sage1_status sage1_geodesic(const double *x, const double *y, size_t n,
                             double lambda, double *out) {
    if (!out || out==x || out==y || !finite_unit(x,n) || !finite_unit(y,n) ||
        !isfinite(lambda) || lambda<0.0 || lambda>1.0) return SAGE1_BAD_INPUT;
    double dot=0.0;
    for (size_t i=0;i<n;i++) dot+=x[i]*y[i];
    dot=fmax(-1.0,fmin(1.0,dot));
    double theta=acos(dot), sin_theta=sin(theta);
    double a,b;
    if (theta<1e-7) {a=1.0-lambda;b=lambda;}
    else {a=sin((1.0-lambda)*theta)/sin_theta;
          b=sin(lambda*theta)/sin_theta;}
    double norm2=0.0;
    for (size_t i=0;i<n;i++) {
        out[i]=a*x[i]+b*y[i];
        if (out[i]<0 && out[i]>-1e-14) out[i]=0;
        if (!isfinite(out[i]) || out[i]<0) return SAGE1_NUMERIC;
        norm2+=out[i]*out[i];
    }
    if (!(norm2>0.0) || !isfinite(norm2)) return SAGE1_NUMERIC;
    double inv=1.0/sqrt(norm2);
    for (size_t i=0;i<n;i++) out[i]*=inv;
    return SAGE1_OK;
}

sage1_status sage1_step(const sage1_csr *m, const double *x, double lambda,
                         double *out, double *q, double *y) {
    if (!m || !out || !q || !y || !x || out==x || out==q || out==y ||
        q==x || y==x || q==y || !finite_unit(x,m->states)) return SAGE1_BAD_INPUT;
    for (size_t i=0;i<m->states;i++) q[i]=x[i]*x[i];
    /* propagate needs distinct buffers; y temporarily carries incoming probabilities */
    sage1_status status=sage1_propagate(m,q,y);
    if (status!=SAGE1_OK) return status;
    for (size_t i=0;i<m->states;i++) y[i]=sqrt(y[i]);
    /* q is now disposable; use it for output-independent scratch only. */
    return sage1_geodesic(x,y,m->states,lambda,out);
}

sage1_status sage1_select(const double *x, const uint8_t *legal, size_t n,
                           size_t *selected, double *probability) {
    if (!legal || !selected || !probability || !finite_unit(x,n)) return SAGE1_BAD_INPUT;
    double total=0.0, best=-1.0;
    size_t winner=0;
    for (size_t i=0;i<n;i++) if (legal[i]) {
        double weight=x[i]*x[i];
        total+=weight;
        if (weight>best) {best=weight;winner=i;}
    }
    if (!(total>0.0)) return SAGE1_DEAD_END;
    *selected=winner;
    *probability=best/total;
    return SAGE1_OK;
}

sage1_status sage1_condition(size_t selected, size_t n, double *x) {
    if (!x || n==0 || selected>=n) return SAGE1_BAD_INPUT;
    memset(x,0,n*sizeof(double));
    x[selected]=1.0;
    return SAGE1_OK;
}
