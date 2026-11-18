#include "sage1.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

int main(void) {
    const size_t rows[]={0,2,3,4};
    const uint32_t cols[]={1,2,2,0};
    const uint64_t counts[]={3,1,1,1};
    const sage1_csr m={3,4,rows,cols,counts};
    const double p[]={1,0,0};
    double q[3],target[3],out[3],x[]={1,0,0};
    assert(sage1_validate(&m)==SAGE1_OK);
    assert(sage1_propagate(&m,p,q)==SAGE1_OK);
    assert(fabs(q[1]-.75)<1e-12 && fabs(q[2]-.25)<1e-12);
    assert(sage1_step(&m,x,1.0,out,q,target)==SAGE1_OK);
    assert(fabs(out[1]*out[1]-.75)<1e-12);
    assert(fabs(out[2]*out[2]-.25)<1e-12);
    assert(sage1_step(&m,x,0.0,out,q,target)==SAGE1_OK);
    assert(fabs(out[0]-1.0)<1e-12);
    assert(sage1_step(&m,x,0.5,out,q,target)==SAGE1_OK);
    assert(fabs(out[0]*out[0]+out[1]*out[1]+out[2]*out[2]-1.0)<1e-12);
    const uint8_t mask[]={0,1,0};
    size_t selected=99; double probability=0;
    assert(sage1_select(target,mask,3,&selected,&probability)==SAGE1_OK);
    assert(selected==1 && fabs(probability-1.0)<1e-12);
    assert(sage1_condition(selected,3,out)==SAGE1_OK);
    assert(out[1]==1.0 && out[0]==0.0 && out[2]==0.0);
    const uint8_t none[]={0,0,0};
    assert(sage1_select(out,none,3,&selected,&probability)==SAGE1_DEAD_END);
    const size_t invalid_rows[]={0,0,3,4};
    const sage1_csr invalid={3,4,invalid_rows,cols,counts};
    assert(sage1_validate(&invalid)==SAGE1_BAD_MATRIX);
    puts("SAGE-1: all checks passed");
}
