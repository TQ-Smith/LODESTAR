
#ifndef _REAL_SYM_EIGEN_H_
#define _REAL_SYM_EIGEN_H_

#include <stdlib.h>

// Thanks to Jochen Voss's matroot.c
// https://github.com/Foadsf/Cmathtuts/blob/a3610505e04c162c48ad8b8eab9745b32f2f0b6d/C_CBLAS/seehuhn/example2/matroot.c#L61

typedef struct {
    double *A, *B, *W, *Z, *WORK;
    int *ISUPPZ, *IWORK;
    int N, M;
} RealSymEigen;

RealSymEigen* init_real_sym_eigen(int N);

int compute_k_eigenpairs(RealSymEigen* eigen, int k);

int compute_k_eigenvalues(RealSymEigen* eigen, int k);

void destroy_real_sym_eigen(RealSymEigen* eigen);

#endif