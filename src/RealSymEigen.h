
#ifndef _REAL_SYM_EIGEN_H_
#define _REAL_SYM_EIGEN_H_

#include <stdlib.h>

// Thanks to Jochen Voss's matroot.c
// https://github.com/Foadsf/Cmathtuts/blob/a3610505e04c162c48ad8b8eab9745b32f2f0b6d/C_CBLAS/seehuhn/example2/matroot.c#L61

static int dsyevr(
    char JOBZ, 
    char RANGE, 
    char UPLO, 
    int N,
    double *A, 
    int LDA, 
    double VL, 
    double VU,
    int IL, 
    int IU, 
    double ABSTOL, 
    int *M,
    double *W, 
    double *Z, 
    int LDZ, 
    int *ISUPPZ,
    double *WORK, 
    int LWORK, 
    int *IWORK, 
    int LIWORK
) {
    
    extern void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, 
                        int *Np, double *A, int *LDAp, double *VLp, 
                        double *VUp, int *ILp, int *IUp, double *ABSTOLp, 
                        int *Mp, double *W, double *Z, int *LDZp, int *ISUPPZ, 
                        double *WORK, int *LWORKp, int *IWORK, int *LIWORKp, int *INFOp);

    int INFO;

    dsyevr_(&JOBZ, &RANGE, &UPLO, 
            &N, A, &LDA, &VL, 
            &VU, &IL, &IU, &ABSTOL,
            M, W, Z, &LDZ, ISUPPZ, 
            WORK, &LWORK, IWORK, &LIWORK, &INFO);
    
    return INFO;

}

static double dlamch(char CMACH) {
  extern double dlamch_(char *CMACHp);
  return dlamch_(&CMACH);
}

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