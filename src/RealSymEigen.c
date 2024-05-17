
// File: RealSymEigen.c
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute the eigenpairs of a real symmetric matrix.
// Reference: Jochen Voss's matroot.c
//            https://github.com/Foadsf/Cmathtuts/blob/a3610505e04c162c48ad8b8eab9745b32f2f0b6d/C_CBLAS/seehuhn/example2/matroot.c#L61
//
//            We are using LAPACK's dsyevr routine.
//            https://netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga2ad9f4a91cddbf67fe41b621bd158f5c.html

#include "RealSymEigen.h"

// A few things to note:
//  - LAPACK requires gfortran which is not a part of CLANG.
//  - Must be compiled with gcc/gfortran.
//  - All matrices in LAPACK are one-dimensional arrays stored in column major
//      form. This does not matter with symmetric matrices.
//  - Matrix A is used to store the matrix we are using in our eigen computations.
//  - Array W will hold the resulting eigenvalues.
//  - Matrix Z will hold the resulting eigenvectors.
//  - Eigenpairs in W and Z are stored from least to greatest.

// Define our LAPACK methods.
extern void dsyevr_(char *JOBZp, char *RANGEp, char *UPLOp, 
                        int *Np, double *A, int *LDAp, double *VLp, 
                        double *VUp, int *ILp, int *IUp, double *ABSTOLp, 
                        int *Mp, double *W, double *Z, int *LDZp, int *ISUPPZ, 
                        double *WORK, int *LWORKp, int *IWORK, int *LIWORKp, int *INFOp);

extern double dlamch_(char *CMACHp);

// We wrap the methods in a function to make them more readable in C.
int dsyevr(
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

    int INFO;

    dsyevr_(&JOBZ, &RANGE, &UPLO, 
            &N, A, &LDA, &VL, 
            &VU, &IL, &IU, &ABSTOL,
            M, W, Z, &LDZ, ISUPPZ, 
            WORK, &LWORK, IWORK, &LIWORK, &INFO);
    
    return INFO;

}

double dlamch(char CMACH) {
  return dlamch_(&CMACH);
}


// Allocate all the necessary memory.
RealSymEigen_t* init_real_sym_eigen(int N) {
    RealSymEigen_t* eigen = malloc(sizeof(RealSymEigen_t));
    eigen -> A = malloc(N * N * sizeof(double));
    eigen -> W = malloc(N * sizeof(double));
    eigen -> Z = malloc(N * N * sizeof(double));
    eigen -> ISUPPZ = malloc(2 * N * sizeof(int));
    eigen -> WORK = malloc(26 * N * sizeof(double));
    eigen -> IWORK = malloc(10 * N * sizeof(int));
    eigen -> N = N;
    // Used for auxilary work in Procrustes.
    eigen -> auxilary = malloc(N * N * sizeof(double));
    return eigen;
}

int compute_k_eigenpairs(RealSymEigen_t* eigen, int k) {
    // If we are computing all eigenpairs, indicate option 'A'.
    if (k == eigen -> N)
        return dsyevr('V', 'A', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
    // If we are computing only k eigenpairs, indicate option 'I'.
    return dsyevr('V', 'I', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, eigen -> N - k + 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
}

int compute_k_eigenvalues(RealSymEigen_t* eigen, int k) {
    // If we are computing all eigenvalues, indicate option 'A'.
    if (k == eigen -> N)
        return dsyevr('N', 'A', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
    // If we are computing only k eigenvalues, indicate option 'I'.
    return dsyevr('N', 'I', 'U', eigen -> N, eigen -> A, eigen -> N, 
                    0, 0, eigen -> N - k + 1, eigen -> N, dlamch('S'), 
                    &(eigen -> M), eigen -> W, eigen -> Z, eigen -> N, 
                    eigen -> ISUPPZ, eigen -> WORK, 26 * eigen -> N, 
                    eigen -> IWORK, 10 * eigen -> N);
}

// Free associated memory.
void destroy_real_sym_eigen(RealSymEigen_t* eigen) {
    if (eigen == NULL)
        return;
    free(eigen -> A);
    free(eigen -> W);
    free(eigen -> Z);
    free(eigen -> ISUPPZ);
    free(eigen -> WORK);
    free(eigen -> IWORK);
    free(eigen -> auxilary);
    free(eigen);
}