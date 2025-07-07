
// File: MatrixOperations.c
// Date: 7 June 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Operations on real valued matrices.
// Reference: Jochen Voss's matroot.c
//            https://github.com/Foadsf/Cmathtuts/blob/a3610505e04c162c48ad8b8eab9745b32f2f0b6d/C_CBLAS/seehuhn/example2/matroot.c#L61
//
//            We are using LAPACK's dsyevr routine.
//            https://netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga2ad9f4a91cddbf67fe41b621bd158f5c.html

#include "MatrixOperations.h"
#include <math.h>

void center_matrix(double** X, double* x0, int n, int k) {
    // Calculate the center.
    for (int i = 0; i < k; i++) {
        x0[i] = 0;
        for (int j = 0; j < n; j++)
            x0[i] += (X[j][i] / n);
    }
    // Center the set of points.
    for (int i = 0; i < n; i++)
        for (int j = 0; j < k; j++) 
            X[i][j] -= x0[j];
}

void normalize_matrix(double** X, int n, int k) {
    // trx is equivalent to trX^TX.
    double trX = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            trX += X[i][j] * X[i][j];
        }
    }
    if (fabs(trX) <= EPS) {
        return 0;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] = X[i][j] / sqrt(trX);
        }
    }
}

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

double compute_classical_mds(RealSymEigen_t* eigen, double* packedDistanceMatrix, int k, double** X) {

    int N = eigen -> N;

    // We use the WORK array to hold row sums.
    for (int i = 0; i < N; i++)
        eigen -> WORK[i] = 0;

    double grand_mean = 0;

    // Square each element of the distance matrix. Compute grand mean and row sums.
    for(int i = 0; i < N; i++) {
        for(int j = i + 1; j < N; j++) {
            // Distance matrices cannot contain negative values.
            if (packedDistanceMatrix[INDEX(i, j, N)] < 0)
                return 1;
            // Square each element.
            packedDistanceMatrix[INDEX(i, j, N)] *= packedDistanceMatrix[INDEX(i, j, N)];
            // Calculate row sums and grand mean.
            eigen -> WORK[i] += packedDistanceMatrix[INDEX(i, j, N)];
            eigen -> WORK[j] += packedDistanceMatrix[INDEX(j, i, N)];
            grand_mean += packedDistanceMatrix[INDEX(i, j, N)];
            grand_mean += packedDistanceMatrix[INDEX(j, i, N)];
        }
        packedDistanceMatrix[INDEX(i, i, N)] *= packedDistanceMatrix[INDEX(i, i, N)];
        eigen -> WORK[i] += packedDistanceMatrix[INDEX(i, i, N)];
        grand_mean += packedDistanceMatrix[INDEX(i, i, N)];
    }

    // Double center matrix and store column major in A (not necessary but it's good to keep it consistent).
    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            eigen -> A[j * N + i] = -0.5 * (packedDistanceMatrix[INDEX(i, j, N)] - (eigen -> WORK[i] / N) - (eigen -> WORK[j] / N) + (grand_mean / (N * N)));
    
    // Compute eigen pairs.
    int INFO = compute_k_eigenpairs(eigen, k);

    // If error code from LAPACK, return code.
    if (INFO != 0)
        return -1;
    
    // Project points down into dimension k.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            // If the eigenvalue is negative, we cannot take square-root. Return error.
            // If one of the eigenvalues are 0, then we project into a dimension less than k,
            //  return error.
            if (eigen -> W[k - j - 1] < 0 || fabs(eigen -> W[k - j - 1]) <= EPS)
                return -1;
            X[i][j] = eigen -> Z[(k - j - 1) * N + i] * sqrt(eigen -> W[k - j - 1]);
        }
    }

    return 0;
}

// A few things to note:
//  - LAPACK requires gfortran which is not a part of CLANG.
//  - Must be compiled with gcc/gfortran.
//  - All matrices in LAPACK are one-dimensional arrays stored in column major
//      form. This does not matter with symmetric matrices.
//  - A is the matrrix to decompose. 
//  - S will contain the singular values.
//  - U will contain the U matrix.
//  - VT will contain V^T. 

// Define out LAPACK SVD routine.
extern void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A,
                int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT,
                double* WORK, int* LWORK, int* INFO);

// We wrap the methods in a function to make them more readable in C.
int dgesvd(
    char JOBU, 
    char JOBVT, 
    int M, 
    int N, 
    double* A,
    int LDA, 
    double* S, 
    double* U, 
    int LDU, 
    double* VT, 
    int LDVT,
    double* WORK, 
    int LWORK 
) {
    int INFO;

    dgesvd_(&JOBU, &JOBVT, &M, &N, A,
            &LDA, S, U, &LDU, VT, &LDVT,
            WORK, &LWORK, &INFO);
    
    return INFO;
}