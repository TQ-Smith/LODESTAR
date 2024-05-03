
// File: RealSymEigen.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute the eigenpairs of a real symmetric matrix.
// Reference: Jochen Voss's matroot.c
//            https://github.com/Foadsf/Cmathtuts/blob/a3610505e04c162c48ad8b8eab9745b32f2f0b6d/C_CBLAS/seehuhn/example2/matroot.c#L61
//
//            We are using LAPACK's dsyevr routine.
//            https://netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga2ad9f4a91cddbf67fe41b621bd158f5c.html

#ifndef _REAL_SYM_EIGEN_H_
#define _REAL_SYM_EIGEN_H_

#include <stdlib.h>

// A structure to hold all the memory used by dsyevr.
typedef struct RealSymEigen {
    double *A, *B, *W, *Z, *WORK;
    int *ISUPPZ, *IWORK;
    int N, M;
} RealSymEigen_t;

// Allocate the necessary memory to perform eigen computations.
// Accepts:
//  int N -> The dimension of the square matrix.
// Returns: RealSymEigen_t*, The created structure.
RealSymEigen_t* init_real_sym_eigen(int N);

// Computes the most significant k eigenpairs of the matrix stored in A.
// Accepts:
//  RealSymEigen_t* eigen -> Our structure with the necessary memory used by dsyevr.
//  int k -> The number of eigenpairs to compute.
// Returns: int, LAPACK success/failure code. 0 for success.
int compute_k_eigenpairs(RealSymEigen_t* eigen, int k);

// Computes the most significant k eigenvalues of the matrix stored in A.
// Accepts:
//  RealSymEigen_t* eigen -> Our structure with the necessary memory used by dsyevr.
//  int k -> The number of eigenvalues to compute.
// Returns: int, LAPACK success/failure code. 0 for success.
int compute_k_eigenvalues(RealSymEigen_t* eigen, int k);

// Destroy structure.
// Accepts:
//  RealSymEigen_t* eigen -> The structure to destroy.
// Returns: void.
void destroy_real_sym_eigen(RealSymEigen_t* eigen);

#endif