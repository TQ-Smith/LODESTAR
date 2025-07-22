
// File: MatrixOperations.h
// Date: 7 June 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Operations on real valued matrices.
// Reference: Jochen Voss's matroot.c
//            https://github.com/Foadsf/Cmathtuts/blob/a3610505e04c162c48ad8b8eab9745b32f2f0b6d/C_CBLAS/seehuhn/example2/matroot.c#L61
//
//            We are using LAPACK's dsyevr routine.
//            https://netlib.org/lapack/explore-html-3.6.1/d2/d8a/group__double_s_yeigen_ga2ad9f4a91cddbf67fe41b621bd158f5c.html

#ifndef _MATRIX_OPERATIONS_H_
#define _MATRIX_OPERATIONS_H_

#include <stdlib.h>

// Index a packed stored matrix. If i > j, then treat the upper triangle as the lower triangle.
#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

// Fortran matrices are stored as vectors in column-major form.
//  Covariance matrices are symmetric so this is not explicitly necessary,
//  but it is good to keep things consistent.
#define COL_MAJOR(i, j, K) (j * K + i)

// Used for comparing a floating point to 0.
#define EPS 1.0e-8

// Create an n-by-k matrix.
double** init_matrix(int n, int k);

// Destroy an n-by-k matrix.
void destroy_matrix(double** X, int n);

// Center a matrix for use in MDS.
// Accepts:
//  double** X -> The matrix to center.
//  double* x0 -> The vector to hold the column means.
//  int n -> The number of rows in X.
//  int k -> The number of columns in X.
// Returns: void.
void center_matrix(double** X, double* x0, int n, int k);

// Normalize a matrix by sqrt(trX^TX).
// Accepts:
//  double** X -> The matrix to normalize.
//  int n -> The number of rows in X.
//  int k -> The number of columns in X.
// Returns: int, 0 for success, -1 if X contains only one point (trX = 0).
int normalize_matrix(double** X, int n, int k);

// A structure to hold all the memory used by dsyevr.
typedef struct {
    double *A, *B, *W, *Z, *WORK, *auxilary;
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

// Perfrom classical MDS on a distance matrix.
// Accepts:
//  RealSymEigen_t* eigen -> Structure with memory to execute LAPACK's dsyevr routine.
//  double* packedDistanceMatrix -> The upper triangle of the distance matrix in packed storage.
//                                  Not perserved.
//  int k -> The dimension to project down into.
//  double** X -> An allocated N-by-k matrix to store the resulting points.
// Returns: double, the effective rank of the matrix. Otherwise, -1, LAPACK failure, k eigenvalues are not all positive,
//              or packedDistanceMatrix contains nan.
double compute_classical_mds(RealSymEigen_t* eigen, double* packedDistanceMatrix, int k, double** X);

// Compute Procrustes statistic between two sets of points.
// Accepts:
//  double** Xc -> The n-by-k mean-centered query point set.
//  double* x0 -> The k-dimensional vector holding the mean of each dimension.
//                  If NULL, x0 is assumed to be the 0 vector.
//  double** Yc -> The n-by-k mean centered target point set.
//  double* y0 -> The k-dimensional vector holding the mean of each dimension.
//                  If NULL, x0 is assumed to be the 0 vector.
//  RealSymEigen_t* eigen -> Used to find eigen-pairs of k-by-k matrix. Allocated
//                  memory is used to store covariance matrix even if Xc is not to be
//                  transformed.
//  int N -> The number of points.
//  int K -> The dimension of each point.
//  bool transform -> If set, x0 will be transformed by Procrustes analysis.
// Returns: double, The Procrustes statistic.
double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen_t* eigen, int N, int K, bool transform);

#endif