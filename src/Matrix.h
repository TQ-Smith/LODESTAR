
// File: Matrix.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define basic matrix management on generic types.

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>
#include <math.h>

// The number of elements in the upper triangle of a symmetric matrix.
#define PACKED_SIZE(N) ((N * (N + 1)) / 2)

// We store the upper triangle of a symmetric matrix in a one-dimensional array,
//  This is known as packed storage. Element Aij -> a_{i + j * (j + 1) / 2}
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

// Generate the create and destroy methods for a matrix of name NAME
//  and type TYPE. Note, the destroy method does not need to be generated,
//  but I kept it for completeness. In the future, we could supply a free method
//  to destroy matrices of pointers to structures.
#define MATRIX_INIT(NAME, TYPE) \
    static TYPE** create_##NAME##_matrix(int n, int m) { \
        TYPE** matrix = (TYPE**) calloc(n, sizeof(TYPE*)); \
        for (int i = 0; i < n; i++) \
            matrix[i] = (TYPE*) calloc(m, sizeof(TYPE)); \
        return matrix; \
    } \
    static void destroy_##NAME##_matrix(TYPE** matrix, int n) { \
        if (matrix == NULL) \
            return; \
        for (int i = 0; i < n; i++) \
            free(matrix[i]); \
        free(matrix); \
    }

// Wrapper methods for creating and destroying a matrix.
#define create_matrix(NAME, N, M) create_##NAME##_matrix(N, M)
#define destroy_matrix(NAME, MATRIX, N) destroy_##NAME##_matrix(MATRIX, N)

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
// Returns: double, tr(X^TX).
double normalize_matrix(double** X, int n, int k);

// Uncenter a matrix.
// Accepts:
//  double** X -> The matrix to uncenter.
//  double* x0 -> The vector to hold the column means.
//  int n -> The number of rows in X.
//  int k -> The number of columns in X.
// Returns: void.
void uncenter_matrix(double** X, double* x0, int n, int k);

#endif