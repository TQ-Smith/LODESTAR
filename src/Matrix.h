
// File: Matrix.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define basic matrix management on generic types.

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>

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
    TYPE** create_##NAME##_matrix(int n, int m) { \
        TYPE** matrix = (TYPE**) calloc(n, sizeof(TYPE*)); \
        for (int i = 0; i < n; i++) \
            matrix[i] = (TYPE*) calloc(m, sizeof(TYPE)); \
        return matrix; \
    } \
    void destroy_##NAME##_matrix(TYPE** matrix, int n) { \
        if (matrix == NULL) \
            return; \
        for (int i = 0; i < n; i++) \
            free(matrix[i]); \
        free(matrix); \
    }

// Wrapper methods for creating and destroying a matrix.
#define create_matrix(NAME, N, M) create_##NAME##_matrix(N, M)
#define destroy_matrix(NAME, MATRIX, N) destroy_##NAME##_matrix(MATRIX, N)

#endif