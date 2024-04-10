
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>

#define PACKED_SIZE(N) ((N * (N + 1)) / 2)
#define PACKED_INDEX(i, j) (i + j * (j + 1) / 2)

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

#define create_matrix(NAME, N, M) create_##NAME##_matrix(N, M)
#define destroy_matrix(NAME, MATRIX, N) destroy_##NAME##_matrix(MATRIX, N)

#endif