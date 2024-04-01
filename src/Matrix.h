
#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>

#include <stdlib.h>

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

#endif