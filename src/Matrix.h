
#ifndef _MATRIX_
#define _MATRIX_

#include <stdio.h>

#include <stdlib.h>

#define MATRIX_INIT(TYPE, PRINT) \
    TYPE** create_##TYPE##_matrix(int n, int m) { \
        TYPE** matrix = (TYPE**) calloc(n, sizeof(TYPE*)); \
        for (int i = 0; i < n; i++) \
            matrix[i] = (TYPE*) calloc(m, sizeof(TYPE)); \
        return matrix; \
    } \
    void print_##TYPE##_matrix(TYPE** matrix, int n, int m) { \
        for (int i = 0; i < n; i++) { \
            for (int j = 0; j < m; j++) { \
                PRINT(matrix[i][j]); \
                printf("\t"); \
            } \
            printf("\n"); \
        } \
    } \
    void destroy_##TYPE##_matrix(TYPE** matrix, int n) { \
        if (matrix == NULL) \
            return; \
        for (int i = 0; i < n; i++) \
            free(matrix[i]); \
        free(matrix); \
    }

#endif