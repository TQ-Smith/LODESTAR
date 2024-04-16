
#include "ProcrustesAnalysis.h"

#include <math.h>

#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

#include <stdio.h>

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen* eigen, int N, int K, bool transform, bool similarity) {

    double trX = 0, trY = 0;

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < N; j++) {
            trX += Xc[j][i] * Xc[j][i];
            trY += Yc[j][i] * Yc[j][i]; 
        }
    }

    double* C = eigen -> Z;
    double* covC = eigen -> A;

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            C[INDEX(i, j, K)] = 0;
            for (int k = 0; k < N; k++)
                C[INDEX(i, j, K)] += Yc[k][i] * Xc[k][j];
        }
    }

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            covC[INDEX(i, j, K)] = 0;
            for (int k = 0; k < K; k++)
                covC[INDEX(i, j, K)] += C[INDEX(k, i, K)] * C[INDEX(k, j, K)];
        }
    }

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            printf("%5f\t", covC[INDEX(i, j, K)]);
        }
        printf("\n");
    }

    return 0;

}