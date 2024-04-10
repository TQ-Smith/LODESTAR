#include "MultidimensionalScaling.h"

#include <math.h>

#include <stdio.h>

#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

int compute_classical_mds(RealSymEigen* eigen, double* packedDistanceMatrix, int k, double** X) {

    int N = eigen -> N;

    for (int i = 0; i < N; i++)
        eigen -> WORK[i] = 0;

    double grand_mean = 0;

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if (i <= j)
                packedDistanceMatrix[INDEX(i, j, N)] *= packedDistanceMatrix[INDEX(i, j, N)];
            eigen -> WORK[i] += packedDistanceMatrix[INDEX(i, j, N)];
            grand_mean += packedDistanceMatrix[INDEX(i, j, N)];
        }
    }

    for(int i = 0; i < N; i++)
        for(int j = 0; j < N; j++)
            eigen -> A[j * N + i] = -0.5 * (packedDistanceMatrix[INDEX(i, j, N)] - (eigen -> WORK[i] / N) - (eigen -> WORK[j] / N) + (grand_mean / (N * N)));
    
    
    int INFO = compute_k_eigenpairs(eigen, k);

    if (INFO != 0)
        return INFO;
    
    for (int i = 0; i < N; i++)
        for (int j = 0; j < k; j++)
            X[i][j] = eigen -> Z[(k - j - 1) * N + i] * sqrt(eigen -> W[k - j - 1]);

    return 0;

}