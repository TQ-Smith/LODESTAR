
// File: MultidimensionalScaling.c
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Perform classical MDS on a distance matrix using LAPACK's dsyevr routine.

#include "MultidimensionalScaling.h"
#include <math.h>
#include <stdio.h>

// Index a packed stored matrix. If i > j, then treat the upper triangle as the lower triangle.
#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

// Used for comparing a floating point to 0.
#define EPS 1.0e-5

int compute_classical_mds(RealSymEigen_t* eigen, double* packedDistanceMatrix, int k, double** X) {

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
        return INFO;
    
    // Project points down into dimension k.
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < k; j++) {
            // If the eigenvalue is negative, we cannot take square-root. Return error.
            // If one of the eigenvalues are 0, then we project into a dimension less than k,
            //  return error.
            if (eigen -> W[k - j - 1] < 0 || fabs(eigen -> W[k - j - 1]) <= EPS)
                return 1;
            X[i][j] = eigen -> Z[(k - j - 1) * N + i] * sqrt(eigen -> W[k - j - 1]);
        }
    }

    return 0;
}

/*
int main () {
    int N = 3;
    int k = 2;
    RealSymEigen_t* eigen = init_real_sym_eigen(N);
    double* D = (double*) calloc(N * (N + 1) / 2, sizeof(double));
    double** X = (double**) calloc(N, sizeof(double*));
    for (int i = 0; i < N; i++)
        X[i] = (double*) calloc(k, sizeof(double));
    int a = 1;
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++)
            D[INDEX(i, j, N)] = a++;
        D[INDEX(i, i, N)] = 0;
    }
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%lf\t", D[INDEX(i, j, N)]);
        }
        printf("\n");
    }
    printf("\n");
    printf("\nINFO: %d\n", compute_classical_mds(eigen, D, k, X));
    for (int i = 0; i < N; i++) {
        for (int j = 0 ; j < k; j++) 
            printf("%lf\t", X[i][j]);
        printf("\n");
    }
    destroy_real_sym_eigen(eigen);
    for (int i = 0; i < N; i++)
        free(X[i]);
    free(X);
    free(D);
}
*/