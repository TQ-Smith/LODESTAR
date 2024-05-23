
// File: Matrix.c
// Date: 17 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Extends header by define basic matrix operations.

#include "Matrix.h"

void center_matrix(double** X, double* x0, int n, int k) {
    // Calculate the center.
    for (int i = 0; i < k; i++) {
        x0[i] = 0;
        for (int j = 0; j < n; j++)
            x0[i] += (X[j][i] / n);
    }
    // Center the set of points.
    for (int i = 0; i < n; i++)
        for (int j = 0; j < k; j++) 
            X[i][j] -= x0[j];
}

void uncenter_matrix(double** X, double* x0, int n, int k) {
    for (int i = 0; i < n; i++) 
        for (int j = 0; j < k; j++)
            X[i][j] += x0[j];
}