//
// File: MultidimensionalScaling.cpp
// Started: 15 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains functions to perform various forms of MDS.
//

// Used for basic matrix operations.
#include "MatrixOperations.h"

// Used to calculate eigenvectors.
#include "NumericalRecipesInC.h"

// Used for the square-root function.
#include <cmath>

#include <iostream>
using namespace std;

void compute_classical_mds(double** X, double* d, double* e, int n, int k) {

    // We double each element and keep track of each row's and
    //  each column's sum. Then, we double center.

    // Note, d and e are used to hold row and column sums before eigen computations.

    // Holds total of each row and each column.
    for(int i = 0; i < n; i++) {
        d[i] = 0;
        e[i] = 0;
    }

    // Holds mean of all of the elements.
    double grand_mean = 0;

    // Double each element and compute totals.
    double twice;
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            twice = X[i][j] * X[i][j];
            d[i] += twice;
            e[j] += twice;
            grand_mean += twice;
        }
    }

    // Double center X.
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            X[i][j] = X[i][j] - (d[i] / n) - (e[j] / n) + (grand_mean / (n * n));
        }
    }

    // Now, calculate the eigenvetors and eigenvalues.
    tred2(X, d, e, n, true);
    tqli(X, d, e, n, true);
    eigsrt(X, d, n);

    // Allocate the matrix to hold the reduced data.
    double** reduced_matrix = create_real_matrix(n, k);

    // Calculate projections in dimension k.
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < k; b++) {
            reduced_matrix[a][b] = 0;
            for (int c = 0; c < n; c++) {
                reduced_matrix[a][b] += X[c][b] * sqrt(d[c]);
            }
        }
    }

    // Reassign X to the reduced_matrix and free eigenvectors matrix.
    double*** temp = &X;
    X = reduced_matrix;
    reduced_matrix = *temp;
    destroy_real_matrix(reduced_matrix, n, n);
    
}

void classical_mds(double** X, int n, int k) {

    // Allocate extra memory.
    double* d = new double[n];
    double* e = new double[n];

    // Compute classical MDS.
    compute_classical_mds(X, d, e, n, k);

    // Free temporary memory.
    delete [] d;
    delete [] e;

}

int main() {

    int n = 7;
    int k = 2;

    double** X = create_real_matrix(n, n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (i == j)
                X[i][j] = 1;
            X[i][j] = (i + 1) * (j + 1);
        }
    }

    print_real_matrix(X, n, n);

    classical_mds(X, n, k);

    print_real_matrix(X, n, k);

}