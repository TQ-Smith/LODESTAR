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

double** compute_classical_mds(double** X, double* d, double* e, int n, int k) {

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
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            X[i][j] = X[i][j] * X[i][j];
            d[i] += X[i][j];
            e[j] += X[i][j];
            grand_mean += X[i][j];
        }
    }

    // Double center X.
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            // Scale by a constant of -0.5
            X[i][j] = -0.5 * (X[i][j] - (d[i] / n) - (e[j] / n) + (grand_mean / (n * n)));
        }
    }

    // Now, calculate the eigenvetors and eigenvalues.
    compute_eigen_pairs(X, d, e, n, true, true);

    // Allocate the matrix to hold the reduced data.
    double** reduced_matrix = create_real_matrix(n, k);

    // Calculate projections in dimension k.
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < k; b++) {
            reduced_matrix[a][b] = 0;
        }
        for (int b = 0; b < k; b++) {
            reduced_matrix[a][b] = sqrt(d[b]) * X[a][b];
        }
    }

    // Destroy X and return the reduced dataset.
    destroy_real_matrix(X, n, n);
    return reduced_matrix;
    
}

double** classical_mds(double** X, int n, int k) {

    // Allocate extra memory.
    double* d = new double[n];
    double* e = new double[n];

    // Compute classical MDS.
    double** reduced = compute_classical_mds(X, d, e, n, k);

    // Free temporary memory.
    delete [] d;
    delete [] e;

    return reduced;

}

int main() {

    int n = 7;
    int k = 2;

    double** X = create_real_matrix(n, n);
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            if (i == j)
                X[i][j] = 1;
            else
                X[i][j] = (i + 1) * (j + 1);
        }
    }

    X = classical_mds(X, n, k);

    print_real_matrix(X, n, k);

}