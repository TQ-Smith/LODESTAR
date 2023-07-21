//
// File: MultidimensionalScaling.hpp
// Started: 21 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains methods for classic multidimensional scaling
//              and the FastMap heuristic.
//

#include "MultidimensionalScaling.hpp"

// Used to calculate eigenvectors.
#include "NumericalRecipesInC.hpp"

// Used for the square-root function.
#include <cmath>

void compute_classical_mds(double** D, double** X, double* d, double* e, bool* doesConverge, int n, int k) {

    // We double each element and keep track of each row's and
    //  each column's sum.

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
            D[i][j] = D[i][j] * D[i][j];
            d[i] += D[i][j];
            e[j] += D[i][j];
            grand_mean += D[i][j];
        }
    }

    // Double center X.
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            // Scale by a constant of -0.5
            D[i][j] = -0.5 * (D[i][j] - (d[i] / n) - (e[j] / n) + (grand_mean / (n * n)));
        }
    }

    // Now, calculate the eigenvetors and eigenvalues.
    compute_eigen_pairs(D, d, e, doesConverge, n, true, true);

    // If eigenpairs could not be computed, exit routine.
    if (!doesConverge) {
        return;
    }

    // Calculate projections in dimension k.
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < k; b++) {
            X[a][b] = 0;
        }
        for (int b = 0; b < k; b++) {
            X[a][b] = sqrt(d[b]) * D[a][b];
        }
    }

}