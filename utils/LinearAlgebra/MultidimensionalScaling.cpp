//
// File: MultidimensionalScaling.hpp
// Date: 21 July 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Contains methods for classic multidimensional scaling
//              and the FastMap heuristic.
//

// NOTE: FastMap functions are adapted from Faloutsos and Lin's implementation found here:
// https://www.cs.cmu.edu/~christos/software.html#:~:text=Dimensionality%20reduction-,FastMap%20tar%20file,-In%20C.%20Paper

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
    if (!(*doesConverge)) {
        return;
    }

    // Calculate projections in dimension k.
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < k; b++) {
            X[a][b] = D[a][b] * sqrt(d[b]);
        }
    }

}

// A helper function for the FastMap algorithm.
//  Calculates the distance between two objects
//  using Equation 4.
// Accepts:
//  double** D -> Our n x n distance matrix.
//  double** X -> An n x k matrix to hold the points in dimension k.
//  int object_a -> The index of the first object.
//  int object_b -> The index of the second object. 
//  int k -> The current projection-dimension.
//      Assume 0 < k < n.
// Returns: double, the distance of the projected objects.
double projection_distance(double** D, double** X, int object_a, int object_b, int k) {
    double temp;
    // Get the distance between the objects.
    double dim0 = D[object_a][object_b];
    // Decrease the distance by the desired number of projections.
    for(int i = 0; i < k; i++) {
        temp = X[object_a][i] - X[object_b][i];
        dim0 -= temp * temp;
    }
    return dim0;
}

// A helper function for the FastMap algorithm.
//  Given an object, return the index of the farthest one.
// Accepts:
//  double** D -> Our n xn distance matrix.
//  double** X -> An n x k matrix to hold the points in dimension k.
//  int object_a -> The index of the initial object.
//  int n -> The dimension of the data.
//  int k -> The current projection-dimension.
//      Assume 0 < k < n.
int furthest(double** D, double** X, int object_a, int n, int k) {
    double temp;
    // The max distance.
    double max = 0;
    // The furthest object.
    int object_b = object_a;
    // Find the max.
    for(int i = 0; i < n; i++) {
        temp = projection_distance(D, X, object_a, i, k);
        if (temp > max) {
            object_b = i;
            max = temp;
        }
    }
    // Return the index of the max element.
    return object_b;
}

void compute_fastmap(double** D, double** X, int* maxDimReached, int n, int k) {

    // The indices of the objects that are furthest 
    //  away from each other.
    int object_a;
    int object_b;

    // Values for distances between objects.
    double distsq_objects_ab, distsq_objects_aj, distsq_objects_bj;

    // For each projection dimension.
    for (int i = 0; i < k; i++) {

        // We start off our search for distant objects with the object at index 0.
        object_a = 0;
        // Pick two distance object.
        object_b = furthest(D, X, object_a, n, k);
        object_a = furthest(D, X, object_b, n, k);

        // If we are still trying to do projections and no objects are
        //  found to be dissimilar, then we set the maximum number of 
        //  projections reached, and we exit the routine.
        distsq_objects_ab = projection_distance(D, X, object_a, object_b, i);
        if (distsq_objects_ab == 0) {
            *maxDimReached = i + 1;
            return;
        }

        // Project the objects on the line.
        for (int j = 0; j < n; j++) {
            distsq_objects_aj = projection_distance(D, X, object_a, j, i);
            distsq_objects_bj = projection_distance(D, X, object_b, j, i);
            X[j][i] = (distsq_objects_aj + distsq_objects_ab - distsq_objects_bj) / (2.0 * sqrt((double) distsq_objects_ab));
        }
    }

    // All k projects have been satisfied.
    *maxDimReached = k;

}