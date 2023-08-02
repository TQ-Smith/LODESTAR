//
// File: Procrustes.cpp
// Date: 27 July 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Unit tests Procrustes analysis and permutation test.
//

#include "../src/Procrustes.hpp"

// Used for creating and deleting matrices.
#include "../utils/LinearAlgebra/MatrixOperations.hpp"

// Used for absolute value.
#include <cmath>

// Used for testing results.
#include <cassert>

#include <iostream>
using namespace std;

// A helper method to add back the centers so we can reuse X.
// Accepts:
//  double** Xc -> The matrix to uncenter.
//  double* x_0 -> The centers for each dimension.
//  int n -> The number of points.
//  int k -> The number of dimensions.
void uncenter(double** Xc, double* x_0, int n , int k) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            // Add back center.
            Xc[i][j] += x_0[j];
        }
    }
}

int main() {

    // Our tests will consist of 3 points with k = 1, 2, 3.
    int n = 3;

    // Allocates our two matrices used for the Procrustes
    //  analysis. We vary k = 1, 2, 3 which denotes the number
    //  of columns we are considering.
    double** X = create_real_matrix(n, n);
    X[0][0] = 1; X[0][1] = 1; X[0][2] = 1;
    X[1][0] = 2; X[1][1] = 4; X[1][2] = 8;
    X[2][0] = 3; X[2][1] = 9; X[2][2] = 27;

    double** Y = create_real_matrix(n, n);
    Y[0][0] = 5; Y[0][1] = 13; Y[0][2] = 21;
    Y[1][0] = 7; Y[1][1] = 17; Y[1][2] = 23;
    Y[2][0] = 11; Y[2][1] = 19; Y[2][2] = 27;

    // Additional memory needed.
    double* x_0 = new double[n];
    double* y_0 = new double[n];
    double** C = create_real_matrix(n, n);
    double** CT_C = create_real_matrix(n, n);

    // Holds the Procrustes statistic.
    double D;

    // Note: All of the asserted values where hand computed.

    // Test when k = 1
    cout << "Testing k = 1" << endl;
    cout << endl;

    int k = 1;

    // Print the X matrix.
    cout << "Our initial matrix X =" << endl;
    print_real_matrix(X, n, k, 1, 1);
    cout << endl;

    // Print the Y matrix.
    cout << "Our initial matrix Y =" << endl;
    print_real_matrix(Y, n, k, 1, 1);
    cout << endl;

    // Print the centered X matrix.
    cout << "Centered matrix X =" << endl;
    center_matrix(X, x_0, n, k);
    print_real_matrix(X, n, k, 1, 4);
    cout << endl;

    // Print the centered Y matrix.
    cout << "Centered matrix Y =" << endl;
    center_matrix(Y, y_0, n, k);
    print_real_matrix(Y, n, k, 1, 4);
    cout << endl;

    // Compute the procrustes statistic.
    D = procrustes_analysis(X, Y, C, CT_C, n, k);
    assert(abs(0.964285 - D) < 0.00001);
    cout << "Procrustes Statistic = " << D << endl;

    // Uncenter the matrices.
    uncenter(X, x_0, n, k);
    uncenter(Y, y_0, n, k);

    cout << endl;

    // Test when k = 2
    cout << "Testing k = 2" << endl;
    cout << endl;

    k = 2;

    cout << "Our initial matrix X =" << endl;
    print_real_matrix(X, n, k, 1, 1);
    cout << endl;

    cout << "Our initial matrix Y =" << endl;
    print_real_matrix(Y, n, k, 1, 1);
    cout << endl;

    cout << "Centered matrix X =" << endl;
    center_matrix(X, x_0, n, k);
    print_real_matrix(X, n, k, 1, 4);
    cout << endl;

    cout << "Centered matrix Y =" << endl;
    center_matrix(Y, y_0, n, k);
    print_real_matrix(Y, n, k, 1, 4);
    cout << endl;

    D = procrustes_analysis(X, Y, C, CT_C, n, k);
    assert(abs(0.958791 - D) < 0.00001);
    cout << "Procrustes Statistic = " << D << endl;

    uncenter(X, x_0, n, k);
    uncenter(Y, y_0, n, k);

    cout << endl;

    // Test when k = 3
    cout << "Testing k = 3" << endl;
    cout << endl;

    k = 3;

    cout << "Our initial matrix X =" << endl;
    print_real_matrix(X, n, k, 1, 1);
    cout << endl;

    cout << "Our initial matrix Y =" << endl;
    print_real_matrix(Y, n, k, 1, 1);
    cout << endl;

    cout << "Centered matrix X =" << endl;
    center_matrix(X, x_0, n, k);
    print_real_matrix(X, n, k, 1, 4);
    cout << endl;

    cout << "Centered matrix Y =" << endl;
    center_matrix(Y, y_0, n, k);
    print_real_matrix(Y, n, k, 1, 4);
    cout << endl;

    D = procrustes_analysis(X, Y, C, CT_C, n, k);
    assert(abs(0.936014 - D) < 0.00001);
    cout << "Procrustes Statistic = " << D << endl;

    uncenter(X, x_0, n, k);
    uncenter(Y, y_0, n, k);

    cout << endl;

    cout << "Lastly, we check the Permutation test procedure." << endl;
    cout << "We use the case when k = 3 and its corresponding statistic." << endl;
    cout << endl;

    // Seed the random number generator with a constant value
    //  to ensure the same sequence each time.
    srand(0);

    // Note, uncomment the print statistic line in Procrustes.cpp
    //  That way, you can see the statistic for each permutation.

    // Create the memory for the copied matrix.
    double** shuffleX = create_real_matrix(n, n);

    center_matrix(X, x_0, n, k);
    center_matrix(Y, y_0, n, k);

    // Execute permutation test.
    double p = permutation_test(3, X, Y, shuffleX, C, CT_C, n, k, D);

    cout << "The p-value using 3 permutations is " << p << endl;
    assert(p == 0.25);

    cout << endl;

    // Destroy all used memory.
    destroy_real_matrix(X, n);
    destroy_real_matrix(shuffleX, n);
    destroy_real_matrix(Y, n);
    destroy_real_matrix(C, k);
    destroy_real_matrix(CT_C, k);
    delete [] x_0;
    delete [] y_0;

}