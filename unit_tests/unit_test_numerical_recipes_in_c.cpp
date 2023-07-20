//
// File: unit_test_numerical_recipes_in_c.cpp
// Started: 20 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Unit tests the eigenvalue and eigenvector computations
//              from the functions in Numerical Recipes in C, Third Edition.
//              Note, this is not for testing the correctness of their algorithm;
//              it is for making sure that I know how to use the functions and that
//              my slight modifications did not break anything.
//

// Used to create and print matrices.
#include "../utils/LinearAlgebra/MatrixOperations.hpp"

#include "../utils/LinearAlgebra/NumericalRecipesInC.hpp"

#include <iostream>

using namespace std;

int main() {

    // The dimension of the matrices.
    int n = 3;

    // Create a matrix with three distinct eigenvectors.
    double** a = create_real_matrix(3, 3);
    a[0][0] = a[1][1] = a[2][2] = 1;
    a[0][1] = a[1][2] = a[1][0] = a[2][1] = 2;
    a[0][2] = a[2][0] = 3;

    // Create a matrix where the tqli routine does not converge.
    double** b = create_real_matrix(3, 3);
    b[0][0] = 4.0/5;
    b[1][1] = 4.0/5;
    b[2][2] = 2;
    b[0][1] = -3.0/5;
    b[1][2] = 0;
    b[1][0] = 3.0/5;
    b[2][1] = 2;
    b[0][2] = 0;
    b[2][0] = 1;

    // Holds eigenvalues.
    double* d = new double[n];

    // Holds off diagonal entries.
    double* e = new double[n];

    cout << "A =" << endl;
    print_real_matrix(a, 3, 3, 1, 1);
    cout << endl;

    // Flag to test if tqli converges.
    bool doesConverge;

    // Perform eigen calculations.
    compute_eigen_pairs(a, d, e, &doesConverge, n, true, true);

    cout << "Eigenvalues of A =" << endl;
    print_real_array(d, n, 1, 1);
    cout << endl;

    cout << "The eigenvectors of A=" << endl;
    print_real_matrix(a, 3, 3, 1, 1);
    cout << endl;

    assert(doesConverge);

    cout << "B =" << endl;
    print_real_matrix(b, 3, 3, 1, 1);
    cout << endl;

    // Perform eigen calculations.
    compute_eigen_pairs(b, d, e, &doesConverge, n, true, true);

    cout << "Eigenvalues of B =" << endl;
    print_real_array(d, n, 1, 1);
    cout << endl;

    cout << "The eigenvectors of B=" << endl;
    print_real_matrix(b, 3, 3, 1, 1);
    cout << endl;

    assert(!doesConverge);

    destroy_real_matrix(a, 3);
    destroy_real_matrix(b, 3);
    delete [] d;
    delete [] e;

}