//
// File: unit_test_matrix_operations.cpp
// Started: 19 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Unit tests the matrix operations.
//

#include "../utils/LinearAlgebra/MatrixOperations.hpp"

#include <iostream>

using namespace std;

int main() {

    // Set random seed for shuffling.
    srand(time(NULL));

    // Test the creation, filling, and deep copying of matrices.
    double** a = create_and_fill_real_matrix(2, 2, 2);
    double** b = create_real_matrix(2, 2);
    fill_real_matrix(b, 3, 2, 2);
    double** c = deep_copy_real_matrix(b, 2, 2);
    double** d = create_and_fill_real_matrix(4, 2, 4);

    // Test the printing of the matrices.
    cout << "Matrix A = " << endl;
    print_real_matrix(a, 2, 2, 1, 3);
    cout << endl;
    cout << "Matrix B = " << endl;
    print_real_matrix(b, 2, 2, 1, 3);
    cout << endl;
    cout << "Matrix C = " << endl;
    print_real_matrix(c, 2, 2, 1, 3);
    cout << endl;
    cout << "Matrix D = " << endl;
    print_real_matrix(d, 2, 4, 1, 3);
    cout << endl;

    // Test addition.
    cout << "A + B = C" << endl;
    cout << "C =" << endl;
    add_matrices(c, a, b, 2, 2);
    print_real_matrix(c, 2, 2, 1, 3);
    cout << endl;

    // Test subtraction.
    cout << "C - A = B" << endl;
    cout << "B =" << endl;
    subtract_matrices(b, c, a, 2, 2);
    print_real_matrix(b, 2, 2, 1, 3);
    cout << endl;

    // Test scaling.
    cout << "2 * A" << endl;
    cout << "A =" << endl;
    scale_matrix(2, a, 2, 2);
    print_real_matrix(a, 2, 2, 1, 3);
    cout << endl;

    // Test multiplication.
    cout << "A * D = " << endl;
    double** e = create_real_matrix(2, 4);
    multiply_matrices(e, a, d, 2, 2, 4);
    print_real_matrix(e, 2, 4, 1, 3);
    cout << endl;

    double** f = create_real_matrix(3, 2);
    f[0][0] = f[0][1] = 0;
    f[1][0] = f[1][1] = 1;
    f[2][0] = f[2][1] = 2;
    cout << "F = " << endl;
    print_real_matrix(f, 3, 2, 1, 3);
    cout << endl;

    // Test shuffling of rows.
    cout << "Shuffle F =" << endl;
    shuffle_real_matrix(f, 3);
    print_real_matrix(f, 3, 2, 1, 3);
    cout << endl;

    // Test deallocation.
    destroy_real_matrix(a, 2);
    destroy_real_matrix(b, 2);
    destroy_real_matrix(c, 2);
    destroy_real_matrix(d, 2);
    destroy_real_matrix(e, 2);
    destroy_real_matrix(f, 3);

}