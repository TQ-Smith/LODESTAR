//
// File: MatrixOperations.hpp
// Started: 19 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains basic operations for manipulating m x n matrices
//              of doubles. Note, we do not account for null pointers,
//              values of m and n, and jagged arrays. We do this for two
//              reasons: first, we would have to abstract the double** to 
//              a matrix structure holding additional information, and second,
//              there would be much more boiler plate.
//

#ifndef _MATRIX_OPERATIONS_HPP_  
#define _MATRIX_OPERATIONS_HPP_

// Used for printing.
#include <iostream>
using namespace std;

// Create a matrix of doubles.
// Accepts:
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, the allocated matrix.
double** create_real_matrix(int m, int n) {
    // Allocate columns.
    double** matrix = new double*[m];
    // Allocate rows.
    for (int i = 0; i < m; i++) {
        matrix[i] = new double[n];
    }
    return matrix;
}

// Create a matrix of doubles.
// Accepts:
//  int value -> The value that fills the matrix.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, the allocated matrix.
double** create_and_fill_real_matrix(int value, int m, int n) {
    double** matrix = new double*[m];
    for (int i = 0; i < m; i++) {
        matrix[i] = new double[n];
        // Fill each row.
        for (int j = 0; j < n; j++) {
            matrix[i][j] = value;
        }
    }
    return matrix;
}

// Adds matrices b and c and stores in a.
// Accepts:
//  double** a -> The resulting matrix.
//  double** b -> The first operand.
//  double** c -> The second operand.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void add_matrices(double** a, double** b, double** c, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            // Add elements.
            a[i][j] = b[i][j] + c[i][j];
        }
    }
}

// Subtracts matrix c from b and stores in a.
// Accepts:
//  double** a -> The resulting matrix.
//  double** b -> The first operand.
//  double** c -> The second operand.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void subtract_matrices(double** a, double** b, double** c, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            // Subtract elements.
            a[i][j] = b[i][j] - c[i][j];
        }
    }
}

// Scale a matrix by a constant.
// Accepts:
//  double c -> The constant.
//  double** a -> The matrix.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void scale_matrix(double c, double** a, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            // Scale the element.
            a[i][j] *= c * a[i][j];
        }
    }
}

// Matrix multiplication. a = b * c.
// Accepts:
//  double** a -> The resulting matrix.
//  double** b -> The first operand.
//  double** c -> The second operand.
//  int m -> The length of the first matrix.
//  int n -> The width of the first matrix.
//  int o -> The width of the second matrix.
// Returns: void.
void multiply_matrices(double** a, double** b, double** c, int m, int n, int o) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < o; j++) {
            a[i][j] = 0;
            for (int k = 0; k < n; k++) {
                a[i][j] += b[i][k] * c[k][j];
            }
        }
    }
}

// Shuffles the rows of a matrix.
//  Used for permutation testing.
// Note: Assumes random number generator has been seeded. 
// Accepts:
//  double** matrix -> The matrix to shuffle.
//  int m -> The length of the matrix.
// Returns: void.
void shuffle_real_matrix(double** matrix, int m) {
    // A pointer used for swapping.
    double* temp = NULL;
    // Integer used to hold random index.
    int j;
    // Use the Fischer-Yates shuffle algorithm.
    for (int i = m - 1; i > 0; i--) {
        // NOTE: Random seed must be set first!
        j = (rand() % (i + 1));
        temp = matrix[j];
        matrix[j] = matrix[i];
        matrix[i] = temp;
    }
}

// Print the contents of an array.
// Accepts:
//  double* array -> The array to print.
//  int n -> The number of elements in the array.
//  int width -> The field witdth.
//  int percision -> The percision of the float.
// Returns: void.
void print_real_array(double* array, int n, int width, int percision) {
    for (int i = 0; i < n; i++) {
        // Width of 7 and precision of 4.
        printf("%*.*lf ", width, percision, array[i]);
    }
}

// Print the contents of a matrix.
// Accepts:
//  double** matrix -> The matrix to print.
//  int n -> The length of the matrix.
//  int m -> The width of the matrix.
//  int width -> The field witdth.
//  int percision -> The percision of the float.
// Returns: void.
void print_real_matrix(double** matrix, int m, int n, int width, int percision) {
     // Iterate through the elements and print row-by-row.
    for (int i = 0; i < m; i++) {
        print_real_array(matrix[i], n, width, percision);
        printf("\n");
    }
}

// Fills matrix with a value.
// Accepts:
//  double** matrix -> The matrix to print.
//  int value -> The value to fill the matrix with.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void fill_real_matrix(double** matrix, int value, int m, int n) {
    // Iterate through the elements and set value.
    for (int i = 0; i < m; i++) {
        for (int j = i; j < n; j++) {
            matrix[i][j] = matrix[j][i] = value;
        }
    }
}

// Deep copy a real matrix.
// Accepts:
//  double** matrix -> The matrix to deep copy.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, n x m matrix of copied values from matrix.
double** deep_copy_real_matrix(double** matrix, int m, int n) {
    // Allocate our copy.
    double** copy = create_real_matrix(m, n);
    // Fill copy with values from matrix.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            copy[i][j] = matrix[i][j]; 
        }
    }
    return copy;
}

// De-allocates a matrix of doubles.
// Accepts:
//  double** matrix -> The matrix to free.
//  int m -> The length of the matrix.
// Returns: void.
void destroy_real_matrix(double** matrix, int m) {
    // Delete all the rows.
    for (int i = 0; i < m; i++) {
        delete [] matrix[i];
    }
    // Delete the matrix.
    delete [] matrix;
}

#endif