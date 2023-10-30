//
// File: MatrixOperations.cpp
// Date: 19 July 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Contains basic operations for manipulating m x n matrices
//              of doubles. Note, we do not account for null pointers,
//              values of m and n, and jagged arrays. We do this for two
//              reasons: first, we would have to abstract the double** to 
//              a matrix structure holding additional information, and second,
//              there would be much more boiler plate.
//

// Used for printing.
#include <iostream>
using namespace std;

double** create_real_matrix(int m, int n) {
    // Allocate columns.
    double** matrix = new double*[m];
    // Allocate rows.
    for (int i = 0; i < m; i++) {
        matrix[i] = new double[n];
    }
    return matrix;
}

double** create_and_fill_real_matrix(int value, int m, int n) {
    double** matrix = create_real_matrix(m, n);
    for (int i = 0; i < m; i++) {
        // Fill each row.
        for (int j = i + 1; j < n; j++) {
            matrix[i][j] = matrix[j][i] = value;
        }
        matrix[i][i] = value;
    }
    return matrix;
}

void add_matrices(double** a, double** b, double** c, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            // Add elements.
            a[i][j] = b[i][j] + c[i][j];
        }
    }
}

void subtract_matrices(double** a, double** b, double** c, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            // Subtract elements.
            a[i][j] = b[i][j] - c[i][j];
        }
    }
}

void scale_matrix(double c, double** a, int m, int n) {
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            // Scale the element.
            a[i][j] *= c * a[i][j];
        }
    }
}

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

void print_real_array(double* array, int n, int width, int percision) {
    for (int i = 0; i < n; i++) {
        // Width of 7 and precision of 4.
        printf("%*.*lf ", width, percision, array[i]);
    }
}

void print_real_matrix(double** matrix, int m, int n, int width, int percision) {
     // Iterate through the elements and print row-by-row.
    for (int i = 0; i < m; i++) {
        print_real_array(matrix[i], n, width, percision);
        printf("\n");
    }
}

void fill_real_matrix(double** matrix, int value, int m, int n) {
    // Iterate through the elements and set value.
    for (int i = 0; i < m; i++) {
        for (int j = i; j < n; j++) {
            matrix[i][j] = matrix[j][i] = value;
        }
    }
}

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

void destroy_real_matrix(double** matrix, int m) {
    // Delete all the rows.
    for (int i = 0; i < m; i++) {
        delete [] matrix[i];
    }
    // Delete the matrix.
    delete [] matrix;
}