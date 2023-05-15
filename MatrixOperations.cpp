//
// File: MatrixOperations.cpp
// Started: 15 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains basic operations for manipulating matrices.
//

#include "MatrixOperations.h"

// Used for input and output.
#include <iostream>

// Unless otherwise specified, our streams come from std.
using namespace std;

double** create_real_matrix(int m, int n) {

    double** matrix = new double*[m];

    // Create each row.
    for (int i = 0; i < m; i++) {
        matrix[i] = new double[n];
    }

    return matrix;

}


void shuffle_real_matrix(double** matrix, int m) {

    // A pointer used for swapping.
    double* temp = NULL;

    // Integer used to hold random index.
    int j;

    // Use the Fischer-Yates shuffle algorithm.
    for (int i = m - 1; i > 0; i--) {
        j = (rand() % (i + 1));
        temp = matrix[j];
        matrix[j] = matrix[i];
        matrix[i] = temp;
    }

} 

void print_real_array(double* array, int n) {
    for (int i = 0; i < n; i++) {
        // Width of 7 and precision of 4.
        printf("%7.4lf ", array[i]);
    }
}

void print_real_matrix(double** matrix, int m, int n) {

    // Iterate through the elements and print row-by-row.
    for (int i = 0; i < m; i++) {
        print_real_array(matrix[i], n);
        printf("\n");
    }

}

void fill_real_matrix(double** matrix, int value, int m, int n) {
    // Iterate through the elements and set value.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = value;
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

void destroy_real_matrix(double** matrix, int m, int n) {
    
    // Delete all the rows.
    for (int i = 0; i < m; i++) {
        delete [] matrix[i];
    }
    // Delete the matrix.
    delete [] matrix;

}