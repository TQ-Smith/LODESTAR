//
// File: MatrixOperations.hpp
// Started: 15 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for MatrixOperations.cpp
//

#ifndef MATRIX_OPERATIONS_HPP_  
#define MATRIX_OPERATIONS_HPP_

// Create a matrix of doubles.
// Accepts:
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, n x m matrix filled with 0s.
double** create_real_matrix(int m, int n);

// Shuffles the rows of a matrix.
//  Used for permutation testing.
// Note: Assumes random number generator has been seeded. 
// Accepts:
//  double** matrix -> The matrix to shuffle.
//  int m -> The length of the matrix.
// Returns: void.
void shuffle_real_matrix(double** matrix, int m);

// Print the contents of an array.
// Accepts:
//  double* array -> The array to print.
//  int n -> The number of elements in the array.
// Returns: void.
void print_real_array(double* array, int n);

// Print the contents of a matrix.
// Accepts:
//  double** matrix -> The matrix to print.
//  int n -> The length of the matrix.
//  int m -> The width of the matrix.
// Returns: void.
void print_real_matrix(double** matrix, int m, int n);

// Fills matrix with a value.
// Accepts:
//  double** matrix -> The matrix to print.
//  int value -> The value to fill the matrix with.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void fill_real_matrix(double** matrix, int value, int m, int n);

// Deep copy a real matrix.
// Accepts:
//  double** matrix -> The matrix to deep copy.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, n x m matrix of copied values from matrix.
double** deep_copy_real_matrix(double** matrix, int m, int n);

// De-allocates a matrix of doubles.
// Accepts:
//  double** matrix -> The matrix to free.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void destroy_real_matrix(double** matrix, int m, int n);

#endif