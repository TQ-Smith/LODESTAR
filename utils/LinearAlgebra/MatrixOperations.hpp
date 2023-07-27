//
// File: MatrixOperations.hpp
// Date: 19 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains basic operations for manipulating m x n matrices
//              of doubles.
//

#ifndef _MATRIX_OPERATIONS_HPP_  
#define _MATRIX_OPERATIONS_HPP_

// Create a matrix of doubles.
// Accepts:
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, the allocated matrix.
double** create_real_matrix(int m, int n);

// Create a matrix of doubles.
// Accepts:
//  int value -> The value that fills the matrix.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, the allocated matrix.
double** create_and_fill_real_matrix(int value, int m, int n);

// Adds matrices b and c and stores in a.
// Accepts:
//  double** a -> The resulting matrix.
//  double** b -> The first operand.
//  double** c -> The second operand.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void add_matrices(double** a, double** b, double** c, int m, int n);

// Subtracts matrix c from b and stores in a.
// Accepts:
//  double** a -> The resulting matrix.
//  double** b -> The first operand.
//  double** c -> The second operand.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void subtract_matrices(double** a, double** b, double** c, int m, int n);

// Scale a matrix by a constant.
// Accepts:
//  double c -> The constant.
//  double** a -> The matrix.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void scale_matrix(double c, double** a, int m, int n);

// Matrix multiplication. a = b * c.
// Accepts:
//  double** a -> The resulting matrix.
//  double** b -> The first operand.
//  double** c -> The second operand.
//  int m -> The length of the first matrix.
//  int n -> The width of the first matrix.
//  int o -> The width of the second matrix.
// Returns: void.
void multiply_matrices(double** a, double** b, double** c, int m, int n, int o);

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
//  int width -> The field witdth.
//  int percision -> The percision of the float.
// Returns: void.
void print_real_array(double* array, int n, int width, int percision);

// Print the contents of a matrix.
// Accepts:
//  double** matrix -> The matrix to print.
//  int n -> The length of the matrix.
//  int m -> The width of the matrix.
//  int width -> The field witdth.
//  int percision -> The percision of the float.
// Returns: void.
void print_real_matrix(double** matrix, int m, int n, int width, int percision);

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
// Returns: void.
void destroy_real_matrix(double** matrix, int m);

#endif