//
// File: NumericalRecipesInC.h
// Started: 15 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for NumericalRecipesInC.cpp
//

#ifndef NUMERICAL_RECIPES_IN_C_H_  
#define NUMERICAL_RECIPES_IN_C_H_

// Taken from Numeric Recipes in C.
// Note: Out of respect, I left as much of this method 
// unchanged as I could.
//  Sorts column eigenvector matrix and correponding
//  eigenvalues in order from greatest to least eigenvalue.
// Accepts:
//  double** v -> The eigenvector matrix.
//  double* d -> The eigenvalues.
// Returns: void.
void eigsrt(double** v, double* d, int n);

// Taken from Numeric Recipes in C.
// Note: Out of respect, I left as much of this method 
// unchanged as I could.
//  Converts a real, symmetric matrix to a tridiagonal matrix
//  with the same eigenvectors, eigenvalues.
// Accepts:
//  double** z -> The real, symmetric martix that will be transformed.
//                  z will be transformed to the orthogonal matrix.
//  double* d -> An array to hold the diagonal elements.
//  double* e -> An array to hold the off diagonal elements.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
// Returns: void.
void tred2(double** z, double* d, double* e, int n, bool yesvecs);

// Taken from Numeric Recipes in C.
// Note: Out of respect, I left as much of this method 
// unchanged as I could.
//  Calculates eigenvalues and/or eigenvectors of a tridiagonal
//  matrix using the QR method. 
// Accepts:
//  double** z -> Column matrix to hold eigenvectors.
//  double* d -> The diagonal elements of the matrix.
//  double* e -> The off diagonal elements of the matrix.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
// Returns: void.
void tqli(double** z, double* d, double* e, int n, bool yesvecs);

#endif