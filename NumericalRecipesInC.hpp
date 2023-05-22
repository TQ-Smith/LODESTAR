//
// File: NumericalRecipesInC.hpp
// Started: 15 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for NumericalRecipesInC.cpp
//

#ifndef NUMERICAL_RECIPES_IN_C_HPP_  
#define NUMERICAL_RECIPES_IN_C_HPP_

// Wrapper function to compute eigenvalues and eigenvectors as described
//  in Numerical Recipes in C. 
// Accepts:
//  double** z -> Column matrix to hold eigenvectors.
//  double* d -> The diagonal elements of the matrix.
//  double* e -> The off diagonal elements of the matrix.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
//  bool sort -> Flag set to sort eigenvectors and eigenvalues.
// Returns: void.
void compute_eigen_pairs(double** z, double* d, double* e, int n, bool yesvecs, bool sort);

#endif