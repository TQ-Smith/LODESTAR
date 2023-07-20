//
// File: NumericalRecipesInC.hpp
// Started: 20 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains a function to compute the eigenvectors and 
//              eigenvalues of a real, symmetric matrix. Taken
//              from Numerical Recipes in C, Third Edition.
//

#ifndef _NUMERICAL_RECIPES_IN_C_HPP_  
#define _NUMERICAL_RECIPES_IN_C_HPP_

// Wrapper function to compute eigenvalues and eigenvectors as described
//  in Numerical Recipes in C. 
// Accepts:
//  double** z -> Column matrix to hold eigenvectors.
//  double* d -> The diagonal elements of the matrix.
//  double* e -> The off diagonal elements of the matrix.
//  double* doesConverge -> Sets boolean if the algorithm converges.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
//  bool sort -> Flag set to sort eigenvectors and eigenvalues.
// Returns: void.
void compute_eigen_pairs(double** z, double* d, double* e, bool* doesConverge, int n, bool yesvecs, bool sort);

#endif