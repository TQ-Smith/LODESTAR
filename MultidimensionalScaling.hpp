//
// File: MultidimensionalScaling.hpp
// Started: 15 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for MultidimensionalScaling.cpp
//

#ifndef MULTIDIMENSIONAL_SCALING_HPP_  
#define MULTIDIMENSIONAL_SCALING_HPP_

// Performs Classical MDS on a real, symmetric matrix.
//  In this application, this is a dissimilarity matrix.
// Accepts:
//  double** X -> Our n x n distance matrix with 1s on the diagonal.
//                  Is set to the mtrix of reduced points.
//  double* d -> Used to hold eigen values.
//  double* e -> Used as temp for eigen computations.
//  int n -> The dimension of n.
//  int k -> The dimension of the reduced data.
//      Assume 0 < k < n.
// Return: void.
void compute_classical_mds(double**& X, double* d, double* e, int n, int k);

// Performs Classical MDS on a real, symmetric matrix.
// This method handles all memory allocations.
//  In this application, this is a dissimilarity matrix.
// Accepts:
//  double** X -> Our n x n distance matrix with 1s on the diagonal.
//                  Is set to the matrix of reduced points.
//  int n -> The dimension of n.
//  int k -> The dimension of the reduced data.
//      Assume 0 < k < n.
// Return: void.
void classical_mds(double** X, int n, int k);

#endif