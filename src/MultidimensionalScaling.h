
// File: MultidimensionalScaling.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Perform classical MDS on a distance matrix using LAPACK's dsyevr routine.

#ifndef _MULTIDIMENSIONAL_SCALING_H_
#define _MULTIDIMENSIONAL_SCALING_H_

#include "RealSymEigen.h"

// Perfrom classical MDS on a distance matrix.
// Accepts:
//  RealSymEigen_t* eigen -> Structure with memory to execute LAPACK's dsyevr routine.
//  double* packedDistanceMatrix -> The upper triangle of the distance matrix in packed storage.
//                                  Not perserved.
//  int k -> The dimension to project down into.
//  double** X -> An allocated N-by-k matrix to store the resulting points.
// Returns: int, 0 for success. Otherwise, LAPACK failure, or k eigenvalues are not all positive.
int compute_classical_mds(RealSymEigen_t* eigen, double* packedDistanceMatrix, int k, double** X);

#endif