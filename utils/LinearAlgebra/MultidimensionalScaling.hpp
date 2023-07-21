//
// File: MultidimensionalScaling.hpp
// Started: 21 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Contains methods for classic multidimensional scaling
//              and the FastMap heuristic.
//

// FastMap Citation:
// Faloutsos, Christos; Lin, King-Ip (1995). FastMap: a fast algorithm for indexing, 
//      data-mining and visualization of traditional and multimedia datasets. 
//      Carnegie Mellon University. Journal contribution. https://doi.org/10.1184/R1/6605624.v1

#ifndef _MULTIDIMENSIONAL_SCALING_HPP_  
#define _MULTIDIMENSIONAL_SCALING_HPP_

// Performs Classical MDS on a real, symmetric matrix.
// NOTE: D will be clobbered at the end of this routine.
// Accepts:
//  double** D -> Our n x n distance matrix.
//  double** X -> An n x k matrix to hold the points in dimension k.
//  double* d -> Used to hold eigen values.
//  double* e -> Used as temp for eigen computations.
//  bool* doesConverge -> Sets bool if eigenpairs were successfully calculated.
//  int n -> The dimension of the data.
//  int k -> The dimension of the reduced data.
//      Assume 0 < k < n.
// Return: void.
void compute_classical_mds(double** D, double** X, double* d, double* e, bool* doesConverge, int n, int k);

// Perfroms FastMap using a distance matrix.
// Accepts:
//  double** D -> Our n x n distance matrix.
//  double** X -> An n x k matrix to hold the points in dimension k.
//  int* maxDimReached -> Sets the maximum number of projections achieved.
//  int n -> The dimension of the data.
//  int k -> The dimension of the reduced data.
//      Assume 0 < k < n.
// Return: void.
void compute_fastmap(double** D, double** X, int* maxDimReached, int n, int k);

#endif