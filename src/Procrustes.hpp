//
// File: Procrustes.hpp
// Date: 27 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Perform Procrustes analysis and permutation test as described in 
//              Wang, Chaolong, Szpiech, Zachary A, Degnan, James H, Jakobsson, Mattias, Pemberton, 
//                  Trevor J, Hardy, John A, Singleton, Andrew B and Rosenberg, Noah A. "Comparing 
//                  Spatial Maps of Human Population-Genetic Variation Using Procrustes Analysis" 
//                  Statistical Applications in Genetics and Molecular Biology, vol. 9, no. 1, 2010. 
//                  https://doi.org/10.2202/1544-6115.1493
//

#ifndef _PROCRUSTES_HPP_
#define _PROCRUSTES_HPP_

// Our function that performs Procrustes analysis.
// Accepts:
//  double** Xc -> The first centered set of points.
//  double** Yc -> The second centered set of points.
//  double** C -> A k x k matrix used to hold the value of YcT_Xc.
//  double** CT_C -> A k x k matrix used to hold the value of CT_C.
//  int n -> The number of points.
//  int k -> The dimension of each point. Assume k = 1, 2, 3.
// Returns: double, The Procrustes statistic.
double procrustes_analysis(double** Xc, double** Yc, double** C, double** CT_C, int n, int k);

// Our function that performs a permutation test with the Procrustes statistic.
// Accepts:
//  int NUM_PERMUTATIONS -> The number of permutations to execute in the test.
//  double** Xc -> The first centered set of points.
//  double** Yc -> The second centered set of points.
//  double** shuffleXc -> A matrix the same size of Xc used to hold the permuted matrix.
//  double** C -> A k x k matrix used to hold the value of YcT_Xc.
//  double** CT_C -> A k x k matrix used to hold the value of CT_C.
//  int n -> The number of points.
//  int k -> The dimension of each point. Assume k = 1, 2, 3.
//  double t_0 -> The initial Procrustes statistic used in the permutation test.
// Returns: double, The p-value from the permutation test.
double permutation_test(int NUM_PERMUTATIONS, double** Xc, double** Yc, double** shuffleXc, double** C, double** CT_C, int n, int k, double t_0);

// A method to mean center a matrix.
//  This could be incorporated in the Procrustes
//  methods; however, seperating the methods saves
//  computation when multiple analyses are executed
//  using the same matrices.
// Accepts:
//  double** X -> The n x k matrix to be centered.
//  double* x_0 -> The 1 x k array used to hold means.
//                  When k is small, we could create x_0
//                  within the method, but this method will
//                  be called alot.
//  int n -> The number of points.
//  int k -> The dimension of each point.
// Returns: void.
void center_matrix(double** X, double* x_0, int n, int k);

#endif