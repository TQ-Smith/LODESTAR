
// File: ProcrustesAnalysis.h
// Date: 5 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute Procrustes statistic between two sets of points and perfrom permutation test.
// References: Wang et al, Procrustes Analysis in Population Genetics, 2010.

#ifndef _PROCRUSTES_ANALYSIS_H_
#define _PROCRUSTES_ANALYSIS_H_

#include "RealSymEigen.h"
#include "Window.h"
#include <stdbool.h>

// Compute Procrustes statistic between two sets of points.
// Accepts:
//  double** Xc -> The n-by-k mean-centered query point set.
//  double* x0 -> The k-dimensional vector holding the mean of each dimension.
//                  If NULL, x0 is assumed to be the 0 vector.
//  double** Yc -> The n-by-k mean centered target point set.
//  double* y0 -> The k-dimensional vector holding the mean of each dimension.
//                  If NULL, x0 is assumed to be the 0 vector.
//  RealSymEigen_t* eigen -> Used to find eigen-pairs of k-by-k matrix. Allocated
//                  memory is used to store covariance matrix even if Xc is not to be
//                  transformed.
//  int N -> The number of points.
//  int K -> The dimension of each point.
//  bool transform -> If set, x0 will be transformed by Procrustes analysis.
//  bool similarity -> If set, statistic represents similarity between points. Otherwise, dissimilarity.
// Returns: double, The Procrustes statistic.
double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen_t* eigen, int N, int K, bool transform, bool similarity);

// Perfrom permutation test to evaluate the significance of the Procrustes statistic t0.
// Accepts:
//  double** Xc -> The n-by-k mean-centered query point set.
//  double** Yc -> The n-by-k mean centered target point set.
//  double** shuffleX -> A copy of Xc used to hold permuted rows of Xc.
//  RealSymEigen_t* eigen -> Used to find eigen-pairs of k-by-k matrix. Allocated
//                  memory is used to store covariance matrix. Points are not transformed 
//                  during permutation test.
//  int N -> The number of points.
//  int K -> The dimension of each point.
//  bool similarity -> If set, statistic represents similarity between points. Otherwise, dissimilarity.
//  double t0 -> The Procrustes statistic used in the permutation test.
//  int NUM_PERMS -> The number of permutations to execute.
// Returns: double, The p-value obtained from the permutation test.
double permutation_test(double** Xc, double** Yc, double** shuffleX, RealSymEigen_t* eigen, int N, int K, bool similarity, double t0, int NUM_PERMS);


// Perform Procrustes Analysis and/or permutation test along a sliding window.
// Accepts:
//  Window_t** windows -> Our array of windows.
//  int numWindows -> The length of windows.
//  double** target -> The target set of points.
//  double* target0 -> The target column means vector.
//  int N -> The number of points.
//  int K -> The dimension of each point.
//  bool similarity -> If set, statistic represents similarity between points. Otherwise, dissimilarity.
//  int NUM_PERMS -> If 0, permutation test is NOT performed.
// Returns: void.
void procrustes_sliding_window(Window_t** windows, int numWindows, double** target, double* target0, int N, int K, bool similarity, int NUM_PERMS);

#endif