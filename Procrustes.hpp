//
// File: Procrustes.hpp
// Started: 16 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for Procrustes.cpp
//

#ifndef PROCRUSTES_HPP_  
#define PROCRUSTES_HPP_

// Perform Procrustes Analysis on two sets of points.
// Each set of points are stored in a n x k matrix where
//  there are n k-dimensional points. We assume there are 
//  no duplicate points.
// Accepts:
//  double** Xc -> The first centered set of points.
//  double** Yc -> The second centered set of points.
//  double** temp1 -> A k x k matrix used as a temp variable.
//  double** temp2 -> A k x k matrix used as a temp variable.
//  int n -> The number of points.
//  int k -> The dimension of each point.
//      Assumes 0 < k < n.
// Returns: double, The Procrustes statistic.
double procrustes_statistic(double** Xc, double** Yc, double** temp1, double** temp2, int n, int k);

// A wrapper function for procrustes_statistic that 
//  preforms procrustes analysis between two sets
//  of points.
// Accepts:
//  double** X -> The first set of points.
//  double** Y -> The second set of points.
//  int n -> The number of points.
//  int k -> The dimension of each point.
//      Assumes 0 < k < n.
// Returns: double, The Procrustes statistic.
double procrustes_analysis(double** X, double** Y, int n, int k);

// Performs a permutation test for the Procrustes statistic.
// Accepts:
//  int NUM_PERMUTATIONS -> The number of permutations to run.
//  double** X -> The first set of points.
//  double** Y -> The second set of points.
//  int n -> The number of points.
//  int k -> The dimension of each point.
//  double* t -> Sets the observed initial procrustes statistic, our t_0.
//  double* p_value -> Sets the p_value given by the permutation test.
//      Assumes 0 < k < n.
// Returns: void.
void permutation_test(int NUM_PERMUTATIONS, double** X, double** Y, int n, int k, double* t, double* p_value);

#endif