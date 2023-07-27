//
// File: Procrustes.cpp
// Date: 27 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Perform Procrustes analysis from two real, symmetric distance
//              matrices in dimension k and a permutation test using
//              the Procrustes statistic.
//
// Note: Our current implementation only handles k = 1, 2, 3 and does
//          not calculate the rotation/scaling factors for the Procrustes
//          analysis.
//

#include "Procrustes.hpp"

// Used for creating, destroying, deep copying, and shuffling matrices.
#include "MatrixOperations.hpp"

// Used for eigenvalue/eigenvector computation when 3 < k < n.
#include "NumericalRecipesInC.hpp"

// Used for the basic math operatons sqrt, cos, and acos.
#include <cmath>

// Our helper function for Procrustes analysis.
// Accepts:
//  double** Xc -> The first centered set of points.
//  double** Yc -> The second centered set of points.
//  double** C -> A k x k matrix used to hold the value of YcT_Xc.
//  double** CT_C -> A k x k matrix used to hold the value of CT_C.
//  int n -> The number of points.
//  int k -> The dimension of each point. Assume k = 1, 2, 3.
// Returns: double, The Procrustes statistic.
double procrustes(double** Xc, double** Yc, double** C, double** CT_C, int n, int k) {
    
}