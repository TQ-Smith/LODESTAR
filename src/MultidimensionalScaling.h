
#ifndef _MULTIDIMENSIONAL_SCALING_H_
#define _MULTIDIMENSIONAL_SCALING_H_

#include "RealSymEigen.h"

int compute_classical_mds(RealSymEigen* eigen, double* packedDistanceMatrix, int k, double** X);

#endif