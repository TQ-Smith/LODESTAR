
#ifndef _PROCRUSTES_ANALYSIS_H_
#define _PROCRUSTES_ANALYSIS_H_

#include "RealSymEigen.h"

#include <stdbool.h>

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen* eigen, int n, int k, bool transform, bool similarity);

double permutation_test(double** Xc, double** Yc, double** shuffleX, RealSymEigen* eigen, int n, int k, bool similarity, double t0, int NUM_PERMS);

#endif