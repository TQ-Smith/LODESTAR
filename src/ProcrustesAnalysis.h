
#ifndef _PROCRUSTES_ANALYSIS_H_
#define _PROCRUSTES_ANALYSIS_H_

#include "RealSymEigen.h"

#include <stdbool.h>

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen* eigen, int n, int k, bool transform, bool similarity);

#endif