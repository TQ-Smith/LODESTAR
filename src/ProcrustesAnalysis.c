
#include "ProcrustesAnalysis.h"

#include <math.h>

#include <stdio.h>

#define COL_MAJOR(i, j, K) (j * K + i)

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen* eigen, int N, int K, bool transform, bool similarity) {

    double trX = 0, trY = 0, trLambda = 0;

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < N; j++) {
            trX += Xc[j][i] * Xc[j][i];
            trY += Yc[j][i] * Yc[j][i]; 
        }
    }

    double* C = eigen -> Z;
    double* covC = eigen -> A;

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            C[COL_MAJOR(i, j, K)] = 0;
            for (int k = 0; k < N; k++)
                C[COL_MAJOR(i, j, K)] += Yc[k][i] * Xc[k][j];
        }
    }

    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            covC[COL_MAJOR(i, j, K)] = 0;
            for (int k = 0; k < K; k++)
                covC[COL_MAJOR(i, j, K)] += C[COL_MAJOR(k, i, K)] * C[COL_MAJOR(k, j, K)];
        }
    }

    if (!transform) {
        double a0, a1, a2, b, c, det, q, r, theta, root1, root2, root3;
        switch (K) {
            case 1:
                trLambda += covC[COL_MAJOR(0, 0, K)];
                break;
            case 2:
                b = -(covC[COL_MAJOR(0, 0, K)] + covC[COL_MAJOR(1, 1, K)]);
                c = (covC[COL_MAJOR(0, 0, K)] * covC[COL_MAJOR(1, 1, K)]) - (covC[COL_MAJOR(0, 1, K)] * covC[COL_MAJOR(1, 0, K)]);
                
                det = sqrt(b * b - 4 * c);
                q = b < 0 ? -0.5 * (b - det) : -0.5 * (b + det);

                trLambda += sqrt(q) + sqrt(c / q);
                break;
            case 3:
                a2 = -(covC[COL_MAJOR(0, 0, K)] + covC[COL_MAJOR(1, 1, K)] + covC[COL_MAJOR(2, 2, K)]);
                a1 = ((covC[COL_MAJOR(1, 1, K)] * covC[COL_MAJOR(2, 2, K)]) - (covC[COL_MAJOR(1, 2, K)] * covC[COL_MAJOR(2, 1, K)]))
                    + ((covC[COL_MAJOR(0, 0, K)] * covC[COL_MAJOR(2, 2, K)]) - (covC[COL_MAJOR(0, 2, K)] * covC[COL_MAJOR(2, 0, K)]))
                    + ((covC[COL_MAJOR(0, 0, K)] * covC[COL_MAJOR(1, 1, K)]) - (covC[COL_MAJOR(0, 1, K)] * covC[COL_MAJOR(1, 0, K)]));
                a0 = -(covC[COL_MAJOR(0, 0, K)] * ((covC[COL_MAJOR(1, 1, K)] * covC[COL_MAJOR(2, 2, K)]) - (covC[COL_MAJOR(2, 1, K)] * covC[COL_MAJOR(1, 2, K)])) 
                    - covC[COL_MAJOR(0, 1, K)] * ((covC[COL_MAJOR(1, 0, K)] * covC[COL_MAJOR(2, 2, K)]) - (covC[COL_MAJOR(1, 2, K)] * covC[COL_MAJOR(2, 0, K)])) 
                    + covC[COL_MAJOR(0, 2, K)] * ((covC[COL_MAJOR(1, 0, K)] * covC[COL_MAJOR(2, 1, K)]) - (covC[COL_MAJOR(1, 1, K)] * covC[COL_MAJOR(2, 0, K)])));
                q = a1 / 3 - (a2 * a2) / 9;
                r = (a1 * a2 - 3 * a0) / 6 - (a2 * a2 * a2) / 27;
                theta = (q == 0) ? 0 : acos(r / sqrt(-q * -q * -q));
                root1 = 2 * sqrt(-q) * cos(theta / 3) - (a2 / 3);
                root2 = 2 * sqrt(-q) * cos((theta / 3) - (2 * M_PI / 3)) - (a2 / 3);
                root3 = 2 * sqrt(-q) * cos((theta / 3) + (2 * M_PI / 3)) - (a2 / 3);

                trLambda += sqrt(abs(root1)) + sqrt(abs(root2)) + sqrt(abs(root3));
                break;
            default:
                break;
        }

    } else {

    }

    return 0;

}