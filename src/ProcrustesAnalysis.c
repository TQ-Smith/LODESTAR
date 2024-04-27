
#include "ProcrustesAnalysis.h"

#include <math.h>

#include <stdio.h>

#include <stdlib.h>

#define COL_MAJOR(i, j, K) (j * K + i)

void shuffle_real_matrix(double** matrix, int m) {
    double* temp = NULL;
    int j;
    for (int i = m - 1; i > 0; i--) {
        j = (rand() % (i + 1));
        temp = matrix[j];
        matrix[j] = matrix[i];
        matrix[i] = temp;
    }
}

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
                trLambda += sqrt(covC[COL_MAJOR(0, 0, K)]);
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
                trLambda += sqrt(root1) + sqrt(root2) + sqrt(root3);
                break;
            default:
                compute_k_eigenvalues(eigen, K);
                for (int i = 0; i < K; i++)
                    trLambda += sqrt(eigen -> W[i]);

                break;
        }

    } else {

        compute_k_eigenpairs(eigen, K);

        for (int i = 0; i < K; i++)
            trLambda += sqrt(eigen -> W[i]);

        double* U = eigen -> A;
        double* V = eigen -> Z;
        double* shift = eigen -> W;

        double rho = trLambda / trX;

        /*
        printf("\nEigenvectors:\n");
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++)
                printf("%lf\t", V[COL_MAJOR(K - i - 1, j, K)]);
            printf("\n");
        }
        printf("\nEigenvalues:\n");
        for (int i = 0; i < K; i++)
            printf("%lf\t", sqrt(eigen -> W[K - i - 1]));
        printf("\n");
        */

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                eigen -> WORK[j] = 0;
                for (int k = 0; k < N; k++)
                    eigen -> WORK[j] += Yc[k][i] * Xc[k][j];
            }
            for (int j = 0; j < K; j++) {
                U[COL_MAJOR(i, j, K)] = 0;
                for (int k = 0; k < K; k++)
                    U[COL_MAJOR(i, j, K)] += eigen -> WORK[k] * V[COL_MAJOR(K - j - 1, k, K)];
                U[COL_MAJOR(i, j, K)] *= (1 / sqrt(eigen -> W[K - j - 1]));
            }
        }

        /*
        printf("\nU:\n");
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++)
                printf("%lf\t", U[COL_MAJOR(i, j, K)]);
            printf("\n");
        }
        */

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                eigen -> WORK[j] = 0;
                for (int k = 0; k < K; k++)
                    eigen -> WORK[j] += U[COL_MAJOR(i, k, K)] * V[COL_MAJOR(K - j - 1, k, K)];
            }
            shift[i] = 0;
            for (int j = 0; j < K; j++) {
                U[COL_MAJOR(i, j, K)] = eigen -> WORK[j];
                if (x0 != NULL && y0 != NULL)
                    shift[i] += U[COL_MAJOR(i, j, K)] * x0[j];
            }
            if (x0 != NULL && y0 != NULL)
                shift[i] = y0[i] - rho * shift[i];
        }

        /*
        printf("\nA^T:\n");
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++)
                printf("%lf\t", U[COL_MAJOR(i, j, K)]);
            printf("\n");
        }

        printf("\nrho: %lf\n", rho);
        printf("\ntrX: %lf\n", trX);
        printf("\ntrY: %lf\n", trY);
        printf("\ntrLambda: %lf\n", trLambda);

        printf("\nShift:\n");
        for (int i = 0; i < K; i++)
            printf("%lf\t", shift[i]);
        printf("\n\n");
        */

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < K; j++) {
                eigen -> WORK[j] = 0;
                for (int k = 0; k < K; k++)
                    eigen -> WORK[j] += U[COL_MAJOR(j, k, K)] * Xc[i][k];
            }
            for (int j = 0; j < K; j++) {
                Xc[i][j] = rho * eigen -> WORK[j];
                if (x0 != NULL && y0 != NULL)
                    Xc[i][j] += shift[j];
            }
        }

    }

    double statistic = 1 - (trLambda * trLambda) / (trX * trY);

    return similarity ? sqrt(1 - statistic) : sqrt(statistic);

}

double permutation_test(double** Xc, double** Yc, double** shuffleX, RealSymEigen* eigen, int n, int k, bool similarity, double t0, int NUM_PERMS) {

    for (int i = 0; i < n; i++)
        for (int j = 0; j < k; j++)
            shuffleX[i][j] = Xc[i][j];

    double t;
    int numSig = 1;

    for (int i = 0; i < NUM_PERMS; i++) {
        shuffle_real_matrix(shuffleX, n);
        t = procrustes_statistic(shuffleX, NULL, Yc, NULL, eigen, n, k, false, similarity);
        if (t > t0)
            numSig++;
    }

    return ((double) numSig) / (NUM_PERMS + 1);

}