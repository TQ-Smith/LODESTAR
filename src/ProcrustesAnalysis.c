
// File: ProcrustesAnalysis.c
// Date: 5 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute Procrustes statistic between two sets of points and perform permutation test.
// References: Wang et al, Procrustes Analysis in Population Genetics, 2010.

#include "ProcrustesAnalysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Fortran matrices are stored as vectors in column-major form.
//  Covariance matrices are symmetric so this is not explicitly necessary,
//  but it is good to keep things consistent.
#define COL_MAJOR(i, j, K) (j * K + i)

// Fischer-Yates shuffle algorithm to shuffle the rows of a matrix.
//  NOTE: Random number generator needs to be seeded before calling function.
// Accepts:
//  double** matrix -> The matrix to be shuffled.
//  int m -> The number of rows in the matrix.
// Returns: void.
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

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen_t* eigen, int N, int K, bool transform, bool similarity) {

    // Trace of (Xc^T)Xc, trace of (Yc^T)Yc.
    double trX = 0, trY = 0;

    // For convenience, we rename Z and A to C and covC, respecively.
    //  We will perfrom SVD on C, which means we will take the eigen-decomposition
    //  of the covariance of C.
    double* C = eigen -> Z;
    double* covC = eigen -> A;

    // Calculate C = (Yc^T)Xc. Simultaneously calculate the trace of (Xc^T)Xc and trace of (Yc^T)Yc.
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            C[COL_MAJOR(i, j, K)] = 0;
            for (int k = 0; k < N; k++) {
                C[COL_MAJOR(i, j, K)] += Yc[k][i] * Xc[k][j];
                if (j == 0) {
                    trX += Xc[k][i] * Xc[k][i];
                    trY += Yc[k][i] * Yc[k][i]; 
                }
            }
        }
    }

    // Calculate covC = (C^T)C.
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            covC[COL_MAJOR(i, j, K)] = 0;
            for (int k = 0; k < K; k++)
                covC[COL_MAJOR(i, j, K)] += C[COL_MAJOR(k, i, K)] * C[COL_MAJOR(k, j, K)];
        }
    }

    // Trace of the lambda matrix.
    double trLambda = 0;

    // In the case we do not transform the X coordinates, we just need the trace of the lambda matrix.
    //  When K = 1, 2, 3, we can calculate the eigenvalues of covC directly as the roots of its characteric equation.
    //  We used Numerical Recipes in C, Third Edition to solve the quadratic and cubic equations.
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
                // Cubic is in the form x^3 + a0x^2 + a1x^1 + a0 = 0.
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
                // When K > 3, we use LAPACK to calculate the eigenvalues.
                compute_k_eigenvalues(eigen, K);
                for (int i = 0; i < K; i++)
                    trLambda += sqrt(eigen -> W[i]);
                break;
        }
    } else {

        // When we transform the coordinates, we need both the eigenvalues and eigenvectors.
        //  For all K, we rely on LAPACK to calculate eigenpairs. We could do this manually 
        //  for K = 1, 2, and 3, but the code would be somewhat lengthy.
        compute_k_eigenpairs(eigen, K);

        // Calculate trLambda.
        for (int i = 0; i < K; i++)
            trLambda += sqrt(eigen -> W[i]);

        // For convenience, we rename A and Z to U and V, respectively. 
        double* U = eigen -> A;
        double* V = eigen -> Z;

        // Caclulate our scaling factor rho.
        double rho = trLambda / trX;

        // Calculate U using the equation Ui = (1 / sigma_i) * C * Vi
        //  We could create another auxilary k-by-k array to store C, which would not require much memory
        //  just more book-keeping. However, since k is usually small, we can recalculate C quickly.
        //  We use eigen -> WORK[j] as a temporary vector.
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

        // Define our shift vector. Known as b in the paper.
        double* shift = eigen -> W;

        // We calculate (A^T) = U(V^T). Again, eigen -> WORK is used as an auxilary vector.
        //  NOTE: U is replaced by (A^T).
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                eigen -> WORK[j] = 0;
                for (int k = 0; k < K; k++)
                    eigen -> WORK[j] += U[COL_MAJOR(i, k, K)] * V[COL_MAJOR(K - k - 1, j, K)];
            }
            // As we form A^T, we also calculate shift.
            shift[i] = 0;
            for (int j = 0; j < K; j++) {
                // Copy over WORK to build A^T.
                U[COL_MAJOR(i, j, K)] = eigen -> WORK[j];
                if (x0 != NULL)
                    shift[i] += U[COL_MAJOR(i, j, K)] * x0[j];
            }
            if (y0 != NULL)
                shift[i] = y0[i] - rho * shift[i];
            else
                shift[i] = -rho * shift[i];
        }

        // Transform each of the samples by the equation rho * (A^T) * X + shift.
        //  We could put this in the previous loop, but since K is usually small,
        //  we keep is seperate for clarity.
        for (int i = 0; i < N; i++) {
            // Sample vector times (A^T).
            for (int j = 0; j < K; j++) {
                eigen -> WORK[j] = 0;
                for (int k = 0; k < K; k++)
                    eigen -> WORK[j] += U[COL_MAJOR(j, k, K)] * Xc[i][k];
            }
            // Scale and shift sample.
            for (int j = 0; j < K; j++)
                Xc[i][j] = rho * eigen -> WORK[j] + shift[j];
        }

    }
    
    // Calculate statistic.
    double statistic = 1 - (trLambda * trLambda) / (trX * trY);

    // Convert to either similarity or dissimilarity.
    return similarity ? sqrt(1 - statistic) : sqrt(statistic);

}

double permutation_test(double** Xc, double** Yc, double** shuffleX, RealSymEigen_t* eigen, int n, int k, bool similarity, double t0, int NUM_PERMS) {

    // Copy X to shuffleX.
    for (int i = 0; i < n; i++)
        for (int j = 0; j < k; j++)
            shuffleX[i][j] = Xc[i][j];

    double t;
    int numSig = 1;

    // Perform NUM_PERMS permutations.
    for (int i = 0; i < NUM_PERMS; i++) {
        // Shuffle the rows of X.
        shuffle_real_matrix(shuffleX, n);
        // Compute statistic between the two sets of points. Do not transform.
        t = procrustes_statistic(shuffleX, NULL, Yc, NULL, eigen, n, k, false, similarity);
        // If our statistic exceeds our threshold, increment counter.
        if (t > t0)
            numSig++;
    }

    // Return p-value.
    return ((double) numSig) / (NUM_PERMS + 1);
}

/*
int main() {
    int N = 4, K = 2;
    double** X = (double**) calloc(N, sizeof(double*));
    double* x0 = (double*) calloc(3, sizeof(double));
    double** Y = (double**) calloc(N, sizeof(double*));
    double* y0 = (double*) calloc(3, sizeof(double));
    for (int i = 0; i < N; i++) {
        X[i] = (double*) calloc(3, sizeof(double));
        Y[i] = (double*) calloc(3, sizeof(double));
    }
    X[0][0] = 1; X[0][1] = 2; X[0][2] = 2;
    X[1][0] =  3; X[1][1] = 4; X[1][2] = 3;
    X[2][0] = 3; X[2][1] = -6; X[2][2] = -4;
    X[3][0] =  7; X[3][1] = -8; X[3][2] = -5;
    x0[0] = -1; x0[1] = 3; x0[2] = -2;

    Y[0][0] = 1; Y[0][1] = 0; Y[0][2] = -1;
    Y[1][0] =  0; Y[1][1] = 2; Y[1][2] = -3;
    Y[2][0] = -5; Y[2][1] = 0; Y[2][2] = 5;
    Y[3][0] =  7; Y[3][1] = 0; Y[3][2] = 7;
    y0[0] = 1; y0[1] = 2; y0[2] = -1;

    RealSymEigen_t* eigen = init_real_sym_eigen(K);

    double t0 = procrustes_statistic(X, x0, Y, y0, eigen, N, K, true, true);

    printf("Statistic: %lf\n", t0);
}
*/