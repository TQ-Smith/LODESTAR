
// File: ProcrustesAnalysis.c
// Date: 5 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Compute Procrustes statistic between two sets of points and perform permutation test.
// References: Wang et al. 2010. and Mardia's Multivariate Analysis.

#include "ProcrustesAnalysis.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <pthread.h>
#include "Logger.h"
#include "Matrix.h"
MATRIX_INIT(double, double)

// Define machine percision for floating point comparison.
#define EPS 1.49e-08

// Fortran matrices are stored as vectors in column-major form.
//  Covariance matrices are symmetric so this is not explicitly necessary,
//  but it is good to keep things consistent.
#define COL_MAJOR(i, j, K) (j * K + i)

// A few things to note:
//  - LAPACK requires gfortran which is not a part of CLANG.
//  - Must be compiled with gcc/gfortran.
//  - All matrices in LAPACK are one-dimensional arrays stored in column major
//      form. This does not matter with symmetric matrices.
//  - A is the matrrix to decompose. 
//  - S will contain the singular values.
//  - U will contain the U matrix.
//  - VT will contain V^T. 

// Define out LAPACK SVD routine.
extern void dgesvd_(char* JOBU, char* JOBVT, int* M, int* N, double* A,
                int* LDA, double* S, double* U, int* LDU, double* VT, int* LDVT,
                double* WORK, int* LWORK, int* INFO);

// We wrap the methods in a function to make them more readable in C.
int dgesvd(
    char JOBU, 
    char JOBVT, 
    int M, 
    int N, 
    double* A,
    int LDA, 
    double* S, 
    double* U, 
    int LDU, 
    double* VT, 
    int LDVT,
    double* WORK, 
    int LWORK 
) {
    int INFO;

    dgesvd_(&JOBU, &JOBVT, &M, &N, A,
            &LDA, S, U, &LDU, VT, &LDVT,
            WORK, &LWORK, &INFO);
    
    return INFO;
}

// Fischer-Yates shuffle algorithm to shuffle the rows of a matrix.
//  NOTE: Random number generator needs to be seeded before calling function.
// Accepts:
//  double** matrix -> The matrix to be shuffled.
//  int m -> The number of rows in the matrix.
// Returns: void.
void shuffle_real_matrix(double** matrix, int m) {
    double* temp = NULL;
    int j;
    for (int i = 0; i < m - 1; i++) {
        j = i + rand() / (RAND_MAX / (m - i) + 1);
        temp = matrix[j];
        matrix[j] = matrix[i];
        matrix[i] = temp;
    }
}

// Calculate Xc^TYc.
// Accepts:
//  double** Xc -> The Xc matrix.
//  double** Yc -> The Yc matrix.
//  int N -> The number of samples.
//  int K -> The dimension to project down into.
//  double* C -> The lower triangle to hold Xc^TYc.
//  double* trX -> Sets tr(Xc^TXc).
//  double* trY -> Sets tr(Yc^TYc).
void crossprod(double** Xc, double** Yc, int N, int K, double* C, double* trX, double* trY) {
    *trX = 0;
    *trY = 0;
    // Calculate C = (Yc^T)Xc. Simultaneously calculate the sum(X^2) and the sum(Y^2).
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            C[COL_MAJOR(i, j, K)] = 0;
            for (int k = 0; k < N; k++) {
                C[COL_MAJOR(i, j, K)] += Yc[k][i] * Xc[k][j];
                if (j == 0) {
                    *trX += Xc[k][i] * Xc[k][i];
                    *trY += Yc[k][i] * Yc[k][i]; 
                }
            }
        }
    }
}

// Calculates the trace of the singular values from the eigenvalue decomposition of C.
// Accepts:
//  double* C -> The matrix to decompose.
//  RealSymEigen_t* eigen -> Used for eigen analysis.
//  int N -> The number of rows.
//  int K -> The number of columns.
//  double* covC -> Holds the covariance matrix of C.
// Returns: double, the trace of singular value matrix.
double trace(double* C, RealSymEigen_t* eigen, int N, int K, double* covC) {
    // Variables to calculate singular values when K = 1, 2, or 3.
    double a0, a1, a2, b, c, det, q, r, theta, root1, root2, root3;

    // Calculate covC = (C^T)C.
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            covC[COL_MAJOR(i, j, K)] = 0;
            for (int k = 0; k < K; k++)
                covC[COL_MAJOR(i, j, K)] += C[COL_MAJOR(k, i, K)] * C[COL_MAJOR(k, j, K)];
        }
    }

    double trLambda = 0;

    // Calculate roots of the characteristic equation.
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
    return trLambda;
}

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen_t* eigen, int N, int K, bool transform, bool similarity) {

    // If either set of points were not given, 
    //  return -1. t-statistic can never be negative.
    if (Xc == NULL || Yc == NULL)
        return -1;

    // Trace of (Xc^T)Xc, trace of (Yc^T)Yc.
    double trX = 0, trY = 0;

    // For convience we rename eigen -> Z to C.
    double* C = eigen -> Z;

    crossprod(Xc, Yc, N, K, C, &trX, &trY);

    /*
    printf("C = \n");
    for (int i = 0; i < K; i++) {
        for (int j = 0; j < K; j++) {
            printf("\t%lf", C[COL_MAJOR(i, j, K)]);
        }
        printf("\n");
    }
    printf("\n");

    printf("ctraceX = %lf\n", trX);
    printf("ctraceY = %lf\n", trY);
    printf("\n");
    */

    // Trace of the lambda matrix.
    double rho, trLambda = 0;

    // In the case we do not transform the X coordinates, we just need the trace of the lambda matrix.
    //  When K = 1, 2, 3, we can calculate the eigenvalues of covC directly as the roots of its characteric equation.
    //  We used Numerical Recipes in C, Third Edition to solve the quadratic and cubic equations.
    if (!transform) {
        // We rename eigen -> A to covC to get the singular values.
        double* covC = eigen -> A;

        trLambda = trace(C, eigen, N, K, covC);

        /*
        printf("covC = \n");
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                printf("\t%lf", covC[COL_MAJOR(i, j, K)]);
            }
            printf("\n");
        }
        printf("\n");
        */

        rho = trLambda / trX;

    } else {

        // Calculate SVD on C.
        //  eigen -> W will contain the singular values.
        //  eigen -> A will contain U.
        //  eigen -> auxilary will contain V^T.
        int INFO = dgesvd('A', 'A', K, K, C, K, eigen -> W, eigen -> A, K, eigen -> auxilary, K, eigen -> WORK, 26 * K);

        // If SVD did not converge, return -1. Procrustes statistic can never have this value.
        if (INFO != 0)
            return -1;

        // Compute A^T which is stored in Z.
        double dot;
        for (int i = 0; i < K; i++) {
            // Set b vector to 0 by default.
            eigen -> WORK[i] = 0;
            for (int j = 0; j < K; j++) {
                dot = 0;
                for (int k = 0; k < K; k++) {
                    dot += eigen -> A[COL_MAJOR(i, k, K)] * eigen -> auxilary[COL_MAJOR(k, j, K)];
                }
                eigen -> Z[COL_MAJOR(i, j, K)] = dot;
            }
            trLambda += eigen -> W[i];
        }

        /*
        printf("A^T = \n");
        for (int i = 0; i < K; i++) {
            for (int j = 0; j < K; j++) {
                printf("\t%lf", eigen -> Z[COL_MAJOR(i, j, K)]);
            }
            printf("\n");
        }
        printf("\n");
        */

        // Calculate rho.
        rho = trLambda / trX;

        // Calculate b and store in WORK.
        //  K is usually small so we can seperate
        //  this calculation from the previous loop.
        if (x0 != NULL || y0 != NULL) {
            for (int i = 0; i < K; i++) {
                if (x0 == NULL) {
                    eigen -> WORK[i] = y0[i];
                } else {
                    dot = 0;
                    for (int j = 0; j < K; j++) {
                        dot += eigen -> Z[COL_MAJOR(i, j, K)] * x0[j];
                    }
                    eigen -> WORK[i] = y0[i] - rho * dot;
                }
            }
        }

        /*
        printf("Translation = \n");
        for (int i = 0; i < K; i++)
            printf("\t%lf", eigen -> WORK[i]);
        printf("\n\n");
        */

        // Transform points.
        //  rho * A^T * x + b.
        //  Use eigen -> W as extra space.
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < K; j++) {
                eigen -> W[j] = 0;
                for (int k = 0; k < K; k++) {
                    eigen -> W[j] += eigen -> Z[COL_MAJOR(j, k, K)] * Xc[i][k];
                }
            }
            for (int j = 0; j < K; j++)
                Xc[i][j] = rho * eigen -> W[j] + eigen -> WORK[j];
        }

    }

    // printf("trLambda = %lf rho = %lf\n", trLambda, rho);
    // printf("\n");
    
    double ss = trY + rho * rho * trX - 2 * rho * trLambda;
    // printf("ss = %lf\n", ss);
    // printf("\n");

    // Convert to either similarity or dissimilarity.
    return similarity ? sqrt(1 - ss) : 1 - sqrt(1 - ss);

}

double correlation(double** Xc, double** Yc, RealSymEigen_t* eigen, int n, int k, bool similarity) {
    // For convience we rename eigen -> Z to C.
    double* C = eigen -> Z;
    // Not used.
    double trX, trY;
    // Calculate crossproduct.
    crossprod(Xc, Yc, n, k, C, &trX, &trY);
    // We rename eigen -> A to covC to get the singular values.
    double* covC = eigen -> A;
    // Calculate trace.
    double trLambda = trace(C, eigen, n, k, covC);
    return similarity ? trLambda : (1 - trLambda);
}

double permutation_test(double** Xc, double** Yc, double** shuffleX, RealSymEigen_t* eigen, int n, int k, bool similarity, double t0, int NUM_PERMS) {

    // If either set of points were not given, 
    //  return 0. p-values should never be 0.
    if (Xc == NULL || Yc == NULL)
        return 0;

    // Copy X to shuffleX.
    for (int i = 0; i < n; i++)
        for (int j = 0; j < k; j++)
            shuffleX[i][j] = Xc[i][j];

    double t;
    int numSig = 0;

    // Perform NUM_PERMS permutations.
    for (int i = 0; i < NUM_PERMS; i++) {
        // Shuffle the rows of X.
        shuffle_real_matrix(shuffleX, n);
        
        // Calculate correlation
        t = correlation(shuffleX, Yc, eigen, n, k, similarity);

        // If our statistic exceeds our threshold, increment counter.
        if (t > t0 - EPS) {
            numSig++;
        }
    }

    // Return p-value.
    return ((double) numSig + 1) / (NUM_PERMS + 1.0);
}

// Our record used for multithreaded Procrustes permutation test.
//  We give each thread a chunk of windows to process.
typedef struct {
    Window_t** windows;
    // The start and end index of the windows a thread must process.
    int startWindow;
    int endWindow;
    double** target;
    double* target0;
    int N;
    int K;
    bool similarity;
    int NUM_PERMS;
} ProcrustesRecord_t;

// When the user specifies multiple threads for Procrustes permutation test, each
//  thread executes this method.
// Accepts:
//  void* arg -> Pointer to the created ProcrustesRecord_t*.
// Returns: void.
void* procrustes_permutation_multi_thread(void* arg) {
    ProcrustesRecord_t* record = (ProcrustesRecord_t*) arg;

    // Current window thread is performing Procrustes Analysis on.
    Window_t* window = NULL;
    RealSymEigen_t* eigen = init_real_sym_eigen(record -> K);
    double** shuffleX = create_matrix(double, record -> N, record -> K);
    double t;

    // Iterate through thread's partition.
    for (int i = record -> startWindow; i <= record -> endWindow; i++) {
        // Get current window.
        window = record -> windows[i];
        // LOG_INFO("Performing Procrustes for window %d ...\n", window -> winNum);
        // Calculate Procrustes statistic for window.
        t = procrustes_statistic(window -> X, NULL, record -> target, NULL, eigen, record -> N, record -> K, false, record -> similarity);
        window -> t = t;
        window -> pval = permutation_test(window -> X, record -> target, shuffleX, eigen, record -> N, record -> K, record -> similarity, t, record -> NUM_PERMS);
        // Transform.
        procrustes_statistic(window -> X, NULL, record -> target, record -> target0, eigen, record -> N, record -> K, true, record -> similarity);
    }
    
    destroy_real_sym_eigen(eigen);
    destroy_matrix(double, shuffleX, record -> N);
    free(record);
    return NULL;
}

void procrustes_sliding_window(Window_t** windows, int numWindows, double** target, double* target0, int N, int K, bool similarity, int NUM_PERMS, int NUM_THREADS) {
    
    // If the user did not enter coordinates, perform Procrustes against genome-wide.
    //  Otherwise, we compare the global to the target, entered by the user.
    int startWindow = 0;
    if (windows[0] -> X == target)
        startWindow = 1;

    // If a single thread or we are not executing a permutation test.
    if (NUM_THREADS == 1) {
        RealSymEigen_t* eigen = init_real_sym_eigen(K);
        double t;
        // For each window, perform Procrustes analysis.
        double** shuffleX = create_matrix(double, N, K);
        for (int i = startWindow; i < numWindows; i++) {
            LOG_INFO("Performing Procrustes for window %d ...\n", windows[i] -> winNum);
            t = procrustes_statistic(windows[i] -> X, NULL, target, NULL, eigen, N, K, false, similarity);
            windows[i] -> t = t;
            // Execute permutation test.
            windows[i] -> pval = permutation_test(windows[i] -> X, target, shuffleX, eigen, N, K, similarity, t, NUM_PERMS);
            // Transform.
            procrustes_statistic(windows[i] -> X, NULL, target, target0, eigen, N, K, true, similarity);
        }
        destroy_matrix(double, shuffleX, N);
        destroy_real_sym_eigen(eigen);
        return;
    }

    // Otherwise, we are multithreading.
    int chunkSize = numWindows / NUM_THREADS;

    // Create threads, giving each a chunk of windows.
    pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS - 1; i++) {
        ProcrustesRecord_t* record = (ProcrustesRecord_t*) calloc(1, sizeof(ProcrustesRecord_t));
        record -> windows = windows;
        record -> startWindow = startWindow;
        record -> endWindow = startWindow + chunkSize;
        record -> target = target;
        record -> target0 = target0;
        record -> N = N;
        record -> K = K;
        record -> similarity = similarity;
        record -> NUM_PERMS = NUM_PERMS;
        pthread_create(&threads[i], NULL, procrustes_permutation_multi_thread, (void*) record);
        startWindow = startWindow + chunkSize + 1;
    }
    ProcrustesRecord_t* record = (ProcrustesRecord_t*) calloc(1, sizeof(ProcrustesRecord_t));
    record -> windows = windows;
    record -> startWindow = startWindow;
    record -> endWindow = numWindows - 1;
    record -> target = target;
    record -> N = N;
    record -> K = K;
    record -> similarity = similarity;
    record -> NUM_PERMS = NUM_PERMS;
    procrustes_permutation_multi_thread((void*) record);

    // Join and destroy threads.
    for (int i = 0; i < NUM_THREADS - 1; i++)
        pthread_join(threads[i], NULL);
    free(threads);
}


/*
int main() {

    double** X = create_matrix(double, 3, 2);
    double* x0 = calloc(2, sizeof(double));
    double** Y = create_matrix(double, 3, 2);
    double* y0 = calloc(2, sizeof(double));
    X[0][0] = 0; X[0][1] = 0;
    X[1][0] = 1; X[1][1] = 0;
    X[2][0] = 0; X[2][1] = 1;
    center_matrix(X, x0, 3, 2);

    Y[0][0] = 0; Y[0][1] = 0;
    Y[1][0] = 1; Y[1][1] = 0;
    Y[2][0] = 0; Y[2][1] = sqrt(3);
    center_matrix(Y, y0, 3, 2);

    RealSymEigen_t* eigen = init_real_sym_eigen(2);

    printf("X = \n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            printf("\t%lf", X[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Y = \n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            printf("\t%lf", Y[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    printf("Similarity statistic: %lf\n\n", procrustes_statistic(X, x0, Y, y0, eigen, 3, 2, true, true));

    printf("Transformed X = \n");
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 2; j++) {
            printf("\t%lf", X[i][j]);
        }
        printf("\n");
    }
    printf("\n");

    destroy_matrix(double, X, 3);
    destroy_matrix(double, Y, 3);
    free(x0);
    free(y0);
    destroy_real_sym_eigen(eigen);

    return 0;

}
*/