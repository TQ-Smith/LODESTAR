
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
#include <pthread.h>
#include "Logger.h"
#include "Matrix.h"
MATRIX_INIT(double, double)

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

double procrustes_statistic(double** Xc, double* x0, double** Yc, double* y0, RealSymEigen_t* eigen, int N, int K, bool transform, bool similarity) {

    // If either set of points were not given, 
    //  return -1. t-statistic can never be negative.
    if (Xc == NULL || Yc == NULL)
        return -1;

    // Trace of (Xc^T)Xc, trace of (Yc^T)Yc.
    double trX = 0, trY = 0;

    // For convience we rename eigen -> Z to C.
    double* C = eigen -> Z;

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

    // Trace of the lambda matrix.
    double trLambda = 0;

    // In the case we do not transform the X coordinates, we just need the trace of the lambda matrix.
    //  When K = 1, 2, 3, we can calculate the eigenvalues of covC directly as the roots of its characteric equation.
    //  We used Numerical Recipes in C, Third Edition to solve the quadratic and cubic equations.
    if (!transform) {
        // We rename eigen -> A to covC to get the singular values.
        double* covC = eigen -> A;

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

        // Calculate rho.
        double rho = trLambda / trX;

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
    
    // Calculate statistic.
    //  NOTE: We take the absolute value to avoid the error of small values close to 0 can be stored as -0.
    //          Will have no effect on the results.
    double statistic = abs(1 - (trLambda * trLambda) / (trX * trY));

    // Convert to either similarity or dissimilarity.
    return similarity ? sqrt(1 - statistic) : sqrt(statistic);

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
    int numSig = 1;

    // Perform NUM_PERMS permutations.
    for (int i = 0; i < NUM_PERMS; i++) {
        // Shuffle the rows of X.
        shuffle_real_matrix(shuffleX, n);
        // Compute statistic between the two sets of points. Do not transform.
        t = procrustes_statistic(shuffleX, NULL, Yc, NULL, eigen, n, k, false, similarity);
        
        // If our statistic exceeds our threshold, increment counter.
        if (t > t0) {
            numSig++;
        }
    }

    // Return p-value.
    return ((double) numSig) / (NUM_PERMS + 1.0);
}

// Threading was slowing Procrustes down.
/*
// Our record used for multithreaded Procrustes analysis.
//  All fields are the same as the arguments to procrustes_sliding_window
//  except the curWinIndex.
typedef struct {
    Window_t** windows;
    // The start and end index of the windows a thread must process.
    int startIndex;
    int endIndex;
    double** target;
    double* target0;
    int N;
    int K;
    bool similarity;
    int NUM_PERMS;
} ProcrustesRecord_t;

// When the user specifies multiple threads for Procrustes Analysis, each
//  thread executes this method.
// Accepts:
//  void* arg -> Pointer to the created ProcrustesRecord_t*.
// Returns: void.
void* procrustes_multi_thread(void* arg) {
    ProcrustesRecord_t* record = (ProcrustesRecord_t*) arg;

    // Current window thread is performing Procrustes Analysis on.
    Window_t* window = NULL;
    RealSymEigen_t* eigen = init_real_sym_eigen(record -> K);
    double** shuffleX = create_matrix(double, record -> N, record -> K);
    double t;

    // Iterate through thread's partition.
    for (int i = record -> startIndex; i <= record -> endIndex; i++) {
        // Get current window.
        window = record -> windows[i];
        // LOG_INFO("Performing Procrustes for window %d ...\n", window -> winNum);
        // Calculate Procrustes statistic for window.
        t = procrustes_statistic(window -> X, NULL, record -> target, NULL, eigen, record -> N, record -> K, false, record -> similarity);
        window -> t = t;
        // If the user wants to do a permutation test, execute permutation test.
        if (record -> NUM_PERMS > 0) {
            window -> pval = permutation_test(window -> X, record -> target, shuffleX, eigen, record -> N, record -> K, record -> similarity, t, record -> NUM_PERMS);
        }
    }

    free(record);
    destroy_real_sym_eigen(eigen);
    destroy_matrix(double, shuffleX, record -> N);
    return NULL;
}
*/

void procrustes_sliding_window(Window_t** windows, int numWindows, double** target, int N, int K, bool similarity, int NUM_PERMS) {
    
    RealSymEigen_t* eigen = init_real_sym_eigen(K);
    int startWindow = 0;
    double t;
    // If the user did not enter coordinates, perform Procrustes against genome-wide.
    //  Otherwise, we compare the global to the target, entered by the user.
    if (windows[0] -> X == target)
        startWindow = 1;
    // For each window, perform Procrustes analysis.
    for (int i = startWindow; i < numWindows; i++) {
        LOG_INFO("Performing Procrustes for window %d ...\n", windows[i] -> winNum);
        t = procrustes_statistic(windows[i] -> X, NULL, target, NULL, eigen, N, K, false, similarity);
        windows[i] -> t = t;
        // If the user wants to do a permutation test, execute permutation test.
        if (NUM_PERMS > 0) {
            double** shuffleX = create_matrix(double, N, K);
            windows[i] -> pval = permutation_test(windows[i] -> X, target, shuffleX, eigen, N, K, similarity, t, NUM_PERMS);
            destroy_matrix(double, shuffleX, N);
        }
    }
    destroy_real_sym_eigen(eigen);
    
    /*
    int startIndex;
    if (windows[0] -> X == target)
        startIndex = 1;
    else
        startIndex = 0;

    int chunkSize = numWindows / NUM_THREADS;
    
    // Create threads.
    pthread_t* threads = (pthread_t*) calloc(NUM_THREADS - 1, sizeof(pthread_t));
    for (int i = 0; i < NUM_THREADS - 1; i++) {
        ProcrustesRecord_t* record = (ProcrustesRecord_t*) calloc(1, sizeof(ProcrustesRecord_t));
        record -> windows = windows;
        record -> startIndex = startIndex;
        record -> endIndex = startIndex + chunkSize;
        record -> target = target;
        record -> N = N;
        record -> K = K;
        record -> similarity = similarity;
        record -> NUM_PERMS = NUM_PERMS;
        pthread_create(&threads[i], NULL, procrustes_multi_thread, (void*) record);
        startIndex = startIndex + chunkSize + 1;
    }
    ProcrustesRecord_t* record = (ProcrustesRecord_t*) calloc(1, sizeof(ProcrustesRecord_t));
    record -> windows = windows;
    record -> startIndex = startIndex;
    record -> endIndex = numWindows - 1;
    record -> target = target;
    record -> N = N;
    record -> K = K;
    record -> similarity = similarity;
    record -> NUM_PERMS = NUM_PERMS;
    procrustes_multi_thread((void*) record);
    
    // Join and destroy threads.
    for (int i = 0; i < NUM_THREADS - 1; i++)
        pthread_join(threads[i], NULL);
    free(threads);
    */
}

/*
int main() {
    // Seed random number generator.
    srand(time(NULL));

    int N = 150, K = 2;
    double** X = (double**) calloc(N, sizeof(double*));
    double** Y = (double**) calloc(N, sizeof(double*));
    for (int i = 0; i < N; i++) {
        X[i] = (double*) calloc(K, sizeof(double));
        Y[i] = (double*) calloc(K, sizeof(double));
    }

    FILE* local = fopen("local.tsv", "r");
    FILE* global = fopen("global.tsv", "r");

    double in;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < K; j++) {
            fscanf(local, "%lf", &in);
            X[i][j] = in;
            fscanf(global, "%lf", &in);
            Y[i][j] = in;
        }
    }

    RealSymEigen_t* eigen = init_real_sym_eigen(K);

    printf("%lf\n", procrustes_statistic(X, NULL, Y, NULL, eigen, N, K, true, false));
    printf("X = \n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < K; j++) {
            printf("%lf\t", X[i][j]);
        }
        printf("\n");
    }
    fclose(local);
    fclose(global);

}
*/