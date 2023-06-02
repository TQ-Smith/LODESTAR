//
// File: Procrustes.cpp
// Started: 16 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Perform Procrustes Analysis from two real, symmetric
//              matrices in dimension k.
//

#include "Procrustes.hpp"

#include "MatrixOperations.hpp"

#include "NumericalRecipesInC.hpp"

// Used for basic math operatons.
#include <cmath>

// Used in QR algorithm to get double EPS.
#include <limits>

double procrustes_statistic(double** Xc, double** Yc, double** temp1, double** temp2, int n, int k) {


    // Next, we calculate the trace of XT_cX_c and YT_cY_c
    //  as described in the paper.
    double trX = 0;
    double trY = 0;

    // The trace is only the diagonal.
    //  We are taking the dot-product of each column
    //  with its opposing column.
    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            trX += (Xc[j][i]) * (Xc[j][i]);
            trY += (Yc[j][i]) * (Yc[j][i]); 
        }
    }

    // The hard part. Calculating tr(Lambda).
    double trLambda = 0;

    // We use our temp matrices in the calculation of 
    //  of C^TC for our SVD. We do not use a matrix 
    //  multiplicaiton routine because this is faster
    //  than generating transposes.

    // Used as a place holder in matrix multiplication.
    double dot;

    // First, calculate C = Y_c^TX_c and store in temp1.
    for (int a = 0; a < k; a++) {
        for (int b = 0; b < k; b++) {
            dot = 0;
            for (int c = 0; c < n; c++) {
                dot += (Yc[c][a]) * (Xc[c][b]);
            }
            temp1[a][b] = dot;
        }
    }

    // Next, we calculate C^TC.
    for (int a = 0; a < k; a++) {
        for (int b = 0; b < k; b++) {
            temp2[a][b] = 0;
            for (int c = 0; c < k; c++) {
                temp2[a][b] += temp1[c][a] * temp1[c][b];
            }
        }
    }

    // Finally, we need to eigenvalues of
    //  C^TC. If k = 1, 2, 3, this is easy
    //  otherwise we use the QR algorithm.

    // If C is a 1 x 1, then that is the singular value.
    if (k == 1) {

        trLambda = sqrt(temp2[0][0]);

    // If C is a 2 x 2, then we directly solve with the quadratic formula.
    } else if (k == 2) {

        // From Numeric Recipes in C.

        // Create the coefficients of the quadratic and the determinant.
        // Our quadratic is in the form of x^2+bx+c.
        double b = -(temp2[0][0] + temp2[1][1]);
        double c = (temp2[0][0] * temp2[1][1]) - (temp2[0][1] * temp2[1][0]);
        double det = sqrt(b * b - 4 * c);

        // Calculate q.
        double q = b < 0 ? -0.5 * (b - det) : -0.5 * (b + det);

        // Calculate trace.
        trLambda += sqrt(q) + sqrt(c / q);
    
    // If C is a 3 x 3, then we directly solve with the cubic formula
    // by the method of Francois Viete.
    } else if (k == 3) {

        // From Numeric Recipes in C.

        // Our cubic is in the form of x^3+ax^2+bx+c.
        // We are using the trace formula for the characteristic polynomial.

        // Calculate our coefficients.
        double a = -(temp2[0][0] + temp2[1][1] + temp2[2][2]);

        // Sum of minors along the diagonal.
        double b = ((temp2[1][1] * temp2[2][2]) - (temp2[1][2] * temp2[2][1]))
                    + ((temp2[0][0] * temp2[2][2]) - (temp2[0][2] * temp2[2][0]))
                    + ((temp2[0][0] * temp2[1][1]) - (temp2[0][1] * temp2[1][0]));

        // -detA
        double c = -(temp2[0][0] * ((temp2[1][1] * temp2[2][2]) - (temp2[2][1] * temp2[1][2])) 
                    - temp2[0][1] * ((temp2[1][0] * temp2[2][2]) - (temp2[1][2] * temp2[2][0])) 
                    + temp2[0][2] * ((temp2[1][0] * temp2[2][1]) - (temp2[1][1] * temp2[2][0])));

        // Now, we solve the cubic.
        double Q = (a * a - 3 * b) / 9;
        double R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
        double theta = acos(R / sqrt(Q * Q * Q));

        // Now calculate the sume of the roots.
        trLambda += sqrt(-2 * sqrt(Q) * cos(theta / 3) - (a/3));
        trLambda += sqrt(-2 * sqrt(Q) * cos((theta + 2 * M_PI) / 3) - (a/3));
        trLambda += sqrt(-2 * sqrt(Q) * cos((theta - 2 * M_PI) / 3) - (a/3));
        //printf("%lf\n", -2 * sqrt(Q) * cos(theta / 3) - (a/3));
        //printf("%lf\n", -2 * sqrt(Q) * cos((theta + 2 * M_PI) / 3) - (a/3));
        //printf("%lf\n", -2 * sqrt(Q) * cos((theta - 2 * M_PI) / 3) - (a/3));
        //printf("%lf\n", trLambda);

    } else {

        // We can clobber temp1 and use the rows of temp2 as d and e.
        //  We do not care about the eigenvectors or sorting.
        compute_eigen_pairs(temp2, temp1[0], temp1[1], k, false, false);

        // Now, we add up the singular values.
        //  They are stored in temp1[0].
        for (int i = 0; i < k; i++) {
            trLambda += sqrt(temp1[0][i]);
        }

    }

    // Finally, we compute D(X, Y).
    return  1 - (trLambda * trLambda) / (trX * trY);
}

double procrustes_analysis(double** X, double** Y, int n, int k) {

    // Create our two temp matrices.
    double** temp1 = create_real_matrix(k, k);
    double** temp2 = create_real_matrix(k, k);

    // We first calculate x_0 and y_0.
    //  Each are of dimension k.
    double* x_0 = new double[k];
    double* y_0 = new double[k];

    // Note, the centering of the values could
    //  be its own function. But, we can center both
    //  sets together, and the code is only used twice.
    //  For those two reasons, we keep everything as is.

    // Set each element equal to 0.
    for (int i = 0; i < k; i++) {
        x_0[i] = 0;
        y_0[i] = 0;
    }
    // Now, calculate the center.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            x_0[j] += (X[i][j] / n);
            y_0[j] += (Y[i][j] / n);
        }
    }
    // Center the sets of points.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] -= x_0[j];
            Y[i][j] -= y_0[j];
        }
    }

    // Calculate our statisitc.
    double stat = procrustes_statistic(X, Y, temp1, temp2, n, k);

    // Uncenter the points.
    //  This is only needed if we want X and Y unchanged.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] += x_0[j];
            Y[i][j] += y_0[j];
        }
    }

    // Deallocate our temp matrices and center points.
    destroy_real_matrix(temp1, k, k);
    destroy_real_matrix(temp2, k, k);

    delete [] x_0;
    delete [] y_0;

    // Return our statistic.
    return stat;

}

void permutation_test(int NUM_PERMUTATIONS, double** X, double** Y, int n, int k, double* t, double* p_value) {

    // First, we allocate all of the necessary memory and variables.

    // Create our temp matrices.
    double** temp1 = create_real_matrix(k, k);
    double** temp2 = create_real_matrix(k, k);

    // Create our arrays to hold the point centers.
    double* x_0 = new double[k];
    double* y_0 = new double[k];

    // Set each element equal to 0.
    for (int i = 0; i < k; i++) {
        x_0[i] = 0;
        y_0[i] = 0;
    }
    // Now, calculate the center.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            x_0[j] += (X[i][j] / n);
            y_0[j] += (Y[i][j] / n);
        }
    }
    // Center the sets of points.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] -= x_0[j];
            Y[i][j] -= y_0[j];
        }
    }

    // We need a copy of X to permute.
    double** permuteX = deep_copy_real_matrix(X, n, k);

    // Calculate our initial t value.
    double D = procrustes_statistic(X, Y, temp1, temp2, n, k);

    // Set our initial t.
    double t_0 = sqrt(1 - D);

    // Keep the count of observations with a value >= t_0.
    int count = 0;

    // Now do the permutation test.
    for (int i = 0; i < NUM_PERMUTATIONS; i++) {

        // Shuffle permuteX.
        shuffle_real_matrix(permuteX, n);
        
        // Do the analysis.
        D = procrustes_statistic(permuteX, Y, temp1, temp2, n, k);

        // See if it was significant.
        if ( t_0 <= sqrt(1 - D) ) {
            count++;
        }

    }

    // Uncenter the points.
    //  This is only needed if we want X and Y unchanged.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] += x_0[j];
            Y[i][j] += y_0[j];
        }
    }

    // Free all the used memory.
    destroy_real_matrix(temp1, k, k);
    destroy_real_matrix(temp2, k, k);
    destroy_real_matrix(permuteX, n, k);
    delete [] x_0;
    delete [] y_0;

    // Set our initial t value.
    *t = t_0;

    // Finally, we set our p_value.
    *p_value = (count + 1.0) / (NUM_PERMUTATIONS + 1);

}

void permutation_test(int NUM_PERMUTATIONS, double** X, double** Y, double** temp1, double** temp2, double* x_0, double* y_0, int n, int k, double* t, double* p_value) {

    // First, we allocate all of the necessary memory and variables.

    // Set each element equal to 0.
    for (int i = 0; i < k; i++) {
        x_0[i] = 0;
        y_0[i] = 0;
    }
    // Now, calculate the center.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            x_0[j] += (X[i][j] / n);
            y_0[j] += (Y[i][j] / n);
        }
    }
    // Center the sets of points.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] -= x_0[j];
            Y[i][j] -= y_0[j];
        }
    }

    // We need a copy of X to permute.
    double** permuteX = deep_copy_real_matrix(X, n, k);

    // Calculate our initial t value.
    double D = procrustes_statistic(X, Y, temp1, temp2, n, k);

    // Set our initial t.
    double t_0 = sqrt(1 - D);

    // Keep the count of observations with a value >= t_0.
    int count = 0;

    // Now do the permutation test.
    for (int i = 0; i < NUM_PERMUTATIONS; i++) {

        // Shuffle permuteX.
        shuffle_real_matrix(permuteX, n);
        
        // Do the analysis.
        D = procrustes_statistic(permuteX, Y, temp1, temp2, n, k);

        // See if it was significant.
        if ( t_0 <= sqrt(1 - D) ) {
            count++;
        }

    }

    // Uncenter the points.
    //  This is only needed if we want X and Y unchanged.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] += x_0[j];
            Y[i][j] += y_0[j];
        }
    }

    // Free all the used memory.
    destroy_real_matrix(permuteX, n, k);

    // Set our initial t value.
    *t = t_0;

    // Finally, we set our p_value.
    *p_value = (count + 1.0) / (NUM_PERMUTATIONS + 1);

}