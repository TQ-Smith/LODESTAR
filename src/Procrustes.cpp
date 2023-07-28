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

// Used for shuffling the rows of a matrix.
//  This does not even need to be imported if we
//  just copy over the one method.
#include "../utils/LinearAlgebra/MatrixOperations.hpp"

// Used for the basic math operatons sqrt, cos, and acos.
#include <cmath>

#include <iostream>
using namespace std;

double procrustes_analysis(double** Xc, double** Yc, double** C, double** CT_C, int n, int k) {

    // Next, we calculate the trace of XTc_Xc and YTc_Yc
    //  We are taking the dot-product of each column
    //  with itself for each matrix. This corresponds
    //  to the diagonal elements.
    double trXTc_Xc = 0;
    double trYTc_Yc = 0;

    for (int i = 0; i < k; i++) {
        for (int j = 0; j < n; j++) {
            trXTc_Xc += (Xc[j][i]) * (Xc[j][i]);
            trYTc_Yc += (Yc[j][i]) * (Yc[j][i]); 
        }
    }

    // Calculate C = YTc_Xc and store in C.
    for (int a = 0; a < k; a++) {
        for (int b = 0; b < k; b++) {
            C[a][b] = 0;
            for (int c = 0; c < n; c++) {
                C[a][b] += (Yc[c][a]) * (Xc[c][b]);
            }
        }
    }

    // Calculate CT_C and store in CT_C.
    for (int a = 0; a < k; a++) {
        for (int b = 0; b < k; b++) {
            CT_C[a][b] = 0;
            for (int c = 0; c < k; c++) {
                CT_C[a][b] += C[c][a] * C[c][b];
            }
        }
    }

    // Calculate the sum of the singular values of C,
    //  or equivalently, the sum of the eigenvalues
    //  for CT_C.
    double trLambda = 0;

    // We solve the characteristic polynomial of the matrix using the quadratic
    //  and cubic formulas for the roots. Code taken from Numerical Recipes
    //  in C, Third Edition.
    switch (k) {
        case 1: {
            trLambda = sqrt(CT_C[0][0]);
            break;
        }
        case 2: {
            // Create the coefficients of the quadratic and the determinant.
            // Our quadratic is in the form of x^2+bx+c.
            double b = -(CT_C[0][0] + CT_C[1][1]);
            double c = (CT_C[0][0] * CT_C[1][1]) - (CT_C[0][1] * CT_C[1][0]);
            double det = sqrt(b * b - 4 * c);

            // Calculate q.
            double q = b < 0 ? -0.5 * (b - det) : -0.5 * (b + det);

            // Calculate trace.
            trLambda += sqrt(q) + sqrt(c / q);
            break;
        }
        case 3: {
            // Our cubic is in the form of ax^3+bx^2+cx+d.
            //  Here, a = 1.

            // Calculate our coefficients.
            double b = -(CT_C[0][0] + CT_C[1][1] + CT_C[2][2]);

            // Sum of minors along the diagonal.
            double c = ((CT_C[1][1] * CT_C[2][2]) - (CT_C[1][2] * CT_C[2][1]))
                        + ((CT_C[0][0] * CT_C[2][2]) - (CT_C[0][2] * CT_C[2][0]))
                        + ((CT_C[0][0] * CT_C[1][1]) - (CT_C[0][1] * CT_C[1][0]));

            // -detA
            double d = -(CT_C[0][0] * ((CT_C[1][1] * CT_C[2][2]) - (CT_C[2][1] * CT_C[1][2])) 
                        - CT_C[0][1] * ((CT_C[1][0] * CT_C[2][2]) - (CT_C[1][2] * CT_C[2][0])) 
                        + CT_C[0][2] * ((CT_C[1][0] * CT_C[2][1]) - (CT_C[1][1] * CT_C[2][0])));

            // Now, we solve the cubic.
            double Q = (3 * c - b * b) / 9.0;
            double R = (9 * b * c - 27 * d - 2 * b * b * b) / 54.0;
            double theta = acos(R / sqrt(-Q * Q * Q));
            cout << sqrt(2 * sqrt(-Q) * cos(theta / 3) - (b / 3)) << endl;
            cout << sqrt(2 * sqrt(-Q) * cos((theta / 3) + ((2 * M_PI) / 3)) - (b / 3)) << endl;
            cout << sqrt(2 * sqrt(-Q) * cos((theta / 3) + ((4 * M_PI) / 3)) - (b / 3)) << endl;
            trLambda += sqrt(2 * sqrt(-Q) * cos(theta / 3) - (b / 3));
            trLambda += sqrt(2 * sqrt(-Q) * cos((theta / 3) + ((2 * M_PI) / 3)) - (b / 3));
            trLambda += sqrt(2 * sqrt(-Q) * cos((theta / 3) + ((4 * M_PI) / 3)) - (b / 3));
            break;
        }
        default:
            // Note: Our default case would use the eigenvalue computations from 
            //          NumericalRecipesInC.hpp.
            break;
    }

    // Return the Procrustes statistic. 
    //  Measures dissimilarity in current state.
    return (trLambda * trLambda) / (trXTc_Xc * trYTc_Yc);

}

double permutation_test(int NUM_PERMUTATIONS, double** Xc, double** Yc, double** shuffleXc, double** C, double** CT_C, int n, int k, double t_0) {

    // Copy contents of Xc to shuffleXc.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            shuffleXc[i][j] = Xc[i][j];
        }
    }

    // The Procrustes statistic for each permutation.
    double D;

    // Keep the count of observations with a value >= t_0.
    int count = 0;

    // Now do the permutation test.
    for (int i = 0; i < NUM_PERMUTATIONS; i++) {

        // Shuffle shuffleXc.
        shuffle_real_matrix(shuffleXc, n);
        
        // Perfrom Procrustes.
        D = procrustes_analysis(shuffleXc, Yc, C, CT_C, n, k);

        // Increment if it was significant.
        if ( t_0 <= sqrt(1 - D) ) {
            count++;
        }

    }

    // Calculate and return p-value.
    return (count + 1.0) / (NUM_PERMUTATIONS + 1);

}

void center_matrix(double** X, double* x_0, int n, int k) {

    // Calculate the center.
    for (int i = 0; i < k; i++) {
        x_0[i] = 0;
        for (int j = 0; j < n; j++) {
            x_0[i] += (X[j][i] / n);
        }
    }

    // Center the set of points.
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            X[i][j] -= x_0[j];
        }
    }

}