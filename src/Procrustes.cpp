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

// Used for creating, destroying, deep copying, and shuffling matrices.
#include "MatrixOperations.hpp"

// Used for the basic math operatons sqrt, cos, and acos.
#include <cmath>

// Our helper function for Procrustes analysis.
// Accepts:
//  double** Xc -> The first centered set of points.
//  double** Yc -> The second centered set of points.
//  double** C -> A k x k matrix used to hold the value of YcT_Xc.
//  double** CT_C -> A k x k matrix used to hold the value of CT_C.
//  int n -> The number of points.
//  int k -> The dimension of each point. Assume k = 1, 2, 3.
// Returns: double, The Procrustes statistic.
double procrustes(double** Xc, double** Yc, double** C, double** CT_C, int n, int k) {

    // Next, we calculate the trace of XTc_Xc and YTc_Yc
    //  We are taking the dot-product of each column
    //  with itself for each matrix. This corresponds
    //  to the diagonal elements.
    double trXTc_Xc = 0;
    double trY = 0;

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
        case 1:
            trLambda = sqrt(CT_C[0][0]);
            break;

        case 2:
            // Create the coefficients of the quadratic and the determinant.
            // Our quadratic is in the form of x^2+bx+c.
            double b = -(temp2[0][0] + temp2[1][1]);
            double c = (temp2[0][0] * temp2[1][1]) - (temp2[0][1] * temp2[1][0]);
            double det = sqrt(b * b - 4 * c);

            // Calculate q.
            double q = b < 0 ? -0.5 * (b - det) : -0.5 * (b + det);

            // Calculate trace.
            trLambda += sqrt(q) + sqrt(c / q);
            break;

        case 3:
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
            break;

        default:
            // Note: Our default case would use the eigenvalue computations from 
            //          NumericalRecipesInC.hpp.
            break;
    }

    // Return the Procrustes statisitc.
    return (trLambda * trLambda) / (trXTc_Xc * trYTc_Yc);

}