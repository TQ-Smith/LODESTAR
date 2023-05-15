//
// File: procrustes.cpp
// Started: 30 Janurary 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Perform Procrustes Analysis from two real, symmetric
//              matrices in dimension k.
//

// Used for input and output.
#include <iostream>

// Used for basic math operatons.
#include <cmath>

// Used for random numbers.
#include <cstdlib>

// Used to seed the random number generator.
#include <ctime>

// Used in QR algorithm to get double EPS.
#include <limits>

// Unless otherwise specified, our streams come from std.
using namespace std;

////////////////////////////// BASIC MATRIX OPERATIONS ///////////////////////////////

// Create a matrix of doubles.
// Accepts:
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, n x m matrix filled with 0s.
double** create_real_matrix(int m, int n) {

    double** matrix = new double*[m];

    // Create each row.
    for (int i = 0; i < m; i++) {
        matrix[i] = new double[n];
    }

    return matrix;

}

// Shuffles the rows of a matrix.
//  Used for permutation testing.
// Note: Assumes random number generator has been seeded. 
// Accepts:
//  double** matrix -> The matrix to shuffle.
//  int m -> The length of the matrix.
// Returns: void.
void shuffle_real_matrix(double** matrix, int m) {

    // A pointer used for swapping.
    double* temp = NULL;

    // Integer used to hold random index.
    int j;

    // Use the Fischer-Yates shuffle algorithm.
    for (int i = m - 1; i > 0; i--) {
        j = (rand() % (i + 1));
        temp = matrix[j];
        matrix[j] = matrix[i];
        matrix[i] = temp;
    }

} 

// Print the contents of an array.
// Accepts:
//  double* array -> The array to print.
//  int n -> The number of elements in the array.
// Returns: void.
void print_real_array(double* array, int n) {
    for (int i = 0; i < n; i++) {
        // Width of 7 and precision of 4.
        printf("%7.4lf ", array[i]);
    }
}

// Print the contents of a matrix.
// Accepts:
//  double** matrix -> The matrix to print.
//  int n -> The length of the matrix.
//  int m -> The width of the matrix.
// Returns: void.
void print_real_matrix(double** matrix, int m, int n) {

    // Iterate through the elements and print row-by-row.
    for (int i = 0; i < m; i++) {
        print_real_array(matrix[i], n);
        printf("\n");
    }

}

// Fills matrix with a value.
// Accepts:
//  double** matrix -> The matrix to print.
//  int value -> The value to fill the matrix with.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void fill_real_matrix(double** matrix, int value, int m, int n) {
    // Iterate through the elements and set value.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = value;
        }
    }
}

// Deep copy a real matrix.
// Accepts:
//  double** matrix -> The matrix to deep copy.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: double**, n x m matrix of copied values from matrix.
double** deep_copy_real_matrix(double** matrix, int m, int n) {

    // Allocate our copy.
    double** copy = create_real_matrix(m, n);

    // Fill copy with values from matrix.
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            copy[i][j] = matrix[i][j]; 
        }
    }

    return copy;

}

// De-allocates a matrix of doubles.
// Accepts:
//  double** matrix -> The matrix to free.
//  int m -> The length of the matrix.
//  int n -> The width of the matrix.
// Returns: void.
void destroy_real_matrix(double** matrix, int m, int n) {
    
    // Delete all the rows.
    for (int i = 0; i < m; i++) {
        delete [] matrix[i];
    }
    // Delete the matrix.
    delete [] matrix;

}

//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////// PRINCIPAL COMPONENT ANALYSIS /////////////////////////////

// Defined as a macro in Numeric Recipes in C.
//  Determine sign of z based on sign of b.
// Accepts:
//  double z -> The first number.
//  double b -> The second number.
// Returns: double, The sign change of z.
double SIGN(double z, double b) {
	return (b >= 0 ? (z >= 0 ? z : -z) : (z >= 0 ? -z : z));
}

// Defined as a macro in Numeric Recipes in C.
//  Squares a number.
// Accepts:
//  double z -> The first number to be squared.
// Returns: double, The square of z.
double SQR(double z) {
	return z * z;
}

// Taken from Numeric Recipes in C.
//  Computes the length of the hypotenuse of 
//  a right angeled triange.
// Accepts:
//  const double a -> The magnitude of the one side.
//  const double b -> The magnitude of the other side.
// Returns: double, The magnitude of the hypotenuse.
double pythag(const double a, const double b) {
	double absa=abs(a), absb=abs(b);
	return (absa > absb ? absa*sqrt(1.0+SQR(absb/absa)) :
		(absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb))));
}

// Taken from Numeric Recipes in C.
// Note: Out of respect, I left as much of this method 
// unchanged as I could.
//  Sorts column eigenvector matrix and correponding
//  eigenvalues in order from greatest to least eigenvalue.
// Accepts:
//  double** v -> The eigenvector matrix.
//  double* d -> The eigenvalues.
// Returns: void.
void eigsrt(double** v, double* d, int n) {
	int k;
	for (int i = 0; i < n-1; i++) {
		double p = d[k = i];
		for (int j = i; j < n; j++) {
			if (d[j] >= p) {
				p = d[k = j];
			}
		}
		if (k != i) {
			d[k]=d[i];
			d[i]=p;
			if (v != NULL) {
				for (int j = 0; j < n; j++) {
					p = v[j][i];
					v[j][i]= v[j][k];
					v[j][k]=p;
				}
			}
		}
	}
}

// Taken from Numeric Recipes in C.
// Note: Out of respect, I left as much of this method 
// unchanged as I could.
//  Converts a real, symmetric matrix to a tridiagonal matrix
//  with the same eigenvectors, eigenvalues.
// Accepts:
//  double** z -> The real, symmetric martix that will be transformed.
//                  z will be transformed to the orthogonal matrix.
//  double* d -> An array to hold the diagonal elements.
//  double* e -> An array to hold the off diagonal elements.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
// Returns: void.
void tred2(double** z, double* d, double* e, int n, bool yesvecs) {
	int l,k,j,i;
	double scale,hh,h,g,f;
	for (i=n-1;i>0;i--) {
		l=i-1;
		h=scale=0.0;
		if (l > 0) {
			for (k=0;k<i;k++)
				scale += abs(z[i][k]);
			if (scale == 0.0)
				e[i]=z[i][l];
			else {
				for (k=0;k<i;k++) {
					z[i][k] /= scale;
					h += z[i][k]*z[i][k];
				}
				f=z[i][l];
				g=(f >= 0.0 ? -sqrt(h) : sqrt(h));
				e[i]=scale*g;
				h -= f*g;
				z[i][l]=f-g;
				f=0.0;
				for (j=0;j<i;j++) {
					if (yesvecs)
						z[j][i]=z[i][j]/h;
					g=0.0;
					for (k=0;k<j+1;k++)
						g += z[j][k]*z[i][k];
					for (k=j+1;k<i;k++)
						g += z[k][j]*z[i][k];
					e[j]=g/h;
					f += e[j]*z[i][j];
				}
				hh=f/(h+h);
				for (j=0;j<i;j++) {
					f=z[i][j];
					e[j]=g=e[j]-hh*f;
					for (k=0;k<j+1;k++)
						z[j][k] -= (f*e[k]+g*z[i][k]);
				}
			}
		} else
			e[i]=z[i][l];
		d[i]=h;
	}
	if (yesvecs) d[0]=0.0;
	e[0]=0.0;
	for (i=0;i<n;i++) {
		if (yesvecs) {
			if (d[i] != 0.0) {
				for (j=0;j<i;j++) {
					g=0.0;
					for (k=0;k<i;k++)
						g += z[i][k]*z[k][j];
					for (k=0;k<i;k++)
						z[k][j] -= g*z[k][i];
				}
			}
			d[i]=z[i][i];
			z[i][i]=1.0;
			for (j=0;j<i;j++) z[j][i]=z[i][j]=0.0;
		} else {
			d[i]=z[i][i];
		}
	}
}

// Taken from Numeric Recipes in C.
// Note: Out of respect, I left as much of this method 
// unchanged as I could.
//  Calculates eigenvalues and/or eigenvectors of a tridiagonal
//  matrix using the QR method. 
// Accepts:
//  double** z -> Column matrix to hold eigenvectors.
//  double* d -> The diagonal elements of the matrix.
//  double* e -> The off diagonal elements of the matrix.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
// Returns: void.
void tqli(double** z, double* d, double* e, int n, bool yesvecs) {
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	const double EPS=numeric_limits<double>::epsilon();
	for (i=1;i<n;i++) e[i-1]=e[i];
	e[n-1]=0.0;
	for (l=0;l<n;l++) {
		iter=0;
		do {
			for (m=l;m<n-1;m++) {
				dd=abs(d[m])+abs(d[m+1]);
				if (abs(e[m]) <= EPS*dd) break;
			}
			if (m != l) {
				if (iter++ == 30) {
					// Change this later. If no convergence, exit.
					printf("Too many iterations in tqli.\nExiting...");
					exit(1);
				}
				g=(d[l+1]-d[l])/(2.0*e[l]);
				r=pythag(g,1.0);
				g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
				s=c=1.0;
				p=0.0;
				for (i=m-1;i>=l;i--) {
					f=s*e[i];
					b=c*e[i];
					e[i+1]=(r=pythag(f,g));
					if (r == 0.0) {
						d[i+1] -= p;
						e[m]=0.0;
						break;
					}
					s=f/r;
					c=g/r;
					g=d[i+1]-p;
					r=(d[i]-g)*s+2.0*c*b;
					d[i+1]=g+(p=s*r);
					g=c*r-b;
					if (yesvecs) {
						for (k=0;k<n;k++) {
							f=z[k][i+1];
							z[k][i+1]=s*z[k][i]+c*f;
							z[k][i]=c*z[k][i]-s*f;
						}
					}
				}
				if (r == 0.0 && i >= l) continue;
				d[l] -= p;
				e[l]=g;
				e[m]=0.0;
			}
		} while (m != l);
	}
}

// Performs PCA on a real, symmetric matrix.
//  In this application, this is a distance matrix.
//  NOTE: This method should be modified to accept the needed temporary memory
//          since allocation for large matrices could be time consuming.
// Accepts:
//  double** X -> Our n x n distance matrix with 1s on the diagonal.
//  int n -> The dimension of n.
//  int k -> The dimension of the reduced data.
//      Assume 0 < k < n.
// Return: double**, A n x k matrix representing the reduction of X
//          into dimension k. 
double** distance_matrix_pca(double** X, int n, int k) {

    // Allocate all of the additional memory.

    // Used to compute Z^TZ, and then, hold the eigenvectors.
    double** eigvecs = create_real_matrix(n, n);

    // d and e used in eigenvector calculation.

    // Holds eigenvalues.
    double* d = new double[n];

    // Used as auxilary vector in tped2 and tqli.
    double* e = new double[n];

    // Used to hold means of columns.
    double* m = new double[n];

    // Allocate our set of points to return.
    double** points = create_real_matrix(n, k);

    // Calculate mean of each column of the matrix.
    for (int i = 0; i < n; i++) {
        m[i] = 0;
        // Calculate mean.
        for (int j = 0; j < n; j++) {
            m[i] += X[j][i];
        }
        m[i] /= n;
    }

    // Calculate Z^TZ.
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < n; b++) {
            eigvecs[a][b] = 0;
            for (int c = 0; c < n; c++) {
                eigvecs[a][b] += (X[c][a] - m[a]) * (X[c][b] - m[b]);
            }
            eigvecs[a][b] /= (n - 1);
        }
    }

    // Now, calculate the eigenvalues.
    tred2(eigvecs, d, e, n, true);
    tqli(eigvecs, d, e, n, true);
    eigsrt(eigvecs, d, n);

    // Calculate projections in dimension k.
    for (int a = 0; a < n; a++) {
        for (int b = 0; b < k; b++) {
            points[a][b] = 0;
            for (int c = 0; c < n; c++) {
                points[a][b] += (X[a][c] - m[c]) * eigvecs[c][b];
            }
        }
    }

    // Free all allocated memory.
    destroy_real_matrix(eigvecs, n, n);
    delete [] d;
    delete [] e;
    delete [] m;

    // Return the points.
    return points;

}

//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////// PROCRUSTES ANALYSIS //////////////////////////////////

// Perform Procrustes Analysis on two sets of points.
// Each set of points are stored in a n x k matrix where
//  there are n k-dimensional points. We assume there are 
//  no duplicate points.
// Accepts:
//  double** Xc -> The first centered set of points.
//  double** Yc -> The second centered set of points.
//  double** temp1 -> A k x k matrix used as a temp variable.
//  double** temp2 -> A k x k matrix used as a temp variable.
//  int n -> The number of points.
//  int k -> The dimension of each point.
//      Assumes 0 < k < n.
// Returns: double, The Procrustes statistic.
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
        tred2(temp2, temp1[0], temp1[1], k, false);
        tqli(temp2, temp1[0], temp1[1], k, false);

        // Now, we add up the singular values.
        //  They are stored in temp1[0].
        for (int i = 0; i < k; i++) {
            trLambda += sqrt(temp1[0][i]);
        }

    }

    // Finally, we compute D(X, Y).
    return  1 - (trLambda * trLambda) / (trX * trY);
}

// A wrapper function for procrustes_statistic that 
//  preforms procrustes analysis between two sets
//  of points.
// Accepts:
//  double** X -> The first set of points.
//  double** Y -> The second set of points.
//  int n -> The number of points.
//  int k -> The dimension of each point.
//      Assumes 0 < k < n.
// Returns: double, The Procrustes statistic.
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

//////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////// PERMUTATION TESTING /////////////////////////////////

// Performs a permutation test for the Procrustes statistic.
// Accepts:
//  int NUM_PERMUTATIONS -> The number of permutations to run.
//  double** X -> The first set of points.
//  double** Y -> The second set of points.
//  int n -> The number of points.
//  int k -> The dimension of each point.
//  double* t -> Sets the observed initial procrustes statistic, our t_0.
//  double* p_value -> Sets the p_value given by the permutation test.
//      Assumes 0 < k < n.
// Returns: void.
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

//////////////////////////////////////////////////////////////////////////////////////


// Our main method unit tests Procrustes Analysis.
int main() {

    // Seed the generator.
    srand((unsigned int) time(NULL));

    /*

    // Test 5 x 5 distance matrices.
    int n = 5;

    // Create matrices.
    double** A = create_real_matrix(n, n);
    double** B = create_real_matrix(n, n);

    // Fill A.
    A[0][0] = 1; A[0][1] = 2; A[0][2] = 3; A[0][3] = 4; A[0][4] = 5;
    A[1][0] = 2; A[1][1] = 1; A[1][2] = 2; A[1][3] = 3; A[1][4] = 4;
    A[2][0] = 3; A[2][1] = 2; A[2][2] = 1; A[2][3] = 2; A[2][4] = 3;
    A[3][0] = 4; A[3][1] = 3; A[3][2] = 2; A[3][3] = 1; A[3][4] = 2;
    A[4][0] = 5; A[4][1] = 4; A[4][2] = 3; A[4][3] = 2; A[4][4] = 1;

    // Fill B.
    B[0][0] = 1; B[0][1] = 2; B[0][2] = 6; B[0][3] = 9; B[0][4] = 11;
    B[1][0] = 2; B[1][1] = 1; B[1][2] = 3; B[1][3] = 7; B[1][4] = 10;
    B[2][0] = 6; B[2][1] = 3; B[2][2] = 1; B[2][3] = 4; B[2][4] = 8;
    B[3][0] = 9; B[3][1] = 7; B[3][2] = 4; B[3][3] = 1; B[3][4] = 5;
    B[4][0] = 11; B[4][1] = 10; B[4][2] = 8; B[4][3] = 5; B[4][4] = 1;

    // Print the contents of our matrices for convience.
    printf("\n");
    printf("Our A matrix:\n");
    print_real_matrix(A, n, n);
    printf("\n");

    printf("\n");
    printf("Our B matrix:\n");
    print_real_matrix(B, n, n);
    printf("\n");

    // Compute PCA for A and B for k = 1, 2, 3, 4
    double** a1 = distance_matrix_pca(A, n, 1);
    double** a2 = distance_matrix_pca(A, n, 2);
    double** a3 = distance_matrix_pca(A, n, 3);
    double** a4 = distance_matrix_pca(A, n, 4);
    double** b1 = distance_matrix_pca(B, n, 1);
    double** b2 = distance_matrix_pca(B, n, 2);
    double** b3 = distance_matrix_pca(B, n, 3);
    double** b4 = distance_matrix_pca(B, n, 4);

    // Print PCA for a_k and b_k for k = 1, 2, 3, 4
    printf("\n");
    printf("PCA on A with k = 1:\n");
    print_real_matrix(a1, n, 1);
    printf("\n");

    
    printf("\n");
    printf("PCA on A with k = 2:\n");
    print_real_matrix(a2, n, 2);
    printf("\n");

    printf("\n");
    printf("PCA on A with k = 3:\n");
    print_real_matrix(a3, n, 3);
    printf("\n");

    printf("\n");
    printf("PCA on A with k = 4:\n");
    print_real_matrix(a4, n, 4);
    printf("\n");

    printf("\n");
    printf("PCA on B with k = 1:\n");
    print_real_matrix(b1, n, 1);
    printf("\n");

    printf("\n");
    printf("PCA on B with k = 2:\n");
    print_real_matrix(b2, n, 2);
    printf("\n");

    printf("\n");
    printf("PCA on B with k = 3:\n");
    print_real_matrix(b3, n, 3);
    printf("\n");

    printf("\n");
    printf("PCA on B with k = 4:\n");
    print_real_matrix(b4, n, 4);
    printf("\n");

    double p1 = procrustes_analysis(a1, b1, n, 1);
    double p2 = procrustes_analysis(a2, b2, n, 2);
    double p3 = procrustes_analysis(a3, b3, n, 3);
    double p4 = procrustes_analysis(a4, b4, n, 4);

    printf("The Procrustes statistic between A and B with k = 1:\n");
    printf("%lf\n", p1);
    printf("\n");

    printf("The Procrustes statistic between A and B with k = 2:\n");
    printf("%lf\n", p2);
    printf("\n");

    printf("The Procrustes statistic between A and B with k = 3:\n");
    printf("%lf\n", p3);
    printf("\n");

    printf("The Procrustes statistic between A and B with k = 4:\n");
    printf("%lf\n", p4);
    printf("\n");
    
    // Destroy all free memory.
    destroy_real_matrix(A, n, n);
    destroy_real_matrix(B, n, n);
    destroy_real_matrix(a1, n, 1);
    destroy_real_matrix(a2, n, 2);
    destroy_real_matrix(a3, n, 3);
    destroy_real_matrix(a4, n, 4);
    destroy_real_matrix(b1, n, 1);
    destroy_real_matrix(b2, n, 2);
    destroy_real_matrix(b3, n, 3);
    destroy_real_matrix(b4, n, 4);

    */

    int n = 500;

    int k = 2;

    double** A = create_real_matrix(n, n);
    double** B = create_real_matrix(n, n);

    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            if (i == j) {
                A[i][j] = 1;
                B[i][j] = 1;
            }
            A[i][j] = A[j][i] = (rand() % 10);
            B[i][j] = B[j][i] = (rand() % 10);
        }
    }

    double** X = distance_matrix_pca(A, n, k);
    double** Y = distance_matrix_pca(B, n, k);

    double t, p;

    permutation_test(1000, X, Y, n, k, &t, &p);

    printf("P(t > t_0 = %lf) = %lf\n", t, p);

    destroy_real_matrix(A, n, n);
    destroy_real_matrix(B, n, n);
    destroy_real_matrix(X, n, k);
    destroy_real_matrix(Y, n, k);

}