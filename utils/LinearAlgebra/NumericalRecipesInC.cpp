//
// File: NumericalRecipesInC.cpp
// Date: 20 July 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Contains operations from Numerical Recipes in C, Third Edition.
//              Functions used to compute the eigenvalues and eigenvectors 
//              of real symmetric matrices are included. I do not unit test 
//				these functions because I made minimal changes.
//

#include "NumericalRecipesInC.hpp"

// Used for basic math operatons.
#include <cmath>

// Used in QR algorithm to get double EPS.
#include <limits>

// Determine sign of z based on sign of b.
// Accepts:
//  double z -> The first number.
//  double b -> The second number.
// Returns: double, The sign change of z.
#define SIGN(z, b) (b >= 0 ? (z >= 0 ? z : -z) : (z >= 0 ? -z : z))

// Squares a number.
// Accepts:
//  double z -> The first number to be squared.
// Returns: double, The square of z.
#define SQR(z) (z * z)

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
//  Calculates eigenvalues and/or eigenvectors of a tridiagonal
//  matrix using the QR method. 
// Accepts:
//  double** z -> Column matrix to hold eigenvectors.
//  double* d -> The diagonal elements of the matrix.
//  double* e -> The off diagonal elements of the matrix.
//  bool* doesConverge -> Sets boolean if the algorithm converges.
//  int n -> The dimension of z.
//  bool yesvecs -> Flag set to calculate eigenvectors.
// Returns: void.
void tqli(double** z, double* d, double* e, bool* doesConverge, int n, bool yesvecs) {
	int m,l,iter,i,k;
	double s,r,p,g,f,dd,c,b;
	const double EPS=std::numeric_limits<double>::epsilon();
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
                    // The algorithm does not converge.
                    //  Set flag and return.
					*doesConverge = false;
                    return;
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
    // Algorithm converges.
    *doesConverge = true;
}

void compute_eigen_pairs(double** z, double* d, double* e, bool* doesConverge, int n, bool yesvecs, bool sort) {
	tred2(z, d, e, n, yesvecs);
    tqli(z, d, e, doesConverge, n, yesvecs);
	if (sort)
		eigsrt(z, d, n);
}