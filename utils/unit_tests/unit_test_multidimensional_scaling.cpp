//
// File: MultidimensionalScaling.hpp
// Date: 21 July 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Tests the classical MDS method and the FastMap
//              heuristic.
//

#include "../LinearAlgebra/MatrixOperations.hpp"

#include "../LinearAlgebra/MultidimensionalScaling.hpp"

#include <iostream>

#include <cassert>

using namespace std;

int main() {

    // Our dimensions.
    int n = 3;
    int k = 2;

    // Create our symmetric matrix.
    double** D = create_real_matrix(3, 3);
    D[0][0] = D[1][1] = D[2][2] = 1;
    D[0][1] = D[1][0] = D[1][2] = D[2][1] = 3;
    D[0][2] = D[2][0] = 5;

    cout << "D =" << endl;
    print_real_matrix(D, n, n, 1, 1);
    cout << endl;

    // Allocate auxilary memory.
    double** X = create_real_matrix(3, 2);
    double* d = new double[n];
    double* e = new double[n];
    bool doesConverge;

    // Compute classical mds.
    compute_classical_mds(D, X, d, e, &doesConverge, n, k);

    assert(doesConverge);

    // Print MDS given by R.
    cout << "Classical MDS given by cmdscale" << endl;
    cout << "2.449490e+00  0.4714045" << endl;
    cout << "-3.845925e-16 -0.9428090" << endl;
    cout << "-2.449490e+00  0.4714045" << endl;
    cout << endl;

    // Classical MDS was tested using the cmdscale R function.
    //  Verifying the results are not always as straight-forward
    //  as in this example, since not every MDS implementation
    //  chooses the signs of the eigenvectors in the same way.
    cout << "Classical MDS" << endl;
    cout << "X =" << endl;
    print_real_matrix(X, n, k, 1, 4);
    cout << endl;

    // Free memory.
    destroy_real_matrix(D, n);
    destroy_real_matrix(X, n);
    delete [] d;
    delete [] e;

    // We now test our FastMap algorithm.
    //  The example is taken from tst.small from 
    //  https://www.cs.cmu.edu/~christos/software.html#:~:text=Dimensionality%20reduction-,FastMap%20tar%20file,-In%20C.%20Paper

    // Print FastMap given by paper's implementation.
    cout << "From FastMap Paper" << endl;
    cout << "*********original matrix" << endl;
    cout << "10	20	30	40" << endl;
    cout << "20	25	30	40" << endl;
    cout << "10	20	30	50" << endl;
    cout << "1	2	3	4" << endl;
    cout << endl;
    cout << "*******FASTMAP" << endl;
    cout << "8.06893	4.25942" << endl;
    cout << "4.91152	14.1731" << endl;
    cout << "0	0" << endl;
    cout << "57.0088	1.00266e-15" << endl;
    cout << endl;

    // Reallocate memory and fill D.
    X = create_real_matrix(4, 4);
    X[0][0] = 10; X[0][1] = 20; X[0][2] = 30; X[0][3] = 40;
    X[1][0] = 20; X[1][1] = 25; X[1][2] = 30; X[1][3] = 40;
    X[2][0] = 10; X[2][1] = 20; X[2][2] = 30; X[2][3] = 50;
    X[3][0] = 1; X[3][1] = 2; X[3][2] = 3; X[3][3] = 4;

    // We must calculate the distance matrix from the original.
    D = create_real_matrix(4, 4);
    double dist;
    for(int i = 0; i < 4; i++) {
        for(int j = 0; j < 4; j++) {
            dist = 0;
            for (int k = 0; k < 4; k++) {
                dist += (X[i][k] - X[j][k]) * (X[i][k] - X[j][k]);
            }
            D[i][j] = dist;
        }
    } 

    int maxDimReached;

    // Print distance matrix.
    cout << "We must create the distance matrix from the original matrix." << endl;
    cout << "Distance matrix =" << endl;
    print_real_matrix(D, 4, 4, 1, 1);
    cout << endl;

    // Perform our implementation.
    cout << "FastMap MDS in current implementation" << endl;
    compute_fastmap(D, X, &maxDimReached, 4, 2);
    assert(maxDimReached == 2);
    print_real_matrix(X, 4, 2, 1, 4);
    cout << endl;

}