//
// File: MultidimensionalScaling.hpp
// Started: 21 July 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Tests the classical MDS method and the FastMap
//              heuristic.
//

#include "../utils/LinearAlgebra/MatrixOperations.hpp"

#include "../utils/LinearAlgebra/MultidimensionalScaling.hpp"

#include <iostream>

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

    // Classical MDS was tested using the cmdscale R function.
    //  Verifying the results are not always as straight-forward
    //  as in this example, since not every MDS implementation
    //  chooses the signs of the eigenvectors in the same way.
    cout << "X =" << endl;
    print_real_matrix(X, n, k, 1, 1);
    cout << endl;

    destroy_real_matrix(D, n);
    destroy_real_matrix(X, n);
    delete [] d;
    delete [] e;

}