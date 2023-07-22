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
    cout << "Classical MDS" << endl;
    cout << "X =" << endl;
    print_real_matrix(X, n, k, 1, 1);
    cout << endl;

    destroy_real_matrix(D, n);
    destroy_real_matrix(X, n);
    delete [] d;
    delete [] e;

    // We now test our FastMap algorithm.
    //  The example is taken from tst.small from 
    //  https://www.cs.cmu.edu/~christos/software.html#:~:text=Dimensionality%20reduction-,FastMap%20tar%20file,-In%20C.%20Paper

    cout << "From FastMap Paper:"
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

    cout << "";

}