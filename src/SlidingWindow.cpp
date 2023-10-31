//
// File: SlidingWindow.cpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines the main algorithm for a sliding window along the genome.
//

#include "SlidingWindow.hpp"

// Use classical MDS and FastMap.
#include "../utils/LinearAlgebra/MultidimensionalScaling.hpp"

// Used to create and destroy matrices.
#include "../utils/LinearAlgebra/MatrixOperations.hpp"

// Used for the ceiling function.
#include <cmath>

#define LOCUS_ASD(i, j) (!!(i ^ j) + !(i & j))

#include <iostream>
using namespace std;

window* create_mds_window(
    string chromosome, int start_position, int end_position, int num_loci, 
    double** D, int* maxDimReached, bool* doesConverge, double* d, double* e, 
    int n, int k, bool useFastMap) {
    
    double** points = create_real_matrix(n, k);
    if (useFastMap) {
        compute_fastmap(D, points, maxDimReached, n, k);
    } else {
        compute_classical_mds(D, points, d, e, doesConverge, n, k);
    }
    window* w = new window;
    w -> chromosome = chromosome;
    w -> start_position = start_position;
    w -> end_position = end_position;
    w -> num_loci = num_loci;
    w -> points = points;
    return w;

}

list<window*>* window_genome(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n, int k, bool useFastMap) {
    
    list<window*>* windows = new list<window*>;

    return windows;

}

void destroy_window(window** w, int n) {

    // Deallocate points matrix.
    destroy_real_matrix((*w) -> points, n);

    // Deallocate structure.
    delete *w;

    // Set pointer to NULL.
    *w = NULL;

}