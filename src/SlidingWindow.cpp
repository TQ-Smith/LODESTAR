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

#include <iostream>
using namespace std;

int HAPLOTYPE_ASD(int asd, int left_a, int right_a, int left_b, int right_b) {
    if (asd == 2) {
        return 2;
    } else if (!(left_a ^ right_a ^ left_b ^ right_b)) {
        return 0;
    } else if (left_a != left_b && left_a != right_b && right_a != left_b && right_a != right_b ) {
        return 2;
    } else {
        return 1;
    }
}

// A helper function used to allocate a window structure, set fields, and perfrom MDS.
//  NOTE: Will change with threading.
// Accepts:
//  string chrom -> The chromosome the window lies on.
//  int start_pos -> The starting position of the window.
//  int end_pos -> The ending position of the window.
//  int num_loci -> The number of loci within the window.
//  int n -> The number of samples.
//  int k -> The dimension to project into.
//  bool useFastMap -> Flag to indicate use of FastMap.
//  double** D -> Our distance matrix between the samples.
//  double* d -> An auxilary array used by classical MDS, holds eigenvalues. NULL if FastMap is used. -| -> If this function is used for threading,
//  double* e -> Another auxilary array used by classical MDS. NULL if FastMap is used. ---------------| -> each thread must have their own d and e. 
// Returns: window*, The pointer to the newly allocated window structure.
window* create_window_with_mds(string chrom, int start_pos, int end_pos, int num_loci, int n, int k, bool useFastMap, double** D, double* d, double* e ) {
    
    // Used for classical MDS to test convergence.
    bool doesConverge;

    // Used by FastMap to get the maximum dimension reached by the algorithm.
    int maxDimReached = k;

    // Allocate matrix to hold dimension reduced samples.
    double** points = create_real_matrix(n, k);

    // Perform MDS.
    if (useFastMap) {
        compute_fastmap(D, points, &maxDimReached, n, k);
    } else {
        compute_classical_mds(D, points, d, e, &doesConverge, n, k);
    }

    // Set fields of window.
    window* w = new window;
    w -> chromosome = chrom;
    w -> start_position = start_pos;
    w -> end_position = end_pos;
    w -> num_loci = num_loci;

    // If divergence of cMDS or k was not reached in FastMap, do not save the points.
    //  Ask Zach if this is okay.
    if (maxDimReached != k || !doesConverge) {
        destroy_real_matrix(points, n);
        w -> points = NULL;
    } else {
        w -> points = points;
    }
    
    return w;

}

list<window*>* window_genome(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n, int k, bool useFastMap) {
    
    list<window*>* windows = new list<window*>;

    
    return windows;

}

void destroy_window(window** w, int n) {

    // Deallocate points matrix, if it exists.
    if ((*w) -> points != NULL) {
        destroy_real_matrix((*w) -> points, n);
    }

    // Deallocate structure.
    delete *w;

    // Set pointer to NULL.
    *w = NULL;

}