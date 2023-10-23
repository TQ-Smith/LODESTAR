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

// We define a macro to calculate the number of different alleles between two diploid individuals
//  at a single locus with the "Genotype" encoding from VCFParser.
#define ASD(i, j) (!!(i ^ j) + !(i & j))

#include <iostream>

list<window*>* window_genome(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n, int k, bool useFastMap) {

    list<window*>* windows = new list<window*>;

    Genotype* genotypes = new Genotype[n];

    double** window_haplotype_counts = create_and_fill_real_matrix(0, n, n);
    double** auxilary_haplotype_counts = create_and_fill_real_matrix(0, n, n);
    double** next_window_haplotype_counts = create_and_fill_real_matrix(0, n, n);
    double** global_haplotype_counts = create_and_fill_real_matrix(0, n, n);
    double** temp;

    string chromosome;
    int position; 
    bool isMonomorphic; 
    bool isComplete;

    double* d = new double[n];
    double* e = new double[n];
    bool doesConverge;
    int maxDimReached;


    delete [] genotypes;

    destroy_real_matrix(window_haplotype_counts, n);
    destroy_real_matrix(auxilary_haplotype_counts, n);
    destroy_real_matrix(next_window_haplotype_counts, n);
    destroy_real_matrix(global_haplotype_counts, n);

    return windows;

}

void destroy_window(window** w) {

    // Deallocate points matrix.
    delete [] (*w) -> points;

    // Deallocate structure.
    delete *w;

    // Set pointer to NULL.
    *w = NULL;

}