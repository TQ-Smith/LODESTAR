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

    Genotype* genotypes = new Genotype[n];

    double** window_haplotype_counts = create_real_matrix(n, n);
    double** auxilary_haplotype_counts = create_real_matrix(n, n);
    double** next_window_haplotype_counts = create_real_matrix(n, n);
    double** global_haplotype_counts = create_real_matrix(n, n);

    string chromosome;
    int position; 
    bool isMonomorphic; 
    bool isComplete;

    double* d = new double[n];
    double* e = new double[n];
    bool doesConverge;
    int maxDimReached;

    string previous_chromosome = "";
    int start_position;
    int start_next_position;
    int previous_position;
    int num_loci_window = 0;
    int num_global_haplotypes = 0;
    int num_global_loci = 0;

    while (parser -> getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes)) {
        
        if (isMonomorphic || !isComplete) {
            continue;
        }

        if (num_loci_window != 0 && (chromosome != previous_chromosome || num_loci_window == (hap_size * window_hap_size))) {
            int num_haps = ceil(num_loci_window / hap_size);
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    next_window_haplotype_counts[i][j] = next_window_haplotype_counts[j][i] = window_haplotype_counts[i][j] - auxilary_haplotype_counts[j][i];
                    window_haplotype_counts[i][j] = window_haplotype_counts[j][i]
                        += (auxilary_haplotype_counts[i][j] > 2 ? 0 : auxilary_haplotype_counts[i][j]);
                    window_haplotype_counts[i][j] = window_haplotype_counts[j][i] /= (2 * num_haps);
                    auxilary_haplotype_counts[j][i] = 0;
                }
            }
            windows -> push_back(create_mds_window(previous_chromosome, start_position, previous_position, num_loci_window, window_haplotype_counts, &maxDimReached, &doesConverge, d, e, n, k, useFastMap));
            double** temp = window_haplotype_counts;
            window_haplotype_counts = next_window_haplotype_counts;
            next_window_haplotype_counts = temp;
            num_loci_window -= hap_size * (num_loci_window / hap_size);
            start_position = window_hap_size == offset_hap_size ? position : start_next_position;
            num_global_haplotypes++;
        }

        if (chromosome != previous_chromosome) {
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    global_haplotype_counts[i][j] = global_haplotype_counts[j][i]
                        += (auxilary_haplotype_counts[i][j] > 2 ? 0 : auxilary_haplotype_counts[i][j]);
                    window_haplotype_counts[i][j] = window_haplotype_counts[j][i] 
                    = auxilary_haplotype_counts[i][j] = auxilary_haplotype_counts[j][i]
                    = next_window_haplotype_counts[i][j] = next_window_haplotype_counts[j][i] = 0;
                }
                window_haplotype_counts[i][i] = auxilary_haplotype_counts[i][i] = next_window_haplotype_counts[i][i] = 0;
            }
            num_loci_window = 0;
            start_next_position = position;
            start_position = position;
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (num_loci_window != 0 && num_loci_window % hap_size == 0) {
                    window_haplotype_counts[i][j] = window_haplotype_counts[j][i] 
                        += (auxilary_haplotype_counts[i][j] > 2 ? 0 : auxilary_haplotype_counts[i][j]);
                    global_haplotype_counts[i][j] = global_haplotype_counts[j][i] 
                        += (auxilary_haplotype_counts[i][j] > 2 ? 0 : auxilary_haplotype_counts[i][j]);
                    if (num_loci_window <= hap_size * offset_hap_size) {
                        auxilary_haplotype_counts[j][i] += window_haplotype_counts[i][j];
                    }
                    auxilary_haplotype_counts[i][j] = 0;
                    num_global_haplotypes++;
                }
                auxilary_haplotype_counts[i][j] += ASD(genotypes[i], genotypes[j]);
            }
        }

        if (num_loci_window == (hap_size * offset_hap_size + 1)) {
            start_next_position = position;
        }
        previous_chromosome = chromosome;
        previous_position = position;
        num_loci_window++;
        num_global_loci++;
    }

    double allele_counts;
    int num_haps = ceil(num_loci_window / hap_size);
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            allele_counts = auxilary_haplotype_counts[i][j] > 2 ? 0 : auxilary_haplotype_counts[i][j];
            window_haplotype_counts[i][j] = window_haplotype_counts[j][i] += allele_counts;
            window_haplotype_counts[i][j] = window_haplotype_counts[j][i] /= num_haps;
            global_haplotype_counts[i][j] = global_haplotype_counts[j][i] += allele_counts;
            global_haplotype_counts[i][j] = global_haplotype_counts[j][i] /= num_global_haplotypes;
        }
    }
    num_global_haplotypes += num_haps;

    windows -> push_back(create_mds_window(previous_chromosome, start_position, previous_position, num_loci_window, window_haplotype_counts, &maxDimReached, &doesConverge, d, e, n, k, useFastMap));
    windows -> push_back(create_mds_window("Global", windows -> front() -> start_position, previous_position, num_global_loci, global_haplotype_counts, &maxDimReached, &doesConverge, d, e, n, k, useFastMap));

    cout << "Number of Global Haplotypes: " << num_global_haplotypes << endl;
    delete [] genotypes;

    destroy_real_matrix(window_haplotype_counts, n);
    destroy_real_matrix(auxilary_haplotype_counts, n);
    destroy_real_matrix(next_window_haplotype_counts, n);
    destroy_real_matrix(global_haplotype_counts, n);

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