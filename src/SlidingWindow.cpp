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

    Genotype* genotypes = new Genotype[n];
    string chrom;
    int pos;
    bool isMonomorphic, isComplete, nextRecord;

    double* d = NULL; double* e = NULL;
    if (!useFastMap) {
        d = new double[n]; e = new double[n];
    }

    string prev_chrom;
    int start_pos, next_start_pos, prev_pos, num_haps, num_loci = 0, num_global_haps = 0;

    double** window_asd = new double*[n]; double** global_asd = new double*[n];  
    double** next_asd = new double*[n]; int** auxilary_asd = new int*[n]; 
    for (int i = 0; i < n; i++) {
        window_asd[i] = new double[n]; next_asd[i] = new double[n];
        global_asd[i] = new double[n]; auxilary_asd[i] = new int[n];
        for (int j = 0; j < n; j++) {
            window_asd[i][j] = next_asd[i][j] = auxilary_asd[i][j] = global_asd[i][j] = 0;
        }
    }

    while ( true ) {

        nextRecord = parser -> getNextLocus(&chrom, &pos, &isMonomorphic, &isComplete, genotypes);

        if (nextRecord && (isMonomorphic || !isComplete)) {
            continue;
        }

        if (num_loci != 0 && (!nextRecord || chrom != prev_chrom || num_loci == (hap_size * window_hap_size))) {
            num_haps = ceil(num_loci / (double) hap_size);
            num_global_haps += offset_hap_size;
            if (chrom != prev_chrom || !nextRecord) {
                num_global_haps += ceil((num_loci - hap_size * offset_hap_size) / (double) hap_size);
            }

            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    next_asd[i][j] = next_asd[j][i] = window_asd[i][j] - auxilary_asd[j][i];
                    window_asd[i][j] = window_asd[j][i] = (window_asd[i][j] + auxilary_asd[i][j]) / (2.0 * num_haps);
                    if (!nextRecord) {
                        if (chrom == prev_chrom && num_loci != (hap_size * window_hap_size)) {
                            global_asd[i][j] = global_asd[j][i] += auxilary_asd[i][j];
                        }
                        global_asd[i][j] = global_asd[j][i] = global_asd[i][j] / (2.0 * num_global_haps);
                    }
                    auxilary_asd[j][i] = 0;
                }
                window_asd[i][i] = next_asd[i][i] = 0;
            }
            
            windows -> push_back(create_window_with_mds(prev_chrom, start_pos, prev_pos, num_loci, n, k, useFastMap, window_asd, d, e));

            num_loci -= (num_loci / hap_size) * hap_size;
            start_pos = (window_hap_size == offset_hap_size) ? pos : next_start_pos;
        }

        if (!nextRecord) {
            break;
        }

        if (prev_chrom != chrom) {
            num_loci = 0;
            next_start_pos = pos;
            start_pos = pos;
        }

        num_loci++;

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (prev_chrom != chrom) {
                    global_asd[i][j] = global_asd[j][i] += auxilary_asd[i][j];
                    window_asd[i][j] = window_asd[j][i] = auxilary_asd[i][j] = auxilary_asd[j][i] = next_asd[i][j] = next_asd[j][i] = 0;
                }
                auxilary_asd[i][j] += !(hap_size ^ 1) * UNORIENTED_ASD(genotypes[i], genotypes[j]) 
                                    + !!(hap_size ^ 1) * HAPLOTYPE_ASD(auxilary_asd[i][j], ORIENTED_ASD(genotypes[i], genotypes[j]), genotypes[i], genotypes[j]);
                if (num_loci % hap_size == 0) {
                    window_asd[i][j] = window_asd[j][i] += auxilary_asd[i][j];
                    global_asd[i][j] = global_asd[j][i] += auxilary_asd[i][j];
                    if (num_loci <= hap_size * offset_hap_size) {
                        auxilary_asd[j][i] += auxilary_asd[i][j];
                    }
                    auxilary_asd[i][j] = 0;
                }
            }
        }

        if (num_loci == (offset_hap_size * hap_size) + 1) {
            next_start_pos = pos;
        }

        prev_chrom = chrom;
        prev_pos = pos;
    
    }

    windows -> push_back(create_window_with_mds("Global", windows -> front() -> start_position, prev_pos, num_global_haps, n, k, useFastMap, global_asd, d, e));

    delete [] genotypes;
    if (!useFastMap) {
        delete [] d;
        delete [] e;
    }
    for (int i = 0; i < n; i++) {
        delete [] window_asd; delete [] global_asd;
        delete [] auxilary_asd; delete [] next_asd;
    }
    delete [] window_asd; delete [] global_asd; delete [] auxilary_asd; delete [] next_asd;

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