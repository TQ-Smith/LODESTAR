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

// Used for debugging.
// #include <iostream>
// using namespace std;

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
// Returns: window*, The pointer to the newly allocated window structure.
window* create_window_with_mds(string chrom, int start_pos, int end_pos, int num_loci, int n, int k, bool useFastMap, double** D) {
    
    // Used for classical MDS to test convergence.
    bool doesConverge;

    // Used by FastMap to get the maximum dimension reached by the algorithm.
    int maxDimReached = k;

    // Used for classical MDS. This is fast because k = 1, 2, or 3.
    double d[k];
    double e[k];

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

list<window*>* window_genome(VCFParser* parser, int hap_size, int win_size, int step_size, int n, int k, bool useFastMap) {
    
    list<window*>* windows = new list<window*>;

    double** allele_counts  = new double*[n];
    double** local_and_global_window = new double*[n];
    Genotype** prev_genotypes  = new Genotype*[n];
    for (int i = 0; i <  n; i++) {
        allele_counts[i] = new double[n];
        local_and_global_window[i] = new double[n];
        prev_genotypes[i] = new Genotype[n];
    }

    string prev_chrom = "";
    int start_pos, next_start_pos, prev_pos, num_haps, num_loci = 0, num_global_haps = 0;
    int left_i, right_i, left_j, right_j;
    bool newWindow = true;

    string chrom;
    int pos;
    bool isMonomorphic, isComplete, nextRecord;
    Genotype* genotypes = new Genotype[n];

    while (true) {

        nextRecord = parser -> getNextLocus(&chrom, &pos, &isMonomorphic, &isComplete, genotypes);

        if (nextRecord && (isMonomorphic || !isComplete)) {
            continue;
        }

        if (num_loci != 0 && (!nextRecord || chrom != prev_chrom || num_loci == (hap_size * win_size))) {

            num_haps = ceil(num_loci / (double) hap_size);
            num_global_haps += step_size;
            if (chrom != prev_chrom || !nextRecord) {
                num_global_haps += ceil((num_loci - hap_size * step_size) / (double) hap_size);
            }

            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    allele_counts[i][i] = local_and_global_window[i][j];
                    if (hap_size > 1 && num_loci == 1) {
                        allele_counts[i][j] = UNORIENTED_ASD(prev_genotypes[i][j], prev_genotypes[j][i]);
                    }
                    local_and_global_window[i][j] -= allele_counts[j][i];
                    if ( hap_size != 1) {
                        local_and_global_window[j][i] += allele_counts[i][j];
                    }
                    if (!nextRecord || prev_chrom != chrom) {
                        local_and_global_window[j][i] += local_and_global_window[i][j];
                    }
                    allele_counts[i][j] = (allele_counts[i][j] + allele_counts[i][i]) / (2 * num_haps);
                    allele_counts[j][i] = allele_counts[i][j];
                }
                allele_counts[i][i] = 0;
            }

            // cout << "Window on " << prev_chrom << " from " << start_pos << " to " << prev_pos << endl;
            // print_real_matrix(allele_counts, n, n, 2, 2);
            // cout << endl;

            windows -> push_back(create_window_with_mds(prev_chrom, start_pos, prev_pos, num_loci, n, k, useFastMap, allele_counts));

            num_loci -= step_size * hap_size;
            start_pos = (win_size == step_size) ? pos : next_start_pos;
            newWindow = true;

        }

        if (!nextRecord) {
            break;
        }

        if (prev_chrom != chrom) {
            num_loci = 0;
            next_start_pos = pos;
            start_pos = pos;
        }

        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                if (chrom != prev_chrom) {
                    local_and_global_window[i][j] = 0;
                }
                if (newWindow) {
                    allele_counts[i][j] = allele_counts[j][i] = 0;
                }
                if (num_loci % hap_size == 0) {
                    if (hap_size == 1) {
                        allele_counts[i][j] = UNORIENTED_ASD(genotypes[i], genotypes[j]);
                    }
                    if (num_loci < hap_size * step_size) {
                        allele_counts[j][i] += allele_counts[i][j];
                    }
                    local_and_global_window[i][j] += allele_counts[i][j];
                    local_and_global_window[j][i] += allele_counts[i][j];
                    prev_genotypes[i][j] = genotypes[i];
                    prev_genotypes[j][i] = genotypes[j];
                    allele_counts[i][j] = 0;
                }
                if (hap_size != 1 && num_loci > 0 && allele_counts[i][j] != 2  && ORIENTED_ASD(genotypes[i], genotypes[j]) != 0) {
                    left_i = LEFT_HAPLOTYPE(prev_genotypes[i][j], genotypes[i]); 
                    right_i = RIGHT_HAPLOTYPE(prev_genotypes[i][j], genotypes[i]);
                    left_j = LEFT_HAPLOTYPE(prev_genotypes[j][i], genotypes[j]); 
                    right_j = RIGHT_HAPLOTYPE(prev_genotypes[j][i], genotypes[j]);
                    allele_counts[i][j] = HAPLOTYPE_ASD(left_i, right_i, left_j, right_j);
                    prev_genotypes[i][j] = genotypes[i];
                    prev_genotypes[j][i] = genotypes[j];
                }
            }
        }

        newWindow = false;
        prev_chrom = chrom;
        prev_pos = pos;
        if (num_loci == step_size * hap_size) {
            next_start_pos = pos;
        }
        num_loci++;

    }

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            local_and_global_window[j][i] = local_and_global_window[j][i] / (2 * num_global_haps);
            local_and_global_window[i][j] = local_and_global_window[j][i];
        }
        local_and_global_window[i][i] = 0;
    }

    // cout << "Global Windows" << endl;
    // print_real_matrix(local_and_global_window, n, n, 2, 2);
    // cout << endl;

    windows -> push_back(create_window_with_mds("Global", windows -> front() -> start_position, prev_pos, num_global_haps, n, k, useFastMap, local_and_global_window));

    for (int i = 0; i < n; i++) {
        delete [] allele_counts[i];
        delete [] local_and_global_window[i];
        delete [] prev_genotypes[i];
    }
    delete [] allele_counts;
    delete [] local_and_global_window;
    delete [] prev_genotypes;

    delete [] genotypes;

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