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
    bool doesConverge = true;

    // Used by FastMap to get the maximum dimension reached by the algorithm.
    int maxDimReached = k;

    // Used for classical MDS. This is fast because k = 1, 2, or 3.
    double* d = new double[n];
    double* e = new double[n];

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

    // Ask Zach about maxDimReached and doesConverge?
    if (maxDimReached == k && doesConverge) {
        w -> points = points;
    } else {
        destroy_real_matrix(points, n);
        w -> points = NULL;
    }

    delete [] d;
    delete [] e;

    return w;

}

// A method that deines the behavior of advancing a window.
//  Not needed as a seperate method but improves the readability of the window_genome function.
// Accepts: As defined in window_genome.
// Returns: void.
void advance_window(double** allele_counts, double** local_and_global_window, Genotype** prev_genotypes, int num_loci, int hap_size, int step_size, int win_size, int n) {
    
    // Calculate number of haplotypes in the window.
    int num_haps = ceil(num_loci / (double) hap_size);

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // Add the remaining haplotype's ASD to the local window.
            local_and_global_window[i][j] += allele_counts[i][j];
            if (num_loci != hap_size * win_size) {
                local_and_global_window[j][i] += (local_and_global_window[i][j] - allele_counts[j][i]);
            } else {
                local_and_global_window[j][i] += allele_counts[i][j];
            }
            // Is this the best solution? Why?
            if (step_size * win_size == 1) {
                local_and_global_window[i][j] -= allele_counts[j][i];
                allele_counts[i][j] = local_and_global_window[i][j] / (2.0 * num_haps);
            } else {
                allele_counts[i][j] = local_and_global_window[i][j] / (2.0 * num_haps);
                local_and_global_window[i][j] -= allele_counts[j][i];
            }
            // The allele_counts matrix becomes the window.
            allele_counts[j][i] = allele_counts[i][j];
        }
        // Set diagonal to 0 for MDS.
        allele_counts[i][i] = 0;
    }

}

// A method that defines the behavior of processing a VCF entry.
//  Note needed as a seperate method but improves the readability of the window_genome function.
// Accepts: As defined in window_genome, except newChrom -> boolean to indicate the start of a new chromosome.
// Returns: void.
void calculate_dissimilarity(double** allele_counts, double** local_and_global_window, Genotype** prev_genotypes, Genotype* genotypes, int num_loci, int hap_size, int step_size, bool newChrom, bool newWindow, int n) {

    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            // If we are on a new chromosome, then we zero our haplotype, step, and window counts.
            if (newChrom) {
                local_and_global_window[i][j] = 0;
                allele_counts[i][j] = allele_counts[j][i] = 0;
            }
            // If we are in a new window, copy over leftovers from previous window
            //  and zero haplotype counts.
            if (newWindow) {
                allele_counts[j][i] = local_and_global_window[i][j];
                allele_counts[i][j] = 0;
            }
            // If we are starting a new haplotype..
            if (num_loci % hap_size == 0) {
                // Add counts to offset if in interval.
                if (num_loci <= hap_size * step_size) {
                    allele_counts[j][i] += allele_counts[i][j];
                }
                // Add haplotype to window.
                local_and_global_window[i][j] += allele_counts[i][j];
                local_and_global_window[j][i] += allele_counts[i][j];
                // Set previous genotypes (only important for mini-haplotypes).
                prev_genotypes[i][j] = genotypes[i];
                prev_genotypes[j][i] = genotypes[j];
                // Reset haplotype counts.
                allele_counts[i][j] = 0;
            }
            // If 1 SNP haplotypes or we are on the first locus.
            if (hap_size == 1 || num_loci == 0) {
                // Calculate unoriented dissimilarity.
                allele_counts[i][j] = UNORIENTED_DISSIMILARITY(genotypes[i], genotypes[j]);
            // Otherewise, we are dealing with mini-haplotypes of at least size 2.
            } else {
                // If dissimilarity is not already two, we calculate the updated haplotype dissimilarity.
                if (allele_counts[i][j] != 2) {
                    allele_counts[i][j] = HAPLOTYPE_DISSIMILARITY(LEFT_HAPLOTYPE(prev_genotypes[i][j], genotypes[i]), RIGHT_HAPLOTYPE(prev_genotypes[i][j], genotypes[i]), LEFT_HAPLOTYPE(prev_genotypes[j][i], genotypes[j]), RIGHT_HAPLOTYPE(prev_genotypes[j][i], genotypes[j]));
                }
                // We only reset previous genotypyes if current ones are not identical.
                if (UNORIENTED_DISSIMILARITY(genotypes[i], genotypes[j]) != 0) {
                    prev_genotypes[i][j] = genotypes[i];
                    prev_genotypes[j][i] = genotypes[j];
                }
            }
        }
    }

}

list<window*>* window_genome(VCFParser* parser, int hap_size, int win_size, int step_size, int n, int k, bool useFastMap) {
    
    // Allocate our list of windows.
    list<window*>* windows = new list<window*>;

    // Allocate and fill our three main matrices.
    // allele_counts:
    //  upper triangle -> The dissimilarity counts between samples in the current haplotype.
    //  lower triangle -> The total number of dissimilar haplotypes between samples in the step regions.
    //                      Subtracted from the window's counts when the window is advanced.
    // local_and_global_window:
    //  upper trangle -> The number of dissimilar haplotypes between samples in the current window.
    //  lower triangle -> The number of dissimilar haplotypes between samples in globally.
    // prev_genotypes:
    //  To keep track of the number of dissimilar haplotypes between samples, we have to keep
    //  track of the previous genotypes between pairs of samples.
    //  prev_genotypes[i][j] -> previous genotype of i between i and j.
    //  prev_genotypes[j][i] -> previous genotype of j between i and j.
    //  When a pair of samples share a the same oriented genotpye, it does not
    //  affect the dissimilarity and must be skipped.
    double** allele_counts  = new double*[n];
    double** local_and_global_window = new double*[n];
    Genotype** prev_genotypes  = new Genotype*[n];
    for (int i = 0; i < n; i++) {
        allele_counts[i] = new double[n]; 
        local_and_global_window[i] = new double[n];
        prev_genotypes[i] = new Genotype[n];
        for (int j = 0; j < n; j++) {
            allele_counts[i][j] = 0;
            local_and_global_window[i][j] = 0;
            prev_genotypes[i][j] = 0;
        }
        allele_counts[i][i] = local_and_global_window[i][i] =  prev_genotypes[i][i] = 0;
    }

    // Used for keeping track of the windows.
    // prev_chrom -> The chromosome of the previous entry.
    // start_pos -> The start position of the current window.
    // next_start_pos -> The start position of the next window. Entry immediately after step.
    // prev_pos -> The position of the previous entry.
    // num_loci -> The number of loci in the current window.
    // num_global_haps -> The number of global haplotypes.
    // newWindow -> Set's flag if current entry defines a new window.
    string prev_chrom = "";
    int start_pos, next_start_pos, prev_pos, num_loci = 0, num_global_haps = 0;
    bool newWindow = true;

    // Used to read in a VCF entry.
    string chrom;
    int pos;
    bool isMonomorphic, isComplete, nextRecord;
    Genotype* genotypes = new Genotype[n];

    // Read in the haplotype.
    //  We do this so we don't have to seperately process the last window after the loop.
    while (true) {

        // Get VCF entry. If end of file, sets nextRecord to false.
        nextRecord = parser -> getNextLocus(&chrom, &pos, &isMonomorphic, &isComplete, genotypes);

        if (nextRecord && (isMonomorphic || !isComplete)) {
            continue;
        }

        // If it is not the first locus and all the conditions to indicate a finished window.
        if (num_loci != 0 && (!nextRecord || chrom != prev_chrom || num_loci == (hap_size * win_size))) {

            // Increment the global number of haplotypes by the step size and the additional
            //  haplotypes if EOF was reached or chromosome finished.
            num_global_haps += step_size;
            if (chrom != prev_chrom || !nextRecord) {
                num_global_haps += ceil((num_loci - hap_size * step_size) / (double) hap_size);
            }

            // Advance the window.
            advance_window(allele_counts, local_and_global_window, prev_genotypes, num_loci, hap_size, step_size, win_size, n);

            // cout << "Local and Global:" << endl;
            // print_real_matrix(local_and_global_window, n, n, 2, 2);
            // cout << endl;
            cout << "Window on " << prev_chrom << " from " << start_pos << " to " << prev_pos << endl;
            // print_real_matrix(allele_counts, n, n, 2, 2);
            // cout << endl;

            // Perfrom MDS and add to list.
            windows -> push_back(create_window_with_mds(prev_chrom, start_pos, prev_pos, num_loci, n, k, useFastMap, allele_counts));
            // Determine number of loci and start position of the next window.
            num_loci -= (step_size * hap_size);
            start_pos = (win_size == step_size) ? pos : next_start_pos;

            // Set flag current entry is in new window.
            newWindow = true;

        }

        // If EOF was reached, break.
        if (!nextRecord) {
            break;
        }

        // If entry is on new chromosome,
        //  reset the number of loci and starting postion of the window.
        if (prev_chrom != chrom) {
            num_loci = 0;
            next_start_pos = pos;
            start_pos = pos;
        }

        // Update dissimilarities.
        calculate_dissimilarity(allele_counts, local_and_global_window, prev_genotypes, genotypes, num_loci, hap_size, step_size, prev_chrom != chrom, newWindow, n);

        // cout << endl;
        // cout << "Counts at " << pos << ":" << endl;
        // print_real_matrix(allele_counts, n, n, 2, 2);
        // cout << endl;

        newWindow = false;
        prev_chrom = chrom;
        prev_pos = pos;
        // For sliding window, next window's start is the entry after the step size.
        if (num_loci == step_size * hap_size) {
            next_start_pos = pos;
        }
        // Next locus.
        num_loci++;

    }

    // local_and_global_window is converted to global ASD matrix.
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            local_and_global_window[j][i] = local_and_global_window[j][i] / (2.0 * num_global_haps);
            local_and_global_window[i][j] = local_and_global_window[j][i];
        }
        local_and_global_window[i][i] = 0;
    }

    // cout << "Global Window" << endl;
    // print_real_matrix(local_and_global_window, n, n, 2, 2);
    // cout << endl;

    // Perform MDS on global and add to list.
    //  NOTE: num_loci for global is the number of haplotypes over the genome.
    windows -> push_back(create_window_with_mds("Global", windows -> front() -> start_position, prev_pos, num_global_haps, n, k, useFastMap, local_and_global_window));

    // Free all dynamically allocated memory.
    for (int i = 0; i < n; i++) {
        delete [] allele_counts[i]; 
        delete [] local_and_global_window[i]; 
        delete [] prev_genotypes[i];
    }
    delete [] allele_counts; 
    delete [] local_and_global_window; 
    delete [] prev_genotypes;
    delete [] genotypes;

    // Return our list of windows.
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