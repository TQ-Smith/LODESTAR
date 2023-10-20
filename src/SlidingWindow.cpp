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

// We define a macro to calculate the number of different alleles between two diploid individuals
//  at a single locus with the "Genotype" encoding from VCFParser.
#define ASD(i, j) (!!(i ^ j) + !(i & j))

#include <iostream>

// Allocates a window with the specified attributes.
// Accepts:
//  string chromosome -> The chromosome the window is on.
//  int start_position -> The base position of the first locus used in the window.
//  int end_position -> The base position of the last locus used in the window.
//  int num_loci -> The total number of loci used in the window.
//  double** points -> The set of points to represent the individuals.
// Returns: window*, The allocated window.
window* create_window(string chromosome, int start_position, int end_position, int num_loci, double** points) {

    window* w = new window;

    // Set fields.
    w -> chromosome = chromosome;
    w -> start_position = start_position;
    w -> end_position = end_position;
    w -> num_loci = num_loci;
    w -> points = points;

    return w;

}

list<window*>* sliding_window(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n, int k, bool useFastMap) {

    // Create our list to hold the windows.
    list<window*>* windows = new list<window*>;

    // Variables needed to read in a locus from the VCF parser.
    string chromosome;
    int position; 
    bool isMonomorphic, isComplete;
    Genotype* genotypes = new Genotype[n];

    // Variables needed to keep track of window information.
    string previous_chromosome = "";
    int previous_position, start_position, num_current_loci = 0;

    // Allocate matrices for ASD calculations.
    //  window_asd holds the ASD matrix over the whole window.
    //  auxilary_asd holds two pieces of information. In the top triangle,
    //      it holds the counts for the current haplotype. In the lower triangle,
    //      it holds the counts for the current offset.
    //  newWindow holds ASD values for the nextWindow
    double** window_asd = create_and_fill_real_matrix(0, n, n);
    double** auxilary_asd = create_and_fill_real_matrix(0, n, n);
    double** newWindow = create_real_matrix(n, n);

    // Allocate memory used for classicalMDS.
    double* d = new double[n];
    double* e = new double[n];
    bool doesConverge = true;
    int maxDimReached;

    // Read in the loci from the VCF.
    while(parser -> getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes)) {
        
        // If the site is monomorphic or incomplete, then skip over the site.
        if (isMonomorphic || !isComplete) {
            continue;
        }

        // If new chromosome or finished window, then preform MDS, save window, and reset values for new window.
        //  Make sure we aren't in the first window of a new chromosome.
        if (num_current_loci != 0 && chromosome != previous_chromosome || (hap_size * window_hap_size) == num_current_loci) {

            int num_haps = (num_current_loci % hap_size == 0) ? (num_current_loci / hap_size) : ((num_current_loci / hap_size) + 1);
            
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    newWindow[i][j] = newWindow[j][i] = window_asd[i][j] - auxilary_asd[j][i];
                    auxilary_asd[j][i] = 0;
                    window_asd[i][j] = window_asd[j][i] += (auxilary_asd[i][j] >= 2 ? 2 : auxilary_asd[i][j]) / 2.0;
                    window_asd[i][j] = window_asd[j][i] = window_asd[i][j] / num_haps;
                }
            }

            double** X = create_real_matrix(n, k);

            if (useFastMap) {
                compute_fastmap(window_asd, X, &maxDimReached, n, k);
            } else {
                compute_classical_mds(window_asd, X, d, e, &doesConverge, n, k);
            }
            
            if (doesConverge) {
                windows -> push_back(create_window(previous_chromosome, start_position, previous_position, num_current_loci, X));
            }

            // New window becomes the current window.
            double** temp = window_asd;
            window_asd = newWindow;
            newWindow = temp; 

            num_current_loci = num_current_loci - hap_size * (num_current_loci / hap_size);

        }

        if (num_current_loci == 0 || num_current_loci == hap_size * offset_hap_size) {
            start_position = position;
        }

        // Here, we handle the IBS counts within the haplotype, window, and offset.
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                // If we are not in the first haplotype of a chromosme and we have ended the current haplotype on the previous locus.
                if (num_current_loci != 0 && num_current_loci % hap_size == 1) {
                    // Convert the different allele counts to ASD.
                    auxilary_asd[i][j] = (auxilary_asd[i][j] >= 2 ? 2 : auxilary_asd[i][j]) / 2.0;
                    // Add haplotype ASD to the offset matrix.
                    if (num_current_loci <= (hap_size * offset_hap_size)) {
                        auxilary_asd[j][i] += auxilary_asd[i][j];
                    }
                    // Add the haplotype ASD to the window matrix.
                    window_asd[i][j] = window_asd[j][i] += auxilary_asd[i][j];
                    // We are starting a new haplotype. Set counter to 0.
                    auxilary_asd[i][j] = 0;
                }
                // Count the number of different alleles at the locus.
                auxilary_asd[i][j] += ASD(genotypes[i], genotypes[j]);
            }
        }

        num_current_loci++;
        previous_chromosome = chromosome;
        previous_position = position;

    }
    windows -> push_back(create_window(previous_chromosome, start_position, previous_position, num_current_loci, NULL));

    // Destroy our ASD matrices.
    destroy_real_matrix(window_asd, n);
    destroy_real_matrix(newWindow, n);
    destroy_real_matrix(auxilary_asd, n);

    // Free allocated genotypes array.
    delete [] genotypes;

    // Delete vectors used for classical MDS.
    delete [] d;
    delete [] e;
    
    // Return the list of windows.
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