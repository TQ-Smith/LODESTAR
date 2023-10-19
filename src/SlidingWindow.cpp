//
// File: SlidingWindow.cpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines the main algorithm for a sliding window along the genome.
//

#include "SlidingWindow.hpp"

// Use classical MDS and FastMap.
#include "MultidimensionalScaling.hpp"

// Used to create and destroy matrices.
#include "MatrixOperations.hpp"

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
    bool isMonomorphic;
    bool isComplete;
    Genotype* genotypes = new Genotype[n];

    // Variables needed to keep track of window information.
    string previous_chromosome = "";
    int previous_position;
    int start_position;
    int num_current_loci = 0;

    // Read in the loci from the VCF.
    while(parser -> getNextLocus(&chromosome, &position, &isMonomorphic, &isComplete, genotypes)) {
        
        // If the site is monomorphic or incomplete, then skip over the site.
        if (isMonomorphic || !isComplete) {
            continue;
        }

        // If new chromosome or finished window, then preform MDS, save window, and reset values for new window.
        if (chromosome != previous_chromosome || (hap_size * window_hap_size) == num_current_loci) {

            // We have to create the first window.
            if (num_current_loci == 0) {
            } else {
                windows -> push_back(create_window(previous_chromosome, start_position, previous_position, num_current_loci, NULL));
                num_current_loci = 0;
            }

            start_position = position;

        }

        num_current_loci++;
        previous_chromosome = chromosome;
        previous_position = position;

    }
    windows -> push_back(create_window(previous_chromosome, start_position, previous_position, num_current_loci, NULL));


    // Free allocated genotypes array.
    delete [] genotypes;
    
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