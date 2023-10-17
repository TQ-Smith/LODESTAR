//
// File: SlidingWindow.cpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines the main algorithm for a sliding window along the genome.
//

#include "SlidingWindow.hpp"

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

list<window*>* sliding_window(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n) {

    // Create our list to hold the windows.
    list<window*>* windows = new list<window*>;

    // Variables needed to read in a locus.
    string chromosome;
    int position; 
    bool isMonomorphic;
    bool isComplete;
    Genotype* genotypes = new Genotype[n];

}

void destroy_window(window** w) {

    // Deallocate points matrix.
    delete [] (*w) -> points;

    // Deallocate structure.
    delete *w;

    // Set pointer to NULL.
    *w = NULL;

}