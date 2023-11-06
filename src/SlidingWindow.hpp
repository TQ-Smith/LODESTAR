//
// File: SlidingWindow.hpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines the attributes for a window of the genome.
//              Gives defintion for a sliding window along a genome.
//

// Used for chromosome name.
#include <string>

// Used to read in VCF files.
#include "VCFParser.hpp"

// Used to hold a list of windows.
#include <list>

// Define the contents of a window.
struct window {

    // The chromosome the window is on.
    string chromosome;

    // The base position of the first locus used in the window.
    int start_position;

    // The base position of the last locus used in the window.
    int end_position;

    // The total number of loci used in the window.
    int num_loci;

    // Finally, the set of points used to represent the
    //  individuals in the window.
    double** points;

};

// Our main windowing algorithm. For each window, calculate ASD
//  matrix and perform MDS.
// Accepts:
//  VCFParser* parser -> The VCF file reader.
//  int hap_size -> The number of SNPs in each haplotype.
//  int win_size -> The number of haplotypes in each window.
//  int step_size -> The number of haplotypes in each step.
//  int n -> The number of samples.
//  int k -> The dimension to project into.
//  bool useFastMap -> If set, use FastMap MDS heuristic.
// Returns: list<windows*>*, the windows along the genome. Last one
//              window* is the global window.
list<window*>* window_genome(VCFParser* parser, int hap_size, int win_size, int step_size, int n, int k, bool useFastMap);

// Destroy the memory allocated to a window.
// Accepts:
//  window** w -> A pointer to a pointer to a window. Sets pointer to NULL.
//  int n -> The number of samples.
// Returns: void.
void destroy_window(window** w, int n);