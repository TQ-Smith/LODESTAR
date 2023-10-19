//
// File: SlidingWindow.hpp
// Date: 17 October 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Defines the attributes for a window of the genome.
//              Gives defintion for a sliding window along a genome.
//

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

// Our sliding window algorithm. Computes ASD matrices and preforms MDS.
// Accepts:
//  VCFParser* vcf_reader -> Read in the VCF file.
//  int hap_size -> The number of SNPs in a haplotype.
//  int window_hap_size -> The number of haplotypes in a window.
//  int offset_hap_size -> The number of haplotypes in the offset.
//  int n -> The number of samples in the VCF file.
//  int k -> The dimension to project the data into. Assumes k = 1, 2, 3.
//  bool useFastMap -> If true, preform MDS using FastMap heuristic instead of classical MDS.
// Returns: list<window*>, The list of windows created in the genome.
list<window*>* sliding_window(VCFParser* parser, int hap_size, int window_hap_size, int offset_hap_size, int n, int k, bool useFastMap);

// Deallocates a window.
// Accepts:
//  window** w -> A reference to the pointer to deallocate.
//                  Sets pointer to NULL.
// Returns: void.
void destroy_window(window** w);