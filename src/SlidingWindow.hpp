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

list<window*>* window_genome(VCFParser* parser, int hap_size, int win_size, int step_size, int n, int k, bool useFastMap);

void destroy_window(window** w, int n);