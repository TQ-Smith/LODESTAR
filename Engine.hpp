//
// File: Engine.hpp
// Started: 17 May 2023
// Author: TQ Smith
// Principle Investigator: Dr. Zachary Szpiech
// Purpose: Header file for Engine.cpp
//

#ifndef ENGINE_H_  
#define ENGINE_H_

// Used for chromosome name and file name.
#include <string>

// Needed for string header.
using namespace std;

// Our structure that defines a window of the genome.
struct Window {
    // The name of the current chromosome.
    string chromosome;

    // The start position of the window.
    int start_position;

    // The end position of the chromosome.
    int end_position;

    // The number of complete loci in the window.
    int num_loci;

    // The set of dimension reduced points.
    double** points;
};

void lodestar_pipeline(string input_file_name, string unit, int window_width, int window_offset, int k);

#endif