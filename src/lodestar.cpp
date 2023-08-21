//
// File: lodestar.cpp
// Date: 21 August 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Main executable of LODESTAR.
//

// Used for parsing user command line arguments.
#include "../utils/CommandLineArgumentParser.hpp"

// Manipulate strings.
#include <string>

// Used for IO and printing.
#include <iostream>
using namespace std;

// Assert used for CommandLineArgumentParser.
#include <assert>

int main(int argc, char *argv[]) {

    // Create command line argument parser.
    CommandLineArgumentParser parser;

    // Used to get the number of arguments given by a flag.
    int num_arguments;

    // Flag used to tell if operation was successful.
    bool successfulOperation;

    // Add all of our options.
    parser.addOption("help", "Flag to print help menu.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("vcf", "<file.vcf> Input is formatted as VCF.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("out", "<file> Name of the output file. Prints to terminal on default.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("unit", "<bp/snp> The units used to determine the sliding window.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("window", "<int> <int> Sets the window width and the offset for the sliding window.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("k", "<1/2/3> The dimension to project into k = 1, 2, or 3.", &successfulOperation);
    assert(successfulOperation);
    parser.addOption("fastmap", "Flag to use FastMap MDS heuristic. Deafult false.", &successfulOperation);
    assert(successfulOperation);
    
}