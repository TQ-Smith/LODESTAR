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
#include <cassert>

int main(int argc, char *argv[]) {

    // Create command line argument parser.
    CommandLineArgumentParser parser;

    // Used to get the number of arguments given by a flag.
    int num_arguments;

    // Flag used to tell if operation was successful.
    bool successfulOperation;

    // Add all of our options.
    //  Creates variables with defaults set.

    bool help = false;
    parser.addOption("help", "Flag to print help menu.", &successfulOperation);
    assert(successfulOperation);

    string vcf = "";
    parser.addOption("vcf", "<file.vcf> Input is formatted as VCF.", &successfulOperation);
    assert(successfulOperation);

    string out = "";
    parser.addOption("out", "<file> Name of the output file. Prints to terminal on default.", &successfulOperation);
    assert(successfulOperation);

    string unit = "snp";
    parser.addOption("unit", "<bp/snp> The units used to determine the sliding window. Default snp.", &successfulOperation);
    assert(successfulOperation);

    int width;
    int offset;
    parser.addOption("window", "<int> <int> Sets the window width and the offset for the sliding window.", &successfulOperation);
    assert(successfulOperation);

    int k = 2;
    parser.addOption("k", "<1/2/3> The dimension to project into k = 1, 2, or 3. Deafult 2.", &successfulOperation);
    assert(successfulOperation);

    bool fastmap = false;
    parser.addOption("fastmap", "Flag to use FastMap MDS heuristic. Default false.", &successfulOperation);
    assert(successfulOperation);

    // Now, we parse the command line arguments.
    parser.parseCommandLine(argc, argv, &successfulOperation);
    assert(successfulOperation);

    // Get arguments supplied by the user.

    // Test if help flag was set. If set, print help dialogue and exit.
    delete[] parser.getOptionArguments<bool>("--help", &num_arguments, &successfulOperation);
    if (num_arguments == 0) {
        // Print dialogue.
        // Exit.
        return 0;
    }
    if (num_arguments > 0 || !successfulOperation) {
        cout << "--help does not take any arguments." << endl;
        return 0;
    }

    // Parse the VCF input file name.
    string* vcf_buffer = parser.getOptionArguments<string>("--vcf", &num_arguments, &successfulOperation);
    if (num_arguments != 1 || !successfulOperation) {
        cout << "--vcf accepts one VCF file name. Exiting..." << endl;
        delete[] vcf_buffer;
        return 0;
    }
    vcf = vcf_buffer[0];
    delete[] vcf_buffer;

    // Parse the output file name.
    string* out_buffer = parser.getOptionArguments<string>("--out", &num_arguments, &successfulOperation);
    if (num_arguments != 1 || !successfulOperation) {
        cout << "--out accepts one file name. Exiting..." << endl;
        delete[] out_buffer;
        return 0;
    }
    out = out_buffer[0];
    delete[] out_buffer;

    // Parse the unit argument.
    string* unit_buffer = parser.getOptionArguments<string>("--unit", &num_arguments, &successfulOperation);
    if (num_arguments == 1 && (unit_buffer[0] == "bp" || unit_buffer[0] == "snp")) {
        unit = unit_buffer[0];
        delete[] unit_buffer;
    }
    if (unit_buffer != NULL || !successfulOperation) {
        cout << "--unit accepts either \"bp\" or \"snp\". Default \"snp\". Exiting..." << endl;
        delete[] unit_buffer;
        return 0;
    }

    // Parse window arguments..
    int* window_buffer = parser.getOptionArguments<int>("--window", &num_arguments, &successfulOperation);
    if (num_arguments != 2 || !successfulOperation || window_buffer[0] <= 0 || window_buffer[1] <= 0 || window_buffer[1] > window_buffer[0]) {
        cout << "--window accepts two positive integers where the second integer is less than or equal to the first." << endl;
        delete[] window_buffer;
        return 0;
    }
    width = window_buffer[0];
    offset = window_buffer[1];
    delete[] window_buffer;

    // Parse the k argument.
    int* k_buffer = parser.getOptionArguments<int>("-k", &num_arguments, &successfulOperation);
    if (num_arguments == 1 && k >= 1 && k <= 3) {
        k = k_buffer[0];
        delete[] k_buffer;
    }
    if (k_buffer != NULL || !successfulOperation) {
        cout << "-k accepts 1, 2, or 3. Default 2. Exiting..." << endl;
        delete[] k_buffer;
        return 0;
    }

    // Parse the fastmap argument.
    delete[] parser.getOptionArguments<bool>("--fastmap", &num_arguments, &successfulOperation);
    if (num_arguments == 0) {
        fastmap = true;
    }
    if (num_arguments > 0 || !successfulOperation) {
        cout << "--fastmap does not take any arguments." << endl;
        return 0;
    }

}