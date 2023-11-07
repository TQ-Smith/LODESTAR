//
// File: lodestar.cpp
// Date: 3 November 2023
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

int main(int argc, char *argv[]) {

    // Create command line parser.
    CommandLineArgumentParser cmd_parser;
    bool successfulOperation;
    int num_arguments;
    string* str_args = NULL;
    int* int_args = NULL;
    double* double_args = NULL;

    // A list of all command line option variables.
    //  Default values are immediately set.
    string input_file;
    string output_file;
    string format = "txt";
    int h = 10;
    int w = 10;
    int s = 10;
    int k = 2;
    bool fastmap = false;
    int num_perms  = 10000;

    // Add the options.
    cmd_parser.addOption("help", "Display command line option descriptions.", &successfulOperation);
    cmd_parser.addOption("i", "The input file in .vcf or .vcf.gz format. Required.",  &successfulOperation);
    cmd_parser.addOption("o", "The output file name.", &successfulOperation);
    cmd_parser.addOption("format", "The output file format \"txt\" or \"json\". Default \"txt\".",  &successfulOperation);
    cmd_parser.addOption("h", "The haplotype size in SNPs.", &successfulOperation);
    cmd_parser.addOption("w", "The window size in number of haplotypes.",  &successfulOperation);
    cmd_parser.addOption("s", "The step size in number of haplotypes.", &successfulOperation);
    cmd_parser.addOption("k", "The dimension to project into. Default 2.", &successfulOperation);
    cmd_parser.addOption("fastmap", "Flag to use the fastmap heuristic.", &successfulOperation);
    cmd_parser.addOption("num_perms", "The number of permutations to execute. Default 10000.", &successfulOperation);

    // Parse the command line.
    cmd_parser.parseCommandLine(argc, argv, &successfulOperation);

}