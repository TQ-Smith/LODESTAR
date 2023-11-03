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
    string*  str_args = NULL;
    int* int_args = NULL;
    double* double_args = NULL;

    // A list of all command line option variables.
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

    // Parser help option.
    //  If help is set, print option descriptions and exit.
    //  If help was give parameters, print error message and exit.
    if ((num_arguments = cmd_parser.getNumberOfArguments("help", &successfulOperation)) == 0) {
        cout << "LODESTAR: LOcal DEcomposition and Similarity To All Regions" << endl << endl;
        cout << "Options:" << endl;
        cout << "--------" << endl << endl;
        cmd_parser.printOptionDescriptions();
        return 0;
    } else if (num_arguments > 0) {
        cout << "--help does not take any paramters!! Exiting ..." << endl;
        return 0;
    }

    // Parse i option.
    str_args = cmd_parser.getOptionArguments<string>("i", &num_arguments, &successfulOperation);
    if (num_arguments == 1) {
        input_file = str_args[0];
    } else {
        cout << "-i has one required argument! Exiting ..." << endl;
        return 0;
    }
    delete [] str_args;

    // Parse o option.
    str_args = cmd_parser.getOptionArguments<string>("o", &num_arguments, &successfulOperation);
    if (num_arguments == 1) {
        output_file = str_args[0];
    } else {
        cout << "-o has one required argument! Exiting ..." << endl;
        return 0;
    }
    delete [] str_args;

    // Parse format option.
    str_args = cmd_parser.getOptionArguments<string>("format", &num_arguments, &successfulOperation);
    if (num_arguments == 1 && (str_args[0] == "txt" || str_args[0] == "json")) {
        format = str_args[0];
        delete [] str_args;
    } else if (num_arguments != 1) {
        cout << "--format was not given a valid argument! Exiting ..." << endl;
        delete [] str_args;
        return 0;
    }

    // Parse h option.
    int_args = cmd_parser.getOptionArguments<int>("h", &num_arguments, &successfulOperation);
    if (num_arguments == 1 && int_args[0] > 0) {
        h = int_args[0];
        delete [] int_args;
    } else if (num_arguments != -1) {
        cout << "-h accepts one positive integer! Exiting ..." << endl;
        delete [] int_args;
        return 0;
    }

    // Parse w option.
    int_args = cmd_parser.getOptionArguments<int>("w", &num_arguments, &successfulOperation);
    if (num_arguments == 1 && int_args[0] > 0) {
        w = int_args[0];
        delete [] int_args;
    } else if (num_arguments != -1) {
        cout << "-w accepts one positive integer! Exiting ..." << endl;
        delete [] int_args;
        return 0;
    }

    // Parse s option.
    int_args = cmd_parser.getOptionArguments<int>("s", &num_arguments, &successfulOperation);
    if (num_arguments == 1 && int_args[0] > 0 && int_args[0] <= w) {
        s = int_args[0];
        delete [] int_args;
    } else if (num_arguments != -1 || int_args[0] > w || s > w) {
        cout << "-s accepts one positive integer less than or equal to the window size! Exiting ..." << endl;
        delete [] int_args;
        return 0;
    }

}