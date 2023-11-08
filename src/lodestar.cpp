//
// File: lodestar.cpp
// Date: 3 November 2023
// Author: TQ Smith
// Principal Investigator: Dr. Zachary Szpiech
// Purpose: Main executable of LODESTAR.
//

// Used for parsing user command line arguments.
#include "../utils/CommandLineArgumentParser.hpp"

// Used for creating and destroying matrices.
#include "../utils/LinearAlgebra/MatrixOperations.hpp"

// Used for sliding window algorithm.
#include "SlidingWindow.hpp"

// Used for Procrustes analysis.
#include "Procrustes.hpp"

// Manipulate strings.
#include <string>

// Used for IO and printing.
#include <iostream>
using namespace std;

// Used for setting width.
#include <iomanip>

// Used to set nan.
#include <limits>

// Print the help menu.
// Accepts:
//  CommandLineArgumentParser* cmd_parser -> The CLI parser.
// Returns: void.
void print_help(CommandLineArgumentParser* cmd_parser) {

    // Print header.
    cout << "LODESTAR: LOcal DEcomposition and Similarity To All Regions" << endl;
    cout << "-----------------------------------------------------------" << endl;
    cout << endl;
    cout << "Options:" << endl;

    bool successfulOperation;

    // Get each option in alphabetical order and print description.
    string* options = cmd_parser -> getOptions();
    for (int i = 0; i < cmd_parser -> getNumberOfOptions(); i++) {
        cout << setw(15) << options[i] << " ";
        cout << setw(60) << cmd_parser -> getDescription(options[i], &successfulOperation) << endl;
    }
    delete [] options;

}

int main(int argc, char *argv[]) {

    ////////////////////////////// FIRST WE PARSE COMMAND LINE ARGUMENTS //////////////////////////

    // Create command line parser.
    CommandLineArgumentParser cmd_parser;
    bool successfulOperation;
    int num_arguments;
    string* str_args = NULL;
    int* int_args = NULL;

    // A list of all command line option variables.
    //  Default values are immediately set.
    string input_file;
    string output_file;
    string format = "txt";
    int h = 1;
    int w = 1;
    int s = 1;
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

    // Parse the --help option.
    switch (cmd_parser.getNumberOfArguments("--help", &successfulOperation)) {
        // If set, print menu and exit.
        case 0:
            print_help(&cmd_parser);
            return 0;
            break;
        // Do nothing if option is not set.
        case -1: break;
        // Print error if arguments were given.
        default:
            cout << "--help does not take any arguments! Exiting ..." << endl;
            break;
    };

    // Parse -i option.
    switch(cmd_parser.getNumberOfArguments("-i", &successfulOperation)) {
        // If one argument was given, set file name.
        case 1:
            str_args = cmd_parser.getOptionArguments<string>("-i", &num_arguments, &successfulOperation);
            input_file = str_args[0];
            delete [] str_args;
            break;
        // Otherwise, print error and exit.
        default:
            cout << "-i must be given one argument! Exiting ..." << endl;
            return 0;
            break;
    };
    /*
    // Parse -o option.
    switch(cmd_parser.getNumberOfArguments("-o", &successfulOperation)) {
        // If one argument was given, set file name.
        case 1:
            str_args = cmd_parser.getOptionArguments<string>("-o", &num_arguments, &successfulOperation);
            output_file = str_args[0];
            delete [] str_args;
            break;
        // Otherwise, print error and exit.
        default:
            cout << "-o must be given one argument! Exiting ..." << endl;
            return 0;
            break;
    };
    */

    // Parse --format option.
    switch (cmd_parser.getNumberOfArguments("--format", &successfulOperation)) {
        // If one argument make sure it it "txt" or "json".
        case 1:
            str_args = cmd_parser.getOptionArguments<string>("--format", &num_arguments, &successfulOperation);
            if (str_args[0] != "txt" || str_args[0] != "json" || !successfulOperation) {
                cout << "--format must be given \"txt\" or \"json\"! Exiting ..." << endl;
                delete [] str_args;
                return 0;
            }
            format = str_args[0];
            delete [] str_args;
            break;
        // Use default
        case -1: break;
        default:
            cout << "--format must be given \"txt\" or \"json\"! Exiting ..." << endl;
            return 0;
            break;
    };

    // Parse -h option.
    switch (cmd_parser.getNumberOfArguments("-h", &successfulOperation)) {
        case 1:
            int_args = cmd_parser.getOptionArguments<int>("-h", &num_arguments, &successfulOperation);
            if (int_args[0] <= 0 || !successfulOperation) {
                cout << "-h must be given a positive integer! Exiting ..." << endl;
                delete [] int_args;
                return 0;
            }
            h = int_args[0];
            delete [] int_args;
            break;
        // Use default
        case -1: break;
        default:
            cout << "-h must be give one positive integer! Exiting ..." << endl;
            return 0;
            break;
    };

    // Parse -w option.
    switch (cmd_parser.getNumberOfArguments("-w", &successfulOperation)) {
        case 1:
            int_args = cmd_parser.getOptionArguments<int>("-w", &num_arguments, &successfulOperation);
            if (int_args[0] <= 0 || !successfulOperation) {
                cout << "-w must be given a positive integer! Exiting ..." << endl;
                delete [] int_args;
                return 0;
            }
            w = int_args[0];
            delete [] int_args;
            break;
        // Use default
        case -1: break;
        default:
            cout << "-w must be give one positive integer! Exiting ..." << endl;
            return 0;
            break;
    };

    // Parse -s option.
    switch (cmd_parser.getNumberOfArguments("-s", &successfulOperation)) {
        case 1:
            int_args = cmd_parser.getOptionArguments<int>("-s", &num_arguments, &successfulOperation);
            if (int_args[0] <= 0 || int_args[0] > w || !successfulOperation) {
                cout << "-s must be given a positive integer less than or equal to window width! Exiting ..." << endl;
                delete [] int_args;
                return 0;
            }
            s = int_args[0];
            delete [] int_args;
            break;
        // Use default
        case -1: break;
        default:
            cout << "-s must be give one positive integer! Exiting ..." << endl;
            return 0;
            break;
    };

    // Parse -k option.
    switch (cmd_parser.getNumberOfArguments("-k", &successfulOperation)) {
        case 1:
            int_args = cmd_parser.getOptionArguments<int>("-k", &num_arguments, &successfulOperation);
            if (int_args[0] < 1 || int_args[0] > 3 || !successfulOperation) {
                cout << "-k must be 1, 2, or 3! Exiting ..." << endl;
                delete [] int_args;
                return 0;
            }
            k = int_args[0];
            delete [] int_args;
            break;
        // Use default
        case -1: break;
        default:
            cout << "-k must be 1, 2, or 3! Exiting ..." << endl;
            return 0;
            break;
    };

    // Parse the --fastmap option.
    switch (cmd_parser.getNumberOfArguments("--fastmap", &successfulOperation)) {
        // If set, print menu and exit.
        case 0:
            fastmap = true;
            break;
        // Do nothing if option is not set.
        case -1: break;
        // Print error if arguments were given.
        default:
            cout << "--fastmap does not take any arguments! Exiting ..." << endl;
            break;
    };

    // Parse --num_perms option.
    switch (cmd_parser.getNumberOfArguments("--num_perms", &successfulOperation)) {
        case 1:
            int_args = cmd_parser.getOptionArguments<int>("--num_perms", &num_arguments, &successfulOperation);
            if (int_args[0] <= 0 || !successfulOperation) {
                cout << "-num_perms must be a positive integer! Exiting ..." << endl;
                delete [] int_args;
                return 0;
            }
            num_perms = int_args[0];
            delete [] int_args;
            break;
        // Use default
        case -1: break;
        default:
            cout << "--num_perms must be a positive integer! Exiting ..." << endl;
            return 0;
            break;
    };

    //////////////////////////// FINISHED PARSING COMMAND LINE ARGUMENTS ////////////////////////////



    ////////////////////////////// NOW WE START THE LODESTAR PIPELINE ///////////////////////////////

    // Open VCF file.
    VCFParser* vcf_parser = new VCFParser(input_file);
    if (!(vcf_parser -> isOpen())) {
        cout << "Input file does not exist!" << endl;
        delete vcf_parser;
        return 0;
    }

    // Get the number of samples.
    int n = vcf_parser -> getNumberOfSamples();

    // Now, we window the genome.
    list<window*>* windows = window_genome(vcf_parser, h, w, s, n, k, fastmap);

    // Get the global window.
    window* global = windows -> back();
    windows -> pop_back();

    // Allocate memory for Procrustes analysis.
    double* x_0 = new double[n];
    double** C = create_real_matrix(n, n);
    double** CT_C = create_real_matrix(n, n);
    double** shuffleX = create_real_matrix(n, n);

    // Hold dissimilarity statistic for Procrustes.
    double statistic;

    // Hold p-value for permutation test.
    double p_value;

    // Center global points.
    center_matrix(global -> points, x_0, n, k);

    // Perform Procrustes analysis for each window.
    window* current_window;
    int index = 0;
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        current_window = *it;

        // Center points.
        center_matrix(current_window -> points, x_0, n, k);

        // Get Procrustes dissimilarity statistic.
        statistic = procrustes_analysis(current_window -> points, global -> points, C, CT_C, n, k);

        // Execute permutation test.
        p_value = permutation_test(num_perms, current_window -> points, global -> points, shuffleX, C, CT_C, n, k, statistic);

        // Set values for window.
        current_window -> statistic = statistic;
        current_window -> p_value = p_value;
        
        current_window -> index = index;

        index++;
    }

    // Free memory used for Procrustes analysis.
    delete [] x_0;
    destroy_real_matrix(C, n);
    destroy_real_matrix(CT_C, n);
    destroy_real_matrix(shuffleX, n);

    // Free our vcf_parser. We are done reading the input file.
    delete vcf_parser;

    ////////////////////////////////// FINISHED LODESTAR PIPELINE ///////////////////////////////////



    /////////////////////////////////// OUTPUT LODESTAR RESULTS /////////////////////////////////////

    // We could output as we iterate but this is more flexible.
    cout << endl;
    for (list<window*>::iterator it = windows -> begin(); it != windows -> end(); it++) {
        current_window = *it;

        cout << "Window: " << (current_window -> index) << endl;
        cout << "Chromosome: " << (current_window -> chromosome) << endl;
        cout << "Start Position: " << (current_window -> start_position) << endl;
        cout << "End Position: " << (current_window -> end_position) << endl;
        cout << "Number of Loci: " << (current_window -> num_loci) << endl;
        cout << "Statistic: " << (current_window -> statistic) << endl;
        cout << "P-Value: " << (current_window -> p_value) << endl;
        cout << "Points:" << endl;
        print_real_matrix(current_window -> points, n, k, 4, 4);
        cout << endl;
        cout << endl;
    }

    cout << "Global" << endl;
    cout << "Number of Haplotypes: " << (global -> num_loci) << endl;
    cout << "Points:" << endl;
    print_real_matrix(global -> points, n, k, 4, 4);
    cout << endl;

    // We are done with the windows.
    delete global;
    delete windows;

    /////////////////////////////////// FINISHED LODESTAR OUTPUT ////////////////////////////////////

}