
// File: Interface.h
// Date: 18 February 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines LODESTAR configuration from command line interface and
//              output formatting of LODESTAR results.

#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "Window.h"
#include "VCFLocusParser.h"
#include "HaplotypeEncoder.h"
#include "../lib/kstring.h"

// Defines all the command line options supplied
//  by user and parameters to run the LODESTAR analysis.
typedef struct {
    // The input VCF file.
    char* inputFileName;
    // The basename for output files.
    char* outputBasename;
    // The haplotype size in number of loci.
    int HAP_SIZE;
    // The windows size in number of haplotypes.
    int WINDOW_SIZE;
    // The step size in number of haplotypes.
    int STEP_SIZE;
    // The dimension of the reduced points.
    int k;
    // The number of threads to use.
    int threads;
    // Flag to indicate if Procrustes analysis is measuring
    //  similarity or dissimilarity.
    bool similarity;
    // Flag to indicate if we are just calculating global window.
    bool global;
    // Name of file for user specified points to perform Procrustes
    //  analysis with.
    char* targetFileName;
    // Number of permutations to execute in permutation test.
    int NUM_PERMS;
    // Minor allele frequency threshhold.
    double maf;
    // Missing allele frequency threshold.
    double afMissing;
    // Max possible gap between adjacent loci.
    int MAX_GAP;
} LodestarConfiguration_t;

// Prints the configuration of the LODESTAR in JSON format.
// Accepts:
//  FILE* output -> The output stream to print to.
//  LodestarConfiguration_t* lodestar_config -> The LODESTAR configuration.
// Returns: void.
void print_configuration(FILE* output, LodestarConfiguration_t* lodestar_config);

// Open a tsv file for target points.
// Accepts:
//  char* targetFileName -> The name of the file containing the target points.
//  int N -> The number of samples.
//  int K -> The reduced dimension.
// Returns: double**, the matrix of points; otherwise, NULL if file has less than N * K entries.
double** open_target_file(char* targetFileName, int N, int K); 

// Print the single entry fields for a window.
// Accepts:
//  FILE* output -> The output stream to print to.
//  Window_t* window -> The window whose fields will be printed.
// Returns: void.
void print_window_summary(FILE* output, Window_t* window);

// Print the attributes of a window in JSON format.
// Accepts:
//  FILE* output -> The output stream to print to.
//  kstring_t** sampleNames -> Array of sample names corresponding to each point in the window.
//              Only used when useLongOutput is set.
//  Window_t* window -> The window whose attributes will be printed.
//  int N -> The number of samples.
//  int K -> The dimension of the projected points.
// Returns: void.
void print_window(FILE* output, kstring_t** sampleNames, Window_t* window, int N, int K);

// Parse commandline options and return LODESTAR configuration.
// Accepts:
//  int argc, char *argv[] -> Our command line arguments.
// Returns: LodestarConfiguration_t*, The configuration or NULL if error while parsing command line options.
LodestarConfiguration_t* init_lodestar_config(int argc, char *argv[]);

#endif