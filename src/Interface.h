
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

// Used to index the lower triangle of a packed matrix.
#define INDEX(i, j, N) (i <= j ? i + j * (j + 1) / 2 : j + i * (i + 1) / 2)

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
    // Name of the file assigning each sample to a group.
    char* groupFileName;
    // Pvalue threshold used to print points.
    double pthresh;
    // Number of permutations to execute in permutation test.
    int NUM_PERMS;
    // t-statistic threshold used to print points.
    double tthresh;
    // Minor allele frequency threshhold.
    double maf;
    // Missing allele frequency threshold.
    double afMissing;
    // Max possible gap between adjacent loci.
    int MAX_GAP;
    // Output file in JSON format instead of text file.
    bool useJsonOutput;
} LodestarConfiguration_t;

// Open a tsv file to get the group of each sample.
// Accepts:
//  char* groupFileName -> The tsv file to open.
//  kstring_t** sampleNames -> The names of the samples.
//  int N -> The number of samples.
//  int* numGroups -> Sets the number of groups.
// Returns: int*, an array assigning each sample a group number, NULL if invalid file.
int* open_group_file(char* groupFileName, kstring_t** sampleNames, int N, int* numGroups);

// Open a tsv file for target points.
// Accepts:
//  char* targetFileName -> The name of the file containing the target points.
//  int N -> The number of samples.
//  int K -> The reduced dimension.
// Returns: double**, the matrix of points; otherwise, NULL if file has less than N * K entries.
double** open_target_file(char* targetFileName, int N, int K); 

// Print verbose LODESTAR configuration.
// Accepts:
//  FILE* output -> The output stream to print the configuration.
//  LodestarConfiguration_t* lodestart_config -> THe configuration to print.
// Returns: void.
void print_configuration(FILE* output, LodestarConfiguration_t* lodestar_config);

// Print the single entry fields for a window.
// Accepts:
//  FILE* output -> The output stream to print to.
//  Window_t* window -> The window whose fields will be printed.
// Returns: void.
void print_window_summary(FILE* output, Window_t* window);

// Print the attributes of a window.
//  Do not like how this method is written, but I do not think there is a way around it.
// Accepts:
//  FILE* output -> The output stream to print to.
//  kstring_t** sampleNames -> Array of sample names corresponding to each point in the window.
//              Only used when useLongOutput is set.
//  Window_t* window -> The window whose attributes will be printed.
//  int N -> The number of samples.
//  int K -> The dimension of the projected points.
//  bool useJsonOutput -> Print the attributes in JSON format. Otherwise, print in plain text.
//  bool printCoords -> Print the coordinates of the window if they exist.
// Returns: void.
void print_window(FILE* output, kstring_t** sampleNames, Window_t* window, int N, int K, bool useJsonOutput, bool printCoords);

// Parse commandline options and return LODESTAR configuration.
// Accepts:
//  int argc, char *argv[] -> Our command line arguments.
// Returns: LodestarConfiguration_t*, The configuration or NULL if error while parsing command line options.
LodestarConfiguration_t* init_lodestar_config(int argc, char *argv[]);

#endif