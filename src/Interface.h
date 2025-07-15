
// File: Interface.h
// Date: 18 February 2025
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Defines LODESTAR configuration from command line interface and
//              output formatting of LODESTAR results.

#ifndef _INTERFACE_H_
#define _INTERFACE_H_

#include "kstring.h"

// Defines all the command line options supplied
//  by user and parameters to run the LODESTAR analysis.
typedef struct {
    // The input VCF file.
    char* inputFileName;
    // The basename for output files.
    char* outputBasename;
    // The haplotype size in number of loci.
    int HAP_SIZE;
    // The block size in base pairs.
    int BLOCK_SIZE;
    // The dimension of the reduced points.
    int k;
    // The number of threads to use.
    int threads;
    // Name of file for user specified points to perform Procrustes
    //  analysis with.
    char* targetFileName;
    // Minor allele frequency threshhold.
    double maf;
    // Missing allele frequency threshold.
    double afMissing;
    // Drop block if there are less than threshold number of haps.
    int dropThreshold;
    // The full command.
    char* cmd;
} LodestarConfig_t;

// Parse commandline options and return LODESTAR configuration.
// Accepts:
//  int argc, char *argv[] -> Our command line arguments.
// Returns: LodestarConfig_t*, The configuration or NULL if error while parsing command line options.
LodestarConfig_t* init_lodestar_config(int argc, char *argv[]);

// Free the memory associated with LODESTAR.
void destroy_lodestar_config(LodestarConfig_t* config);

#endif