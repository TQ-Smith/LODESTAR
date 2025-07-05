
// File: Main.c
// Date: 
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Run LODESTAR analysis.

#include "Interface.h"

double** open_target_file(char* targetFileName, int N, int K) {
    // Create target matrix.
    double** targetPoints = (double*) calloc(N, sizeof(double));
    for (int i = 0; i < N; i++)
        targetPoints[i] = calloc(K, sizeof(double));

    // Read in values from target.
    int numLines = 0;
    double inVal = 0;
    size_t length = 0;
    char* line = calloc(1024, sizeof(char));
    FILE* targetFile = fopen(targetFileName, "r");
    while (getline(&line, &length, targetFile) != -1) {
        int numFields = 0;

        // If invalid formatting, free all memory and exit.
        if (numFields != K || numLines == N) {
            fclose(targetFile);
            free(line);
            for (int i = 0; i < N; i++)
                free(targetPoints[i]);
            free(targetPoints);
            return NULL;
        }
        numLines++;
    }
    fclose(targetFile);
    free(line);
    for (int i = 0; i < N; i++)
        free(targetPoints[i]);
    free(targetPoints);
    return targetPoints;
}


int main (int argc, char *argv[]) {

    // Get configuration.
    LodestarConfig_t* lodestarConfig = init_lodestar_config(argc, argv);
    if (lodestarConfig == NULL)
        return -1;

    // We create the VCF parser.
    VCFLocusParser_t* parser = init_vcf_locus_parser(lodestar_config -> inputFileName, lodestar_config -> maf, lodestar_config -> afMissing, true);

    // Check to make sure k is valid.
    if (lodestar_config -> k >= parser -> numSamples) {
        fprintf(stderr, "k = %d >= numSamples %d. Exiting.\n", lodestar_config -> k, parser -> numSamples);
        free(lodestar_config);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // User defined points to perform Procrustes against.
    double** userPoints = NULL;

    // If target file supplied, make sure it is valid, at least numSamples-by-k
    if (lodestar_config -> targetFileName != NULL) {
        userPoints = open_target_file(lodestar_config -> targetFileName, parser -> numSamples, lodestar_config -> k);
        if (userPoints == NULL) {
            fprintf(stderr, "-y %s is not an %d-by-%d matrix. Exiting!\n", lodestar_config -> targetFileName, parser -> numSamples, lodestar_config -> k);
            free(lodestar_config);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
    }

}