
// File: Main.c
// Date: 
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Run LODESTAR analysis.

#include "Interface.h"
#include "VCFLocusParser.h"
#include <stdio.h>
#include <stdlib.h>

double** open_target_file(char* targetFileName, int N, int K) {
    // Create target matrix.
    double** targetPoints = (double**) calloc(N, sizeof(double*));
    for (int i = 0; i < N; i++)
        targetPoints[i] = (double*) calloc(K, sizeof(double));

    // Read in values from target.
    int numLines = 0;
    size_t length = 0;
    char* line = calloc(1024, sizeof(char));
    FILE* targetFile = fopen(targetFileName, "r");
    while (getline(&line, &length, targetFile) != -1) {
        char* tok = strtok(line, "\t");
        for (int i = 0; tok != NULL && i < K; i++) {
            targetPoints[numLines][i] = atof(tok);
            tok = strtok(NULL, "\t");
        }
        // If invalid formatting, free all memory and exit.
        if (tok != NULL || numLines == N) {
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
    return targetPoints;
}


int main (int argc, char *argv[]) {

    // Get configuration.
    LodestarConfig_t* lodestarConfig = init_lodestar_config(argc, argv);
    if (lodestarConfig == NULL)
        return -1;

    // We create the VCF parser.
    VCFLocusParser_t* parser = init_vcf_locus_parser(lodestarConfig -> inputFileName, lodestarConfig -> maf, lodestarConfig -> afMissing, true);

    // Check to make sure k is valid.
    if (lodestarConfig -> k >= parser -> numSamples) {
        fprintf(stderr, "k = %d >= numSamples %d. Exiting!\n", lodestarConfig -> k, parser -> numSamples);
        free(lodestarConfig);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // User defined points to perform Procrustes against.
    double** userPoints = NULL;

    // If target file supplied, make sure it is valid, at least numSamples-by-k
    if (lodestarConfig -> targetFileName != NULL) {
        userPoints = open_target_file(lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
        if (userPoints == NULL) {
            fprintf(stderr, "-y %s is not an %d-by-%d matrix. Exiting!\n", lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
            free(lodestarConfig);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
    }

}