
// File: Main.c
// Date: 
// Version 1: 
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Run LODESTAR analysis.

#include "Interface.h"
#include "LODESTAR.h"
#include <stdio.h>
#include <stdlib.h>

double** open_target_file(char* targetFileName, int N, int K) {
    // Create target matrix.
    double** targetPoints = init_matrix(N, K);

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
            destroy_matrix(targetPoints, N);
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
   
    // We create the VCF parser and haplotype encoder.
    VCFLocusParser_t* parser = init_vcf_locus_parser(lodestarConfig -> inputFileName, lodestarConfig -> maf, lodestarConfig -> afMissing, true);
    HaplotypeEncoder_t* encoder = init_haplotype_encoder(parser -> numSamples);
    
    // Check to make sure k is valid.
    if (lodestarConfig -> k >= parser -> numSamples) {
        fprintf(stderr, "k = %d >= numSamples %d. Exiting!\n", lodestarConfig -> k, parser -> numSamples);
        free(lodestarConfig);
        destroy_vcf_locus_parser(parser);
        return -1;
    }

    // User defined points to perform Procrustes against.
    double** userPoints = NULL;
    double** y = NULL;
    double* y0 = NULL;

    // If target file supplied, make sure it is valid, at least numSamples-by-k
    if (lodestarConfig -> targetFileName != NULL) {
        userPoints = open_target_file(lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
        if (userPoints == NULL) {
            fprintf(stderr, "-y %s is not an %d-by-%d matrix. Exiting!\n", lodestarConfig -> targetFileName, parser -> numSamples, lodestarConfig -> k);
            free(lodestarConfig);
            destroy_vcf_locus_parser(parser);
            return -1;
        }
        // Center and normalize.
        y = init_matrix(parser -> numSamples, lodestarConfig -> k);
        y0 = calloc(lodestarConfig -> k, sizeof(double));
    }
    
    // Partition genome into blocks and calculate IBS within the blocks.
    fprintf(stderr, "Beginning blocking algorithm ...\n");
    BlockList_t* globalList = block_allele_sharing(parser, encoder, encoder -> numSamples, lodestarConfig -> BLOCK_SIZE, lodestarConfig -> HAP_SIZE, lodestarConfig -> dropThreshold, lodestarConfig -> threads);

    // Convert IBS to ASD and compute MDS on global.
    fprintf(stderr, "\nFinished blocking algorithm. Starting genome-wide MDS calculations ...\n");
    double* asd = calloc(PACKED_SIZE(encoder -> numSamples), sizeof(double));
    RealSymEigen_t* eigen = init_real_sym_eigen(encoder -> numSamples);
    for (int i = 0; i < encoder -> numSamples; i++)
        for (int j = i + 1; j < encoder -> numSamples; j++)
            asd[PACKED_INDEX(i, j)] = ibs_to_asd(globalList -> alleleCounts[PACKED_INDEX(i, j)]);
    globalList -> X = init_matrix(encoder -> numSamples, lodestarConfig -> k);
    globalList -> effectRank = compute_classical_mds(eigen, asd, lodestarConfig -> k, globalList -> X);
    globalList -> procrustesT = procrustes_statistic(globalList -> X, NULL, y, y0, eigen, eigen -> N, lodestarConfig -> k, true);
    destroy_real_sym_eigen(eigen);
    free(asd);

    // Convert IBS to ASD and calculate jackknifed procrustes statistic.
    fprintf(stderr, "Finished genome-wide MDS calulations. Starting Procrustes ...\n");
    procrustes(globalList, y, y0, lodestarConfig -> k, lodestarConfig -> threads);

    fprintf(stderr, "\nFinished Procrustes. Writing results to output files...\n");


    // Free used memory.
    if (y0 != NULL)
        free(y0);
    destroy_matrix(y, encoder -> numSamples);
    destroy_block_list(globalList);
    destroy_vcf_locus_parser(parser);
    destroy_haplotype_encoder(encoder);
    destroy_lodestar_config(lodestarConfig);
    fprintf(stderr, "Done!\n");
    return 0;
}